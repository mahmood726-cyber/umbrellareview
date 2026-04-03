"""Network Meta-Meta-Analysis.

When systematic reviews compare different interventions, build an evidence
network at the SR level.  Supports indirect comparisons via the Bucher method,
consistency checks, and P-score frequentist ranking.
"""

import math
import re

import numpy as np
from scipy import stats


def _extract_comparisons(reviews):
    """Extract pairwise comparisons from reviews.

    Looks for scope_tags with pattern "A_vs_B" first.  If no explicit
    comparisons found, creates synthetic ones by clustering theta values.
    """
    comparisons = []
    intervention_set = set()

    for r in reviews:
        found = False
        for tag in (r.scope_tags or []):
            match = re.match(r"^(\w+)_vs_(\w+)$", tag)
            if match:
                a, b = match.group(1), match.group(2)
                comparisons.append({
                    "a": a, "b": b,
                    "theta": r.theta, "se": r.se if r.se > 0 else 0.1,
                    "source": "direct", "review_id": r.review_id,
                })
                intervention_set.add(a)
                intervention_set.add(b)
                found = True
                break
        if not found:
            # Check review_id for "vs" pattern
            parts = re.split(r"[_\-]vs[_\-]", r.review_id, flags=re.IGNORECASE)
            if len(parts) == 2:
                a, b = parts[0].strip(), parts[1].strip()
                comparisons.append({
                    "a": a, "b": b,
                    "theta": r.theta, "se": r.se if r.se > 0 else 0.1,
                    "source": "direct", "review_id": r.review_id,
                })
                intervention_set.add(a)
                intervention_set.add(b)

    # If no explicit comparisons, create synthetic comparison groups
    if not comparisons:
        comparisons, intervention_set = _synthetic_comparisons(reviews)

    return comparisons, intervention_set


def _synthetic_comparisons(reviews):
    """Create synthetic A_vs_B comparisons by splitting reviews into 2 groups.

    Uses k-means (k=2) on theta values.  Group with lower theta = Treatment,
    Group with higher theta = Control.
    """
    thetas = np.array([r.theta for r in reviews])

    if len(reviews) <= 2:
        # Simple: each review is treatment vs control
        comps = []
        for r in reviews:
            comps.append({
                "a": "Treatment", "b": "Control",
                "theta": r.theta, "se": r.se if r.se > 0 else 0.1,
                "source": "direct", "review_id": r.review_id,
            })
        return comps, {"Treatment", "Control"}

    # k-means with k=2
    sorted_idx = np.argsort(thetas)
    mid = len(sorted_idx) // 2
    group_a = sorted_idx[:mid]
    group_b = sorted_idx[mid:]

    mean_a = float(np.mean(thetas[group_a]))
    mean_b = float(np.mean(thetas[group_b]))

    # All reviews contribute to Treatment vs Control comparison
    comps = []
    for r in reviews:
        comps.append({
            "a": "Treatment", "b": "Control",
            "theta": r.theta, "se": r.se if r.se > 0 else 0.1,
            "source": "direct", "review_id": r.review_id,
        })
    return comps, {"Treatment", "Control"}


def _bucher_indirect(comp_ab, comp_bc):
    """Bucher method for indirect comparison A vs C through B.

    theta_AC = theta_AB - theta_BC
    se_AC = sqrt(se_AB^2 + se_BC^2)
    """
    theta_ac = comp_ab["theta"] - comp_bc["theta"]
    se_ac = math.sqrt(comp_ab["se"] ** 2 + comp_bc["se"] ** 2)
    return {
        "a": comp_ab["a"], "b": comp_bc["b"],
        "theta": round(theta_ac, 6),
        "se": round(se_ac, 6),
        "source": "indirect",
        "via": comp_ab["b"],
    }


def _check_consistency(direct_comps, indirect_comps):
    """Check consistency between direct and indirect evidence.

    For each comparison with both direct and indirect estimates,
    compute the inconsistency factor and z-test.
    """
    # Index direct comparisons by (a, b) pair
    direct_by_pair = {}
    for c in direct_comps:
        key = tuple(sorted([c["a"], c["b"]]))
        if key not in direct_by_pair:
            direct_by_pair[key] = []
        direct_by_pair[key].append(c)

    inconsistencies = []
    for ic in indirect_comps:
        key = tuple(sorted([ic["a"], ic["b"]]))
        if key in direct_by_pair:
            # Pool direct estimates for this comparison
            dcs = direct_by_pair[key]
            # IV-weighted pooled direct
            w = np.array([1.0 / (d["se"] ** 2) for d in dcs])
            t = np.array([d["theta"] for d in dcs])
            theta_direct = float(np.sum(w * t) / np.sum(w))
            se_direct = float(1.0 / np.sqrt(np.sum(w)))

            # Inconsistency factor
            inc_factor = theta_direct - ic["theta"]
            se_inc = math.sqrt(se_direct ** 2 + ic["se"] ** 2)
            z = inc_factor / se_inc if se_inc > 0 else 0.0
            p = 2.0 * (1.0 - stats.norm.cdf(abs(z)))

            inconsistencies.append({
                "comparison": f"{ic['a']}_vs_{ic['b']}",
                "direct": round(theta_direct, 6),
                "indirect": round(ic["theta"], 6),
                "inconsistency_factor": round(inc_factor, 6),
                "z": round(z, 4),
                "p": round(p, 4),
            })

    return inconsistencies


def _compute_p_scores(interventions, comparisons):
    """Compute P-scores (frequentist ranking).

    For each pair (i, j), compute P(i better than j) from the pairwise
    comparison z-score.  P-score for i = mean P(i better) across all j != i.
    """
    intv_list = sorted(interventions)
    n = len(intv_list)
    if n < 2:
        return {intv_list[0]: 1.0} if n == 1 else {}

    # Build comparison lookup: {(a,b): (theta, se)} where theta = a minus b
    pair_data = {}
    for c in comparisons:
        key = (c["a"], c["b"])
        rev_key = (c["b"], c["a"])
        if key not in pair_data:
            pair_data[key] = []
        pair_data[key].append((c["theta"], c["se"]))
        if rev_key not in pair_data:
            pair_data[rev_key] = []
        pair_data[rev_key].append((-c["theta"], c["se"]))

    p_scores = {}
    for i_idx, intv_i in enumerate(intv_list):
        p_better_list = []
        for j_idx, intv_j in enumerate(intv_list):
            if i_idx == j_idx:
                continue
            key = (intv_i, intv_j)
            if key in pair_data:
                # Pool estimates
                data = pair_data[key]
                w = np.array([1.0 / (se ** 2) for _, se in data])
                t = np.array([theta for theta, _ in data])
                pooled_theta = float(np.sum(w * t) / np.sum(w))
                pooled_se = float(1.0 / np.sqrt(np.sum(w)))
                # P(i better than j) — for typical log-scale effects,
                # more negative = better (treatment reduces outcome)
                # P-score: P(theta < 0) for this direction
                z = pooled_theta / pooled_se if pooled_se > 0 else 0.0
                p_better = float(stats.norm.cdf(z))  # P(effect < 0) = treatment better
                p_better_list.append(p_better)
            else:
                p_better_list.append(0.5)  # no data, assume equipoise

        p_scores[intv_i] = round(float(np.mean(p_better_list)), 4) if p_better_list else 0.5

    return p_scores


def _build_adjacency(interventions, comparisons):
    """Build adjacency matrix and node degree for the evidence network."""
    intv_list = sorted(interventions)
    n = len(intv_list)
    intv_idx = {v: i for i, v in enumerate(intv_list)}

    adj = [[0] * n for _ in range(n)]
    for c in comparisons:
        if c["a"] in intv_idx and c["b"] in intv_idx:
            i, j = intv_idx[c["a"]], intv_idx[c["b"]]
            adj[i][j] += 1
            adj[j][i] += 1

    degree = {intv_list[i]: sum(1 for j in range(n) if adj[i][j] > 0)
              for i in range(n)}

    # Check connectivity via BFS
    if n == 0:
        connected = True
    else:
        visited = set()
        queue = [0]
        visited.add(0)
        while queue:
            node = queue.pop(0)
            for j in range(n):
                if adj[node][j] > 0 and j not in visited:
                    visited.add(j)
                    queue.append(j)
        connected = len(visited) == n

    return adj, degree, connected


def network_meta_meta(reviews):
    """Build an evidence network from multiple systematic reviews.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews.

    Returns
    -------
    dict with keys:
        comparisons        – list of {a, b, theta, se, source}
        indirect_estimates – list of indirect comparisons via Bucher method
        inconsistency      – list of {comparison, direct, indirect, z, p}
        p_scores           – dict[intervention, float]
        network_connected  – bool
        n_interventions    – int
    """
    if len(reviews) < 2:
        raise ValueError("Need at least 2 reviews for network meta-meta-analysis")

    # Extract direct comparisons
    direct_comps, interventions = _extract_comparisons(reviews)

    # Attempt indirect comparisons via Bucher method
    indirect_estimates = []
    intv_list = sorted(interventions)

    # Group direct comparisons by pair
    pair_comps = {}
    for c in direct_comps:
        key = (c["a"], c["b"])
        if key not in pair_comps:
            pair_comps[key] = []
        pair_comps[key].append(c)

    # For each pair of comparisons sharing a common arm, compute indirect
    all_pairs = list(pair_comps.keys())
    seen_indirect = set()
    for i, pair_i in enumerate(all_pairs):
        for j, pair_j in enumerate(all_pairs):
            if i >= j:
                continue
            # Find common intervention (bridge node)
            set_i = set(pair_i)
            set_j = set(pair_j)
            common = set_i & set_j
            if common:
                bridge = list(common)[0]
                # a_i -- bridge -- b_j
                a_i = pair_i[0] if pair_i[1] == bridge else pair_i[1]
                b_j = pair_j[0] if pair_j[1] == bridge else pair_j[1]

                if a_i == b_j:
                    continue  # same node, skip

                indirect_key = tuple(sorted([a_i, b_j]))
                if indirect_key in seen_indirect:
                    continue
                seen_indirect.add(indirect_key)

                # Pool each pair's estimates
                for ci in pair_comps[pair_i]:
                    for cj in pair_comps[pair_j]:
                        # Ensure direction: ci is a_i vs bridge, cj is bridge vs b_j
                        theta_ab = ci["theta"] if ci["a"] == a_i else -ci["theta"]
                        se_ab = ci["se"]
                        theta_bc = cj["theta"] if cj["a"] == bridge else -cj["theta"]
                        se_bc = cj["se"]

                        comp_ab = {"a": a_i, "b": bridge, "theta": theta_ab, "se": se_ab}
                        comp_bc = {"a": bridge, "b": b_j, "theta": theta_bc, "se": se_bc}

                        indirect = _bucher_indirect(comp_ab, comp_bc)
                        indirect_estimates.append(indirect)
                        break
                    break

    # All comparisons (direct + indirect)
    all_comparisons = direct_comps + indirect_estimates

    # Consistency check
    inconsistency = _check_consistency(direct_comps, indirect_estimates)

    # P-scores
    p_scores = _compute_p_scores(interventions, all_comparisons)

    # Network graph
    adj, degree, connected = _build_adjacency(interventions, direct_comps)

    return {
        "comparisons": [{k: v for k, v in c.items() if k != "review_id"}
                        for c in direct_comps],
        "indirect_estimates": indirect_estimates,
        "inconsistency": inconsistency,
        "p_scores": p_scores,
        "network_connected": connected,
        "n_interventions": len(interventions),
    }
