"""Hierarchical Clustering with Cophenetic Analysis for Umbrella Reviews.

Build a dendrogram of reviews based on Jaccard distance of study overlap,
compute cophenetic correlation, optimal cluster count (gap statistic),
and ultrametric violation score.
"""

from __future__ import annotations
import math
import random
from itertools import combinations


# ── Distance helpers ────────────────────────────────────────────────

def _jaccard_distance(a: set, b: set) -> float:
    """1 - Jaccard similarity.  Returns 1.0 if both empty."""
    if not a and not b:
        return 1.0
    inter = len(a & b)
    union = len(a | b)
    return 1.0 - inter / union


def _build_distance_matrix(reviews) -> tuple[list[list[float]], list[str]]:
    """Pairwise Jaccard distance matrix and ordered list of review IDs."""
    n = len(reviews)
    ids = [r.review_id for r in reviews]
    sets = [set(r.study_ids) for r in reviews]
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            d = _jaccard_distance(sets[i], sets[j])
            dist[i][j] = d
            dist[j][i] = d
    return dist, ids


# ── UPGMA agglomerative clustering ─────────────────────────────────

def _upgma(dist: list[list[float]], ids: list[str]):
    """Average-linkage (UPGMA) agglomerative clustering.

    Returns:
        dendrogram  — list of merge dicts
        cophenetic  — cophenetic distance matrix (n x n)
    """
    n = len(ids)
    # Active cluster -> list of original indices
    clusters: dict[int, list[int]] = {i: [i] for i in range(n)}
    # Working distance: keyed by (min_id, max_id)
    d: dict[tuple[int, int], float] = {}
    for i in range(n):
        for j in range(i + 1, n):
            d[(i, j)] = dist[i][j]

    next_id = n
    dendrogram = []
    cophenetic = [[0.0] * n for _ in range(n)]

    active = set(range(n))

    for _ in range(n - 1):
        # Find closest pair among active clusters
        best_pair = None
        best_dist = float("inf")
        for ci, cj in combinations(sorted(active), 2):
            key = (min(ci, cj), max(ci, cj))
            dd = d.get(key, float("inf"))
            if dd < best_dist:
                best_dist = dd
                best_pair = (ci, cj)

        if best_pair is None:
            break

        ci, cj = best_pair
        members_i = clusters[ci]
        members_j = clusters[cj]

        # Record merge
        dendrogram.append({
            "cluster_a": ci,
            "cluster_b": cj,
            "distance": best_dist,
            "new_cluster": next_id,
        })

        # Fill cophenetic distances for all pairs across merging clusters
        for a in members_i:
            for b in members_j:
                cophenetic[a][b] = best_dist
                cophenetic[b][a] = best_dist

        # Create new cluster
        new_members = members_i + members_j
        clusters[next_id] = new_members

        # Compute distances from new cluster to all remaining active clusters
        active.discard(ci)
        active.discard(cj)
        for ck in sorted(active):
            # Average linkage: mean of all pairwise distances
            total = 0.0
            count = 0
            for a in new_members:
                for b in clusters[ck]:
                    total += dist[a][b]
                    count += 1
            avg = total / count if count > 0 else 0.0
            key = (min(next_id, ck), max(next_id, ck))
            d[key] = avg

        active.add(next_id)
        next_id += 1

    return dendrogram, cophenetic


# ── Cophenetic correlation ──────────────────────────────────────────

def _pearson(x: list[float], y: list[float]) -> float:
    """Pearson correlation between two lists."""
    n = len(x)
    if n < 2:
        return 0.0
    mx = sum(x) / n
    my = sum(y) / n
    sxx = sum((xi - mx) ** 2 for xi in x)
    syy = sum((yi - my) ** 2 for yi in y)
    sxy = sum((xi - mx) * (yi - my) for xi, yi in zip(x, y))
    denom = math.sqrt(sxx * syy)
    if denom < 1e-15:
        return 0.0
    return sxy / denom


def _cophenetic_correlation(dist: list[list[float]],
                            coph: list[list[float]]) -> float:
    """Pearson correlation of upper-triangle entries."""
    n = len(dist)
    orig = []
    coph_vals = []
    for i in range(n):
        for j in range(i + 1, n):
            orig.append(dist[i][j])
            coph_vals.append(coph[i][j])
    return _pearson(orig, coph_vals)


# ── Cluster cutting ────────────────────────────────────────────────

def _cut_dendrogram(dendrogram, n_reviews: int, k: int) -> dict[int, int]:
    """Cut dendrogram to produce exactly k clusters.

    Returns mapping: original_index -> cluster_label (0..k-1).
    """
    if k >= n_reviews:
        return {i: i for i in range(n_reviews)}

    # Sort merges by distance ascending
    merges = sorted(dendrogram, key=lambda m: m["distance"])

    # Union-Find
    parent: dict[int, int] = {i: i for i in range(n_reviews)}

    # Also track which original indices are in each cluster
    members: dict[int, list[int]] = {i: [i] for i in range(n_reviews)}

    n_clusters = n_reviews
    for m in merges:
        if n_clusters <= k:
            break
        a, b, new = m["cluster_a"], m["cluster_b"], m["new_cluster"]

        # Find roots
        def root(x):
            while parent.get(x, x) != x:
                x = parent[x]
            return x

        ra, rb = root(a), root(b)
        if ra == rb:
            continue
        # Merge into new cluster id
        parent[ra] = new
        parent[rb] = new
        parent[new] = new
        members[new] = members.pop(ra, [ra]) + members.pop(rb, [rb])
        n_clusters -= 1

    # Assign labels
    assignments: dict[int, int] = {}
    label = 0
    for cid, mems in members.items():
        if parent.get(cid, cid) == cid:  # root cluster
            for idx in mems:
                assignments[idx] = label
            label += 1

    return assignments


# ── Gap statistic ───────────────────────────────────────────────────

def _wcss(dist: list[list[float]], assignments: dict[int, int]) -> float:
    """Within-cluster sum of squared distances."""
    clusters: dict[int, list[int]] = {}
    for idx, lab in assignments.items():
        clusters.setdefault(lab, []).append(idx)
    total = 0.0
    for members in clusters.values():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                total += dist[members[i]][members[j]] ** 2
    return total


def _optimal_k(dist: list[list[float]], dendrogram, n_reviews: int,
               max_k: int = 5, n_ref: int = 20, seed: int = 42) -> int:
    """Gap statistic for optimal k."""
    rng = random.Random(seed)
    max_k = min(max_k, n_reviews - 1)
    if max_k < 2:
        return 1

    # Observed WCSS for each k
    ks = list(range(2, max_k + 1))
    obs_log_w = []
    for k in ks:
        asgn = _cut_dendrogram(dendrogram, n_reviews, k)
        w = _wcss(dist, asgn)
        obs_log_w.append(math.log(w + 1e-15))

    # Reference: uniform random distances in [min_d, max_d]
    flat = []
    for i in range(n_reviews):
        for j in range(i + 1, n_reviews):
            flat.append(dist[i][j])
    lo = min(flat) if flat else 0.0
    hi = max(flat) if flat else 1.0

    ref_log_w = [0.0] * len(ks)
    for _ in range(n_ref):
        # Generate random distance matrix
        rdist = [[0.0] * n_reviews for _ in range(n_reviews)]
        for i in range(n_reviews):
            for j in range(i + 1, n_reviews):
                d = rng.uniform(lo, hi)
                rdist[i][j] = d
                rdist[j][i] = d
        # Build dendrogram for reference
        ref_ids = [str(i) for i in range(n_reviews)]
        ref_dend, _ = _upgma(rdist, ref_ids)
        for ki, k in enumerate(ks):
            asgn = _cut_dendrogram(ref_dend, n_reviews, k)
            w = _wcss(rdist, asgn)
            ref_log_w[ki] += math.log(w + 1e-15) / n_ref

    # Gap = E[log(W_ref)] - log(W_obs)
    gaps = [ref_log_w[i] - obs_log_w[i] for i in range(len(ks))]
    best_idx = max(range(len(gaps)), key=lambda i: gaps[i])
    return ks[best_idx]


# ── Ultrametric violation score ─────────────────────────────────────

def _ultrametric_violations(coph: list[list[float]]) -> float:
    """Proportion of triplets violating the ultrametric inequality.

    Ultrametric: d(i,k) <= max(d(i,j), d(j,k)) for all i,j,k.
    For a valid cophenetic matrix from UPGMA, violations should be 0.
    """
    n = len(coph)
    if n < 3:
        return 0.0
    total = 0
    violations = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                total += 3  # 3 orientations
                dij = coph[i][j]
                dik = coph[i][k]
                djk = coph[j][k]
                if dik > max(dij, djk) + 1e-12:
                    violations += 1
                if dij > max(dik, djk) + 1e-12:
                    violations += 1
                if djk > max(dij, dik) + 1e-12:
                    violations += 1
    return violations / total if total > 0 else 0.0


# ── Public API ──────────────────────────────────────────────────────

def cophenetic_analysis(reviews) -> dict:
    """Hierarchical clustering with cophenetic analysis.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with:
        dendrogram             — list of merge dicts
        cophenetic_correlation — float in [-1, 1]
        optimal_k              — int
        cluster_assignments    — dict[review_id, int]
        cluster_summaries      — list of {cluster, mean_theta, n_reviews}
        ultrametric_score      — float in [0, 1] (proportion of violations)
    """
    n = len(reviews)
    if n < 2:
        return {
            "dendrogram": [],
            "cophenetic_correlation": 1.0,
            "optimal_k": 1,
            "cluster_assignments": {reviews[0].review_id: 0} if reviews else {},
            "cluster_summaries": [{
                "cluster": 0,
                "mean_theta": reviews[0].theta if reviews else 0.0,
                "n_reviews": n,
            }] if reviews else [],
            "ultrametric_score": 0.0,
        }

    dist, ids = _build_distance_matrix(reviews)
    dendrogram, cophenetic_mat = _upgma(dist, ids)

    coph_corr = _cophenetic_correlation(dist, cophenetic_mat)

    opt_k = _optimal_k(dist, dendrogram, n)

    # Get cluster assignments at optimal k
    idx_assignments = _cut_dendrogram(dendrogram, n, opt_k)

    # Map back to review IDs
    id_assignments = {ids[i]: lab for i, lab in idx_assignments.items()}

    # Cluster summaries
    cluster_members: dict[int, list] = {}
    for i, lab in idx_assignments.items():
        cluster_members.setdefault(lab, []).append(i)

    summaries = []
    for lab in sorted(cluster_members):
        members = cluster_members[lab]
        thetas = [reviews[i].theta for i in members]
        summaries.append({
            "cluster": lab,
            "mean_theta": sum(thetas) / len(thetas),
            "n_reviews": len(members),
        })

    ultrametric = _ultrametric_violations(cophenetic_mat)

    return {
        "dendrogram": dendrogram,
        "cophenetic_correlation": coph_corr,
        "optimal_k": opt_k,
        "cluster_assignments": id_assignments,
        "cluster_summaries": summaries,
        "ultrametric_score": ultrametric,
    }
