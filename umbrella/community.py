"""Louvain community detection on evidence networks (review overlap)."""

import math
import numpy as np


def _jaccard(set_a, set_b):
    """Jaccard similarity between two sets."""
    if not set_a and not set_b:
        return 0.0
    intersection = len(set_a & set_b)
    union = len(set_a | set_b)
    return intersection / union if union > 0 else 0.0


def _build_graph(reviews):
    """Build weighted adjacency from Jaccard similarity of study_ids.

    Returns
    -------
    ids : list[str]
        Review IDs in order.
    adj : dict[int, dict[int, float]]
        Adjacency list: adj[i][j] = weight.
    m : float
        Total edge weight (sum of all weights / 2).
    """
    n = len(reviews)
    ids = [r.review_id for r in reviews]
    sets = [set(r.study_ids) for r in reviews]

    adj = {i: {} for i in range(n)}
    total_weight = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            w = _jaccard(sets[i], sets[j])
            if w > 0:
                adj[i][j] = w
                adj[j][i] = w
                total_weight += w

    return ids, adj, total_weight


def _degree(adj, node):
    """Sum of weights incident to node."""
    return sum(adj[node].values())


def _louvain_phase1(adj, n, m, community, gamma=1.0):
    """Local move phase: greedily move nodes to best neighbor community."""
    if m <= 0:
        return community, False

    improved = True
    any_move = False
    while improved:
        improved = False
        for i in range(n):
            current_comm = community[i]
            k_i = _degree(adj, i)

            # Compute weights from i to each neighboring community
            comm_weights = {}
            for j, w in adj[i].items():
                c = community[j]
                comm_weights[c] = comm_weights.get(c, 0.0) + w

            # Compute sum_tot for current community (excluding i)
            # and sum_in for current community (excluding i's contribution)
            best_comm = current_comm
            best_delta = 0.0

            # Remove i from its community: loss
            k_i_in_current = comm_weights.get(current_comm, 0.0)
            sum_tot_current = sum(
                _degree(adj, node) for node in range(n)
                if community[node] == current_comm and node != i
            )

            for target_comm, k_i_in in comm_weights.items():
                if target_comm == current_comm:
                    continue
                sum_tot_target = sum(
                    _degree(adj, node) for node in range(n)
                    if community[node] == target_comm
                )

                # delta_Q for moving i from current_comm to target_comm
                delta = (
                    (k_i_in - k_i_in_current) / m
                    - gamma * k_i * (sum_tot_target - sum_tot_current) / (2.0 * m * m)
                )

                if delta > best_delta:
                    best_delta = delta
                    best_comm = target_comm

            if best_comm != current_comm:
                community[i] = best_comm
                improved = True
                any_move = True

    return community, any_move


def _compute_modularity(adj, n, community, m, gamma=1.0):
    """Compute modularity Q_gamma = sum_c [e_cc - gamma * a_c^2]."""
    if m <= 0:
        return 0.0

    comms = set(community.values()) if isinstance(community, dict) else set(community)
    Q = 0.0
    for c in comms:
        members = [i for i in range(n) if community[i] == c]
        # e_cc: fraction of total weight within community
        e_cc = 0.0
        for i in members:
            for j in adj[i]:
                if community[j] == c and j > i:
                    e_cc += adj[i][j]
        e_cc /= m  # fraction

        # a_c: fraction of total degree in community
        a_c = sum(_degree(adj, i) for i in members) / (2.0 * m)
        Q += e_cc - gamma * a_c * a_c

    return float(Q)


def _louvain(adj, n, m, gamma=1.0):
    """Full Louvain: phase1 local moves + phase2 aggregation, repeated."""
    community = list(range(n))  # each node in its own community

    community, moved = _louvain_phase1(adj, n, m, community, gamma)
    if not moved:
        return community

    # Phase 2: aggregate into super-nodes and repeat (one level of aggregation)
    comm_ids = sorted(set(community))
    comm_map = {c: idx for idx, c in enumerate(comm_ids)}
    n_super = len(comm_ids)

    if n_super == n:
        return community  # no improvement

    # Build super-graph
    super_adj = {i: {} for i in range(n_super)}
    for i in range(n):
        ci = comm_map[community[i]]
        for j, w in adj[i].items():
            cj = comm_map[community[j]]
            if ci != cj:
                super_adj[ci][cj] = super_adj[ci].get(cj, 0.0) + w

    # Normalize (edges were double-counted)
    for i in range(n_super):
        for j in super_adj[i]:
            super_adj[i][j] /= 2.0
    # Rebuild symmetric
    clean_adj = {i: {} for i in range(n_super)}
    for i in range(n_super):
        for j, w in super_adj[i].items():
            if w > 0:
                clean_adj[i][j] = w
                clean_adj[j][i] = w

    super_m = sum(sum(v for v in clean_adj[i].values()) for i in range(n_super)) / 2.0
    if super_m > 0 and n_super > 1:
        super_community = list(range(n_super))
        super_community, _ = _louvain_phase1(
            clean_adj, n_super, super_m, super_community, gamma
        )
        # Map back to original nodes
        final = [0] * n
        for i in range(n):
            final[i] = super_community[comm_map[community[i]]]
        return final

    return community


def _nmi(labels_a, labels_b):
    """Normalized Mutual Information between two label assignments."""
    n = len(labels_a)
    if n == 0:
        return 0.0

    # Contingency table
    classes_a = sorted(set(labels_a))
    classes_b = sorted(set(labels_b))
    map_a = {c: i for i, c in enumerate(classes_a)}
    map_b = {c: i for i, c in enumerate(classes_b)}
    na = len(classes_a)
    nb = len(classes_b)

    contingency = np.zeros((na, nb))
    for i in range(n):
        contingency[map_a[labels_a[i]]][map_b[labels_b[i]]] += 1

    # Marginals
    row_sums = contingency.sum(axis=1)
    col_sums = contingency.sum(axis=0)

    # Entropies
    def entropy(counts, total):
        h = 0.0
        for c in counts:
            if c > 0:
                p = c / total
                h -= p * math.log(p)
        return h

    H_A = entropy(row_sums, n)
    H_B = entropy(col_sums, n)

    # Mutual information
    mi = 0.0
    for i in range(na):
        for j in range(nb):
            if contingency[i][j] > 0:
                p_ij = contingency[i][j] / n
                p_i = row_sums[i] / n
                p_j = col_sums[j] / n
                mi += p_ij * math.log(p_ij / (p_i * p_j))

    denom = H_A + H_B
    if denom < 1e-12:
        return 1.0  # both are trivial partitions
    return 2.0 * mi / denom


def community_detection(reviews, resolutions=None):
    """Detect communities in the evidence overlap network using Louvain.

    Parameters
    ----------
    reviews : list[ReviewInput]
        Reviews with study_ids for overlap computation.
    resolutions : list[float] or None
        Resolution parameters gamma to scan. Default: [0.5, 1.0, 1.5, 2.0].

    Returns
    -------
    dict with keys: communities, modularity, n_communities, community_meta,
                    resolution_scan, nmi_stability
    """
    if resolutions is None:
        resolutions = [0.5, 1.0, 1.5, 2.0]

    n = len(reviews)
    if n < 2:
        single_comm = [[r.review_id for r in reviews]] if reviews else []
        return {
            "communities": single_comm,
            "modularity": 0.0,
            "n_communities": len(single_comm),
            "community_meta": {},
            "resolution_scan": [],
            "nmi_stability": 1.0,
        }

    ids, adj, m = _build_graph(reviews)

    # Run Louvain at multiple resolutions
    resolution_results = []
    all_labels = []
    for gamma in resolutions:
        if m > 0:
            labels = _louvain(adj, n, m, gamma)
        else:
            labels = list(range(n))  # no edges, each node is its own community
        Q = _compute_modularity(adj, n, labels, m, gamma) if m > 0 else 0.0
        n_comms = len(set(labels))
        resolution_results.append({
            "gamma": gamma,
            "n_communities": n_comms,
            "modularity": round(Q, 4),
        })
        all_labels.append(labels)

    # Default result: gamma=1.0 (index 1 if available, else first)
    default_idx = next((i for i, g in enumerate(resolutions) if g == 1.0), 0)
    default_labels = all_labels[default_idx]
    default_Q = resolution_results[default_idx]["modularity"]

    # Build community lists
    comm_groups = {}
    for i, lbl in enumerate(default_labels):
        comm_groups.setdefault(lbl, []).append(ids[i])
    communities = list(comm_groups.values())

    # Community-level meta-analysis (IV pooling within each community)
    community_meta = {}
    for comm_id, (lbl, members) in enumerate(comm_groups.items()):
        member_indices = [i for i in range(n) if default_labels[i] == lbl]
        thetas = np.array([reviews[i].theta for i in member_indices])
        ses = np.array([reviews[i].se for i in member_indices])
        ses = np.where(ses > 0, ses, 0.1)

        w = 1.0 / (ses ** 2)
        pooled_theta = float(np.sum(w * thetas) / np.sum(w))
        pooled_se = float(1.0 / np.sqrt(np.sum(w)))

        community_meta[comm_id] = {
            "theta": round(pooled_theta, 4),
            "se": round(pooled_se, 4),
            "n_reviews": len(member_indices),
        }

    # NMI stability: average NMI between adjacent resolutions
    nmi_values = []
    for i in range(len(all_labels) - 1):
        nmi_val = _nmi(all_labels[i], all_labels[i + 1])
        nmi_values.append(nmi_val)
    nmi_stability = float(np.mean(nmi_values)) if nmi_values else 1.0

    return {
        "communities": communities,
        "modularity": default_Q,
        "n_communities": len(communities),
        "community_meta": community_meta,
        "resolution_scan": resolution_results,
        "nmi_stability": round(nmi_stability, 4),
    }
