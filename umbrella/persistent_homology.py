"""Topological Data Analysis of Evidence Structure.

Apply persistent homology to the review similarity space defined
by Jaccard distance on study overlap.  Computes Betti curves,
persistence diagrams, and persistence entropy for the Vietoris-Rips
filtration built from the pairwise Jaccard distance matrix.
"""

import math

# ---------------------------------------------------------------------------
# Union-Find (disjoint-set) with path compression + union by rank
# ---------------------------------------------------------------------------

class _UnionFind:
    """Minimal union-find for tracking connected components."""

    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n
        self.n_components = n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # path halving
            x = self.parent[x]
        return x

    def union(self, x, y):
        """Unite sets containing x and y.  Return True if a merge happened."""
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False
        if self.rank[rx] < self.rank[ry]:
            rx, ry = ry, rx
        self.parent[ry] = rx
        if self.rank[rx] == self.rank[ry]:
            self.rank[rx] += 1
        self.n_components -= 1
        return True


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def persistent_homology(reviews, n_steps=50):
    """Persistent homology of the Vietoris-Rips filtration on Jaccard distance.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews with study_ids.
    n_steps : int
        Number of filtration steps from 0 to 1 (inclusive endpoints).

    Returns
    -------
    dict with keys:
        distance_matrix   – list[list[float]], Jaccard distance
        betti_curve_0     – list[int], beta_0 at each epsilon
        betti_curve_1     – list[int], beta_1 at each epsilon
        persistence_diagram – list[dict] with {dim, birth, death}
        persistence_entropy – float
        total_persistence   – float
        filtration_values   – list[float]
    """
    n = len(reviews)
    if n < 2:
        raise ValueError("Need at least 2 reviews for persistent homology")

    # ------------------------------------------------------------------
    # 1. Jaccard distance matrix
    # ------------------------------------------------------------------
    sets = [set(r.study_ids) for r in reviews]
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            inter = len(sets[i] & sets[j])
            union = len(sets[i] | sets[j])
            jacc = inter / union if union > 0 else 0.0
            d = 1.0 - jacc
            dist[i][j] = d
            dist[j][i] = d

    # ------------------------------------------------------------------
    # 2. Sort all edges by distance (for filtration)
    # ------------------------------------------------------------------
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((dist[i][j], i, j))
    edges.sort()

    # Filtration thresholds
    filtration = [i / (n_steps - 1) for i in range(n_steps)]

    # ------------------------------------------------------------------
    # 3. Walk through filtration and record Betti numbers + persistence
    # ------------------------------------------------------------------
    betti_0 = []
    betti_1 = []
    persistence_diagram = []

    # Track component birth/merge for beta_0 persistence
    # Every vertex is born at eps=0.
    # A component "dies" when it merges into another at some eps.
    uf = _UnionFind(n)
    edge_ptr = 0  # pointer into sorted edge list

    # Adjacency bookkeeping for cycle (beta_1) tracking
    adj = [set() for _ in range(n)]  # neighbours of each vertex
    triangles_filled = set()  # frozenset of triangles already counted
    active_edges = 0
    active_triangles = 0

    # beta_0 persistence: each vertex starts as its own component at birth=0
    # When two components merge at eps, the younger one dies.
    # "younger" = smaller component index by convention (arbitrary but consistent).
    component_birth = [0.0] * n  # all born at 0

    for step_idx, eps in enumerate(filtration):
        # Add all edges with weight <= eps
        while edge_ptr < len(edges) and edges[edge_ptr][0] <= eps + 1e-12:
            w, u, v = edges[edge_ptr]
            edge_ptr += 1
            active_edges += 1
            adj[u].add(v)
            adj[v].add(u)

            # Check for new triangles (u, v, common_neighbor)
            common = adj[u] & adj[v]
            for c in common:
                tri = frozenset((u, v, c))
                if tri not in triangles_filled:
                    triangles_filled.add(tri)
                    active_triangles += 1

            # Union-find merge for beta_0
            merged = uf.union(u, v)
            if merged:
                # Record death of one component
                persistence_diagram.append({
                    "dim": 0,
                    "birth": 0.0,
                    "death": round(w, 6),
                })

        # Betti numbers at this epsilon
        beta0 = uf.n_components
        betti_0.append(beta0)

        # beta_1 per connected component: sum of (E_c - V_c + 1)
        # Global formula: beta_1 = E - V + components - triangles_that_kill_cycles
        # Correct via Euler characteristic of 2-complex:
        #   chi = V - E + T  (vertices - edges + triangles)
        #   beta_0 - beta_1 + beta_2 = chi
        #   For a VR complex embedded in low dim, beta_2 ~ 0 in practice
        #   beta_1 = beta_0 - chi = beta_0 - V + E - T
        beta1 = max(0, beta0 - n + active_edges - active_triangles)
        betti_1.append(beta1)

    # ------------------------------------------------------------------
    # 4. Add persistence for cycles (beta_1)
    #    A cycle is born when an edge closes a loop (E increases without
    #    decreasing beta_0, and no triangle kills it).
    #    A cycle dies when a triangle fills it.
    #    We approximate by tracking edges and triangles in order.
    # ------------------------------------------------------------------
    # Re-walk edges to track cycle births/deaths
    cycle_persistence = _compute_cycle_persistence(edges, n)
    persistence_diagram.extend(cycle_persistence)

    # ------------------------------------------------------------------
    # 5. Persistence entropy and total persistence
    # ------------------------------------------------------------------
    lifetimes = []
    for feat in persistence_diagram:
        lt = feat["death"] - feat["birth"]
        if lt > 1e-12:
            lifetimes.append(lt)

    total_persistence = sum(lifetimes)

    if total_persistence > 0 and len(lifetimes) > 0:
        probs = [lt / total_persistence for lt in lifetimes]
        persistence_entropy = -sum(p * math.log(p) for p in probs if p > 0)
    else:
        persistence_entropy = 0.0

    # Summary stats
    max_persistence = max(lifetimes) if lifetimes else 0.0
    n_significant = sum(1 for lt in lifetimes if lt > 0.1)

    return {
        "distance_matrix": [[round(d, 6) for d in row] for row in dist],
        "betti_curve_0": betti_0,
        "betti_curve_1": betti_1,
        "persistence_diagram": persistence_diagram,
        "persistence_entropy": round(persistence_entropy, 6),
        "total_persistence": round(total_persistence, 6),
        "max_persistence": round(max_persistence, 6),
        "n_significant_features": n_significant,
        "filtration_values": [round(f, 6) for f in filtration],
    }


def _compute_cycle_persistence(edges, n):
    """Track 1-dimensional persistence (cycles born/killed).

    A cycle is born when adding an edge does NOT merge two components
    (i.e., both endpoints are already connected).  It dies when a
    triangle fills it.  We use union-find for component tracking and
    adjacency for triangle detection.
    """
    uf = _UnionFind(n)
    adj = [set() for _ in range(n)]
    pending_cycles = []  # list of (birth_eps,) — FIFO
    diagram = []

    for w, u, v in edges:
        # Does this edge merge two components?
        ru, rv = uf.find(u), uf.find(v)
        if ru == rv:
            # Both already connected -> this edge creates a cycle
            pending_cycles.append(w)

        uf.union(u, v)

        # Check if this edge completes any triangles
        common = adj[u] & adj[v]
        for _ in common:
            if pending_cycles:
                birth = pending_cycles.pop(0)
                diagram.append({
                    "dim": 1,
                    "birth": round(birth, 6),
                    "death": round(w, 6),
                })

        adj[u].add(v)
        adj[v].add(u)

    # Remaining cycles that never die (persist to infinity ~ 1.0)
    for birth in pending_cycles:
        diagram.append({
            "dim": 1,
            "birth": round(birth, 6),
            "death": 1.0,
        })

    return diagram
