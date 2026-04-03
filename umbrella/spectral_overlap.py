"""Spectral Analysis of Study Overlap.

Models the review-study relationship as a bipartite graph and uses
spectral methods (SVD, Laplacian) to quantify overlap structure,
identify review clusters, and measure algebraic connectivity.
"""

import numpy as np
from scipy import linalg


def spectral_overlap(reviews):
    """Spectral analysis of the review-study bipartite graph.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews with study_ids.

    Returns
    -------
    dict with keys: spectral_concentration, fiedler_value,
        svd_explained (list), review_clusters (list of int),
        co_citation_matrix (list of lists)
    """
    n_reviews = len(reviews)
    if n_reviews < 2:
        raise ValueError("Need at least 2 reviews for spectral analysis")

    # Build the universe of unique studies
    all_studies = sorted(set(sid for r in reviews for sid in r.study_ids))
    n_studies = len(all_studies)
    study_idx = {sid: j for j, sid in enumerate(all_studies)}

    # Incidence matrix A: reviews x studies
    A = np.zeros((n_reviews, n_studies), dtype=float)
    for i, r in enumerate(reviews):
        for sid in r.study_ids:
            A[i, study_idx[sid]] = 1.0

    # SVD: A = U Sigma V^T
    U, sigma, Vt = linalg.svd(A, full_matrices=False)

    # Spectral concentration: ratio of first singular value to sum
    sigma_sum = float(np.sum(sigma))
    spectral_concentration = float(sigma[0] / sigma_sum) if sigma_sum > 0 else 0.0

    # SVD explained variance (proportion of each singular value)
    svd_explained = [round(float(s / sigma_sum), 6) for s in sigma] if sigma_sum > 0 else []

    # Co-citation matrix: C = A @ A^T (reviews x reviews)
    C = A @ A.T
    co_citation_matrix = [[int(C[i][j]) for j in range(n_reviews)] for i in range(n_reviews)]

    # Normalized Laplacian and Fiedler value
    degrees = np.diag(C).copy()
    # Degree matrix
    D = np.diag(degrees)
    L = D - C  # unnormalized Laplacian

    # Normalized Laplacian: D^{-1/2} L D^{-1/2}
    d_inv_sqrt = np.zeros(n_reviews)
    for i in range(n_reviews):
        if degrees[i] > 0:
            d_inv_sqrt[i] = 1.0 / np.sqrt(degrees[i])
    D_inv_sqrt = np.diag(d_inv_sqrt)
    L_norm = D_inv_sqrt @ L @ D_inv_sqrt

    # Eigenvalues of normalized Laplacian
    eigenvalues = np.sort(np.real(linalg.eigvalsh(L_norm)))

    # Fiedler value = 2nd smallest eigenvalue (algebraic connectivity)
    if len(eigenvalues) >= 2:
        fiedler_value = float(max(0.0, eigenvalues[1]))
    else:
        fiedler_value = 0.0

    # Cluster reviews using left singular vectors (k=2 clusters via sign of 2nd left SV)
    if U.shape[1] >= 2:
        # Use sign of 2nd left singular vector as cluster assignment
        cluster_vec = U[:, 1]
        review_clusters = [int(0 if v >= 0 else 1) for v in cluster_vec]
    else:
        review_clusters = [0] * n_reviews

    return {
        "spectral_concentration": round(spectral_concentration, 6),
        "fiedler_value": round(fiedler_value, 6),
        "svd_explained": svd_explained,
        "review_clusters": review_clusters,
        "co_citation_matrix": co_citation_matrix,
    }
