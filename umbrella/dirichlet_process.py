"""Dirichlet Process Clustering of Systematic Reviews.

Nonparametric Bayesian clustering to discover latent review types using
a stick-breaking construction and collapsed Gibbs sampling via the
Chinese Restaurant Process (CRP).
"""

import math
import numpy as np


def _log_normal_pdf(x, mu, var):
    """Log of N(x | mu, var)."""
    if var <= 0:
        return -1e300
    return -0.5 * math.log(2 * math.pi * var) - 0.5 * (x - mu) ** 2 / var


def _update_cluster_params(thetas, ses, assignments, k):
    """Compute posterior mean for cluster k given assigned reviews.

    Model: theta_i | z_i=k ~ N(mu_k, se_i^2 + sigma_k^2)
    Prior: mu_k ~ N(0, 10), sigma_k^2 ~ InvGamma(1, 1)

    Returns (mu_k, sigma_k^2) — point estimates.
    """
    idx = [i for i, z in enumerate(assignments) if z == k]
    if len(idx) == 0:
        return 0.0, 1.0

    # Estimate sigma_k^2 from residual variance (MAP of InvGamma posterior)
    t = thetas[idx]
    s = ses[idx]
    mu_hat = float(np.mean(t))
    resid_var = float(np.mean((t - mu_hat) ** 2))
    # InvGamma(1,1) MAP: (beta)/(alpha+1) with posterior update
    # Posterior: InvGamma(alpha + n/2, beta + 0.5*sum(resid^2))
    alpha_post = 1.0 + len(idx) / 2.0
    beta_post = 1.0 + 0.5 * np.sum((t - mu_hat) ** 2)
    sigma_k2 = float(beta_post / (alpha_post + 1.0))

    # Posterior mean for mu_k: weighted combination of prior and data
    # Prior: N(0, 10), Likelihood: product of N(theta_i, se_i^2 + sigma_k^2)
    prior_prec = 1.0 / 10.0
    data_prec = float(np.sum(1.0 / (s ** 2 + sigma_k2)))
    post_prec = prior_prec + data_prec
    post_mean = float(np.sum(t / (s ** 2 + sigma_k2))) / post_prec

    return post_mean, sigma_k2


def dirichlet_process_cluster(reviews, alpha=1.0, n_iter=500, burn_in=200,
                               k_max=10, seed=42):
    """Cluster reviews using a Dirichlet Process mixture model.

    Uses collapsed Gibbs sampling via the Chinese Restaurant Process.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews.
    alpha : float
        Concentration parameter (higher = more clusters expected).
    n_iter : int
        Total Gibbs iterations.
    burn_in : int
        Iterations to discard before collecting samples.
    k_max : int
        Maximum number of components (truncation level).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict with keys:
        n_clusters             – int, posterior mode of number of clusters
        cluster_assignments    – dict[review_id, int], modal cluster per review
        cluster_means          – dict[int, float], posterior mean per cluster
        cluster_sizes          – dict[int, int], number of reviews per cluster
        concentration_alpha    – float, the alpha used
        assignment_certainty   – dict[review_id, float], proportion of post-
                                 burn-in iterations in modal cluster
    """
    n = len(reviews)
    if n < 2:
        raise ValueError("Need at least 2 reviews for clustering")

    rng = np.random.RandomState(seed)
    thetas = np.array([r.theta for r in reviews], dtype=float)
    ses = np.array([r.se for r in reviews], dtype=float)
    ses = np.where(ses > 0, ses, 0.1)

    # Initialize: all reviews in one cluster
    assignments = [0] * n
    # Cluster parameters: {cluster_id: (mu, sigma2)}
    cluster_params = {0: _update_cluster_params(thetas, ses, assignments, 0)}
    next_cluster_id = 1

    # Storage for post-burn-in samples
    n_post = max(1, n_iter - burn_in)
    assignment_samples = []  # list of assignment vectors

    for iteration in range(n_iter):
        for i in range(n):
            # Remove review i from its current cluster
            old_k = assignments[i]

            # Count members in each cluster excluding i
            cluster_counts = {}
            for j, z in enumerate(assignments):
                if j != i:
                    cluster_counts[z] = cluster_counts.get(z, 0) + 1

            # If old cluster is now empty, remove it
            if cluster_counts.get(old_k, 0) == 0:
                cluster_counts.pop(old_k, None)
                cluster_params.pop(old_k, None)

            # Compute CRP probabilities for each existing cluster + new
            log_probs = []
            cluster_ids = []

            for k, count in cluster_counts.items():
                mu_k, sigma_k2 = cluster_params.get(k, (0.0, 1.0))
                var_i = ses[i] ** 2 + sigma_k2
                log_p = math.log(count) + _log_normal_pdf(thetas[i], mu_k, var_i)
                log_probs.append(log_p)
                cluster_ids.append(k)

            # New cluster probability
            new_var = ses[i] ** 2 + 10.0  # prior variance for mu
            log_p_new = math.log(alpha) + _log_normal_pdf(thetas[i], 0.0, new_var)
            log_probs.append(log_p_new)
            cluster_ids.append(-1)  # sentinel for new cluster

            # Normalize in log-space
            max_lp = max(log_probs)
            probs = [math.exp(lp - max_lp) for lp in log_probs]
            total = sum(probs)
            probs = [p / total for p in probs]

            # Sample cluster assignment
            chosen_idx = _categorical_sample(probs, rng)
            chosen_k = cluster_ids[chosen_idx]

            if chosen_k == -1:
                # New cluster
                chosen_k = next_cluster_id
                next_cluster_id += 1
                cluster_params[chosen_k] = (float(thetas[i]), 1.0)
                cluster_counts[chosen_k] = 0

            assignments[i] = chosen_k

        # Update cluster parameters
        active_clusters = set(assignments)
        for k in active_clusters:
            cluster_params[k] = _update_cluster_params(thetas, ses, assignments, k)

        # Remove empty clusters from params
        for k in list(cluster_params.keys()):
            if k not in active_clusters:
                del cluster_params[k]

        # Enforce k_max truncation: if too many clusters, merge smallest
        if len(active_clusters) > k_max:
            counts = {}
            for z in assignments:
                counts[z] = counts.get(z, 0) + 1
            sorted_clusters = sorted(counts.items(), key=lambda x: x[1])
            # Merge smallest into second-smallest until at k_max
            while len(set(assignments)) > k_max:
                smallest_k = sorted_clusters[0][0]
                target_k = sorted_clusters[1][0]
                for j in range(n):
                    if assignments[j] == smallest_k:
                        assignments[j] = target_k
                sorted_clusters = sorted_clusters[1:]
                sorted_clusters[0] = (target_k, sorted_clusters[0][1] + 1)

        # Collect post-burn-in samples
        if iteration >= burn_in:
            assignment_samples.append(list(assignments))

    # ------------------------------------------------------------------
    # Summarize posterior
    # ------------------------------------------------------------------
    # Relabel clusters to 0-based consecutive integers based on final sample
    final_assignments = assignment_samples[-1] if assignment_samples else assignments
    unique_labels = sorted(set(final_assignments))
    label_map = {old: new for new, old in enumerate(unique_labels)}

    # Modal assignment per review across post-burn-in samples
    modal_assignments = {}
    assignment_certainty = {}
    for i in range(n):
        counts = {}
        for sample in assignment_samples:
            lab = sample[i]
            counts[lab] = counts.get(lab, 0) + 1
        if counts:
            mode_lab = max(counts, key=counts.get)
            modal_assignments[reviews[i].review_id] = label_map.get(mode_lab, 0)
            assignment_certainty[reviews[i].review_id] = counts[mode_lab] / len(assignment_samples)
        else:
            modal_assignments[reviews[i].review_id] = 0
            assignment_certainty[reviews[i].review_id] = 1.0

    # Cluster sizes and means from modal assignments
    cluster_sizes = {}
    cluster_members = {}
    for rid, k in modal_assignments.items():
        cluster_sizes[k] = cluster_sizes.get(k, 0) + 1
        if k not in cluster_members:
            cluster_members[k] = []
        # Find the review
        for r in reviews:
            if r.review_id == rid:
                cluster_members[k].append(r.theta)
                break

    cluster_means = {}
    for k, members in cluster_members.items():
        cluster_means[k] = round(float(np.mean(members)), 6)

    # Posterior mode of number of clusters
    n_clusters_samples = [len(set(s)) for s in assignment_samples] if assignment_samples else [1]
    n_clusters_counts = {}
    for nc in n_clusters_samples:
        n_clusters_counts[nc] = n_clusters_counts.get(nc, 0) + 1
    n_clusters = max(n_clusters_counts, key=n_clusters_counts.get)

    return {
        "n_clusters": n_clusters,
        "cluster_assignments": modal_assignments,
        "cluster_means": cluster_means,
        "cluster_sizes": cluster_sizes,
        "concentration_alpha": alpha,
        "assignment_certainty": assignment_certainty,
    }


def _categorical_sample(probs, rng):
    """Sample from a categorical distribution given probabilities."""
    u = rng.random()
    cumsum = 0.0
    for i, p in enumerate(probs):
        cumsum += p
        if u < cumsum:
            return i
    return len(probs) - 1
