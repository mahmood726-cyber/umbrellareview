"""Finite Mixture Meta-Meta-Analysis.

Model the set of systematic reviews as a mixture of latent populations.
EM algorithm fits K=1..3 components; BIC selects the best model.
Each review gets posterior probability of belonging to each component.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.stats import norm


def _log_likelihood(thetas, ses, mus, sigmas, pis):
    """Compute log-likelihood of the mixture model."""
    n = len(thetas)
    K = len(mus)
    ll = 0.0
    for i in range(n):
        density = 0.0
        for j in range(K):
            v_ij = ses[i] ** 2 + sigmas[j] ** 2
            density += pis[j] * norm.pdf(thetas[i], loc=mus[j], scale=math.sqrt(v_ij))
        if density > 0:
            ll += math.log(density)
        else:
            ll += -1e10  # avoid log(0)
    return ll


def _fit_mixture(thetas, ses, K, max_iter=100, tol=1e-6, seed=42):
    """Fit a K-component Gaussian mixture via EM.

    Returns dict with component parameters, posterior probs, BIC, etc.
    """
    n = len(thetas)
    rng = np.random.RandomState(seed)

    # Initialise: spread mus across the range of thetas
    if K == 1:
        mus = np.array([np.mean(thetas)])
    else:
        quantiles = np.linspace(0, 1, K + 2)[1:-1]
        mus = np.quantile(thetas, quantiles).copy()
        # Add small jitter for stability
        mus += rng.normal(0, 0.01, size=K)

    sigmas_sq = np.full(K, np.var(thetas) / K)
    pis = np.full(K, 1.0 / K)

    # Ensure sigmas_sq positive
    sigmas_sq = np.maximum(sigmas_sq, 1e-8)

    prev_ll = -np.inf
    converged = False

    gamma = np.zeros((n, K))

    for iteration in range(max_iter):
        # ── E-step ──
        for i in range(n):
            for j in range(K):
                v_ij = ses[i] ** 2 + sigmas_sq[j]
                gamma[i, j] = pis[j] * norm.pdf(thetas[i], loc=mus[j], scale=math.sqrt(v_ij))
            row_sum = gamma[i, :].sum()
            if row_sum > 0:
                gamma[i, :] /= row_sum
            else:
                gamma[i, :] = 1.0 / K

        # ── M-step ──
        for j in range(K):
            gamma_j = gamma[:, j]
            n_j = gamma_j.sum()
            if n_j < 1e-10:
                continue

            # Update pi
            pis[j] = n_j / n

            # Update mu (weighted by 1/v_ij)
            v_j = ses ** 2 + sigmas_sq[j]
            w_j = gamma_j / v_j
            w_sum = w_j.sum()
            if w_sum > 0:
                mus[j] = (w_j * thetas).sum() / w_sum

            # Update sigma_sq via profile: maximise over sigma_j^2
            # Use moment estimator: sigma_j^2 = weighted variance - weighted mean(se^2)
            resid_sq = (thetas - mus[j]) ** 2
            sigma_new = (gamma_j * resid_sq).sum() / n_j - (gamma_j * ses ** 2).sum() / n_j
            sigmas_sq[j] = max(sigma_new, 1e-8)

        # ── Check convergence ──
        ll = _log_likelihood(thetas, ses, mus, np.sqrt(sigmas_sq), pis)
        if abs(ll - prev_ll) < tol:
            converged = True
            prev_ll = ll
            break
        prev_ll = ll

    # Final log-likelihood
    ll = prev_ll

    # BIC: k_params = K means + K variances + (K-1) mixing weights
    k_params = K + K + (K - 1)
    bic = -2.0 * ll + k_params * math.log(n) if n > 0 else float("inf")

    # Build posterior probs dict keyed by review index
    posterior = gamma.tolist()

    return {
        "K": K,
        "mus": mus.tolist(),
        "sigmas_sq": sigmas_sq.tolist(),
        "pis": pis.tolist(),
        "posterior": posterior,
        "log_likelihood": ll,
        "bic": bic,
        "converged": converged,
    }


def compute_mixture_meta(reviews, max_k=3, seed=42):
    """Run finite mixture meta-meta-analysis.

    Parameters
    ----------
    reviews : list[ReviewInput]
    max_k : int
        Maximum number of components to try (1..max_k).
    seed : int
        Random seed for initialisation.

    Returns
    -------
    dict with keys:
        n_components, component_means, component_variances, weights,
        posterior_probs (dict review_id -> list), bic_by_k, log_likelihood,
        converged
    """
    if len(reviews) < 2:
        raise ValueError("Need at least 2 reviews for mixture model")

    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)

    # Fit K=1..max_k
    fits = {}
    bic_by_k = {}
    for k in range(1, max_k + 1):
        if k > len(reviews):
            break
        fit = _fit_mixture(thetas, ses, k, seed=seed)
        fits[k] = fit
        bic_by_k[k] = fit["bic"]

    # Select best K by BIC (lowest)
    best_k = min(bic_by_k, key=lambda x: bic_by_k[x])
    best = fits[best_k]

    # Build posterior_probs keyed by review_id
    posterior_probs = {}
    for i, r in enumerate(reviews):
        posterior_probs[r.review_id] = best["posterior"][i]

    return {
        "n_components": best["K"],
        "component_means": best["mus"],
        "component_variances": best["sigmas_sq"],
        "weights": best["pis"],
        "posterior_probs": posterior_probs,
        "bic_by_k": bic_by_k,
        "log_likelihood": best["log_likelihood"],
        "converged": best["converged"],
    }
