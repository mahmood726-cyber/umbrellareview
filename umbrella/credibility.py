"""Credibility Theory for Review Weighting.

Actuarial Buhlmann credibility applied to evidence synthesis:
each SR's estimate is shrunk toward the grand mean based on its
number of primary studies relative to within/between variance ratio.
"""

from __future__ import annotations

import math
import numpy as np


def compute_credibility(reviews):
    """Run Buhlmann credibility analysis on a list of ReviewInput objects.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with keys:
        credibility_factors  — dict review_id -> Z_i
        credibility_estimates — dict review_id -> theta_cred_i
        k_cred               — float (credibility constant)
        full_credibility_threshold — float (n_full)
        experience_ranking   — list of review_ids sorted by cred estimate
        grand_mean           — float
    """
    if len(reviews) < 2:
        raise ValueError("Need at least 2 reviews for credibility analysis")

    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)
    ks = np.array([max(r.k, 1) for r in reviews])

    # Grand mean (simple average)
    theta_grand = float(np.mean(thetas))

    # Within-review variance: mean of se_i^2
    sigma_within_sq = float(np.mean(ses ** 2))

    # Between-review variance: method of moments (DerSimonian-Laird style)
    # tau2 = max(0, var(thetas) - mean(se^2))
    tau2_mm = float(np.var(thetas, ddof=1) - sigma_within_sq)
    sigma_between_sq = max(tau2_mm, 1e-8)

    # Credibility constant
    k_cred = sigma_within_sq / sigma_between_sq

    # Credibility factors: Z_i = n_i / (n_i + k_cred)
    credibility_factors = {}
    credibility_estimates = {}
    for i, r in enumerate(reviews):
        n_i = float(ks[i])
        z_i = n_i / (n_i + k_cred)
        theta_cred_i = z_i * thetas[i] + (1.0 - z_i) * theta_grand
        credibility_factors[r.review_id] = round(z_i, 6)
        credibility_estimates[r.review_id] = round(theta_cred_i, 6)

    # Full credibility threshold: n_full = (z/r)^2 * (sigma_within^2 / sigma_between^2)
    z_val = 1.96  # 95%
    r_val = 0.05  # 5% precision
    full_credibility_threshold = (z_val / r_val) ** 2 * (sigma_within_sq / sigma_between_sq)

    # Experience ranking: sort by credibility-adjusted estimate (descending by absolute effect)
    ranking = sorted(
        credibility_estimates.keys(),
        key=lambda rid: credibility_estimates[rid],
    )

    return {
        "credibility_factors": credibility_factors,
        "credibility_estimates": credibility_estimates,
        "k_cred": round(k_cred, 6),
        "full_credibility_threshold": round(full_credibility_threshold, 4),
        "experience_ranking": ranking,
        "grand_mean": round(theta_grand, 6),
    }
