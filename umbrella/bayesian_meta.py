"""Bayesian Hierarchical Meta-Meta-Analysis.

Random-effects meta-meta-analysis treating each SR's pooled estimate
as a "study".  Estimates between-review heterogeneity (tau_mm^2) via
DerSimonian-Laird, then produces posterior summaries, prediction
intervals, shrinkage estimates, and BIC model comparison.
"""

import math
import numpy as np
from scipy import stats


def bayesian_meta_meta(reviews):
    """Run Bayesian hierarchical meta-meta-analysis.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews with theta and se > 0.

    Returns
    -------
    dict with keys: mu, mu_ci, tau_mm2, i2_mm, prediction_interval,
        shrinkage_estimates, bic_fixed, bic_random
    """
    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)  # fallback for missing SE
    n = len(reviews)

    if n < 2:
        raise ValueError("Need at least 2 reviews for meta-meta-analysis")

    variances = ses ** 2

    # --- Fixed-effect model (for BIC comparison) ---
    w_fixed = 1.0 / variances
    mu_fixed = float(np.sum(w_fixed * thetas) / np.sum(w_fixed))
    se_fixed = float(1.0 / np.sqrt(np.sum(w_fixed)))
    q_fixed = float(np.sum(w_fixed * (thetas - mu_fixed) ** 2))
    ll_fixed = float(np.sum(stats.norm.logpdf(thetas, loc=mu_fixed, scale=ses)))
    bic_fixed = -2.0 * ll_fixed + 1.0 * math.log(n)  # 1 parameter (mu)

    # --- DerSimonian-Laird for tau_mm^2 ---
    c = float(np.sum(w_fixed) - np.sum(w_fixed ** 2) / np.sum(w_fixed))
    tau_mm2 = max(0.0, (q_fixed - (n - 1)) / c) if c > 0 else 0.0

    # --- Random-effects model ---
    w_re = 1.0 / (variances + tau_mm2)
    mu_re = float(np.sum(w_re * thetas) / np.sum(w_re))
    se_mu = float(1.0 / np.sqrt(np.sum(w_re)))
    mu_ci = (mu_re - 1.96 * se_mu, mu_re + 1.96 * se_mu)

    # I^2 at meta-meta level
    if q_fixed > 0 and n > 1:
        i2_mm = max(0.0, (q_fixed - (n - 1)) / q_fixed * 100.0)
    else:
        i2_mm = 0.0

    # Prediction interval: mu +/- t_{k-1, 0.975} * sqrt(tau_mm2 + se_mu^2)
    if n > 2:
        t_crit = float(stats.t.ppf(0.975, df=n - 1))
    else:
        t_crit = float(stats.t.ppf(0.975, df=1))
    pi_half = t_crit * math.sqrt(tau_mm2 + se_mu ** 2)
    prediction_interval = (mu_re - pi_half, mu_re + pi_half)

    # Shrinkage estimates: each review pulled toward grand mean
    shrinkage_estimates = {}
    for i, r in enumerate(reviews):
        B_i = tau_mm2 / (variances[i] + tau_mm2) if (variances[i] + tau_mm2) > 0 else 0.0
        shrunk = mu_re + B_i * (thetas[i] - mu_re)
        shrinkage_estimates[r.review_id] = round(float(shrunk), 6)

    # BIC for random-effects model (2 parameters: mu, tau^2)
    re_scales = np.sqrt(variances + tau_mm2)
    ll_random = float(np.sum(stats.norm.logpdf(thetas, loc=mu_re, scale=re_scales)))
    bic_random = -2.0 * ll_random + 2.0 * math.log(n)

    return {
        "mu": round(mu_re, 6),
        "mu_ci": (round(mu_ci[0], 6), round(mu_ci[1], 6)),
        "tau_mm2": round(tau_mm2, 6),
        "i2_mm": round(i2_mm, 1),
        "prediction_interval": (round(prediction_interval[0], 6), round(prediction_interval[1], 6)),
        "shrinkage_estimates": shrinkage_estimates,
        "bic_fixed": round(bic_fixed, 4),
        "bic_random": round(bic_random, 4),
    }
