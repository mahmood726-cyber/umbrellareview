"""Advanced Resampling Methods for Umbrella Reviews."""

import numpy as np
from scipy import stats


def _iv_pool(thetas, ses):
    """Inverse-variance fixed-effect pooling. Returns (theta, se, Q)."""
    w = 1.0 / (ses ** 2)
    theta_hat = float(np.sum(w * thetas) / np.sum(w))
    se_hat = float(1.0 / np.sqrt(np.sum(w)))
    q = float(np.sum(w * (thetas - theta_hat) ** 2))
    return theta_hat, se_hat, q


def _re_pool(thetas, ses):
    """DerSimonian-Laird random-effects pooling. Returns theta_re."""
    n = len(thetas)
    w = 1.0 / (ses ** 2)
    theta_fe = float(np.sum(w * thetas) / np.sum(w))
    q = float(np.sum(w * (thetas - theta_fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (n - 1)) / c) if c > 0 and n > 1 else 0.0
    w_re = 1.0 / (ses ** 2 + tau2)
    theta_re = float(np.sum(w_re * thetas) / np.sum(w_re))
    return theta_re


def compute_resampling(reviews, seed=42, B=1000, n_perm=999):
    """Compute jackknife, bootstrap, permutation, and subsampling diagnostics.

    Parameters
    ----------
    reviews : list[ReviewInput]
    seed : int
        Random seed for reproducibility.
    B : int
        Number of bootstrap replicates.
    n_perm : int
        Number of permutation replicates.

    Returns
    -------
    dict with keys: jackknife_estimate, jackknife_se, bootstrap_ci_percentile,
                    bootstrap_ci_bca, permutation_p, subsampling_stability
    """
    n = len(reviews)
    if n < 3:
        raise ValueError("Need at least 3 reviews for resampling analysis")

    rng = np.random.default_rng(seed)

    thetas = np.array([r.theta for r in reviews], dtype=float)
    ses = np.array([r.se for r in reviews], dtype=float)
    ses = np.where(ses > 0, ses, 0.1)

    # Full RE pooled estimate
    theta_full = _re_pool(thetas, ses)

    # --- Jackknife ---
    jack_estimates = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        jack_estimates[i] = _re_pool(thetas[mask], ses[mask])

    jack_mean = float(np.mean(jack_estimates))
    # Jackknife bias
    jack_bias = (n - 1) * (jack_mean - theta_full)
    # Bias-corrected estimate
    jackknife_estimate = float(theta_full - jack_bias)
    # Jackknife variance and SE
    jack_var = (n - 1) / n * np.sum((jack_estimates - jack_mean) ** 2)
    jackknife_se = float(np.sqrt(jack_var))

    # --- Bootstrap ---
    boot_estimates = np.empty(B)
    for b in range(B):
        idx = rng.integers(0, n, size=n)
        boot_estimates[b] = _re_pool(thetas[idx], ses[idx])

    # Percentile CI
    alpha = 0.05
    lo_pct = float(np.percentile(boot_estimates, 100 * alpha / 2))
    hi_pct = float(np.percentile(boot_estimates, 100 * (1 - alpha / 2)))
    bootstrap_ci_percentile = (round(lo_pct, 6), round(hi_pct, 6))

    # BCa CI
    # z0: bias correction
    prop_below = np.mean(boot_estimates < theta_full)
    # Clip to avoid infinite z0
    prop_below = np.clip(prop_below, 1e-10, 1 - 1e-10)
    z0 = float(stats.norm.ppf(prop_below))

    # a: acceleration (from jackknife influence values)
    d = jack_mean - jack_estimates  # influence values
    sum_d2 = np.sum(d ** 2)
    sum_d3 = np.sum(d ** 3)
    if sum_d2 > 0:
        a = float(sum_d3 / (6.0 * (sum_d2 ** 1.5)))
    else:
        a = 0.0

    # Adjusted percentiles
    z_alpha_lo = stats.norm.ppf(alpha / 2)
    z_alpha_hi = stats.norm.ppf(1 - alpha / 2)

    denom_lo = 1 - a * (z0 + z_alpha_lo)
    denom_hi = 1 - a * (z0 + z_alpha_hi)
    if abs(denom_lo) < 1e-15:
        denom_lo = 1e-15
    if abs(denom_hi) < 1e-15:
        denom_hi = 1e-15

    adj_lo = stats.norm.cdf(z0 + (z0 + z_alpha_lo) / denom_lo)
    adj_hi = stats.norm.cdf(z0 + (z0 + z_alpha_hi) / denom_hi)

    # Clip to valid percentile range
    adj_lo = np.clip(adj_lo, 0.001, 0.999)
    adj_hi = np.clip(adj_hi, 0.001, 0.999)

    bca_lo = float(np.percentile(boot_estimates, 100 * adj_lo))
    bca_hi = float(np.percentile(boot_estimates, 100 * adj_hi))
    bootstrap_ci_bca = (round(bca_lo, 6), round(bca_hi, 6))

    # --- Permutation test for heterogeneity ---
    _, _, q_obs = _iv_pool(thetas, ses)
    perm_count = 0
    for _ in range(n_perm):
        perm_thetas = rng.permutation(thetas)
        _, _, q_perm = _iv_pool(perm_thetas, ses)
        if q_perm >= q_obs:
            perm_count += 1
    permutation_p = round(float((perm_count + 1) / (n_perm + 1)), 6)

    # --- Subsampling stability ---
    # Pool for subsets of size m = 2..n-1, measure max deviation
    max_delta = 0.0
    for m in range(2, n):
        # Use first m reviews (deterministic for reproducibility)
        sub_theta = _re_pool(thetas[:m], ses[:m])
        delta = abs(sub_theta - theta_full)
        if delta > max_delta:
            max_delta = delta

    if abs(theta_full) > 1e-15:
        subsampling_stability = round(float(1 - max_delta / abs(theta_full)), 6)
    else:
        subsampling_stability = 1.0

    return {
        "jackknife_estimate": round(jackknife_estimate, 6),
        "jackknife_se": round(jackknife_se, 6),
        "bootstrap_ci_percentile": bootstrap_ci_percentile,
        "bootstrap_ci_bca": bootstrap_ci_bca,
        "permutation_p": permutation_p,
        "subsampling_stability": subsampling_stability,
    }
