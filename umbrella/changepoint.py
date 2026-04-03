"""Temporal Change-Point Detection in Systematic Review Estimates.

Detect whether the body of SR evidence has shifted over time using
CUSUM, Bayesian single change-point model, and effect drift analysis.
"""

import math

import numpy as np
from scipy import stats


def _cusum(thetas, ses, theta_bar):
    """Compute CUSUM statistic: cumulative standardized deviations.

    S_t = sum_{i=1}^{t} (theta_i - theta_bar) / se_i
    """
    n = len(thetas)
    cusum = []
    s = 0.0
    for i in range(n):
        se_i = ses[i] if ses[i] > 0 else 0.1
        s += (thetas[i] - theta_bar) / se_i
        cusum.append(round(s, 6))
    return cusum


def _cusum_bootstrap_p(thetas, ses, theta_bar, observed_max, n_perm=999, seed=42):
    """Bootstrap p-value for CUSUM: permute review order and compare max|S|."""
    rng = np.random.RandomState(seed)
    n = len(thetas)
    n_exceed = 0

    for _ in range(n_perm):
        perm = rng.permutation(n)
        t_perm = thetas[perm]
        s_perm = ses[perm]
        cusum_perm = _cusum(t_perm, s_perm, theta_bar)
        max_perm = max(abs(c) for c in cusum_perm)
        if max_perm >= observed_max:
            n_exceed += 1

    return round((n_exceed + 1) / (n_perm + 1), 4)


def _bayesian_changepoint(thetas, ses, years):
    """Bayesian single change-point detection.

    Model: theta_i ~ N(mu_1, se_i^2) for i <= tau
           theta_i ~ N(mu_2, se_i^2) for i > tau

    With uniform prior on tau in {1, ..., k-1}.
    """
    n = len(thetas)
    if n < 3:
        # Need at least 3 points for a meaningful change-point
        return {
            "changepoint_year": int(years[0]) if len(years) > 0 else 0,
            "posterior": {int(years[0]): 1.0} if len(years) > 0 else {},
        }

    log_liks = []
    for tau in range(1, n):
        # Before change-point (indices 0..tau-1)
        t_before = thetas[:tau]
        s_before = ses[:tau]
        w_before = 1.0 / (s_before ** 2)
        mu1_hat = float(np.sum(w_before * t_before) / np.sum(w_before))

        # After change-point (indices tau..n-1)
        t_after = thetas[tau:]
        s_after = ses[tau:]
        w_after = 1.0 / (s_after ** 2)
        mu2_hat = float(np.sum(w_after * t_after) / np.sum(w_after))

        # Log-likelihood
        ll = 0.0
        for i in range(tau):
            var_i = s_before[i] ** 2
            ll += -0.5 * math.log(2 * math.pi * var_i) - 0.5 * (t_before[i] - mu1_hat) ** 2 / var_i
        for i in range(len(t_after)):
            var_i = s_after[i] ** 2
            ll += -0.5 * math.log(2 * math.pi * var_i) - 0.5 * (t_after[i] - mu2_hat) ** 2 / var_i

        log_liks.append(ll)

    # Posterior: proportional to exp(LL) with uniform prior
    max_ll = max(log_liks)
    unnorm = [math.exp(ll - max_ll) for ll in log_liks]
    total = sum(unnorm)
    posterior_probs = [p / total for p in unnorm]

    # Map to years: tau=1 means change after index 0 = after year[0]
    # Report the year AFTER which the change occurs
    posterior = {}
    for idx, prob in enumerate(posterior_probs):
        tau = idx + 1  # tau ranges 1..n-1
        # Change-point year: between years[tau-1] and years[tau]
        cp_year = int(years[tau]) if tau < n else int(years[-1])
        posterior[cp_year] = posterior.get(cp_year, 0.0) + round(prob, 6)

    # MAP estimate
    map_idx = int(np.argmax(posterior_probs))
    map_tau = map_idx + 1
    map_year = int(years[map_tau]) if map_tau < n else int(years[-1])

    return {
        "changepoint_year": map_year,
        "posterior": posterior,
    }


def _effect_drift(thetas, ses, years):
    """Weighted linear regression of theta on year.

    Weight = 1/se^2.  Returns slope, CI, and p-value.
    """
    n = len(thetas)
    w = 1.0 / (ses ** 2)

    # Weighted linear regression: theta = a + b * year
    x = years.astype(float)
    y = thetas

    sum_w = float(np.sum(w))
    sum_wx = float(np.sum(w * x))
    sum_wy = float(np.sum(w * y))
    sum_wxx = float(np.sum(w * x ** 2))
    sum_wxy = float(np.sum(w * x * y))

    denom = sum_w * sum_wxx - sum_wx ** 2
    if abs(denom) < 1e-15:
        return {"slope": 0.0, "slope_ci": (0.0, 0.0), "drift_p": 1.0}

    b = (sum_w * sum_wxy - sum_wx * sum_wy) / denom
    a = (sum_wy - b * sum_wx) / sum_w

    # Residuals and SE of slope
    resid = y - (a + b * x)
    # Weighted residual variance
    if n > 2:
        s2 = float(np.sum(w * resid ** 2) / (n - 2))
    else:
        s2 = 1.0

    se_b = math.sqrt(s2 * sum_w / denom) if denom > 0 else 0.0

    # t-test for slope
    if se_b > 0:
        t_stat = b / se_b
        df = n - 2
        p_val = 2.0 * (1.0 - stats.t.cdf(abs(t_stat), df)) if df > 0 else 1.0
    else:
        p_val = 1.0

    ci_lo = b - 1.96 * se_b
    ci_hi = b + 1.96 * se_b

    return {
        "slope": round(float(b), 6),
        "slope_ci": (round(float(ci_lo), 6), round(float(ci_hi), 6)),
        "drift_p": round(float(p_val), 4),
    }


def _before_after_comparison(thetas, ses, years, cp_year):
    """Compare pooled estimates before and after change-point."""
    before_mask = years <= cp_year
    after_mask = years > cp_year

    # Handle edge cases: if all on one side, shift to include at least 1
    if not np.any(before_mask):
        before_mask[0] = True
        after_mask[0] = False
    if not np.any(after_mask):
        after_mask[-1] = True
        before_mask[-1] = False

    t_before = thetas[before_mask]
    s_before = ses[before_mask]
    t_after = thetas[after_mask]
    s_after = ses[after_mask]

    # IV-weighted pooled estimates
    w_before = 1.0 / (s_before ** 2)
    theta_before = float(np.sum(w_before * t_before) / np.sum(w_before))

    w_after = 1.0 / (s_after ** 2)
    theta_after = float(np.sum(w_after * t_after) / np.sum(w_after))

    se_before = float(1.0 / np.sqrt(np.sum(w_before)))
    se_after = float(1.0 / np.sqrt(np.sum(w_after)))

    delta = theta_after - theta_before
    se_delta = math.sqrt(se_before ** 2 + se_after ** 2)
    ci_lo = delta - 1.96 * se_delta
    ci_hi = delta + 1.96 * se_delta

    return {
        "before_theta": round(theta_before, 6),
        "after_theta": round(theta_after, 6),
        "delta_theta": round(delta, 6),
        "delta_ci": (round(ci_lo, 6), round(ci_hi, 6)),
    }


def changepoint_detection(reviews, n_perm=999, seed=42):
    """Detect temporal change-points in SR estimates.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 3 reviews with year field.
    n_perm : int
        Number of permutations for CUSUM bootstrap p-value.
    seed : int
        Random seed.

    Returns
    -------
    dict with keys:
        cusum_statistic          – list[float]
        cusum_changepoint_year   – int
        cusum_p_value            – float
        bayesian_changepoint_year – int
        bayesian_posterior       – dict[int, float] (year -> probability)
        drift_slope              – float
        drift_p                  – float
        before_theta             – float
        after_theta              – float
        delta_theta              – float
        delta_ci                 – tuple[float, float]
    """
    if len(reviews) < 3:
        raise ValueError("Need at least 3 reviews for change-point detection")

    # Sort by year
    sorted_reviews = sorted(reviews, key=lambda r: (r.year, r.review_id))

    thetas = np.array([r.theta for r in sorted_reviews], dtype=float)
    ses = np.array([r.se for r in sorted_reviews], dtype=float)
    ses = np.where(ses > 0, ses, 0.1)
    years = np.array([r.year for r in sorted_reviews], dtype=float)

    # Weighted overall mean
    w = 1.0 / (ses ** 2)
    theta_bar = float(np.sum(w * thetas) / np.sum(w))

    # ----- CUSUM -----
    cusum = _cusum(thetas, ses, theta_bar)
    abs_cusum = [abs(c) for c in cusum]
    max_idx = int(np.argmax(abs_cusum))
    cusum_cp_year = int(years[max_idx])
    observed_max = abs_cusum[max_idx]

    cusum_p = _cusum_bootstrap_p(thetas, ses, theta_bar, observed_max,
                                  n_perm=n_perm, seed=seed)

    # ----- Bayesian change-point -----
    bayes = _bayesian_changepoint(thetas, ses, years)

    # ----- Effect drift -----
    drift = _effect_drift(thetas, ses, years)

    # ----- Before/after comparison at CUSUM change-point -----
    ba = _before_after_comparison(thetas, ses, years, float(cusum_cp_year))

    return {
        "cusum_statistic": cusum,
        "cusum_changepoint_year": cusum_cp_year,
        "cusum_p_value": cusum_p,
        "bayesian_changepoint_year": bayes["changepoint_year"],
        "bayesian_posterior": bayes["posterior"],
        "drift_slope": drift["slope"],
        "drift_p": drift["drift_p"],
        "before_theta": ba["before_theta"],
        "after_theta": ba["after_theta"],
        "delta_theta": ba["delta_theta"],
        "delta_ci": ba["delta_ci"],
    }
