"""Robust Meta-Meta-Analysis: Outlier-resistant pooling methods."""

import numpy as np


def _iv_pool(thetas, ses):
    """Standard inverse-variance fixed-effect pooling."""
    w = 1.0 / (ses ** 2)
    theta_hat = float(np.sum(w * thetas) / np.sum(w))
    se_hat = float(1.0 / np.sqrt(np.sum(w)))
    ci_lo = theta_hat - 1.96 * se_hat
    ci_hi = theta_hat + 1.96 * se_hat
    return {"theta": round(theta_hat, 6), "se": round(se_hat, 6),
            "ci_lo": round(ci_lo, 6), "ci_hi": round(ci_hi, 6)}


def _weighted_median(thetas, weights):
    """Compute weighted median."""
    order = np.argsort(thetas)
    sorted_t = thetas[order]
    sorted_w = weights[order]
    cum_w = np.cumsum(sorted_w)
    half = cum_w[-1] / 2.0
    idx = np.searchsorted(cum_w, half)
    idx = min(idx, len(sorted_t) - 1)
    return float(sorted_t[idx])


def _median_pool(thetas, ses, n_boot=1000, seed=42):
    """Weighted median with bootstrap CI."""
    weights = 1.0 / (ses ** 2)
    med = _weighted_median(thetas, weights)

    rng = np.random.RandomState(seed)
    boot_meds = []
    n = len(thetas)
    for _ in range(n_boot):
        idx = rng.choice(n, size=n, replace=True)
        boot_med = _weighted_median(thetas[idx], weights[idx])
        boot_meds.append(boot_med)

    boot_meds = np.array(boot_meds)
    se_boot = float(np.std(boot_meds, ddof=1))
    ci_lo = float(np.percentile(boot_meds, 2.5))
    ci_hi = float(np.percentile(boot_meds, 97.5))

    return {"theta": round(med, 6), "se": round(se_boot, 6),
            "ci_lo": round(ci_lo, 6), "ci_hi": round(ci_hi, 6)}


def _huber_pool(thetas, ses, c=1.345, max_iter=50, tol=1e-6):
    """Huber M-estimation via iteratively reweighted least squares.

    Parameters
    ----------
    c : float
        Huber tuning constant (1.345 gives 95% efficiency under normality).
    """
    weights_base = 1.0 / (ses ** 2)

    # Initialize with weighted median
    theta_hat = _weighted_median(thetas, weights_base)

    converged = False
    outlier_indices = set()

    for iteration in range(1, max_iter + 1):
        # Standardized residuals
        resid = (thetas - theta_hat) / ses
        abs_resid = np.abs(resid)
        # Huber weight function (guard division by zero when residual is exactly 0)
        huber_w = np.where(abs_resid <= c, 1.0, c / np.where(abs_resid > 0, abs_resid, 1.0))
        # Combined weights
        w_combined = huber_w * weights_base

        # Update estimate
        theta_new = float(np.sum(w_combined * thetas) / np.sum(w_combined))

        if abs(theta_new - theta_hat) < tol:
            theta_hat = theta_new
            converged = True
            break
        theta_hat = theta_new

    # Final residuals for outlier flagging
    final_resid = np.abs((thetas - theta_hat) / ses)
    outlier_indices = set(np.where(final_resid > c)[0])

    # SE estimate: sandwich-type
    # se = 1 / sqrt(sum(huber_w * weights_base))
    final_huber_w = np.where(final_resid <= c, 1.0, c / np.where(final_resid > 0, final_resid, 1.0))
    w_final = final_huber_w * weights_base
    se_hat = float(1.0 / np.sqrt(np.sum(w_final))) if np.sum(w_final) > 0 else 0.0

    ci_lo = theta_hat - 1.96 * se_hat
    ci_hi = theta_hat + 1.96 * se_hat

    return {
        "theta": round(theta_hat, 6),
        "se": round(se_hat, 6),
        "ci_lo": round(ci_lo, 6),
        "ci_hi": round(ci_hi, 6),
        "iterations": iteration,
        "converged": converged,
    }, outlier_indices


def _winsorized_pool(thetas, ses):
    """Winsorized mean: replace extremes beyond +/- 2*IQR, then IV pool."""
    q1, q3 = np.percentile(thetas, [25, 75])
    iqr = q3 - q1
    lower = q1 - 2.0 * iqr
    upper = q3 + 2.0 * iqr

    thetas_w = np.clip(thetas, lower, upper)

    w = 1.0 / (ses ** 2)
    theta_hat = float(np.sum(w * thetas_w) / np.sum(w))
    se_hat = float(1.0 / np.sqrt(np.sum(w)))
    ci_lo = theta_hat - 1.96 * se_hat
    ci_hi = theta_hat + 1.96 * se_hat

    return {"theta": round(theta_hat, 6), "se": round(se_hat, 6),
            "ci_lo": round(ci_lo, 6), "ci_hi": round(ci_hi, 6)}


def _compute_breakdown_point(thetas, ses, theta_standard):
    """Compute breakdown point: proportion of reviews that can be replaced
    before the estimate changes by >50%.

    Uses a simulation approach: replace reviews one by one with extreme values
    and check when the estimate shifts by more than 50% of the original.
    """
    n = len(thetas)
    if abs(theta_standard) < 1e-10:
        threshold = 0.5  # absolute threshold when estimate is near zero
    else:
        threshold = abs(theta_standard) * 0.5

    # Sort reviews by influence (largest weight first = most influential)
    weights = 1.0 / (ses ** 2)
    order = np.argsort(-weights)  # descending by weight

    for count in range(1, n):
        # Replace 'count' most influential reviews with extreme value
        thetas_mod = thetas.copy()
        extreme = 100.0 * np.sign(theta_standard) if theta_standard != 0 else 100.0
        for i in range(count):
            thetas_mod[order[i]] = extreme

        w = 1.0 / (ses ** 2)
        theta_mod = float(np.sum(w * thetas_mod) / np.sum(w))

        if abs(theta_mod - theta_standard) > threshold:
            return round((count - 1) / n, 4) if count > 1 else 0.0

    return round((n - 1) / n, 4)


def compute_robust_pooling(reviews, seed=42):
    """Compare standard IV, median, Huber, and winsorized pooling.

    Parameters
    ----------
    reviews : list[ReviewInput]
    seed : int
        Random seed for bootstrap.

    Returns
    -------
    dict with keys: standard, median, huber, winsorized, breakdown_point,
                    outlier_flags
    """
    if len(reviews) < 2:
        raise ValueError("Need at least 2 reviews for pooling")

    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)

    # Standard IV
    standard = _iv_pool(thetas, ses)

    # Weighted median with bootstrap CI
    median_result = _median_pool(thetas, ses, seed=seed)

    # Huber M-estimation
    huber_result, outlier_idx = _huber_pool(thetas, ses)

    # Winsorized mean
    winsorized = _winsorized_pool(thetas, ses)

    # Outlier flags
    outlier_flags = [reviews[i].review_id for i in sorted(outlier_idx)]

    # Breakdown point
    breakdown_point = _compute_breakdown_point(thetas, ses, standard["theta"])

    return {
        "standard": standard,
        "median": median_result,
        "huber": huber_result,
        "winsorized": winsorized,
        "breakdown_point": breakdown_point,
        "outlier_flags": outlier_flags,
    }
