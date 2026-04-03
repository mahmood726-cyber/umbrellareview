"""Publication Bias at the Systematic Review Level.

Detect whether the set of SRs itself is subject to publication/reporting
bias, using adapted funnel-plot methods: Egger's test, trim-and-fill,
Peters' test, and Begg's rank correlation (Kendall tau).
"""

import math
import numpy as np


def _kendall_tau(x, y):
    """Kendall's tau-b with z-test for significance.

    Returns (tau, z, p_two_sided).
    """
    n = len(x)
    if n < 2:
        return 0.0, 0.0, 1.0

    concordant = 0
    discordant = 0
    for i in range(n):
        for j in range(i + 1, n):
            dx = x[i] - x[j]
            dy = y[i] - y[j]
            prod = dx * dy
            if prod > 0:
                concordant += 1
            elif prod < 0:
                discordant += 1
            # ties -> neither

    n_pairs = n * (n - 1) // 2
    if n_pairs == 0:
        return 0.0, 0.0, 1.0

    tau = (concordant - discordant) / n_pairs

    # Variance under H0 (accounting for no ties in simplified form)
    var_tau = (2.0 * (2.0 * n + 5.0)) / (9.0 * n * (n - 1.0))
    se_tau = math.sqrt(var_tau) if var_tau > 0 else 1e-15

    z = tau / se_tau
    # Two-sided p via normal approximation
    p = 2.0 * _normal_sf(abs(z))

    return tau, z, p


def _normal_sf(z):
    """Survival function (1 - CDF) of standard normal, via error function."""
    return 0.5 * math.erfc(z / math.sqrt(2.0))


def _ols_with_t(X, y):
    """OLS returning beta, se, t-stats, p-values (two-sided)."""
    n, p = X.shape
    XtX = X.T @ X
    try:
        XtX_inv = np.linalg.inv(XtX)
    except np.linalg.LinAlgError:
        XtX_inv = np.linalg.pinv(XtX)

    beta = XtX_inv @ (X.T @ y)
    residuals = y - X @ beta
    df = max(1, n - p)
    mse = float(np.sum(residuals ** 2)) / df
    cov_beta = mse * XtX_inv
    se = np.sqrt(np.maximum(np.diag(cov_beta), 0.0))

    # t-statistics and p-values
    t_stats = np.zeros(p)
    p_vals = np.zeros(p)
    for j in range(p):
        if se[j] > 1e-15:
            t_stats[j] = beta[j] / se[j]
            try:
                from scipy.stats import t as t_dist
                p_vals[j] = 2.0 * t_dist.sf(abs(t_stats[j]), df)
            except ImportError:
                p_vals[j] = 2.0 * _normal_sf(abs(t_stats[j]))
        else:
            t_stats[j] = 0.0
            p_vals[j] = 1.0

    return beta, se, t_stats, p_vals, df


def funnel_meta(reviews):
    """Publication bias assessment at the systematic review level.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 3 reviews with theta and se.

    Returns
    -------
    dict with keys:
        funnel_x – list[float] (theta values)
        funnel_y – list[float] (precision = 1/se)
        egger_intercept – float
        egger_slope – float
        egger_p – float
        n_missing_trimfill – int
        adjusted_theta – float
        adjusted_se – float
        peters_intercept – float
        peters_p – float
        begg_tau – float
        begg_p – float
        bias_detected – bool (any test p < 0.1)
    """
    n = len(reviews)
    if n < 3:
        raise ValueError("Need at least 3 reviews for funnel meta-analysis")

    theta = np.array([r.theta for r in reviews], dtype=float)
    se = np.array([r.se for r in reviews], dtype=float)

    # Ensure no zero SE
    se = np.where(se > 1e-15, se, 1e-6)

    # ------------------------------------------------------------------
    # 1. Funnel plot data
    # ------------------------------------------------------------------
    precision = 1.0 / se
    funnel_x = [round(float(t), 6) for t in theta]
    funnel_y = [round(float(p), 6) for p in precision]

    # ------------------------------------------------------------------
    # 2. Egger's test: regress (theta / se) on (1 / se)
    #    z_i = theta_i/se_i = intercept + slope * (1/se_i)
    #    Test H0: intercept = 0
    # ------------------------------------------------------------------
    z_scores = theta / se
    X_egger = np.column_stack([np.ones(n), precision])
    beta_eg, se_eg, t_eg, p_eg, df_eg = _ols_with_t(X_egger, z_scores)

    egger_intercept = float(beta_eg[0])
    egger_slope = float(beta_eg[1])
    egger_p = float(p_eg[0])

    # ------------------------------------------------------------------
    # 3. Trim and fill (L0 estimator, simplified)
    # ------------------------------------------------------------------
    # Pooled estimate (inverse-variance weighted)
    w = 1.0 / (se ** 2)
    theta_pooled = float(np.sum(w * theta) / np.sum(w))

    # Count reviews on each side of pooled estimate
    deviations = theta - theta_pooled
    n_left = int(np.sum(deviations < 0))
    n_right = int(np.sum(deviations >= 0))
    n_missing = max(0, abs(n_right - n_left))

    # Determine which side is sparse
    if n_right > n_left:
        # Missing on the left (lower theta side)
        sparse_side = "left"
    else:
        # Missing on the right
        sparse_side = "right"

    # Fill missing reviews by reflecting across pooled estimate
    filled_theta = list(theta)
    filled_se = list(se)

    if n_missing > 0:
        # Sort reviews by distance from pooled on the heavy side
        if sparse_side == "left":
            # Reflect rightmost extreme reviews to the left
            right_idx = np.where(deviations >= 0)[0]
            dists = [(abs(deviations[i]), i) for i in right_idx]
            dists.sort(reverse=True)
            for count, (_, idx) in enumerate(dists):
                if count >= n_missing:
                    break
                mirror_theta = 2.0 * theta_pooled - theta[idx]
                filled_theta.append(mirror_theta)
                filled_se.append(se[idx])
        else:
            # Reflect leftmost extreme reviews to the right
            left_idx = np.where(deviations < 0)[0]
            dists = [(abs(deviations[i]), i) for i in left_idx]
            dists.sort(reverse=True)
            for count, (_, idx) in enumerate(dists):
                if count >= n_missing:
                    break
                mirror_theta = 2.0 * theta_pooled - theta[idx]
                filled_theta.append(mirror_theta)
                filled_se.append(se[idx])

    # Re-pool with filled reviews
    filled_theta = np.array(filled_theta)
    filled_se = np.array(filled_se)
    w_filled = 1.0 / (filled_se ** 2)
    adjusted_theta = float(np.sum(w_filled * filled_theta) / np.sum(w_filled))
    adjusted_se = float(1.0 / np.sqrt(np.sum(w_filled)))

    # ------------------------------------------------------------------
    # 4. Peters' test: regress theta on 1/N_total
    #    Approximate N_total = k * 100 (median study size)
    # ------------------------------------------------------------------
    k_studies = np.array([r.k if r.k > 0 else 1 for r in reviews], dtype=float)
    n_total = k_studies * 100.0  # approximate total sample size
    inv_n = 1.0 / n_total

    X_peters = np.column_stack([np.ones(n), inv_n])
    beta_pt, se_pt, t_pt, p_pt, df_pt = _ols_with_t(X_peters, theta)

    peters_intercept = float(beta_pt[0])
    peters_p = float(p_pt[1])  # test on the slope (1/N coefficient)

    # ------------------------------------------------------------------
    # 5. Begg's rank correlation (Kendall's tau between theta and se)
    # ------------------------------------------------------------------
    begg_tau, begg_z, begg_p = _kendall_tau(
        [float(t) for t in theta],
        [float(s) for s in se],
    )

    # ------------------------------------------------------------------
    # 6. Overall bias detection
    # ------------------------------------------------------------------
    bias_detected = (egger_p < 0.1) or (begg_p < 0.1) or (peters_p < 0.1)

    return {
        "funnel_x": funnel_x,
        "funnel_y": funnel_y,
        "egger_intercept": round(egger_intercept, 6),
        "egger_slope": round(egger_slope, 6),
        "egger_p": round(egger_p, 6),
        "n_missing_trimfill": n_missing,
        "adjusted_theta": round(adjusted_theta, 6),
        "adjusted_se": round(adjusted_se, 6),
        "peters_intercept": round(peters_intercept, 6),
        "peters_p": round(peters_p, 6),
        "begg_tau": round(begg_tau, 6),
        "begg_p": round(begg_p, 6),
        "bias_detected": bias_detected,
    }
