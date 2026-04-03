"""Causal Modeling of Review Disagreement.

Model WHY systematic reviews disagree using structural equations,
counterfactual analysis, and mediation decomposition (Baron-Kenny).
The causal graph assumes:

    methodology -> effect_estimate
    quality     -> effect_estimate
    overlap     -> effect_estimate
    scope       -> inclusion -> overlap

Each review is encoded as a feature vector and theta is regressed
on the causal predictors via OLS.
"""

import math
import numpy as np


def _quality_score(amstar_items):
    """Convert AMSTAR-2 item dict to 1-4 quality score.

    High=4 (>=12 yes), Moderate=3 (>=8), Low=2 (>=4), CriticallyLow=1.
    """
    n_yes = sum(1 for v in amstar_items.values() if v == "yes")
    if n_yes >= 12:
        return 4
    elif n_yes >= 8:
        return 3
    elif n_yes >= 4:
        return 2
    else:
        return 1


def _mean_jaccard(i, reviews, sets):
    """Mean Jaccard similarity of review i with all others."""
    n = len(reviews)
    if n <= 1:
        return 0.0
    total = 0.0
    for j in range(n):
        if j == i:
            continue
        inter = len(sets[i] & sets[j])
        union = len(sets[i] | sets[j])
        total += (inter / union) if union > 0 else 0.0
    return total / (n - 1)


def _ols(X, y):
    """OLS fit: beta = (X^T X)^{-1} X^T y.

    Returns beta, residuals, and R-squared.
    """
    XtX = X.T @ X
    Xty = X.T @ y
    # Use pseudo-inverse for numerical safety
    beta = np.linalg.lstsq(XtX, Xty, rcond=None)[0]
    y_hat = X @ beta
    residuals = y - y_hat
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r_sq = 1.0 - (ss_res / ss_tot) if ss_tot > 1e-15 else 0.0
    r_sq = max(0.0, min(1.0, r_sq))
    return beta, residuals, r_sq, y_hat


def _t_test_p(beta_j, se_j, df):
    """Two-sided p-value for t = beta_j / se_j with df degrees of freedom.

    Uses the incomplete beta function approximation via scipy if available,
    otherwise a normal approximation for df > 30.
    """
    if se_j < 1e-15:
        return 0.0
    t_stat = abs(beta_j / se_j)
    try:
        from scipy.stats import t as t_dist
        return float(2.0 * t_dist.sf(t_stat, df))
    except ImportError:
        # Normal approximation
        from scipy.stats import norm
        return float(2.0 * norm.sf(t_stat))


def causal_discordance(reviews):
    """Causal structural equation model of review disagreement.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 3 reviews (need degrees of freedom for regression).

    Returns
    -------
    dict with keys:
        structural_coefficients – dict variable -> {beta, se, p}
        counterfactual_high_quality – dict {mean_theta, delta}
        mediation – dict {total, direct, indirect, prop_mediated}
        r_squared – float
        residual_discordance – float (RMSE of residuals)
    """
    n = len(reviews)
    if n < 3:
        raise ValueError("Need at least 3 reviews for causal discordance")

    sets = [set(r.study_ids) for r in reviews]

    # ------------------------------------------------------------------
    # 1. Feature encoding
    # ------------------------------------------------------------------
    years = [r.year for r in reviews]
    mean_year = sum(years) / n

    quality = np.array([_quality_score(r.amstar_items) for r in reviews], dtype=float)
    year_centered = np.array([r.year - mean_year for r in reviews], dtype=float)
    k_studies = np.array([r.k for r in reviews], dtype=float)
    scope_breadth = np.array([len(r.scope_tags) for r in reviews], dtype=float)
    overlap = np.array([_mean_jaccard(i, reviews, sets) for i in range(n)], dtype=float)
    theta = np.array([r.theta for r in reviews], dtype=float)

    # ------------------------------------------------------------------
    # 2. Structural equation: theta = b0 + b1*quality + b2*k + b3*scope + b4*overlap + eps
    # ------------------------------------------------------------------
    # Design matrix with intercept
    X = np.column_stack([
        np.ones(n),
        quality,
        k_studies,
        scope_breadth,
        overlap,
    ])
    var_names = ["intercept", "quality", "k", "scope", "overlap"]

    beta, residuals, r_sq, y_hat = _ols(X, theta)

    # Standard errors of coefficients
    df = max(1, n - X.shape[1])
    mse = float(np.sum(residuals ** 2)) / df
    try:
        cov_beta = mse * np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        cov_beta = mse * np.linalg.pinv(X.T @ X)

    se_beta = np.sqrt(np.maximum(np.diag(cov_beta), 0.0))

    structural_coefficients = {}
    for idx, name in enumerate(var_names):
        b = float(beta[idx])
        se = float(se_beta[idx])
        p = _t_test_p(b, se, df)
        structural_coefficients[name] = {
            "beta": round(b, 6),
            "se": round(se, 6),
            "p": round(p, 6),
        }

    # ------------------------------------------------------------------
    # 3. Counterfactual: if all reviews had High quality (=4)
    # ------------------------------------------------------------------
    X_cf = X.copy()
    X_cf[:, 1] = 4.0  # set quality column to 4
    theta_cf = X_cf @ beta
    mean_theta_cf = float(np.mean(theta_cf))
    mean_theta_obs = float(np.mean(theta))
    delta_cf = mean_theta_cf - mean_theta_obs

    counterfactual = {
        "mean_theta": round(mean_theta_cf, 6),
        "delta": round(delta_cf, 6),
    }

    # ------------------------------------------------------------------
    # 4. Mediation analysis (Baron-Kenny)
    #    Path: quality -> overlap -> theta
    #    a: quality predicts overlap
    #    b: overlap predicts theta (controlling quality)
    #    c: total effect of quality on theta
    #    c': direct effect (controlling overlap)
    #    indirect = a * b
    #    prop_mediated = a*b / c
    # ------------------------------------------------------------------
    # Path c: total effect (quality -> theta, simple regression)
    X_c = np.column_stack([np.ones(n), quality])
    beta_c, _, _, _ = _ols(X_c, theta)
    c_total = float(beta_c[1])

    # Path a: quality -> overlap
    X_a = np.column_stack([np.ones(n), quality])
    beta_a, _, _, _ = _ols(X_a, overlap)
    a_coeff = float(beta_a[1])

    # Path b and c': theta ~ quality + overlap
    X_bc = np.column_stack([np.ones(n), quality, overlap])
    beta_bc, _, _, _ = _ols(X_bc, theta)
    c_prime = float(beta_bc[1])  # direct effect
    b_coeff = float(beta_bc[2])  # overlap -> theta controlling quality

    indirect = a_coeff * b_coeff
    prop_med = (indirect / c_total) if abs(c_total) > 1e-12 else 0.0

    mediation = {
        "total": round(c_total, 6),
        "direct": round(c_prime, 6),
        "indirect": round(indirect, 6),
        "prop_mediated": round(prop_med, 6),
    }

    # ------------------------------------------------------------------
    # 5. Residual discordance (RMSE)
    # ------------------------------------------------------------------
    rmse = float(np.sqrt(np.mean(residuals ** 2)))

    return {
        "structural_coefficients": structural_coefficients,
        "counterfactual_high_quality": counterfactual,
        "mediation": mediation,
        "r_squared": round(r_sq, 6),
        "residual_discordance": round(rmse, 6),
    }
