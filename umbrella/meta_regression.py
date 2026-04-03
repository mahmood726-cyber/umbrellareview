"""Meta-Meta-Regression: Explain between-review heterogeneity with covariates."""

import numpy as np
from umbrella.amstar import score_amstar

AMSTAR_QUALITY_MAP = {
    "High": 4,
    "Moderate": 3,
    "Low": 2,
    "Critically Low": 1,
}


def _build_covariates(reviews):
    """Build covariate matrix and labels from reviews.

    Covariates: year, AMSTAR quality score, k, scope breadth.
    Returns (X, covariate_names, valid_mask) where X includes intercept column.
    """
    n = len(reviews)
    years = []
    quality = []
    ks = []
    breadth = []
    valid = []

    for r in reviews:
        yr = r.year if r.year > 0 else None
        # AMSTAR quality
        if r.amstar_items:
            res = score_amstar(r.review_id, r.amstar_items)
            q = AMSTAR_QUALITY_MAP.get(res.confidence, None)
        else:
            q = None
        k_val = r.k if r.k > 0 else None
        scope = len(set(r.scope_tags)) if r.scope_tags else None

        if yr is not None and q is not None and k_val is not None and scope is not None:
            years.append(yr)
            quality.append(q)
            ks.append(k_val)
            breadth.append(scope)
            valid.append(True)
        else:
            valid.append(False)

    covariate_names = ["intercept", "year", "amstar_quality", "k", "scope_breadth"]
    if sum(valid) < 3:
        return None, covariate_names, valid

    X = np.column_stack([
        np.ones(sum(valid)),
        np.array(years, dtype=float),
        np.array(quality, dtype=float),
        np.array(ks, dtype=float),
        np.array(breadth, dtype=float),
    ])
    return X, covariate_names, valid


def compute_meta_regression(reviews, n_perm=999, seed=42):
    """Run weighted least-squares meta-meta-regression.

    Parameters
    ----------
    reviews : list[ReviewInput]
    n_perm : int
        Number of permutations for permutation test.
    seed : int
        Random seed for permutation test.

    Returns
    -------
    dict with keys: coefficients, r_squared, tau2_residual, permutation_p,
                    n_reviews_used, residuals
    """
    X, cov_names, valid_mask = _build_covariates(reviews)

    valid_reviews = [r for r, v in zip(reviews, valid_mask) if v]
    n_used = len(valid_reviews)

    if X is None or n_used < 3:
        raise ValueError(f"Need at least 3 reviews with complete covariates, got {n_used}")

    thetas = np.array([r.theta for r in valid_reviews])
    ses = np.array([r.se for r in valid_reviews])
    ses = np.where(ses > 0, ses, 0.1)

    weights = 1.0 / (ses ** 2)
    W = np.diag(weights)

    # WLS: beta = (X^T W X)^{-1} X^T W Y
    XtWX = X.T @ W @ X
    try:
        XtWX_inv = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        XtWX_inv = np.linalg.pinv(XtWX)

    beta = XtWX_inv @ X.T @ W @ thetas
    se_beta = np.sqrt(np.diag(XtWX_inv))

    # Fitted and residuals
    fitted = X @ beta
    residuals = thetas - fitted

    # Q_total and Q_residual
    theta_bar = float(np.sum(weights * thetas) / np.sum(weights))
    q_total = float(np.sum(weights * (thetas - theta_bar) ** 2))
    q_residual = float(np.sum(weights * residuals ** 2))

    # R^2 analogue
    r_squared = max(0.0, 1.0 - q_residual / q_total) if q_total > 0 else 0.0

    # Residual tau^2 via DerSimonian-Laird on residuals
    n = n_used
    p = X.shape[1]
    df_resid = n - p
    w_sum = float(np.sum(weights))
    w2_sum = float(np.sum(weights ** 2))
    c = w_sum - w2_sum / w_sum if w_sum > 0 else 1.0
    tau2_resid = max(0.0, (q_residual - df_resid) / c) if c > 0 and df_resid > 0 else 0.0

    # Coefficients dict
    coefficients = {}
    for i, name in enumerate(cov_names):
        z = beta[i] / se_beta[i] if se_beta[i] > 0 else 0.0
        p_val = float(2 * (1 - _normal_cdf(abs(z))))
        coefficients[name] = {
            "beta": round(float(beta[i]), 6),
            "se": round(float(se_beta[i]), 6),
            "z": round(float(z), 4),
            "p": round(p_val, 4),
        }

    # Permutation test for R^2
    rng = np.random.RandomState(seed)
    count_ge = 0
    for _ in range(n_perm):
        y_perm = rng.permutation(thetas)
        beta_perm = XtWX_inv @ X.T @ W @ y_perm
        fitted_perm = X @ beta_perm
        resid_perm = y_perm - fitted_perm
        q_resid_perm = float(np.sum(weights * resid_perm ** 2))
        # Recompute Q_total for permuted Y
        ybar_perm = float(np.sum(weights * y_perm) / np.sum(weights))
        q_total_perm = float(np.sum(weights * (y_perm - ybar_perm) ** 2))
        r2_perm = max(0.0, 1.0 - q_resid_perm / q_total_perm) if q_total_perm > 0 else 0.0
        if r2_perm >= r_squared:
            count_ge += 1

    permutation_p = (count_ge + 1) / (n_perm + 1)

    return {
        "coefficients": coefficients,
        "r_squared": round(r_squared, 4),
        "tau2_residual": round(tau2_resid, 6),
        "permutation_p": round(permutation_p, 4),
        "n_reviews_used": n_used,
        "residuals": [round(float(r), 6) for r in residuals],
    }


def _normal_cdf(x):
    """Standard normal CDF via error function."""
    import math
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2)))
