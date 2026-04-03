"""Prediction Intervals and Probability Statements.

Computes probabilistic statements about future reviews using
the meta-meta posterior distribution: direction agreement probability,
threshold exceedance, predictive intervals, sign reversal probability,
and a fragility index.
"""

import math
import numpy as np
from scipy import stats


def compute_predictions(reviews, threshold=0.0):
    """Compute prediction intervals and probability statements.

    Parameters
    ----------
    reviews : list[ReviewInput]
        At least 2 reviews with theta and se > 0.
    threshold : float
        Effect threshold for P(true effect > threshold). Default 0.

    Returns
    -------
    dict with keys: p_agree_direction, p_exceeds_threshold,
        predictive_interval, p_sign_reversal, fragility_index
    """
    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)
    n = len(reviews)

    if n < 2:
        raise ValueError("Need at least 2 reviews for prediction")

    variances = ses ** 2

    # --- DerSimonian-Laird random-effects ---
    w = 1.0 / variances
    mu_fe = float(np.sum(w * thetas) / np.sum(w))
    q = float(np.sum(w * (thetas - mu_fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau_mm2 = max(0.0, (q - (n - 1)) / c) if c > 0 else 0.0

    w_re = 1.0 / (variances + tau_mm2)
    mu = float(np.sum(w_re * thetas) / np.sum(w_re))
    se_mu = float(1.0 / np.sqrt(np.sum(w_re)))

    # Majority direction
    majority_sign = np.sign(np.median(thetas))
    if majority_sign == 0:
        majority_sign = -1.0

    # --- Predictive distribution for next review ---
    # Var(next) = tau_mm2 + se_mu^2 + median(se^2) for a "typical" new review
    typical_var = float(np.median(variances))
    pred_var = tau_mm2 + se_mu ** 2
    pred_se = math.sqrt(pred_var) if pred_var > 0 else 1e-6

    # Use t-distribution with n-1 df for prediction
    df = max(n - 1, 1)
    t_crit = float(stats.t.ppf(0.975, df=df))
    predictive_interval = (
        round(mu - t_crit * pred_se, 6),
        round(mu + t_crit * pred_se, 6),
    )

    # P(next review agrees with majority direction)
    # Using the predictive t-distribution
    if majority_sign > 0:
        # P(next > 0)
        t_val = (0.0 - mu) / pred_se
        p_agree_direction = float(1.0 - stats.t.cdf(t_val, df=df))
    else:
        # P(next < 0)
        t_val = (0.0 - mu) / pred_se
        p_agree_direction = float(stats.t.cdf(t_val, df=df))

    # P(true effect > threshold) via posterior Normal approximation
    if se_mu > 0:
        z_thresh = (threshold - mu) / se_mu
        p_exceeds_threshold = float(1.0 - stats.norm.cdf(z_thresh))
    else:
        p_exceeds_threshold = 1.0 if mu > threshold else 0.0

    # P(sign reversal): probability the next review has opposite sign to majority
    p_sign_reversal = round(1.0 - p_agree_direction, 6)

    # Fragility index: minimum number of reviews that would need to change
    # direction to flip the majority conclusion
    n_agree = int(np.sum(np.sign(thetas) == majority_sign))
    n_disagree = n - n_agree
    # Need to flip enough to make disagree >= agree
    # After flipping f reviews: agree-f vs disagree+f
    # Need agree-f <= disagree+f => f >= (agree - disagree) / 2
    fragility_index = max(1, math.ceil((n_agree - n_disagree + 1) / 2))
    # But cannot exceed n_agree (can't flip more than exist)
    fragility_index = min(fragility_index, n_agree)

    return {
        "p_agree_direction": round(float(p_agree_direction), 6),
        "p_exceeds_threshold": round(float(p_exceeds_threshold), 6),
        "predictive_interval": predictive_interval,
        "p_sign_reversal": round(float(p_sign_reversal), 6),
        "fragility_index": fragility_index,
    }
