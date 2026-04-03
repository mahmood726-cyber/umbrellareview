"""Leave-One-Review-Out Influence Diagnostics for Umbrella Reviews."""

import numpy as np
from itertools import combinations


def _iv_pool(thetas, ses):
    """Inverse-variance fixed-effect pooling. Returns (theta, se, Q, I2)."""
    w = 1.0 / (ses ** 2)
    theta_hat = float(np.sum(w * thetas) / np.sum(w))
    se_hat = float(1.0 / np.sqrt(np.sum(w)))
    n = len(thetas)
    q = float(np.sum(w * (thetas - theta_hat) ** 2))
    i2 = max(0.0, (q - (n - 1)) / q * 100) if q > 0 and n > 1 else 0.0
    return theta_hat, se_hat, q, i2


def _re_pool(thetas, ses):
    """DerSimonian-Laird random-effects pooling. Returns (theta, se)."""
    w = 1.0 / (ses ** 2)
    n = len(thetas)
    theta_fe = float(np.sum(w * thetas) / np.sum(w))
    q = float(np.sum(w * (thetas - theta_fe) ** 2))
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    tau2 = max(0.0, (q - (n - 1)) / c) if c > 0 and n > 1 else 0.0
    w_re = 1.0 / (ses ** 2 + tau2)
    theta_re = float(np.sum(w_re * thetas) / np.sum(w_re))
    se_re = float(1.0 / np.sqrt(np.sum(w_re)))
    return theta_re, se_re


def compute_influence(reviews):
    """Compute leave-one-review-out influence diagnostics.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with keys: loo_estimates, cooks_d, dfbetas, i2_influence,
                    most_influential, tipping_point
    """
    n = len(reviews)
    if n < 2:
        raise ValueError("Need at least 2 reviews for influence analysis")

    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)

    # Full pooling
    theta_full_fe, se_full_fe, _, i2_full = _iv_pool(thetas, ses)
    theta_full_re, se_full_re = _re_pool(thetas, ses)

    loo_estimates = []
    cooks_d = {}
    dfbetas_dict = {}
    i2_influence = {}

    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        t_loo = thetas[mask]
        s_loo = ses[mask]

        theta_fe_i, se_fe_i, _, i2_i = _iv_pool(t_loo, s_loo)
        theta_re_i, se_re_i = _re_pool(t_loo, s_loo)

        rid = reviews[i].review_id

        loo_estimates.append({
            "review_id": rid,
            "theta_without_fe": round(theta_fe_i, 6),
            "se_without_fe": round(se_fe_i, 6),
            "theta_without_re": round(theta_re_i, 6),
            "se_without_re": round(se_re_i, 6),
        })

        # Cook's distance analogue (FE-based)
        var_full = se_full_fe ** 2
        if var_full > 0:
            d_i = (theta_full_fe - theta_fe_i) ** 2 / var_full
        else:
            d_i = 0.0
        cooks_d[rid] = round(d_i, 6)

        # DFBETAS (FE-based)
        if se_full_fe > 0:
            dfb = (theta_full_fe - theta_fe_i) / se_full_fe
        else:
            dfb = 0.0
        dfbetas_dict[rid] = round(dfb, 6)

        # Influence on I^2
        i2_influence[rid] = round(i2_full - i2_i, 4)

    # Most influential (largest Cook's D)
    most_influential = max(cooks_d, key=cooks_d.get)

    # Tipping point: minimum removals to change sign of pooled FE estimate
    sign_full = np.sign(theta_full_fe)
    if sign_full == 0:
        tipping = 0
    else:
        tipping = _find_tipping_point(thetas, ses, sign_full, reviews, max_exact=3)

    return {
        "loo_estimates": loo_estimates,
        "cooks_d": cooks_d,
        "dfbetas": dfbetas_dict,
        "i2_influence": i2_influence,
        "most_influential": most_influential,
        "tipping_point": tipping,
    }


def _find_tipping_point(thetas, ses, sign_full, reviews, max_exact=3):
    """Find minimum number of reviews to remove to change sign of pooled estimate.

    Try all combinations up to size max_exact, then greedy after that.
    """
    n = len(thetas)

    # Try exact subsets of size 1..max_exact
    for size in range(1, min(max_exact + 1, n)):
        for combo in combinations(range(n), size):
            mask = np.ones(n, dtype=bool)
            for idx in combo:
                mask[idx] = False
            if np.sum(mask) < 1:
                continue
            t_rem = thetas[mask]
            s_rem = ses[mask]
            theta_rem, _, _, _ = _iv_pool(t_rem, s_rem)
            if np.sign(theta_rem) != sign_full or theta_rem == 0:
                return size

    # Greedy: iteratively remove the review that shifts estimate most
    # toward the opposite sign
    remaining = list(range(n))
    for removed_count in range(max_exact + 1, n):
        best_idx = None
        best_theta = None
        for i_pos, i in enumerate(remaining):
            trial = [j for j in remaining if j != i]
            if len(trial) < 1:
                continue
            t_trial = thetas[np.array(trial)]
            s_trial = ses[np.array(trial)]
            theta_trial, _, _, _ = _iv_pool(t_trial, s_trial)
            # Pick the removal that pushes estimate furthest from original sign
            if best_theta is None:
                best_theta = theta_trial
                best_idx = i
            else:
                # If sign_full < 0, we want the largest theta_trial
                if sign_full < 0 and theta_trial > best_theta:
                    best_theta = theta_trial
                    best_idx = i
                elif sign_full > 0 and theta_trial < best_theta:
                    best_theta = theta_trial
                    best_idx = i

        if best_idx is not None:
            remaining.remove(best_idx)
            if np.sign(best_theta) != sign_full or best_theta == 0:
                return removed_count

    # Cannot change sign even removing all but 1
    return n
