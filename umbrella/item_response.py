"""Item Response Theory (Rasch model) for AMSTAR-2 quality assessment."""

import math
import numpy as np


def _score_item(value):
    """Convert AMSTAR item value to numeric: Yes=1, Partial Yes=0.5, No=0."""
    v = str(value).strip().lower()
    if v == "yes":
        return 1.0
    if v in ("partial yes", "partial"):
        return 0.5
    return 0.0


def _logit(p, eps=1e-6):
    p = max(eps, min(1.0 - eps, p))
    return math.log(p / (1.0 - p))


def _expit(x):
    # Numerically stable sigmoid: avoid overflow for large |x|
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    else:
        ex = math.exp(x)
        return ex / (1.0 + ex)


def rasch_irt(reviews, max_iter=50, tol=0.001):
    """Fit a Rasch IRT model to AMSTAR-2 item responses.

    Parameters
    ----------
    reviews : list[ReviewInput]
        Reviews with amstar_items dicts (keys = item IDs, values = Yes/No/Partial Yes).
    max_iter : int
        Maximum JMLE iterations.
    tol : float
        Convergence tolerance on max absolute parameter change.

    Returns
    -------
    dict with keys: theta_estimates, beta_estimates, reliability, outfit, infit,
                    item_difficulty_order, converged
    """
    # Extract response matrix
    reviews_with_items = [r for r in reviews if r.amstar_items]
    if len(reviews_with_items) < 2:
        return {
            "theta_estimates": {},
            "beta_estimates": {},
            "reliability": 0.0,
            "outfit": {},
            "infit": {},
            "item_difficulty_order": [],
            "converged": False,
        }

    # Collect all item IDs across reviews
    all_items = set()
    for r in reviews_with_items:
        all_items.update(r.amstar_items.keys())
    all_items = sorted(all_items)

    n_persons = len(reviews_with_items)
    n_items = len(all_items)
    item_idx = {item: j for j, item in enumerate(all_items)}

    # Build response matrix (n_persons x n_items), NaN for missing
    X = np.full((n_persons, n_items), np.nan)
    for i, r in enumerate(reviews_with_items):
        for item_key, item_val in r.amstar_items.items():
            if item_key in item_idx:
                X[i, item_idx[item_key]] = _score_item(item_val)

    # Initialize theta (person ability) and beta (item difficulty)
    person_means = np.nanmean(X, axis=1)
    person_means = np.where(np.isfinite(person_means), person_means, 0.5)
    theta = np.array([_logit(pm) for pm in person_means])

    item_pass_rates = np.nanmean(X, axis=0)
    item_pass_rates = np.where(np.isfinite(item_pass_rates), item_pass_rates, 0.5)
    beta = np.array([-_logit(pr) for pr in item_pass_rates])

    # JMLE iteration
    converged = False
    for iteration in range(max_iter):
        old_theta = theta.copy()
        old_beta = beta.copy()

        # Update theta given beta
        for i in range(n_persons):
            observed = []
            expected = []
            info = []
            for j in range(n_items):
                if np.isfinite(X[i, j]):
                    p = _expit(theta[i] - beta[j])
                    observed.append(X[i, j])
                    expected.append(p)
                    info.append(p * (1.0 - p))
            if info:
                total_info = sum(info)
                if total_info > 1e-8:
                    r_i = sum(observed)
                    e_i = sum(expected)
                    theta[i] += (r_i - e_i) / total_info

        # Clamp theta to prevent divergence with near-perfect/zero scores
        theta = np.clip(theta, -10.0, 10.0)

        # Center theta (identification constraint)
        theta -= np.mean(theta)

        # Update beta given theta
        for j in range(n_items):
            observed = []
            expected = []
            info = []
            for i in range(n_persons):
                if np.isfinite(X[i, j]):
                    p = _expit(theta[i] - beta[j])
                    observed.append(X[i, j])
                    expected.append(p)
                    info.append(p * (1.0 - p))
            if info:
                total_info = sum(info)
                if total_info > 1e-8:
                    s_j = sum(observed)
                    e_j = sum(expected)
                    # For beta, the Newton-Raphson update has opposite sign
                    beta[j] -= (s_j - e_j) / total_info

        # Clamp beta to prevent divergence
        beta = np.clip(beta, -10.0, 10.0)

        # Check convergence
        max_change = max(
            np.max(np.abs(theta - old_theta)),
            np.max(np.abs(beta - old_beta)),
        )
        if max_change < tol:
            converged = True
            break

    # Compute fit statistics (outfit and infit)
    outfit = {}
    infit = {}
    for j in range(n_items):
        item_name = all_items[j]
        z_sq_list = []
        info_list = []
        for i in range(n_persons):
            if np.isfinite(X[i, j]):
                p = _expit(theta[i] - beta[j])
                variance = p * (1.0 - p)
                if variance > 1e-10:
                    residual = X[i, j] - p
                    z_sq = (residual ** 2) / variance
                    z_sq_list.append(z_sq)
                    info_list.append(variance)
        if z_sq_list:
            outfit[item_name] = float(np.mean(z_sq_list))
            total_info = sum(info_list)
            if total_info > 1e-10:
                infit[item_name] = float(
                    sum(z * w for z, w in zip(z_sq_list, info_list)) / total_info
                )
            else:
                infit[item_name] = 1.0
        else:
            outfit[item_name] = 1.0
            infit[item_name] = 1.0

    # Reliability: (var(theta) - mean(SE^2)) / var(theta)
    se_sq_list = []
    for i in range(n_persons):
        info_sum = 0.0
        for j in range(n_items):
            if np.isfinite(X[i, j]):
                p = _expit(theta[i] - beta[j])
                info_sum += p * (1.0 - p)
        if info_sum > 1e-10:
            se_sq_list.append(1.0 / info_sum)
        else:
            se_sq_list.append(0.0)
    var_theta = float(np.var(theta))
    mean_se_sq = float(np.mean(se_sq_list)) if se_sq_list else 0.0
    if var_theta > 1e-10:
        reliability = max(0.0, (var_theta - mean_se_sq) / var_theta)
    else:
        reliability = 0.0

    # Theta estimates dict
    theta_estimates = {
        reviews_with_items[i].review_id: float(theta[i]) for i in range(n_persons)
    }

    # Beta estimates dict
    beta_estimates = {all_items[j]: float(beta[j]) for j in range(n_items)}

    # Item difficulty order (easiest to hardest = lowest beta to highest)
    item_difficulty_order = sorted(all_items, key=lambda it: beta_estimates[it])

    return {
        "theta_estimates": theta_estimates,
        "beta_estimates": beta_estimates,
        "reliability": round(reliability, 4),
        "outfit": outfit,
        "infit": infit,
        "item_difficulty_order": item_difficulty_order,
        "converged": converged,
    }
