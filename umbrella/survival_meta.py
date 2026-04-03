"""Time-to-event meta-meta-analysis for umbrella reviews of survival outcomes."""

import math
import numpy as np


def _dl_tau2(thetas, ses):
    """DerSimonian-Laird tau-squared estimator."""
    w = 1.0 / (ses ** 2)
    k = len(thetas)
    if k < 2:
        return 0.0
    mu_fe = np.sum(w * thetas) / np.sum(w)
    Q = float(np.sum(w * (thetas - mu_fe) ** 2))
    df = k - 1
    c = float(np.sum(w) - np.sum(w ** 2) / np.sum(w))
    if c <= 0:
        return 0.0
    tau2 = max(0.0, (Q - df) / c)
    return tau2


def _iv_pool(thetas, ses, tau2=0.0):
    """Inverse-variance pooling with optional between-study variance."""
    w = 1.0 / (ses ** 2 + tau2)
    mu = float(np.sum(w * thetas) / np.sum(w))
    se_mu = float(1.0 / np.sqrt(np.sum(w)))
    return mu, se_mu


def _compute_i2(thetas, ses):
    """Compute I-squared statistic."""
    k = len(thetas)
    if k < 2:
        return 0.0
    w = 1.0 / (ses ** 2)
    mu_fe = np.sum(w * thetas) / np.sum(w)
    Q = float(np.sum(w * (thetas - mu_fe) ** 2))
    df = k - 1
    if Q <= 0:
        return 0.0
    return max(0.0, (Q - df) / Q * 100)


def survival_meta(reviews, median_survivals=None, rmst_differences=None,
                   median_ages=None):
    """Pool hazard ratios from SRs reporting survival/time-to-event outcomes.

    Parameters
    ----------
    reviews : list[ReviewInput]
        Each review's theta is interpreted as log(HR), se as its SE.
    median_survivals : list[tuple[float, float]] or None
        Optional list of (treatment_median, control_median) per review.
    rmst_differences : list[float] or None
        Optional RMST differences per review.
    median_ages : list[float] or None
        Optional median ages of patient populations per review (for heterogeneity decomposition).

    Returns
    -------
    dict with keys: pooled_log_hr, pooled_hr, hr_ci, tau2_log_hr, frailty_adjusted_hr,
                    median_survival_ratio, i2, heterogeneity_decomposition
    """
    if not reviews:
        return {
            "pooled_log_hr": 0.0,
            "pooled_hr": 1.0,
            "hr_ci": (1.0, 1.0),
            "tau2_log_hr": 0.0,
            "frailty_adjusted_hr": 1.0,
            "median_survival_ratio": None,
            "i2": 0.0,
            "heterogeneity_decomposition": {"clinical": 0.0, "methodological": 0.0},
        }

    thetas = np.array([r.theta for r in reviews])
    ses = np.array([r.se for r in reviews])
    ses = np.where(ses > 0, ses, 0.1)

    # DL tau2 on log-HR scale
    tau2 = _dl_tau2(thetas, ses)

    # Random-effects pooling
    mu_re, se_re = _iv_pool(thetas, ses, tau2)

    # Fixed-effect pooling (for reference)
    mu_fe, se_fe = _iv_pool(thetas, ses, 0.0)

    pooled_log_hr = mu_re
    pooled_hr = math.exp(mu_re)
    hr_ci_lo = math.exp(mu_re - 1.96 * se_re)
    hr_ci_hi = math.exp(mu_re + 1.96 * se_re)

    # I-squared
    i2 = _compute_i2(thetas, ses)

    # Frailty-adjusted HR: exp(mu + sigma_frailty^2 / 2)
    # sigma_frailty^2 approx tau2 (between-review heterogeneity IS the frailty)
    frailty_adjusted_hr = math.exp(mu_re + tau2 / 2.0)

    # Median survival ratio
    msr = None
    if median_survivals is not None and len(median_survivals) > 0:
        log_msrs = []
        msr_ses = []
        for trt_med, ctrl_med in median_survivals:
            if ctrl_med > 0 and trt_med > 0:
                log_msrs.append(math.log(trt_med / ctrl_med))
                # Approximate SE from delta method; use a rough estimate
                msr_ses.append(0.2)  # placeholder SE when not provided
        if log_msrs:
            log_msrs = np.array(log_msrs)
            msr_ses = np.array(msr_ses)
            tau2_msr = _dl_tau2(log_msrs, msr_ses)
            mu_msr, _ = _iv_pool(log_msrs, msr_ses, tau2_msr)
            msr = math.exp(mu_msr)
    if msr is None:
        # Estimate from HR under exponential assumption: median_ratio approx HR^(-1)
        msr = math.exp(-mu_re)

    # RMST pooling
    rmst_result = None
    if rmst_differences is not None and len(rmst_differences) > 0:
        rmst_arr = np.array(rmst_differences, dtype=float)
        # Approximate SEs as 1/sqrt(k) of each review
        rmst_ses = np.full(len(rmst_arr), 0.5)  # placeholder
        tau2_rmst = _dl_tau2(rmst_arr, rmst_ses)
        mu_rmst, se_rmst = _iv_pool(rmst_arr, rmst_ses, tau2_rmst)
        rmst_result = {"pooled_rmst_diff": float(mu_rmst), "se": float(se_rmst)}

    # Heterogeneity decomposition
    het_decomp = {"clinical": 0.0, "methodological": 0.0}
    if tau2 > 0:
        if median_ages is not None and len(median_ages) >= 2:
            ages = np.array(median_ages, dtype=float)
            valid = np.isfinite(ages)
            if np.sum(valid) >= 2:
                clinical_var = float(np.var(ages[valid]))
                # Normalize: clinical component is proportional to age variance
                # relative to the total heterogeneity
                max_age_var = max(clinical_var, 1.0)
                clinical_fraction = min(1.0, clinical_var / (max_age_var + tau2 * 100))
                het_decomp["clinical"] = round(tau2 * clinical_fraction, 6)
                het_decomp["methodological"] = round(
                    tau2 * (1.0 - clinical_fraction), 6
                )
            else:
                het_decomp["methodological"] = round(tau2, 6)
        else:
            het_decomp["methodological"] = round(tau2, 6)

    result = {
        "pooled_log_hr": round(pooled_log_hr, 6),
        "pooled_hr": round(pooled_hr, 6),
        "hr_ci": (round(hr_ci_lo, 6), round(hr_ci_hi, 6)),
        "tau2_log_hr": round(tau2, 6),
        "frailty_adjusted_hr": round(frailty_adjusted_hr, 6),
        "median_survival_ratio": round(msr, 6) if msr is not None else None,
        "i2": round(i2, 1),
        "heterogeneity_decomposition": het_decomp,
    }
    if rmst_result is not None:
        result["rmst_pooled"] = rmst_result

    return result
