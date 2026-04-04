"""Goodness-of-Fit Tests for Systematic Review Effect Sizes."""

import numpy as np
from scipy import stats


def compute_gof(reviews):
    """Goodness-of-fit tests for the distribution of review effect sizes.

    Tests whether the observed thetas follow a Normal distribution by
    computing the Kolmogorov-Smirnov, Anderson-Darling, and
    Cramer-von Mises statistics, plus QQ-plot data.

    Parameters
    ----------
    reviews : list[ReviewInput]

    Returns
    -------
    dict with keys: ks_stat, ks_p, ad_stat, cvm_stat,
                    qq_theoretical, qq_empirical, is_normal
    """
    if len(reviews) < 3:
        raise ValueError("Need at least 3 reviews for goodness-of-fit tests")

    thetas = np.array([r.theta for r in reviews], dtype=float)
    n = len(thetas)
    mu = float(np.mean(thetas))
    sigma = float(np.std(thetas, ddof=1))
    if sigma < 1e-15:
        sigma = 1e-15

    # --- Kolmogorov-Smirnov test against fitted Normal ---
    ks_stat, ks_p = stats.kstest(thetas, "norm", args=(mu, sigma))

    # --- Anderson-Darling statistic ---
    x_sorted = np.sort(thetas)
    F_vals = stats.norm.cdf(x_sorted, loc=mu, scale=sigma)
    # Clip to avoid log(0) or log(1)
    F_vals = np.clip(F_vals, 1e-15, 1 - 1e-15)
    i_arr = np.arange(1, n + 1)
    ad_stat = -n - np.sum(
        (2 * i_arr - 1) / n * (np.log(F_vals) + np.log(1 - F_vals[::-1]))
    )
    ad_stat = float(ad_stat)

    # --- Cramer-von Mises statistic ---
    cvm_stat = float(
        np.sum((F_vals - (2 * i_arr - 1) / (2 * n)) ** 2) + 1 / (12 * n)
    )

    # --- QQ-plot data ---
    # Theoretical quantiles from standard normal, then scale
    pp = (i_arr - 0.5) / n
    qq_theoretical = [round(float(x), 6) for x in stats.norm.ppf(pp, loc=mu, scale=sigma)]
    qq_empirical = [round(float(x), 6) for x in x_sorted]

    # --- Normality decision ---
    # Anderson-Darling and CvM don't have simple p-values from our manual
    # formulas, so we rely on KS p-value as primary gate.
    # Use scipy's Anderson-Darling test for a secondary check.
    ad_result = stats.anderson(thetas, dist="norm")
    # 5% significance level is index 2 in ad_result.significance_level
    ad_reject = ad_result.statistic > ad_result.critical_values[2]

    is_normal = bool(ks_p > 0.05 and not ad_reject)

    return {
        "ks_stat": round(float(ks_stat), 6),
        "ks_p": round(float(ks_p), 6),
        "ad_stat": round(ad_stat, 6),
        "cvm_stat": round(cvm_stat, 6),
        "qq_theoretical": qq_theoretical,
        "qq_empirical": qq_empirical,
        "is_normal": is_normal,
    }
