"""Profile Likelihood for Between-Review Heterogeneity.

Higher-order inference for tau^2 in umbrella reviews:
profile log-likelihood, REML, Bartlett correction, Q-profile method,
and saddlepoint-approximation SE.
"""

from __future__ import annotations
import math


# ── DerSimonian-Laird tau^2 ─────────────────────────────────────────

def _tau2_dl(thetas: list[float], ses: list[float]) -> float:
    """DerSimonian-Laird estimator of tau^2."""
    k = len(thetas)
    if k < 2:
        return 0.0
    w = [1.0 / (s * s) for s in ses]
    sw = sum(w)
    mu = sum(wi * ti for wi, ti in zip(w, thetas)) / sw
    q = sum(wi * (ti - mu) ** 2 for wi, ti in zip(w, thetas))
    c = sw - sum(wi * wi for wi in w) / sw
    if c <= 0:
        return 0.0
    tau2 = max(0.0, (q - (k - 1)) / c)
    return tau2


# ── Weighted mean given tau^2 ───────────────────────────────────────

def _mu_hat(thetas: list[float], ses: list[float], tau2: float) -> float:
    w = [1.0 / (s * s + tau2) for s in ses]
    sw = sum(w)
    if sw < 1e-30:
        return 0.0
    return sum(wi * ti for wi, ti in zip(w, thetas)) / sw


# ── Profile log-likelihood ──────────────────────────────────────────

def _profile_ll(thetas: list[float], ses: list[float], tau2: float) -> float:
    """Profile log-likelihood l_P(tau2)."""
    mu = _mu_hat(thetas, ses, tau2)
    ll = 0.0
    for ti, si in zip(thetas, ses):
        v = si * si + tau2
        ll -= 0.5 * math.log(v)
        ll -= 0.5 * (ti - mu) ** 2 / v
    return ll


# ── REML log-likelihood ────────────────────────────────────────────

def _reml_ll(thetas: list[float], ses: list[float], tau2: float) -> float:
    """Restricted maximum likelihood log-likelihood l_R(tau2)."""
    lp = _profile_ll(thetas, ses, tau2)
    sw = sum(1.0 / (si * si + tau2) for si in ses)
    if sw < 1e-30:
        return lp
    lr = lp - 0.5 * math.log(sw)
    return lr


# ── Grid construction ───────────────────────────────────────────────

def _make_grid(tau2_dl: float, n_steps: int = 100) -> list[float]:
    """Grid from 0 to 5*tau2_DL (minimum upper bound 1.0)."""
    upper = max(5.0 * tau2_dl, 1.0)
    return [upper * i / n_steps for i in range(n_steps + 1)]


# ── Q-profile method (Viechtbauer 2007) ────────────────────────────

def _q_stat(thetas: list[float], ses: list[float], tau2: float) -> float:
    """Generalised Cochran Q at given tau2."""
    mu = _mu_hat(thetas, ses, tau2)
    return sum((ti - mu) ** 2 / (si * si + tau2) for ti, si in zip(thetas, ses))


def _chi2_quantile(p: float, df: int) -> float:
    """Approximate chi-squared quantile using Wilson-Hilferty."""
    if df <= 0:
        return 0.0
    z = _normal_quantile(p)
    # Wilson-Hilferty
    val = df * (1.0 - 2.0 / (9.0 * df) + z * math.sqrt(2.0 / (9.0 * df))) ** 3
    return max(0.0, val)


def _normal_quantile(p: float) -> float:
    """Approximate inverse normal CDF (Abramowitz & Stegun 26.2.23)."""
    if p <= 0:
        return -8.0
    if p >= 1:
        return 8.0
    if p < 0.5:
        return -_normal_quantile(1.0 - p)
    t = math.sqrt(-2.0 * math.log(1.0 - p))
    c0, c1, c2 = 2.515517, 0.802853, 0.010328
    d1, d2, d3 = 1.432788, 0.189269, 0.001308
    return t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t)


def _chi2_cdf(x: float, df: int) -> float:
    """Approximate chi-squared CDF using regularised incomplete gamma."""
    if x <= 0 or df <= 0:
        return 0.0
    return _regularised_gamma_p(df / 2.0, x / 2.0)


def _regularised_gamma_p(a: float, x: float) -> float:
    """Lower regularised incomplete gamma P(a, x) via series expansion."""
    if x < 0:
        return 0.0
    if x == 0:
        return 0.0
    # Series: P(a,x) = e^{-x} x^a sum_{n=0}^{inf} x^n / Gamma(a+n+1)
    # = e^{-x} x^a / Gamma(a) * sum_{n=0} x^n / (a*(a+1)*...*(a+n))
    ln_prefix = a * math.log(x) - x - math.lgamma(a)
    if ln_prefix < -700:
        return 0.0 if x < a else 1.0

    total = 0.0
    term = 1.0 / a
    total = term
    for n in range(1, 300):
        term *= x / (a + n)
        total += term
        if abs(term) < 1e-15 * abs(total):
            break

    result = math.exp(ln_prefix) * total
    return max(0.0, min(1.0, result))


def _q_profile_ci(thetas: list[float], ses: list[float],
                  grid: list[float], alpha: float = 0.05) -> tuple[float, float]:
    """CI for tau^2 by inverting Q-profile test."""
    k = len(thetas)
    df = k - 1
    if df < 1:
        return (0.0, 0.0)

    crit_lo = _chi2_quantile(1.0 - alpha / 2.0, df)
    crit_hi = _chi2_quantile(alpha / 2.0, df)

    # Find tau2 values where Q crosses critical values
    lo = 0.0
    hi = grid[-1]

    # Scan grid for lower bound: Q(tau2) = crit_lo
    # Q decreases as tau2 increases, so lower CI bound is where Q = crit_lo
    for tau2 in grid:
        q = _q_stat(thetas, ses, tau2)
        if q <= crit_lo:
            lo = tau2
            break

    # Upper bound: Q(tau2) = crit_hi
    for tau2 in grid:
        q = _q_stat(thetas, ses, tau2)
        if q <= crit_hi:
            hi = tau2
            break
    else:
        hi = grid[-1]

    return (lo, hi)


# ── Bartlett correction ─────────────────────────────────────────────

def _bartlett_correction(thetas: list[float], ses: list[float],
                         tau2_hat: float) -> tuple[float, float]:
    """Bartlett-corrected LR test for H0: tau^2 = 0.

    Returns (corrected LR statistic, p-value).
    """
    k = len(thetas)
    if k < 2:
        return (0.0, 1.0)

    lr = 2.0 * (_reml_ll(thetas, ses, tau2_hat) - _reml_ll(thetas, ses, 0.0))
    lr = max(0.0, lr)

    # Approximate Bartlett factor
    # a = sum of (v_i - v_bar)^2 / v_bar^2 terms — simplified approach
    # Use third-derivative approximation
    w0 = [1.0 / (si * si) for si in ses]
    sw0 = sum(w0)
    # Under H0: tau2 = 0
    a = sum(wi * wi for wi in w0) / (sw0 * sw0)
    # Bartlett factor: b = 1 + a / (6 * df)
    df = 1
    b = 1.0 + a / (6.0 * df)
    if b < 0.5:
        b = 0.5  # safety clamp

    lr_corrected = lr / b

    # P-value from chi2(1)
    p = 1.0 - _chi2_cdf(lr_corrected, 1)

    return lr_corrected, p


# ── Saddlepoint SE ──────────────────────────────────────────────────

def _saddlepoint_se(thetas: list[float], ses: list[float],
                    tau2_hat: float) -> float:
    """Approximate SE of tau^2_hat via observed Fisher information.

    Uses the second derivative of the REML log-likelihood.
    """
    h = max(tau2_hat * 0.01, 1e-6)
    ll_plus = _reml_ll(thetas, ses, tau2_hat + h)
    ll_mid = _reml_ll(thetas, ses, tau2_hat)
    ll_minus = _reml_ll(thetas, ses, max(0.0, tau2_hat - h))

    d2 = (ll_plus - 2.0 * ll_mid + ll_minus) / (h * h)

    if d2 >= -1e-15:
        # Flat or convex — fall back to DL-based variance approximation
        k = len(thetas)
        w = [1.0 / (si * si + tau2_hat) for si in ses]
        info = 0.5 * sum(wi * wi for wi in w)
        return 1.0 / math.sqrt(info) if info > 0 else 0.0

    info = -d2
    return 1.0 / math.sqrt(info)


# ── Public API ──────────────────────────────────────────────────────

def profile_likelihood(reviews, n_steps: int = 100, alpha: float = 0.05) -> dict:
    """Profile likelihood analysis for between-review heterogeneity.

    Parameters
    ----------
    reviews : list[ReviewInput]
    n_steps : grid resolution (default 100)
    alpha   : significance level (default 0.05)

    Returns
    -------
    dict with:
        profile_ll    — list[float]
        reml_ll       — list[float]
        tau2_grid     — list[float]
        tau2_reml     — float
        tau2_reml_ci  — (float, float)
        bartlett_lr   — float
        bartlett_p    — float
        q_profile_ci  — (float, float)
        saddlepoint_se — float
    """
    thetas = [r.theta for r in reviews]
    ses = [r.se for r in reviews]

    # Ensure positive SE
    for i, s in enumerate(ses):
        if s <= 0:
            ses[i] = (reviews[i].ci_hi - reviews[i].ci_lo) / (2 * 1.96)
            if ses[i] <= 0:
                ses[i] = 0.1  # safety fallback

    k = len(thetas)
    if k < 2:
        return {
            "profile_ll": [],
            "reml_ll": [],
            "tau2_grid": [],
            "tau2_reml": 0.0,
            "tau2_reml_ci": (0.0, 0.0),
            "bartlett_lr": 0.0,
            "bartlett_p": 1.0,
            "q_profile_ci": (0.0, 0.0),
            "saddlepoint_se": 0.0,
        }

    # DL starting estimate
    tau2_dl = _tau2_dl(thetas, ses)
    grid = _make_grid(tau2_dl, n_steps)

    # Compute profile and REML likelihoods on grid
    pll = [_profile_ll(thetas, ses, t) for t in grid]
    rll = [_reml_ll(thetas, ses, t) for t in grid]

    # REML estimate: maximise REML LL
    best_idx = max(range(len(rll)), key=lambda i: rll[i])
    tau2_reml = grid[best_idx]

    # REML CI: invert profile LR
    crit = _chi2_quantile(1.0 - alpha, 1)  # 3.84 for alpha=0.05
    max_rll = rll[best_idx]
    ci_indices = [i for i in range(len(grid))
                  if 2.0 * (max_rll - rll[i]) <= crit]
    if ci_indices:
        tau2_reml_ci = (grid[ci_indices[0]], grid[ci_indices[-1]])
    else:
        tau2_reml_ci = (tau2_reml, tau2_reml)

    # Bartlett correction
    bartlett_lr, bartlett_p = _bartlett_correction(thetas, ses, tau2_reml)

    # Q-profile CI
    qp_ci = _q_profile_ci(thetas, ses, grid, alpha)

    # Saddlepoint SE
    sp_se = _saddlepoint_se(thetas, ses, tau2_reml)

    return {
        "profile_ll": pll,
        "reml_ll": rll,
        "tau2_grid": [float(t) for t in grid],
        "tau2_reml": tau2_reml,
        "tau2_reml_ci": tau2_reml_ci,
        "bartlett_lr": bartlett_lr,
        "bartlett_p": bartlett_p,
        "q_profile_ci": qp_ci,
        "saddlepoint_se": sp_se,
    }
