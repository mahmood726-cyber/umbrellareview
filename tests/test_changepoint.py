import pytest
from umbrella.changepoint import changepoint_detection


def test_cusum_p_value_in_range(statin_reviews):
    """CUSUM p-value must be in [0, 1]."""
    result = changepoint_detection(statin_reviews, n_perm=199, seed=42)
    assert 0.0 <= result["cusum_p_value"] <= 1.0


def test_bayesian_posterior_sums_to_one(ivermectin_reviews):
    """Bayesian posterior over change-point years must sum to ~1."""
    result = changepoint_detection(ivermectin_reviews, n_perm=99, seed=42)
    total = sum(result["bayesian_posterior"].values())
    assert abs(total - 1.0) < 0.01, f"Posterior sums to {total}"


def test_cusum_statistic_length(statin_reviews):
    """CUSUM statistic list must have same length as number of reviews."""
    result = changepoint_detection(statin_reviews, n_perm=99, seed=42)
    assert len(result["cusum_statistic"]) == len(statin_reviews)


def test_delta_ci_contains_delta(statin_reviews):
    """Delta CI should contain delta_theta (basic sanity)."""
    result = changepoint_detection(statin_reviews, n_perm=99, seed=42)
    lo, hi = result["delta_ci"]
    assert lo <= result["delta_theta"] <= hi


def test_drift_slope_is_finite(sglt2i_reviews):
    """Drift slope must be a finite number."""
    result = changepoint_detection(sglt2i_reviews, n_perm=99, seed=42)
    import math
    assert math.isfinite(result["drift_slope"])
    assert math.isfinite(result["drift_p"])
    assert 0.0 <= result["drift_p"] <= 1.0
