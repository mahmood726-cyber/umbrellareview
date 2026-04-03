import pytest
from umbrella.bayesian_meta import bayesian_meta_meta


def test_tau_mm2_non_negative(statin_reviews):
    result = bayesian_meta_meta(statin_reviews)
    assert result["tau_mm2"] >= 0.0


def test_prediction_interval_contains_mu(sglt2i_reviews):
    result = bayesian_meta_meta(sglt2i_reviews)
    pi_lo, pi_hi = result["prediction_interval"]
    assert pi_lo <= result["mu"] <= pi_hi


def test_shrinkage_all_reviews_present(statin_reviews):
    result = bayesian_meta_meta(statin_reviews)
    ids = {r.review_id for r in statin_reviews}
    assert set(result["shrinkage_estimates"].keys()) == ids


def test_i2_mm_in_range(ivermectin_reviews):
    result = bayesian_meta_meta(ivermectin_reviews)
    assert 0.0 <= result["i2_mm"] <= 100.0


def test_bic_values_finite(statin_reviews):
    result = bayesian_meta_meta(statin_reviews)
    assert result["bic_fixed"] != float("inf")
    assert result["bic_random"] != float("inf")
    # Both should be real numbers
    assert result["bic_fixed"] == result["bic_fixed"]  # not NaN
    assert result["bic_random"] == result["bic_random"]  # not NaN
