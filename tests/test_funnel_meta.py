import pytest
from umbrella.funnel_meta import funnel_meta


def test_egger_intercept_is_finite(statin_reviews):
    """Egger's regression intercept must be a finite number."""
    result = funnel_meta(statin_reviews)
    assert isinstance(result["egger_intercept"], float)
    assert result["egger_intercept"] == result["egger_intercept"]  # not NaN
    assert abs(result["egger_intercept"]) < 1e6  # reasonable magnitude


def test_begg_tau_in_range(ivermectin_reviews):
    """Begg's Kendall tau must be in [-1, 1]."""
    result = funnel_meta(ivermectin_reviews)
    assert -1.0 <= result["begg_tau"] <= 1.0


def test_trimfill_n_missing_non_negative(statin_reviews):
    """Number of missing reviews from trim-and-fill must be >= 0."""
    result = funnel_meta(statin_reviews)
    assert result["n_missing_trimfill"] >= 0
    assert isinstance(result["n_missing_trimfill"], int)


def test_funnel_data_correct_length(sglt2i_reviews):
    """Funnel x and y vectors must match number of reviews."""
    result = funnel_meta(sglt2i_reviews)
    n = len(sglt2i_reviews)
    assert len(result["funnel_x"]) == n
    assert len(result["funnel_y"]) == n
    # All precisions (y) must be positive
    for p in result["funnel_y"]:
        assert p > 0


def test_bias_detected_is_bool(ivermectin_reviews):
    """bias_detected must be boolean and consistent with p-values."""
    result = funnel_meta(ivermectin_reviews)
    assert isinstance(result["bias_detected"], bool)
    # If any p < 0.1, bias_detected should be True
    any_sig = (result["egger_p"] < 0.1 or
               result["begg_p"] < 0.1 or
               result["peters_p"] < 0.1)
    assert result["bias_detected"] == any_sig
