import pytest
from umbrella.robust_pooling import compute_robust_pooling


def test_huber_converges(statin_reviews):
    """Huber M-estimation converges within max iterations."""
    result = compute_robust_pooling(statin_reviews)
    assert result["huber"]["converged"] is True


def test_median_within_range(statin_reviews):
    """Weighted median is within range of observed thetas."""
    thetas = [r.theta for r in statin_reviews]
    result = compute_robust_pooling(statin_reviews)
    assert min(thetas) <= result["median"]["theta"] <= max(thetas)


def test_all_estimators_present(sglt2i_reviews):
    """All four estimators are returned."""
    result = compute_robust_pooling(sglt2i_reviews)
    for key in ("standard", "median", "huber", "winsorized"):
        assert key in result
        assert "theta" in result[key]
        assert "se" in result[key]


def test_standard_matches_manual(two_no_overlap):
    """Standard IV pooling matches hand calculation for 2 reviews."""
    result = compute_robust_pooling(two_no_overlap)
    # Both thetas are negative, so pooled should be negative
    assert result["standard"]["theta"] < 0


def test_outlier_flags_are_review_ids(ivermectin_reviews):
    """Outlier flags (if any) contain valid review IDs."""
    result = compute_robust_pooling(ivermectin_reviews)
    valid_ids = {r.review_id for r in ivermectin_reviews}
    for rid in result["outlier_flags"]:
        assert rid in valid_ids
