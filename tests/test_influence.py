import pytest
from umbrella.influence import compute_influence


def test_loo_length(statin_reviews):
    """LOO estimates list has one entry per review."""
    result = compute_influence(statin_reviews)
    assert len(result["loo_estimates"]) == len(statin_reviews)


def test_cooks_d_nonnegative(statin_reviews):
    """Cook's distance is always >= 0."""
    result = compute_influence(statin_reviews)
    for rid, d in result["cooks_d"].items():
        assert d >= 0.0, f"Negative Cook's D for {rid}: {d}"


def test_dfbetas_finite(sglt2i_reviews):
    """DFBETAS values are finite floats."""
    result = compute_influence(sglt2i_reviews)
    for rid, dfb in result["dfbetas"].items():
        assert abs(dfb) < 1e6, f"DFBETAS not finite for {rid}: {dfb}"


def test_most_influential_in_reviews(statin_reviews):
    """Most influential review ID is one of the input review IDs."""
    result = compute_influence(statin_reviews)
    ids = {r.review_id for r in statin_reviews}
    assert result["most_influential"] in ids


def test_tipping_point_positive(ivermectin_reviews):
    """Tipping point is at least 1 (must remove at least one review)."""
    result = compute_influence(ivermectin_reviews)
    assert result["tipping_point"] >= 1
