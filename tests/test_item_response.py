import math
import pytest
from umbrella.item_response import rasch_irt


def test_theta_estimates_finite(statin_reviews):
    """All theta estimates must be finite numbers."""
    result = rasch_irt(statin_reviews)
    assert len(result["theta_estimates"]) > 0
    for rid, t in result["theta_estimates"].items():
        assert math.isfinite(t), f"theta for {rid} is not finite: {t}"


def test_converges_within_iterations(statin_reviews):
    """JMLE should converge on well-behaved AMSTAR data."""
    result = rasch_irt(statin_reviews)
    assert result["converged"] is True


def test_reliability_bounded(statin_reviews):
    """Reliability must be in [0, 1]."""
    result = rasch_irt(statin_reviews)
    assert 0.0 <= result["reliability"] <= 1.0


def test_outfit_infit_positive(statin_reviews):
    """Outfit and infit statistics must be positive for every item."""
    result = rasch_irt(statin_reviews)
    for item, val in result["outfit"].items():
        assert val > 0, f"outfit for {item} is non-positive: {val}"
    for item, val in result["infit"].items():
        assert val > 0, f"infit for {item} is non-positive: {val}"


def test_insufficient_reviews_returns_empty():
    """With fewer than 2 reviews, should return empty/unconverged result."""
    from umbrella.models import ReviewInput
    single = [ReviewInput(review_id="X", theta=0.1, ci_lo=-0.1, ci_hi=0.3,
                          amstar_items={"1": "yes", "2": "no"})]
    result = rasch_irt(single)
    assert result["converged"] is False
    assert result["theta_estimates"] == {}
    assert result["reliability"] == 0.0
