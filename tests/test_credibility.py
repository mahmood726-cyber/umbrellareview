"""Tests for credibility theory review weighting."""

import pytest
from umbrella.models import ReviewInput
from umbrella.credibility import compute_credibility


def _make_review(rid, theta, se, k):
    return ReviewInput(
        review_id=rid, theta=theta,
        ci_lo=theta - 1.96 * se, ci_hi=theta + 1.96 * se,
        se=se, k=k, study_ids=[f"s{i}" for i in range(k)], measure="logOR",
    )


def test_z_in_zero_one(statin_reviews):
    """All credibility factors Z must be in [0, 1]."""
    result = compute_credibility(statin_reviews)
    for rid, z in result["credibility_factors"].items():
        assert 0.0 <= z <= 1.0, f"Z for {rid} = {z} out of [0,1]"


def test_more_studies_higher_credibility():
    """A review with more primary studies should get higher Z."""
    r_small = _make_review("small", -0.3, 0.15, k=3)
    r_large = _make_review("large", -0.4, 0.10, k=30)
    result = compute_credibility([r_small, r_large])
    assert result["credibility_factors"]["large"] > result["credibility_factors"]["small"]


def test_grand_mean_within_range(sglt2i_reviews):
    """Grand mean should be within the range of observed thetas."""
    thetas = [r.theta for r in sglt2i_reviews]
    result = compute_credibility(sglt2i_reviews)
    assert min(thetas) <= result["grand_mean"] <= max(thetas)


def test_k_cred_positive(statin_reviews):
    """Credibility constant k_cred must be positive."""
    result = compute_credibility(statin_reviews)
    assert result["k_cred"] > 0


def test_full_credibility_threshold_positive(statin_reviews):
    """Full credibility threshold must be positive and finite."""
    result = compute_credibility(statin_reviews)
    assert result["full_credibility_threshold"] > 0
    assert result["full_credibility_threshold"] < float("inf")
