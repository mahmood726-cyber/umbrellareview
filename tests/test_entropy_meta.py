"""Tests for information-theoretic meta-meta-analysis."""

import math
import pytest
from umbrella.models import ReviewInput
from umbrella.entropy_meta import compute_entropy_meta


def _make_review(rid, theta, ci_lo, ci_hi, amstar=None):
    return ReviewInput(
        review_id=rid, theta=theta, ci_lo=ci_lo, ci_hi=ci_hi,
        k=5, study_ids=[f"s{rid}_{i}" for i in range(5)], measure="logOR",
        amstar_items=amstar or {},
    )


def test_entropy_non_negative(statin_reviews):
    """Shannon entropy must be >= 0."""
    result = compute_entropy_meta(statin_reviews)
    assert result["conclusion_entropy"] >= 0.0


def test_mutual_info_non_negative(sglt2i_reviews):
    """Mutual information I(Q;C) must be >= 0."""
    result = compute_entropy_meta(sglt2i_reviews)
    assert result["mutual_info_quality_conclusion"] >= 0.0


def test_normalized_entropy_bounded():
    """Normalized entropy must be in [0, 1]."""
    reviews = [
        _make_review("R1", 0.5, 0.2, 0.8),    # beneficial
        _make_review("R2", -0.5, -0.8, -0.2),  # harmful
        _make_review("R3", 0.1, -0.2, 0.4),    # inconclusive
    ]
    result = compute_entropy_meta(reviews)
    assert 0.0 <= result["normalized_entropy"] <= 1.0 + 1e-10


def test_redundancy_plus_norm_entropy_eq_one():
    """Redundancy + normalized_entropy should equal 1."""
    reviews = [
        _make_review("R1", 0.5, 0.2, 0.8),
        _make_review("R2", 0.3, 0.1, 0.5),
        _make_review("R3", -0.1, -0.4, 0.2),
    ]
    result = compute_entropy_meta(reviews)
    total = result["redundancy"] + result["normalized_entropy"]
    assert abs(total - 1.0) < 1e-6, f"redundancy + norm_entropy = {total}"


def test_simpson_diversity_bounded(statin_reviews):
    """Simpson's diversity must be in [0, 1]."""
    result = compute_entropy_meta(statin_reviews)
    assert 0.0 <= result["simpson_diversity"] <= 1.0
