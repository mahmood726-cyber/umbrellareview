"""Tests for hierarchical clustering with cophenetic analysis."""

import pytest
from umbrella.models import ReviewInput
from umbrella.cophenetic import cophenetic_analysis


# ── Fixtures ─────────────────────────────────────────────────────

def _make_reviews(specs):
    """specs: list of (id, theta, study_ids)."""
    return [
        ReviewInput(
            review_id=rid, theta=theta, ci_lo=theta - 0.3, ci_hi=theta + 0.3,
            k=len(sids), study_ids=sids, measure="logOR",
        )
        for rid, theta, sids in specs
    ]


@pytest.fixture
def three_reviews():
    """Three reviews with varying overlap."""
    return _make_reviews([
        ("A", -0.5, ["s1", "s2", "s3", "s4"]),
        ("B", -0.4, ["s2", "s3", "s4", "s5"]),  # overlaps with A
        ("C", -0.1, ["s6", "s7", "s8"]),          # no overlap with A/B
    ])


@pytest.fixture
def identical_reviews():
    """Two reviews with identical study sets — distance=0."""
    return _make_reviews([
        ("X", -0.3, ["s1", "s2", "s3"]),
        ("Y", -0.4, ["s1", "s2", "s3"]),
    ])


@pytest.fixture
def four_clustered():
    """Four reviews forming two natural clusters."""
    return _make_reviews([
        ("A", -0.5, ["s1", "s2", "s3"]),
        ("B", -0.4, ["s1", "s2", "s4"]),   # cluster with A
        ("C",  0.3, ["s10", "s11", "s12"]),
        ("D",  0.2, ["s10", "s11", "s13"]),  # cluster with C
    ])


# ── Tests ────────────────────────────────────────────────────────

def test_cophenetic_correlation_in_range(three_reviews):
    """Cophenetic correlation must be in [-1, 1]."""
    result = cophenetic_analysis(three_reviews)
    cc = result["cophenetic_correlation"]
    assert -1.0 <= cc <= 1.0, f"Cophenetic correlation {cc} out of range"


def test_all_reviews_assigned(three_reviews):
    """Every review gets a cluster assignment."""
    result = cophenetic_analysis(three_reviews)
    assert set(result["cluster_assignments"].keys()) == {"A", "B", "C"}


def test_dendrogram_length(three_reviews):
    """n-1 merges for n reviews."""
    result = cophenetic_analysis(three_reviews)
    assert len(result["dendrogram"]) == 2  # 3 reviews -> 2 merges


def test_ultrametric_score_zero_for_upgma(four_clustered):
    """UPGMA produces an ultrametric — violation score should be 0."""
    result = cophenetic_analysis(four_clustered)
    assert result["ultrametric_score"] < 1e-10, (
        f"Expected 0 violations, got {result['ultrametric_score']}"
    )


def test_cluster_summaries_cover_all(four_clustered):
    """Sum of n_reviews across clusters equals total reviews."""
    result = cophenetic_analysis(four_clustered)
    total = sum(s["n_reviews"] for s in result["cluster_summaries"])
    assert total == 4
