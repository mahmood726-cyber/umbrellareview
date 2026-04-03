"""Tests for profile likelihood heterogeneity module."""

import pytest
from umbrella.models import ReviewInput
from umbrella.profile_likelihood import profile_likelihood


# ── Fixtures ─────────────────────────────────────────────────────

def _make_reviews(specs):
    """specs: list of (id, theta, se)."""
    return [
        ReviewInput(
            review_id=rid, theta=theta, ci_lo=theta - 1.96 * se,
            ci_hi=theta + 1.96 * se, se=se, k=10,
            study_ids=[f"s{i}" for i in range(10)], measure="logOR",
        )
        for rid, theta, se in specs
    ]


@pytest.fixture
def homogeneous():
    """Five reviews with very similar thetas — low tau^2."""
    return _make_reviews([
        ("R1", -0.50, 0.10),
        ("R2", -0.48, 0.12),
        ("R3", -0.51, 0.11),
        ("R4", -0.49, 0.10),
        ("R5", -0.50, 0.13),
    ])


@pytest.fixture
def heterogeneous():
    """Five reviews with spread thetas — substantial tau^2."""
    return _make_reviews([
        ("R1", -0.80, 0.10),
        ("R2", -0.20, 0.12),
        ("R3",  0.10, 0.15),
        ("R4", -0.50, 0.11),
        ("R5", -1.00, 0.09),
    ])


@pytest.fixture
def two_reviews():
    """Minimal case: two reviews."""
    return _make_reviews([
        ("A", -0.3, 0.15),
        ("B", -0.6, 0.20),
    ])


# ── Tests ────────────────────────────────────────────────────────

def test_tau2_reml_nonnegative(homogeneous):
    """REML tau^2 must be >= 0."""
    result = profile_likelihood(homogeneous)
    assert result["tau2_reml"] >= 0.0


def test_profile_ci_contains_reml(heterogeneous):
    """REML CI must contain the REML point estimate."""
    result = profile_likelihood(heterogeneous)
    lo, hi = result["tau2_reml_ci"]
    tau2 = result["tau2_reml"]
    assert lo <= tau2 + 1e-10, f"CI lower {lo} > tau2 {tau2}"
    assert tau2 <= hi + 1e-10, f"tau2 {tau2} > CI upper {hi}"


def test_grid_length(homogeneous):
    """Grid should have n_steps + 1 points (default 101)."""
    result = profile_likelihood(homogeneous)
    assert len(result["tau2_grid"]) == 101
    assert len(result["profile_ll"]) == 101
    assert len(result["reml_ll"]) == 101


def test_bartlett_p_in_range(heterogeneous):
    """Bartlett-corrected p-value must be in [0, 1]."""
    result = profile_likelihood(heterogeneous)
    assert 0.0 <= result["bartlett_p"] <= 1.0
    assert result["bartlett_lr"] >= 0.0


def test_saddlepoint_se_positive(heterogeneous):
    """Saddlepoint SE must be positive for heterogeneous data."""
    result = profile_likelihood(heterogeneous)
    assert result["saddlepoint_se"] > 0.0
