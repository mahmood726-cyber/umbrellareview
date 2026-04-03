"""Tests for Dempster-Shafer evidence theory module."""

import pytest
from umbrella.models import ReviewInput
from umbrella.dempster_shafer import dempster_shafer, _mass_for_review, OMEGA


# ── Fixtures ─────────────────────────────────────────────────────

def _make_review(rid, theta, ci_lo, ci_hi, amstar=None):
    return ReviewInput(
        review_id=rid, theta=theta, ci_lo=ci_lo, ci_hi=ci_hi,
        k=5, study_ids=[f"s{i}" for i in range(5)], measure="logOR",
        amstar_items=amstar or {},
    )


@pytest.fixture
def concordant_beneficial():
    """Three reviews all significantly beneficial."""
    return [
        _make_review("R1", 0.5, 0.2, 0.8),
        _make_review("R2", 0.6, 0.3, 0.9),
        _make_review("R3", 0.4, 0.1, 0.7),
    ]


@pytest.fixture
def conflicting_reviews():
    """One beneficial, one harmful — high conflict expected."""
    return [
        _make_review("R1", 0.5, 0.2, 0.8),   # beneficial
        _make_review("R2", -0.5, -0.8, -0.2),  # harmful
    ]


@pytest.fixture
def inconclusive_reviews():
    """Two reviews with CIs spanning zero."""
    return [
        _make_review("R1", 0.1, -0.2, 0.4),
        _make_review("R2", -0.05, -0.3, 0.2),
    ]


# ── Tests ────────────────────────────────────────────────────────

def test_masses_sum_to_one(concordant_beneficial):
    """Each review's mass function must sum to 1."""
    for r in concordant_beneficial:
        m = _mass_for_review(r)
        total = sum(m.values())
        assert abs(total - 1.0) < 1e-10, f"Mass sums to {total}"


def test_combined_masses_sum_to_one(concordant_beneficial):
    """Combined mass after Dempster's rule must sum to 1."""
    result = dempster_shafer(concordant_beneficial)
    total = sum(result["combined_mass"].values())
    assert abs(total - 1.0) < 1e-10, f"Combined mass sums to {total}"


def test_belief_leq_plausibility(concordant_beneficial):
    """Bel(A) <= Pl(A) for all hypotheses."""
    result = dempster_shafer(concordant_beneficial)
    for h in ["beneficial", "harmful", "inconclusive"]:
        bel = result["belief"][h]
        pl = result["plausibility"][h]
        assert bel <= pl + 1e-10, f"Bel({h})={bel} > Pl({h})={pl}"


def test_conflicting_high_conflict(conflicting_reviews):
    """Opposing beneficial/harmful reviews produce non-zero conflict."""
    result = dempster_shafer(conflicting_reviews)
    assert result["conflict_total"] > 0.0
    # Pairwise conflict should exist
    assert len(result["conflict_pairwise"]) == 1
    k = list(result["conflict_pairwise"].values())[0]
    assert k > 0.1, f"Expected substantial conflict, got K={k}"


def test_pignistic_sums_to_one(concordant_beneficial):
    """Pignistic probabilities must sum to 1."""
    result = dempster_shafer(concordant_beneficial)
    total = sum(result["pignistic"].values())
    assert abs(total - 1.0) < 1e-10, f"Pignistic sums to {total}"
    # Verdict should be beneficial for all-positive reviews
    assert result["verdict"] == "beneficial"
