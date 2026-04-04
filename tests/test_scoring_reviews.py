import pytest
from umbrella.scoring_reviews import compute_scoring


def test_brier_in_unit_range(statin_reviews):
    """Brier scores must be in [0, 1]."""
    result = compute_scoring(statin_reviews)
    for rid, score in result["brier_scores"].items():
        assert 0.0 <= score <= 1.0, f"Brier out of range for {rid}: {score}"


def test_log_scores_finite(sglt2i_reviews):
    """Log scores must be finite (not NaN or inf)."""
    result = compute_scoring(sglt2i_reviews)
    for rid, score in result["log_scores"].items():
        assert score == score, f"Log score is NaN for {rid}"  # NaN != NaN
        assert abs(score) < 1e10, f"Log score too large for {rid}: {score}"


def test_ranking_length(statin_reviews):
    """Ranking list has one entry per review."""
    result = compute_scoring(statin_reviews)
    assert len(result["ranking"]) == len(statin_reviews)


def test_best_review_in_ids(ivermectin_reviews):
    """Best review ID is one of the input review IDs."""
    result = compute_scoring(ivermectin_reviews)
    ids = {r.review_id for r in ivermectin_reviews}
    assert result["best_review"] in ids


def test_skill_scores_bounded(statin_reviews):
    """Skill scores should be <= 1 (perfect skill relative to reference)."""
    result = compute_scoring(statin_reviews)
    for rid, score in result["skill_scores"].items():
        assert score <= 1.0, f"Skill score > 1 for {rid}: {score}"
