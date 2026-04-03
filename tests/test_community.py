import pytest
from umbrella.community import community_detection


def test_all_reviews_assigned(statin_reviews):
    """Every review must appear in exactly one community."""
    result = community_detection(statin_reviews)
    assigned_ids = []
    for comm in result["communities"]:
        assigned_ids.extend(comm)
    review_ids = {r.review_id for r in statin_reviews}
    assert set(assigned_ids) == review_ids
    assert len(assigned_ids) == len(review_ids)  # no duplicates


def test_modularity_bounded(statin_reviews):
    """Modularity Q must be in [-0.5, 1]."""
    result = community_detection(statin_reviews)
    assert -0.5 <= result["modularity"] <= 1.0


def test_n_communities_matches_list(statin_reviews):
    """n_communities must equal length of communities list."""
    result = community_detection(statin_reviews)
    assert result["n_communities"] == len(result["communities"])
    assert result["n_communities"] >= 1


def test_resolution_scan_has_entries(statin_reviews):
    """Resolution scan should return one entry per gamma value."""
    result = community_detection(statin_reviews)
    assert len(result["resolution_scan"]) == 4  # default [0.5, 1.0, 1.5, 2.0]
    for entry in result["resolution_scan"]:
        assert "gamma" in entry
        assert "n_communities" in entry
        assert "modularity" in entry


def test_no_overlap_all_singletons(two_no_overlap):
    """With zero study overlap, each review should be its own community."""
    result = community_detection(two_no_overlap)
    # No edges => every node is its own community
    assert result["n_communities"] == 2
    assert result["modularity"] == 0.0
