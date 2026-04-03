import pytest
from umbrella.dirichlet_process import dirichlet_process_cluster


def test_n_clusters_at_least_one(statin_reviews):
    """DP clustering must find at least 1 cluster."""
    result = dirichlet_process_cluster(statin_reviews, seed=42)
    assert result["n_clusters"] >= 1


def test_all_reviews_assigned(statin_reviews):
    """Every review must appear in cluster_assignments."""
    result = dirichlet_process_cluster(statin_reviews, seed=42)
    ids = {r.review_id for r in statin_reviews}
    assert set(result["cluster_assignments"].keys()) == ids


def test_cluster_sizes_sum_to_n(ivermectin_reviews):
    """Sum of cluster sizes must equal number of reviews."""
    result = dirichlet_process_cluster(ivermectin_reviews, seed=123)
    total = sum(result["cluster_sizes"].values())
    assert total == len(ivermectin_reviews)


def test_assignment_certainty_valid_range(sglt2i_reviews):
    """Assignment certainty for each review must be in (0, 1]."""
    result = dirichlet_process_cluster(sglt2i_reviews, seed=7)
    for rid, cert in result["assignment_certainty"].items():
        assert 0.0 < cert <= 1.0, f"{rid} certainty {cert} out of range"


def test_deterministic_with_seed(statin_reviews):
    """Same seed must produce identical results."""
    r1 = dirichlet_process_cluster(statin_reviews, seed=99)
    r2 = dirichlet_process_cluster(statin_reviews, seed=99)
    assert r1["n_clusters"] == r2["n_clusters"]
    assert r1["cluster_assignments"] == r2["cluster_assignments"]
