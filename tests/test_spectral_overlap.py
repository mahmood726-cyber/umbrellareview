import pytest
from umbrella.spectral_overlap import spectral_overlap


def test_spectral_concentration_in_range(statin_reviews):
    result = spectral_overlap(statin_reviews)
    assert 0.0 <= result["spectral_concentration"] <= 1.0


def test_fiedler_value_non_negative(sglt2i_reviews):
    result = spectral_overlap(sglt2i_reviews)
    assert result["fiedler_value"] >= 0.0


def test_co_citation_matrix_symmetric(statin_reviews):
    result = spectral_overlap(statin_reviews)
    n = len(statin_reviews)
    mat = result["co_citation_matrix"]
    assert len(mat) == n
    for i in range(n):
        for j in range(n):
            assert mat[i][j] == mat[j][i]


def test_review_clusters_correct_length(ivermectin_reviews):
    result = spectral_overlap(ivermectin_reviews)
    assert len(result["review_clusters"]) == len(ivermectin_reviews)
    # Clusters should be 0 or 1
    for c in result["review_clusters"]:
        assert c in (0, 1)


def test_svd_explained_sums_to_one(statin_reviews):
    result = spectral_overlap(statin_reviews)
    total = sum(result["svd_explained"])
    assert abs(total - 1.0) < 0.01
