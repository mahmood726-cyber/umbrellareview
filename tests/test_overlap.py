import pytest
from umbrella.overlap import compute_overlap

def test_cca_zero_no_overlap(two_no_overlap):
    result = compute_overlap(two_no_overlap)
    assert result.cca == 0.0
    assert result.groove == "Slight"

def test_cca_one_full_overlap(two_full_overlap):
    result = compute_overlap(two_full_overlap)
    assert abs(result.cca - 1.0) < 0.01
    assert result.groove == "Very High"

def test_overlap_matrix_symmetric(statin_reviews):
    result = compute_overlap(statin_reviews)
    n = len(statin_reviews)
    for i in range(n):
        for j in range(n):
            assert result.overlap_matrix[i][j] == result.overlap_matrix[j][i]

def test_unique_studies_count(two_no_overlap):
    result = compute_overlap(two_no_overlap)
    assert result.n_unique_studies == 6
    assert result.n_total_citations == 6

def test_study_frequency(two_full_overlap):
    result = compute_overlap(two_full_overlap)
    for sid in ["s1", "s2", "s3"]:
        assert result.study_frequency[sid] == 2

def test_groove_boundaries():
    from umbrella.overlap import _classify_groove
    assert _classify_groove(0.03) == "Slight"
    assert _classify_groove(0.07) == "Moderate"
    assert _classify_groove(0.12) == "High"
    assert _classify_groove(0.20) == "Very High"
