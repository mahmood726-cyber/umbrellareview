import pytest
from umbrella.persistent_homology import persistent_homology


def test_betti0_starts_at_k_ends_at_1(statin_reviews):
    """beta_0 starts at k (all disconnected) and ends at 1 (fully connected)."""
    result = persistent_homology(statin_reviews)
    k = len(statin_reviews)
    assert result["betti_curve_0"][0] == k
    assert result["betti_curve_0"][-1] == 1


def test_persistence_entropy_non_negative(ivermectin_reviews):
    """Persistence entropy must be >= 0 (Shannon entropy is non-negative)."""
    result = persistent_homology(ivermectin_reviews)
    assert result["persistence_entropy"] >= 0.0


def test_distance_matrix_symmetric_unit_diagonal(sglt2i_reviews):
    """Distance matrix: diagonal is 0, off-diagonal in (0,1], symmetric."""
    result = persistent_homology(sglt2i_reviews)
    dm = result["distance_matrix"]
    n = len(sglt2i_reviews)
    assert len(dm) == n
    for i in range(n):
        assert abs(dm[i][i]) < 1e-9  # diagonal = 0
        for j in range(n):
            assert dm[i][j] == pytest.approx(dm[j][i], abs=1e-9)  # symmetric
            assert 0.0 <= dm[i][j] <= 1.0 + 1e-9


def test_filtration_values_correct_length(statin_reviews):
    """Filtration should have n_steps values from 0 to 1."""
    result = persistent_homology(statin_reviews, n_steps=50)
    filt = result["filtration_values"]
    assert len(filt) == 50
    assert abs(filt[0]) < 1e-9
    assert abs(filt[-1] - 1.0) < 1e-9


def test_persistence_diagram_has_features(ivermectin_reviews):
    """Persistence diagram should contain at least k-1 dim-0 features
    (one for each merge event) and all births <= deaths."""
    result = persistent_homology(ivermectin_reviews)
    pd = result["persistence_diagram"]
    dim0 = [f for f in pd if f["dim"] == 0]
    # At least k-1 merges for k components -> 1
    assert len(dim0) >= len(ivermectin_reviews) - 1
    # All births <= deaths
    for f in pd:
        assert f["birth"] <= f["death"] + 1e-9
