import pytest
from umbrella.meta_regression import compute_meta_regression


def test_r_squared_in_range(statin_reviews):
    """R-squared is between 0 and 1."""
    result = compute_meta_regression(statin_reviews)
    assert 0.0 <= result["r_squared"] <= 1.0


def test_coefficients_present(statin_reviews):
    """All expected covariates appear in coefficients."""
    result = compute_meta_regression(statin_reviews)
    expected = {"intercept", "year", "amstar_quality", "k", "scope_breadth"}
    assert set(result["coefficients"].keys()) == expected


def test_residuals_length(statin_reviews):
    """Residuals list matches number of reviews used."""
    result = compute_meta_regression(statin_reviews)
    assert len(result["residuals"]) == result["n_reviews_used"]


def test_permutation_p_in_range(statin_reviews):
    """Permutation p-value is between 0 and 1."""
    result = compute_meta_regression(statin_reviews, n_perm=99)
    assert 0.0 <= result["permutation_p"] <= 1.0


def test_too_few_reviews_raises(two_no_overlap):
    """Regression with <3 reviews raises ValueError."""
    with pytest.raises(ValueError, match="at least 3"):
        compute_meta_regression(two_no_overlap)
