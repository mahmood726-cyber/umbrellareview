import pytest
from umbrella.causal_discordance import causal_discordance


def test_r_squared_in_unit_interval(statin_reviews):
    """R-squared must be in [0, 1]."""
    result = causal_discordance(statin_reviews)
    assert 0.0 <= result["r_squared"] <= 1.0


def test_structural_coefficients_keys(statin_reviews):
    """All expected variables present with beta, se, p."""
    result = causal_discordance(statin_reviews)
    coefs = result["structural_coefficients"]
    for var in ["intercept", "quality", "k", "scope", "overlap"]:
        assert var in coefs
        assert "beta" in coefs[var]
        assert "se" in coefs[var]
        assert "p" in coefs[var]
        assert coefs[var]["se"] >= 0.0  # SE non-negative


def test_counterfactual_returns_finite(ivermectin_reviews):
    """Counterfactual high-quality estimate must be finite."""
    result = causal_discordance(ivermectin_reviews)
    cf = result["counterfactual_high_quality"]
    assert isinstance(cf["mean_theta"], float)
    assert isinstance(cf["delta"], float)
    assert not (cf["mean_theta"] != cf["mean_theta"])  # NaN check


def test_mediation_decomposition_sums(statin_reviews):
    """Indirect + direct should approximately equal total (Baron-Kenny)."""
    result = causal_discordance(statin_reviews)
    med = result["mediation"]
    # total ~= direct + indirect (exact in linear OLS)
    assert med["total"] == pytest.approx(
        med["direct"] + med["indirect"], abs=0.01
    )


def test_residual_discordance_non_negative(sglt2i_reviews):
    """RMSE must be >= 0."""
    result = causal_discordance(sglt2i_reviews)
    assert result["residual_discordance"] >= 0.0
