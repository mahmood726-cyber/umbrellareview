import pytest
from umbrella.resampling import compute_resampling


def test_jackknife_se_positive(statin_reviews):
    """Jackknife SE must be > 0."""
    result = compute_resampling(statin_reviews)
    assert result["jackknife_se"] > 0.0


def test_bootstrap_percentile_contains_estimate(statin_reviews):
    """Percentile bootstrap CI should contain the jackknife estimate."""
    result = compute_resampling(statin_reviews)
    lo, hi = result["bootstrap_ci_percentile"]
    assert lo <= result["jackknife_estimate"] <= hi


def test_bca_ci_contains_estimate(sglt2i_reviews):
    """BCa bootstrap CI should contain the jackknife estimate."""
    result = compute_resampling(sglt2i_reviews)
    lo, hi = result["bootstrap_ci_bca"]
    assert lo <= result["jackknife_estimate"] <= hi


def test_permutation_p_in_range(ivermectin_reviews):
    """Permutation p-value must be in (0, 1]."""
    result = compute_resampling(ivermectin_reviews)
    assert 0.0 < result["permutation_p"] <= 1.0


def test_subsampling_stability_bounded(statin_reviews):
    """Subsampling stability should be <= 1."""
    result = compute_resampling(statin_reviews)
    assert result["subsampling_stability"] <= 1.0
