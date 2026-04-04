import pytest
from umbrella.empirical_process import compute_gof


def test_ks_p_in_range(statin_reviews):
    """KS p-value must be in [0, 1]."""
    result = compute_gof(statin_reviews)
    assert 0.0 <= result["ks_p"] <= 1.0


def test_ks_stat_nonnegative(sglt2i_reviews):
    """KS statistic is always >= 0."""
    result = compute_gof(sglt2i_reviews)
    assert result["ks_stat"] >= 0.0


def test_qq_lengths_match(statin_reviews):
    """QQ theoretical and empirical arrays have same length as input."""
    result = compute_gof(statin_reviews)
    n = len(statin_reviews)
    assert len(result["qq_theoretical"]) == n
    assert len(result["qq_empirical"]) == n


def test_ad_cvm_nonnegative(ivermectin_reviews):
    """Anderson-Darling and Cramer-von Mises stats are >= 0."""
    result = compute_gof(ivermectin_reviews)
    assert result["ad_stat"] >= 0.0
    assert result["cvm_stat"] >= 0.0


def test_is_normal_boolean(statin_reviews):
    """is_normal flag is a boolean."""
    result = compute_gof(statin_reviews)
    assert isinstance(result["is_normal"], bool)
