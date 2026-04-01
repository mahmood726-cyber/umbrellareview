import pytest
from umbrella.concordance import compute_concordance

def test_all_same_direction(sglt2i_reviews):
    result = compute_concordance(sglt2i_reviews)
    assert result.direction_agreement == 1.0

def test_meta_meta_i2_low_for_concordant(sglt2i_reviews):
    result = compute_concordance(sglt2i_reviews)
    assert result.meta_meta_i2 < 50.0

def test_ci_overlap_high_for_concordant(sglt2i_reviews):
    result = compute_concordance(sglt2i_reviews)
    assert result.ci_overlap_fraction > 0.5

def test_range_theta_correct(two_full_overlap):
    result = compute_concordance(two_full_overlap)
    assert result.range_theta[0] <= result.range_theta[1]

def test_quality_weighted_theta_finite(statin_reviews):
    result = compute_concordance(statin_reviews)
    assert -5.0 < result.quality_weighted_theta < 5.0
