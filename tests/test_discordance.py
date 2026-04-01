import pytest
from umbrella.discordance import compute_discordance
from umbrella.overlap import compute_overlap
from umbrella.concordance import compute_concordance

def test_concordant_when_all_agree(sglt2i_reviews):
    overlap = compute_overlap(sglt2i_reviews)
    conc = compute_concordance(sglt2i_reviews)
    result = compute_discordance(sglt2i_reviews, overlap, conc)
    assert result.overall_discordance == "Concordant"

def test_discordant_for_ivermectin(ivermectin_reviews):
    overlap = compute_overlap(ivermectin_reviews)
    conc = compute_concordance(ivermectin_reviews)
    result = compute_discordance(ivermectin_reviews, overlap, conc)
    assert result.overall_discordance in ("Major", "Contradictory")

def test_factors_nonempty_when_discordant(ivermectin_reviews):
    overlap = compute_overlap(ivermectin_reviews)
    conc = compute_concordance(ivermectin_reviews)
    result = compute_discordance(ivermectin_reviews, overlap, conc)
    assert len(result.factors) > 0

def test_factor_contributions_bounded(ivermectin_reviews):
    overlap = compute_overlap(ivermectin_reviews)
    conc = compute_concordance(ivermectin_reviews)
    result = compute_discordance(ivermectin_reviews, overlap, conc)
    for f in result.factors:
        assert 0 <= f.contribution <= 100
