import pytest
from umbrella.pipeline import run_umbrella
from umbrella.models import UmbrellaVerdict

def test_pipeline_returns_verdict(statin_reviews):
    result = run_umbrella(statin_reviews)
    assert isinstance(result, UmbrellaVerdict)

def test_pipeline_overlap_populated(statin_reviews):
    result = run_umbrella(statin_reviews)
    assert result.overlap.n_reviews == 5
    assert result.overlap.cca >= 0

def test_pipeline_amstar_all_scored(statin_reviews):
    result = run_umbrella(statin_reviews)
    assert len(result.amstar_results) == 5

def test_pipeline_concordance_populated(sglt2i_reviews):
    result = run_umbrella(sglt2i_reviews)
    assert result.concordance.direction_agreement == 1.0

def test_pipeline_discordance_for_ivermectin(ivermectin_reviews):
    result = run_umbrella(ivermectin_reviews)
    assert result.discordance.overall_discordance in ("Major", "Contradictory")

def test_pipeline_certification(statin_reviews):
    result = run_umbrella(statin_reviews)
    assert result.certification in ("PASS", "WARN", "REJECT")

def test_pipeline_recommendation_nonempty(statin_reviews):
    result = run_umbrella(statin_reviews)
    assert len(result.recommendation) > 20
