import json
import pytest
from umbrella.models import ReviewInput

def _load_reviews(path):
    with open(path) as f:
        data = json.load(f)
    return [ReviewInput(**r) for r in data["reviews"]]

@pytest.fixture
def statin_reviews():
    return _load_reviews("data/statins.json")

@pytest.fixture
def sglt2i_reviews():
    return _load_reviews("data/sglt2i.json")

@pytest.fixture
def ivermectin_reviews():
    return _load_reviews("data/ivermectin.json")

@pytest.fixture
def two_no_overlap():
    """Two reviews with zero study overlap."""
    return [
        ReviewInput(review_id="A", theta=-0.5, ci_lo=-0.8, ci_hi=-0.2, k=3,
                    study_ids=["s1", "s2", "s3"], measure="logOR"),
        ReviewInput(review_id="B", theta=-0.4, ci_lo=-0.7, ci_hi=-0.1, k=3,
                    study_ids=["s4", "s5", "s6"], measure="logOR"),
    ]

@pytest.fixture
def two_full_overlap():
    """Two reviews with identical study lists."""
    return [
        ReviewInput(review_id="A", theta=-0.5, ci_lo=-0.8, ci_hi=-0.2, k=3,
                    study_ids=["s1", "s2", "s3"], measure="logOR"),
        ReviewInput(review_id="B", theta=-0.6, ci_lo=-0.9, ci_hi=-0.3, k=3,
                    study_ids=["s1", "s2", "s3"], measure="logOR"),
    ]
