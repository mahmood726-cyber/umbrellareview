import math
import pytest
from umbrella.models import ReviewInput
from umbrella.survival_meta import survival_meta


def _make_hr_reviews():
    """Create reviews with log(HR) as theta."""
    return [
        ReviewInput(review_id="SR1", theta=math.log(0.75), ci_lo=math.log(0.65),
                    ci_hi=math.log(0.87), k=5, measure="logHR", year=2018,
                    study_ids=["t1", "t2", "t3"]),
        ReviewInput(review_id="SR2", theta=math.log(0.80), ci_lo=math.log(0.68),
                    ci_hi=math.log(0.94), k=8, measure="logHR", year=2020,
                    study_ids=["t2", "t4", "t5"]),
        ReviewInput(review_id="SR3", theta=math.log(0.70), ci_lo=math.log(0.55),
                    ci_hi=math.log(0.89), k=3, measure="logHR", year=2022,
                    study_ids=["t1", "t5", "t6"]),
    ]


def test_pooled_hr_positive():
    """Pooled HR must always be > 0."""
    result = survival_meta(_make_hr_reviews())
    assert result["pooled_hr"] > 0


def test_hr_ci_brackets_point():
    """HR CI should bracket the pooled HR."""
    result = survival_meta(_make_hr_reviews())
    lo, hi = result["hr_ci"]
    assert lo <= result["pooled_hr"] <= hi


def test_frailty_adjusted_ge_pooled():
    """Frailty-adjusted HR >= pooled HR when tau2 >= 0 and pooled log-HR < 0.
       Under Jensen's correction exp(mu + tau2/2) vs exp(mu)."""
    result = survival_meta(_make_hr_reviews())
    # When mu < 0, exp(mu + tau2/2) > exp(mu) iff tau2 > 0
    if result["tau2_log_hr"] > 0:
        assert result["frailty_adjusted_hr"] >= result["pooled_hr"]


def test_median_survival_ratio_positive():
    """Median survival ratio must be > 0."""
    result = survival_meta(_make_hr_reviews())
    assert result["median_survival_ratio"] is not None
    assert result["median_survival_ratio"] > 0


def test_i2_bounded():
    """I-squared should be in [0, 100]."""
    result = survival_meta(_make_hr_reviews())
    assert 0.0 <= result["i2"] <= 100.0
