from umbrella.models import AMSTARResult

CRITICAL_ITEMS = {"1", "2", "4", "7", "9", "11", "13"}
ALL_ITEMS = {str(i) for i in range(1, 17)}

def score_amstar(review_id, items):
    """Score AMSTAR-2 and determine confidence rating."""
    n_yes = sum(1 for v in items.values() if v == "yes")
    n_critical_yes = sum(1 for k, v in items.items() if k in CRITICAL_ITEMS and v == "yes")
    n_critical_flaw = sum(1 for k, v in items.items() if k in CRITICAL_ITEMS and v == "no")
    n_noncritical_weakness = sum(1 for k, v in items.items()
                                  if k not in CRITICAL_ITEMS and v != "yes")

    if n_critical_flaw == 0:
        if n_noncritical_weakness <= 1:
            confidence = "High"
        else:
            confidence = "Moderate"
    elif n_critical_flaw == 1:
        confidence = "Low"
    else:
        confidence = "Critically Low"

    return AMSTARResult(
        review_id=review_id,
        item_scores=dict(items),
        n_yes=n_yes,
        n_critical_yes=n_critical_yes,
        n_critical_flaw=n_critical_flaw,
        confidence=confidence,
    )
