"""Dempster-Shafer Evidence Theory for Umbrella Reviews.

Model uncertain SR conclusions using belief functions over the frame
of discernment Omega = {beneficial, harmful, inconclusive}.  Each
review produces a mass function; Dempster's rule of combination
aggregates them.  Outputs belief, plausibility, pignistic probability,
and pairwise + overall conflict.
"""

from __future__ import annotations
from itertools import combinations

# ── Frame of discernment ────────────────────────────────────────────

OMEGA = frozenset({"beneficial", "harmful", "inconclusive"})

# All non-empty subsets of Omega (keys for mass functions)
_SUBSETS = [
    frozenset({"beneficial"}),
    frozenset({"harmful"}),
    frozenset({"inconclusive"}),
    frozenset({"beneficial", "harmful"}),
    frozenset({"beneficial", "inconclusive"}),
    frozenset({"harmful", "inconclusive"}),
    OMEGA,
]

_QUALITY_WEIGHT = {
    "High": 0.9,
    "Moderate": 0.7,
    "Low": 0.5,
    "Critically Low": 0.3,
}


# ── Helpers ─────────────────────────────────────────────────────────

def _amstar_confidence(review) -> str:
    """Derive AMSTAR-2 confidence label from amstar_items dict."""
    items = getattr(review, "amstar_items", None) or {}
    if not items:
        return ""
    critical = {"1", "2", "4", "7", "9", "11", "13"}
    n_crit_flaw = sum(1 for k, v in items.items() if k in critical and v == "no")
    n_noncrit_weak = sum(1 for k, v in items.items() if k not in critical and v != "yes")
    if n_crit_flaw == 0:
        return "High" if n_noncrit_weak <= 1 else "Moderate"
    if n_crit_flaw == 1:
        return "Low"
    return "Critically Low"


def _confidence_weight(review) -> float:
    conf = _amstar_confidence(review)
    return _QUALITY_WEIGHT.get(conf, 0.5)


def _mass_for_review(review) -> dict[frozenset, float]:
    """Build a mass function for a single review."""
    cw = _confidence_weight(review)
    theta = review.theta
    ci_lo = review.ci_lo
    ci_hi = review.ci_hi
    m: dict[frozenset, float] = {}

    if theta > 0 and ci_lo > 0:
        # CI excludes 0, beneficial
        m[frozenset({"beneficial"})] = cw
        m[OMEGA] = 1.0 - cw
    elif theta < 0 and ci_hi < 0:
        # CI excludes 0, harmful
        m[frozenset({"harmful"})] = cw
        m[OMEGA] = 1.0 - cw
    else:
        # CI includes 0 — inconclusive
        m[frozenset({"inconclusive"})] = 0.5 * cw
        m[OMEGA] = 1.0 - 0.5 * cw

    return m


# ── Dempster's rule of combination ─────────────────────────────────

def _combine_two(m1: dict[frozenset, float],
                 m2: dict[frozenset, float]) -> tuple[dict[frozenset, float], float]:
    """Combine two mass functions.  Return (combined mass, conflict K)."""
    # Accumulate unnormalised masses
    raw: dict[frozenset, float] = {}
    conflict = 0.0

    for a, ma in m1.items():
        for b, mb in m2.items():
            intersection = a & b
            prod = ma * mb
            if len(intersection) == 0:
                conflict += prod
            else:
                raw[intersection] = raw.get(intersection, 0.0) + prod

    # Normalise
    denom = 1.0 - conflict
    if denom <= 0:
        # Total conflict — degenerate case
        denom = 1e-15

    combined: dict[frozenset, float] = {}
    for s, v in raw.items():
        combined[s] = v / denom

    return combined, conflict


def _belief(m: dict[frozenset, float], a: frozenset) -> float:
    """Bel(A) = sum of m(B) for all B subset of A."""
    return sum(v for s, v in m.items() if s <= a)


def _plausibility(m: dict[frozenset, float], a: frozenset) -> float:
    """Pl(A) = sum of m(B) for all B that intersect A."""
    return sum(v for s, v in m.items() if len(s & a) > 0)


def _pignistic(m: dict[frozenset, float]) -> dict[str, float]:
    """Pignistic probability BetP for each singleton in Omega."""
    betp: dict[str, float] = {h: 0.0 for h in OMEGA}
    for s, ms in m.items():
        card = len(s)
        share = ms / card
        for h in s:
            betp[h] += share
    return betp


# ── Public API ──────────────────────────────────────────────────────

def _subset_label(s: frozenset) -> str:
    return ",".join(sorted(s))


def dempster_shafer(reviews) -> dict:
    """Run Dempster-Shafer analysis on a list of ReviewInput objects.

    Returns dict with:
        combined_mass   — dict[str, float]  (subset label -> mass)
        belief          — dict[str, float]  (singleton label -> Bel)
        plausibility    — dict[str, float]  (singleton label -> Pl)
        pignistic       — dict[str, float]  (singleton label -> BetP)
        conflict_total  — float
        conflict_pairwise — dict[str, float]  (pair label -> K)
        verdict         — str (highest pignistic)
    """
    if not reviews:
        return {
            "combined_mass": {},
            "belief": {},
            "plausibility": {},
            "pignistic": {},
            "conflict_total": 0.0,
            "conflict_pairwise": {},
            "verdict": "inconclusive",
        }

    # Build individual mass functions
    masses = [(r.review_id, _mass_for_review(r)) for r in reviews]

    # Pairwise conflict
    conflict_pairwise: dict[str, float] = {}
    for (id_a, m_a), (id_b, m_b) in combinations(masses, 2):
        _, k = _combine_two(m_a, m_b)
        conflict_pairwise[f"{id_a} vs {id_b}"] = k

    # Sequential combination
    combined = masses[0][1]
    total_conflict = 0.0
    for i in range(1, len(masses)):
        combined, k = _combine_two(combined, masses[i][1])
        total_conflict = 1.0 - (1.0 - total_conflict) * (1.0 - k)

    # Belief, plausibility for each singleton
    singletons = [frozenset({h}) for h in sorted(OMEGA)]
    belief_dict = {h: _belief(combined, frozenset({h})) for h in sorted(OMEGA)}
    pl_dict = {h: _plausibility(combined, frozenset({h})) for h in sorted(OMEGA)}

    # Pignistic
    pig = _pignistic(combined)

    # Verdict
    verdict = max(pig, key=lambda h: pig[h])

    # Serialise combined mass with string keys
    combined_str = {_subset_label(s): v for s, v in combined.items()}

    return {
        "combined_mass": combined_str,
        "belief": belief_dict,
        "plausibility": pl_dict,
        "pignistic": pig,
        "conflict_total": total_conflict,
        "conflict_pairwise": conflict_pairwise,
        "verdict": verdict,
    }
