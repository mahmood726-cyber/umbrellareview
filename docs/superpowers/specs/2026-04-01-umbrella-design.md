# UmbrellaEngine — Systematic Review of Systematic Reviews

**Date**: 2026-04-01
**Status**: Approved
**Target Journal**: Journal of Clinical Epidemiology
**Location**: `C:\Models\UmbrellaReview\`

## Summary

Browser-based tool for conducting umbrella reviews (systematic reviews of systematic reviews). Automates the three hardest parts: study overlap quantification (CCA/GROOVE), quality assessment (AMSTAR-2), and discordance analysis (why do reviews disagree?). No existing browser tool handles all three.

## Architecture

- **Python engine** (`umbrella/`): Pure-function computations, numpy/scipy only
- **Browser app** (`app/umbrella.html`): Single-file HTML with Plotly.js
- **Test suite** (`tests/`): pytest, 25+ tests
- **TruthCert**: Hash-linked provenance

## Data Model

### ReviewInput
```python
@dataclass
class ReviewInput:
    review_id: str          # e.g. "Smith2022"
    theta: float            # pooled effect estimate
    ci_lo: float            # CI lower bound
    ci_hi: float            # CI upper bound
    se: float               # SE (derived from CI if not given)
    k: int                  # number of included studies
    study_ids: list[str]    # list of study identifiers (author-year or PMID)
    measure: str            # "OR", "RR", "HR", "SMD", "MD", "RD"
    amstar_items: dict[str, str]  # item_id -> "yes"/"partial"/"no"/"NA"
    scope_tags: list[str]   # e.g. ["adults", "RCTs_only", "mortality"]
    year: int               # publication year
    label: str              # short display label
```

### OverlapResult
```python
@dataclass
class OverlapResult:
    cca: float                              # corrected covered area [0,1]
    groove: str                             # "Slight"/"Moderate"/"High"/"Very High"
    n_total_citations: int                  # N (total study appearances)
    n_unique_studies: int                   # r (unique studies)
    n_reviews: int                          # c
    overlap_matrix: list[list[int]]         # c×c pairwise shared-study counts
    study_frequency: dict[str, int]         # study_id -> how many reviews include it
```

### AMSTARResult
```python
@dataclass
class AMSTARResult:
    review_id: str
    item_scores: dict[str, str]             # 16 items
    n_yes: int
    n_critical_yes: int                     # out of 7 critical domains
    n_critical_flaw: int                    # critical domains rated "no"
    confidence: str                         # "High"/"Moderate"/"Low"/"Critically Low"
```

### ConcordanceResult
```python
@dataclass
class ConcordanceResult:
    direction_agreement: float              # fraction of reviews agreeing on direction
    all_significant: bool                   # do all reviews find p<0.05?
    ci_overlap_fraction: float              # fraction of review pairs with overlapping CIs
    meta_meta_theta: float                  # pooled effect across reviews
    meta_meta_se: float
    meta_meta_i2: float                     # heterogeneity across reviews
    range_theta: tuple[float, float]        # min/max pooled effects
```

### DiscordanceResult
```python
@dataclass
class DiscordanceResult:
    overall_discordance: str                # "Concordant"/"Minor"/"Major"/"Contradictory"
    factors: list[DiscordanceFactor]        # ordered by contribution
```

### UmbrellaVerdict
```python
@dataclass
class UmbrellaVerdict:
    overlap: OverlapResult
    amstar_results: list[AMSTARResult]
    concordance: ConcordanceResult
    discordance: DiscordanceResult
    quality_weighted_theta: float
    quality_weighted_ci: tuple[float, float]
    overall_confidence: str
    recommendation: str
    certification: dict
```

## Statistical Methods

### 1. Overlap (CCA / GROOVE)

**Corrected Covered Area** (Pieper et al. 2014):
```
CCA = (N - r) / (r * c - r)
```
Where:
- N = total number of study appearances across all reviews (sum of all k)
- r = number of unique studies (union of all study_ids)
- c = number of reviews

**GROOVE classification**:
| CCA Range | Label | Interpretation |
|-----------|-------|----------------|
| < 0.05 | Slight | Reviews cover largely independent evidence |
| 0.05-0.10 | Moderate | Some shared studies |
| 0.10-0.15 | High | Substantial overlap — findings not independent |
| > 0.15 | Very High | Reviews are essentially re-analysing the same studies |

**Pairwise overlap matrix**: For each pair of reviews (i, j), count `|study_ids_i ∩ study_ids_j|`.

**Study frequency table**: For each unique study, count how many reviews include it. Flag studies appearing in all reviews ("ubiquitous") vs only one ("unique").

### 2. AMSTAR-2 Quality Assessment

**16 items** (Shea et al. 2017), 7 critical domains marked with *:

| # | Domain | Critical? |
|---|--------|-----------|
| 1 | PICO components in research questions | * |
| 2 | Protocol registered before review | * |
| 3 | Study design selection explained | |
| 4 | Comprehensive search strategy | * |
| 5 | Duplicate study selection | |
| 6 | Duplicate data extraction | |
| 7 | List of excluded studies with reasons | * |
| 8 | Adequate detail of included studies | |
| 9 | RoB assessment with satisfactory tool | * |
| 10 | Funding sources reported | |
| 11 | Appropriate statistical methods | * |
| 12 | RoB impact on results assessed | |
| 13 | RoB discussed in interpreting results | * |
| 14 | Heterogeneity discussed | |
| 15 | Publication bias investigated | |
| 16 | COI declared | |

**Confidence rating**:
- **High**: No critical flaws, max 1 non-critical weakness
- **Moderate**: No critical flaws, >1 non-critical weakness
- **Low**: 1 critical flaw ± non-critical weaknesses
- **Critically Low**: >1 critical flaw ± non-critical weaknesses

### 3. Effect Concordance

**Direction agreement**: Fraction of reviews where sign(theta) matches the majority direction.

**CI overlap**: For each pair of reviews, check if CIs overlap. Fraction of pairs with overlap.

**Meta-meta-analysis**: Pool the review-level effects using inverse-variance weighting:
```
w_i = 1 / se_i^2
theta_meta = sum(w_i * theta_i) / sum(w_i)
se_meta = 1 / sqrt(sum(w_i))
```
With between-review heterogeneity I²:
```
Q = sum(w_i * (theta_i - theta_meta)^2)
I² = max(0, (Q - (c-1)) / Q)
```

**Quality-weighted pooling**: Weight = AMSTAR-2 weight × precision:
- High = 1.0, Moderate = 0.75, Low = 0.5, Critically Low = 0.25
- `w_i_adj = amstar_weight * (1/se_i^2)`

### 4. Discordance Analysis

Classify overall discordance:
| Level | Criteria |
|-------|----------|
| Concordant | Direction agreement = 100%, all CIs overlap |
| Minor | Direction agreement = 100%, some CIs don't overlap |
| Major | Direction agreement < 100%, or meta-meta I² > 50% |
| Contradictory | Some reviews significant in opposite directions |

**Factor decomposition** — identify WHY reviews disagree:

1. **Scope difference**: Do reviews differ in population, intervention, comparator, or outcome? Compare `scope_tags` pairwise. Jaccard similarity < 0.5 = scope-driven.
2. **Inclusion criteria**: Different k implies different inclusion. Compute k ratio (max/min). Ratio > 2.0 = criteria-driven.
3. **Study overlap**: If CCA < 0.05 and reviews disagree, it's evidence-base-driven (different studies).
4. **Quality handling**: If removing low-AMSTAR reviews changes the conclusion, quality is a driver.
5. **Statistical methods**: If reviews use different measures (OR vs RR vs HR), method is a driver.

Each factor gets a contribution score (0-100) based on how much it explains the observed effect variation.

### 5. Verdict

Synthesize all components:
- **Quality-weighted pooled effect** with CI
- **Overall confidence**: min(GROOVE-implied independence, median AMSTAR-2)
- **Recommendation text**: "X of Y reviews agree on [direction]. Quality-weighted effect is [theta] [CI]. [Discordance level] discordance, primarily driven by [top factor]."

## Browser App Tabs (6)

### Tab 1: Input
- Add reviews via form (effect, CI, k, study IDs, AMSTAR items)
- Paste CSV/JSON batch import
- Built-in examples (3)
- Review table showing all entered reviews

### Tab 2: Overlap
- CCA value + GROOVE badge
- Heatmap: c×c overlap matrix (pairwise shared studies)
- Study frequency bar chart (how many reviews include each study)
- "Ubiquitous studies" list (in all reviews) and "unique studies" (in only one)

### Tab 3: Quality
- AMSTAR-2 traffic light grid: reviews (rows) × 16 items (columns), colored cells
- Critical domains highlighted with border
- Per-review confidence badge (High/Moderate/Low/Critically Low)
- Summary: median confidence, number of critically low reviews

### Tab 4: Effects
- Forest-of-forests plot: each review as a row with diamond + CI
- Quality-weighted pooled effect as summary diamond
- Meta-meta I² displayed
- Direction agreement badge
- CI overlap fraction

### Tab 5: Discordance
- Overall discordance badge (Concordant/Minor/Major/Contradictory)
- Waterfall chart: horizontal bars showing each factor's contribution (scope, inclusion, overlap, quality, methods)
- Pairwise discordance table: for each disagreeing pair, which factors differ

### Tab 6: Report & Certify
- Structured plain-text report
- TruthCert JSON bundle
- Methods paragraph for manuscript supplement
- PRISMA-O (overview) checklist items addressed

## Visualizations (5 Plotly charts)

1. **Overlap heatmap**: c×c matrix, color intensity = shared study count
2. **AMSTAR-2 traffic light**: Grid of colored cells (green/amber/red) with critical domain borders
3. **Forest-of-forests**: Review-level forest plot with quality-weighted summary
4. **Study frequency bar chart**: Studies sorted by frequency of inclusion
5. **Discordance waterfall**: Horizontal bars per factor, sorted by contribution

## Built-in Examples

### 1. Statins for Primary CV Prevention (5 SRs)
- Cholesterol Treatment Trialists 2012 — OR 0.73 (0.67-0.80), k=27
- Ray et al. 2010 — RR 0.91 (0.83-1.01), k=11
- Taylor et al. 2013 — RR 0.86 (0.79-0.94), k=18
- Tonelli et al. 2011 — RR 0.85 (0.77-0.95), k=9
- Chou et al. 2022 — RR 0.82 (0.74-0.91), k=23
- Scope: adults without CVD, all-cause mortality + CV events, RCTs only
- Expected: Minor discordance (Ray non-significant), high overlap

### 2. SGLT2i in Heart Failure (3 SRs)
- Zannad et al. 2020 — HR 0.67 (0.58-0.78), k=4
- Vaduganathan et al. 2022 — HR 0.70 (0.63-0.78), k=6
- McGuire et al. 2021 — HR 0.68 (0.60-0.77), k=5
- Expected: Concordant, high overlap

### 3. Ivermectin for COVID-19 (4 SRs)
- Hill et al. 2021 — OR 0.38 (0.19-0.73), k=15
- Roman et al. 2022 — RR 0.90 (0.57-1.42), k=10
- Popp et al. 2022 — RR 0.67 (0.40-1.12), k=14
- Bryant et al. 2021 — OR 0.31 (0.14-0.65), k=24
- Expected: Major discordance, quality-driven (AMSTAR-2 range: Critically Low to High)

## Test Coverage (25+ tests)

### overlap.py (6 tests)
- CCA for 2 reviews with no overlap = 0
- CCA for 2 reviews with full overlap = 1
- CCA for known example (3 reviews, published value)
- GROOVE classification at boundaries
- Overlap matrix symmetry
- Study frequency counts

### amstar.py (5 tests)
- All "yes" → High confidence
- 1 critical "no" → Low
- 2+ critical "no" → Critically Low
- No critical flaws, 2 non-critical → Moderate
- Item count validation (16 items)

### concordance.py (5 tests)
- All same direction → agreement = 1.0
- Mixed directions → agreement < 1.0
- Meta-meta I² = 0 when all effects identical
- Quality-weighted effect closer to high-quality reviews
- CI overlap detection

### discordance.py (4 tests)
- Concordant when all agree
- Contradictory when opposite significant results
- Scope-driven when scope_tags differ
- Factor contributions sum to ~100

### pipeline.py (5+ tests)
- End-to-end on each example
- Certification PASS/WARN
- Verdict contains all components

## File Structure

```
C:\Models\UmbrellaReview\
  umbrella/
    __init__.py
    models.py           # All dataclasses
    overlap.py           # CCA, GROOVE, intersection matrix
    amstar.py            # AMSTAR-2 scoring
    concordance.py       # Effect agreement, meta-meta-analysis
    discordance.py       # Disagreement decomposition
    verdict.py           # Quality-weighted synthesis
    pipeline.py          # run_umbrella() orchestrator
    certifier.py         # TruthCert
  tests/
    conftest.py          # Fixtures (statin, SGLT2i, ivermectin examples)
    test_overlap.py
    test_amstar.py
    test_concordance.py
    test_discordance.py
    test_pipeline.py
  app/
    umbrella.html        # Single-file browser app
  data/
    statins.json
    sglt2i.json
    ivermectin.json
  setup.py
  README.md
  LICENSE
```

## Out of Scope (v1)

- Automatic search/retrieval of SRs from PubMed/Cochrane
- Network of overlapping SRs (graph visualization)
- ROBIS (quality of SR methodology) — AMSTAR-2 only for v1
- Grading of Recommendations (GRADE-CERQual) — quantitative GRADE only
- Bayesian quality-weighted synthesis — frequentist only for v1
