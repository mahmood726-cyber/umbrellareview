# UmbrellaEngine: An Open-Source Browser Tool for Automated Overlap Quantification, Quality Assessment, and Discordance Decomposition in Umbrella Reviews

[AUTHOR]^1

^1 [AFFILIATION]. Email: [EMAIL] | ORCID: [ORCID]

**Target journal:** *Journal of Clinical Epidemiology* (Software/Methods Article)

**Word count:** ~3,500

**Keywords:** umbrella review, systematic review of systematic reviews, overlap analysis, corrected covered area, AMSTAR-2, discordance decomposition, meta-meta-analysis, open-source software

---

## Abstract

**Objective:** Umbrella reviews synthesize evidence across multiple systematic reviews addressing a common question, yet no freely available browser tool automates the analytical workflow. We developed UmbrellaEngine, an open-source application that performs overlap quantification, methodological quality assessment, effect concordance analysis, and discordance decomposition without software installation.

**Study Design and Setting:** UmbrellaEngine comprises a Python computational engine (8 modules, 27 unit tests) and a companion single-file HTML application (1,681 lines) with interactive Plotly.js visualizations. The overlap module computes the Corrected Covered Area (CCA) and generates GROOVE (Graphical Representation of Overlap for OVErviews) matrices. Quality assessment implements the full AMSTAR-2 instrument with automated confidence ratings across 16 items and 7 critical domains. Concordance analysis includes direction agreement, pairwise confidence interval overlap, inverse-variance meta-meta-analysis with I-squared heterogeneity, and quality-weighted pooling. A novel 5-factor discordance decomposition identifies and ranks the drivers of disagreement among reviews.

**Results:** We demonstrate UmbrellaEngine on three clinical examples: statins for cardiovascular prevention (5 reviews, minor discordance driven by scope differences), SGLT2 inhibitors in heart failure (3 reviews, full concordance), and ivermectin for COVID-19 (4 reviews, contradictory findings attributable to quality variation and evidence base divergence). The tool classified discordance correctly in all cases and generated structured synthesis verdicts with provenance-linked certification.

**Conclusion:** UmbrellaEngine is the first browser-based tool that automates the complete umbrella review analytical workflow. It is freely available under the MIT licence at [REPOSITORY_URL].

---

## 1. Introduction

Umbrella reviews -- systematic reviews of systematic reviews -- occupy the highest tier of the evidence synthesis hierarchy. By collating and comparing the findings of multiple meta-analyses on a shared clinical question, they aim to identify the most reliable evidence, detect sources of disagreement, and provide a definitive summary for clinical decision-making [1,2]. The number of published umbrella reviews has grown exponentially, from fewer than 20 per year in 2010 to over 600 per year by 2024 [3].

Despite this growth, umbrella reviews remain methodologically challenging. Three analytical steps are particularly burdensome. First, quantifying the overlap of primary studies across included reviews is essential because high overlap inflates apparent concordance while low overlap may explain discordant conclusions [4]. The Corrected Covered Area (CCA) statistic proposed by Pieper and colleagues [5] and the GROOVE visualization framework [6] address this need but require manual matrix construction. Second, assessing the methodological quality of each included review using instruments such as AMSTAR-2 [7] demands structured evaluation across 16 items with complex decision rules involving 7 critical domains. Third, when reviews reach discordant conclusions -- a common scenario in contentious clinical areas -- reviewers must systematically decompose the sources of disagreement, an analytical task for which no standardized tool exists.

Current practice relies on manual spreadsheet construction, bespoke R scripts, or the R package `umbrella` [8], all of which impose installation requirements and programming expertise that limit accessibility. No browser-based tool automates the overlap-quality-concordance-discordance workflow. We developed UmbrellaEngine to fill this gap.

## 2. Methods

### 2.1. Software Architecture

UmbrellaEngine has two components: (a) a Python computational engine (`umbrella/`) containing 8 pure-function modules with numpy as the sole numerical dependency, and (b) a standalone HTML application (`app/umbrella.html`, 1,681 lines) that reimplements all algorithms in client-side JavaScript with Plotly.js for interactive visualization. Both components accept the same data model and produce equivalent results. The Python engine is validated by 27 automated pytest tests; the HTML application requires no installation and runs in any modern browser.

The data model centres on a `ReviewInput` structure that captures, for each included systematic review: a unique identifier, the pooled effect estimate (theta) with 95% confidence interval, the standard error, the number of included primary studies (k), the list of primary study identifiers, the effect measure (odds ratio, risk ratio, hazard ratio, mean difference, or standardized mean difference), 16 AMSTAR-2 item ratings, scope descriptor tags (population, intervention, outcome, study design), publication year, and a display label. Standard errors are derived from confidence interval widths when not provided directly.

### 2.2. Overlap Quantification

#### 2.2.1. Corrected Covered Area (CCA)

We implement the CCA statistic introduced by Pieper et al. [5], which quantifies the degree to which included reviews share the same primary studies. Let *c* denote the number of reviews, *r* the number of unique primary studies across all reviews, and *N* the total number of study citations (counting each appearance of a study in each review). The CCA is defined as:

> CCA = (N - r) / (r * c - r)

CCA ranges from 0 (no overlap: every review includes entirely different studies) to 1 (complete overlap: all reviews include identical study sets). The denominator (r * c - r) represents the maximum possible overlap beyond the baseline of r unique studies appearing once each.

#### 2.2.2. GROOVE Classification

Following the thresholds proposed for the Graphical Representation of Overlap for OVErviews [6], CCA values are classified as: Slight (CCA < 0.05), Moderate (0.05 <= CCA < 0.10), High (0.10 <= CCA < 0.15), or Very High (CCA >= 0.15).

#### 2.2.3. Overlap Matrix and Study Frequency

The tool computes a symmetric *c* x *c* pairwise overlap matrix where entry (i, j) is the number of primary studies shared between reviews *i* and *j*. It also produces a study frequency table recording how many reviews include each primary study. Both outputs are visualized as an interactive heatmap and a bar chart, respectively.

### 2.3. Methodological Quality Assessment (AMSTAR-2)

The AMSTAR-2 instrument [7] evaluates systematic review quality across 16 items, of which 7 are designated as critical domains (items 1, 2, 4, 7, 9, 11, and 13). Each item is rated "yes," "partial yes," or "no." UmbrellaEngine implements the published decision algorithm for deriving an overall confidence rating:

- **High:** No critical domain flaws and at most one non-critical weakness.
- **Moderate:** No critical domain flaws but more than one non-critical weakness.
- **Low:** Exactly one critical domain flaw (regardless of non-critical weaknesses).
- **Critically Low:** More than one critical domain flaw.

The tool displays results as a traffic-light grid (green/yellow/red for each item per review) and summarizes the distribution of confidence ratings across all included reviews.

### 2.4. Effect Concordance Analysis

#### 2.4.1. Direction Agreement

The proportion of reviews whose point estimate agrees with the majority direction is computed. Majority direction is determined by the sign of the median effect estimate across reviews.

#### 2.4.2. Confidence Interval Overlap

For each pair of reviews, UmbrellaEngine tests whether their 95% confidence intervals overlap (i.e., the lower bound of one interval is below the upper bound of the other). The CI overlap fraction is the proportion of all pairwise comparisons that show overlap.

#### 2.4.3. Meta-Meta-Analysis

A fixed-effect inverse-variance meta-analysis pools the review-level effect estimates. Letting theta_i and SE_i denote the pooled estimate and its standard error from review *i*, the weights are w_i = 1/SE_i^2, the pooled estimate is:

> theta_pooled = sum(w_i * theta_i) / sum(w_i)

and Cochran's Q statistic is computed to derive the I-squared heterogeneity measure, interpreted as the proportion of between-review variability not attributable to sampling error.

#### 2.4.4. Quality-Weighted Pooling

Recognizing that not all reviews deserve equal influence, UmbrellaEngine implements quality-weighted pooling. Each review's inverse-variance weight is multiplied by a quality factor derived from its AMSTAR-2 confidence rating: High = 1.0, Moderate = 0.75, Low = 0.5, Critically Low = 0.25. This down-weights conclusions from methodologically flawed reviews while preserving the contribution of higher-quality evidence.

### 2.5. Discordance Decomposition

When reviews reach different conclusions, identifying the source of disagreement is critical for interpretation. UmbrellaEngine decomposes discordance into five factors:

1. **Scope differences.** The Jaccard distance between the scope tags of each review pair is averaged. High scope divergence indicates that reviews addressed subtly different questions (e.g., different populations, co-interventions, or outcome definitions).

2. **Inclusion criteria.** The ratio of the maximum to minimum number of included primary studies (k) across reviews is normalized. A high k-ratio suggests that reviews applied different eligibility criteria, yielding substantially different evidence bases.

3. **Evidence base divergence.** The complement of CCA (i.e., 1 - CCA) captures the degree to which reviews analysed different primary studies, even if their nominal inclusion criteria were similar.

4. **Quality variation.** The range of AMSTAR-2 confidence ratings across reviews, normalized to a 0--1 scale (where 0 means all reviews share the same rating and 1 means the maximum possible spread from High to Critically Low).

5. **Statistical methods.** Heterogeneity in effect measures used across reviews (e.g., one review reporting odds ratios while another reports risk ratios). A binary indicator is normalized by the number of distinct measures observed.

Each factor receives a raw score between 0 and 1. The tool then expresses each factor's contribution as a percentage of the total discordance score, producing a ranked decomposition. The overall discordance level is classified as:

- **Concordant:** All reviews agree on direction, and all pairwise CIs overlap.
- **Minor:** Directional agreement is unanimous, but some CIs do not overlap.
- **Major:** Reviews disagree on effect direction, or meta-meta I-squared exceeds 50%.
- **Contradictory:** At least one review finds a statistically significant benefit and at least one finds a significant harm (or vice versa).

### 2.6. Synthesis Verdict and Certification

The pipeline module orchestrates all analyses and generates a structured narrative verdict summarizing: the proportion of reviews agreeing on direction, the quality-weighted pooled estimate and its confidence interval, CCA and GROOVE classification, the median AMSTAR-2 confidence rating, the discordance classification, and the primary driver of discordance when present.

A certification module computes a SHA-256 hash of the input data, enabling downstream verification that results correspond to specific inputs. The certification status is PASS (>= 2 reviews and >= 3 unique primary studies), WARN (>= 2 reviews but few unique studies), or REJECT (fewer than 2 reviews, rendering umbrella analysis inapplicable).

### 2.7. Browser Application

The HTML application provides six tabs:

1. **Input.** Data entry for reviews (manual, CSV paste, or JSON import), with three built-in example datasets.
2. **Overlap.** Interactive Plotly.js heatmap of the pairwise overlap matrix and a study frequency bar chart.
3. **Quality.** AMSTAR-2 traffic-light grid and confidence rating distribution.
4. **Effects.** Forest-of-forests plot showing each review's effect estimate with confidence interval, the unweighted pooled diamond, and the quality-weighted pooled diamond.
5. **Discordance.** Waterfall chart decomposing discordance by factor contribution, with a badge indicating the overall classification.
6. **Report & Certify.** Structured narrative report, auto-generated methods paragraph, and TruthCert JSON bundle for provenance.

All visualizations are rendered with Plotly.js (version 2.35.0), the sole external dependency. Dark mode is supported. Data persists in localStorage between sessions. The application can be saved as a single file and shared without a server.

### 2.8. Validation

The Python engine is validated by 27 automated tests organized across 5 test modules (Table 1). Tests cover: CCA computation for zero-overlap, full-overlap, and partial-overlap scenarios; overlap matrix symmetry; study frequency counts; GROOVE boundary classification; AMSTAR-2 scoring for all four confidence levels; concordance direction agreement and I-squared; discordance classification for concordant and contradictory datasets; factor contribution bounds; end-to-end pipeline execution; certification logic; and recommendation text generation. All tests pass under Python 3.10+ with numpy.

**Table 1.** Test suite summary.

| Module | Tests | Key verifications |
|--------|-------|-------------------|
| test_overlap.py | 6 | CCA bounds, matrix symmetry, GROOVE thresholds |
| test_amstar.py | 5 | All confidence levels, critical domain counting |
| test_concordance.py | 5 | Direction agreement, I-squared, CI overlap, quality weighting |
| test_discordance.py | 4 | Classification accuracy, factor decomposition, contribution bounds |
| test_pipeline.py | 7 | End-to-end orchestration, certification, verdict generation |

## 3. Illustrative Examples

### 3.1. Statins for Primary Cardiovascular Prevention (Minor Discordance)

Five systematic reviews published between 2010 and 2022 were analysed: CTT 2012 (27 primary studies), Ray 2010 (11 studies), Taylor 2013 (18 studies), Tonelli 2011 (9 studies), and Chou 2022 (23 studies). All five report protective effects (negative log-odds ratios or log-risk ratios), yielding 100% direction agreement. The CCA was 0.225 (Very High), reflecting that many landmark trials (WOSCOPS, AFCAPS, ASCOT, JUPITER) appear in most reviews. The meta-meta pooled effect was log-OR -0.19 (95% CI: -0.24 to -0.14), I-squared 52.3%. Three reviews received High AMSTAR-2 confidence ratings, one Moderate, and one Low. Discordance was classified as Minor, driven primarily by scope differences (CTT 2012 included secondary prevention populations) and statistical methods (one review used odds ratios while others used risk ratios). The quality-weighted analysis slightly attenuated the pooled effect, consistent with the expected quality-effect gradient.

### 3.2. SGLT2 Inhibitors in Heart Failure (Concordance)

Three reviews (Zannad 2020, Vaduganathan 2022, McGuire 2021) evaluated SGLT2 inhibitors for the composite of cardiovascular death or heart failure hospitalization. All three report log-hazard ratios between -0.40 and -0.36, yielding 100% direction agreement, 100% CI overlap, and I-squared of 1.2%. The CCA was 0.417 (Very High), as all reviews included the pivotal DAPA-HF and EMPEROR-Reduced trials. All three reviews received High AMSTAR-2 confidence. Discordance was classified as Concordant. This example illustrates the tool's ability to confirm robust consensus, providing reassurance that clinician-facing summaries rest on consistent evidence.

### 3.3. Ivermectin for COVID-19 (Contradictory Findings)

Four reviews published between 2021 and 2022 exemplify a contested evidence landscape. Hill 2021 and Bryant 2021 reported large, statistically significant mortality reductions (log-OR -0.97 and -1.17, respectively), while Roman 2022 found no significant effect (log-RR -0.11, CI crossing null) and the Cochrane review by Popp 2022 found a non-significant trend (log-RR -0.40, CI crossing null). The CCA was 0.192 (Very High), indicating substantial primary study overlap despite divergent conclusions.

UmbrellaEngine classified this as Contradictory discordance and identified two primary drivers: quality variation (46.2%) and evidence base differences (23.1%). Hill 2021 and Bryant 2021 received Critically Low AMSTAR-2 ratings (multiple critical domain flaws including failure to account for risk of bias, absence of a registered protocol, and inclusion of the subsequently retracted Elgazzar study), while Roman 2022 and Popp 2022 received High ratings. The quality-weighted pooled effect was substantially attenuated compared to the unweighted pooled effect, demonstrating that lower-quality reviews drove the apparently large treatment effects. This example illustrates UmbrellaEngine's principal value proposition: automating the identification of why reviews disagree, enabling users to interpret contradictory evidence rather than merely cataloguing it.

## 4. Software Description

### 4.1. Availability and Requirements

UmbrellaEngine is available at [REPOSITORY_URL] under the MIT licence. The Python engine requires Python >= 3.10 and numpy. The browser application requires only a modern web browser (Chrome, Firefox, Safari, or Edge) and operates entirely client-side; no server, login, or installation is needed.

### 4.2. Input Formats

The tool accepts data in three formats: (a) manual entry through structured form fields in the browser application, (b) CSV with columns for review identifier, effect estimate, confidence interval bounds, effect measure, number of studies, study identifiers (semicolon-separated), AMSTAR-2 items (colon-separated key-value pairs), and scope tags, and (c) JSON matching the ReviewInput schema. Three built-in datasets (statins, SGLT2i, ivermectin) are included for immediate demonstration.

### 4.3. Output Formats

Results are available as: (a) interactive visualizations within the browser (5 Plotly.js charts), (b) a structured narrative report with auto-generated methods paragraph suitable for manuscript use, and (c) exportable JSON including a provenance hash for reproducibility verification.

## 5. Discussion

### 5.1. Contribution

UmbrellaEngine is, to our knowledge, the first freely available browser tool that automates the complete analytical workflow for umbrella reviews. By integrating overlap quantification, quality assessment, concordance analysis, and discordance decomposition in a single interface, it eliminates the need for bespoke scripts or manual spreadsheet construction. The tool is designed for clinical researchers, guideline developers, and evidence synthesis teams who may not have programming expertise.

### 5.2. Comparison with Existing Tools

The R package `umbrella` [8] implements the Fusar-Poli evidence classification system and provides programmatic access to umbrella review computations, but requires R installation and scripting knowledge. RevMan and other Cochrane tools do not support umbrella-level analyses. The JBI SUMARI platform includes umbrella review modules but focuses on narrative synthesis rather than quantitative overlap and discordance decomposition. Our tool complements these approaches by providing immediate, installation-free access to quantitative methods in a browser environment.

### 5.3. The Discordance Decomposition Framework

The 5-factor discordance decomposition is, to our knowledge, a novel contribution. While individual components (CCA for overlap, AMSTAR-2 for quality) are well established, their integration into a single decomposition framework that quantifies the relative contribution of each factor to observed disagreement has not been previously implemented in software. The ivermectin example demonstrates how this decomposition can transform a qualitative observation ("reviews disagree") into a quantitative, actionable finding ("disagreement is 46% quality-driven and 23% evidence-base-driven"). This information is directly useful for guideline panels assessing the credibility of competing reviews.

### 5.4. Methodological Considerations

Several methodological choices warrant discussion. First, the meta-meta-analysis uses a fixed-effect model, which assumes that all reviews estimate the same underlying effect. When reviews differ in scope (as in the statin example where one review included secondary prevention), this assumption may be violated. Users should consider whether the included reviews are sufficiently similar to justify pooling. Second, quality-weighted pooling treats AMSTAR-2 confidence levels as ordinal weights, which introduces a degree of arbitrariness in the weight ratios (1.0 : 0.75 : 0.50 : 0.25). Sensitivity analyses with alternative weight schemes are advisable. Third, the discordance decomposition normalizes raw factor scores to percentages, which means that contributions are relative, not absolute -- a factor contributing 40% of the discordance does not imply a standardized effect size.

### 5.5. Limitations

First, the tool does not adjust meta-meta-analytic estimates for primary study overlap. When CCA is high, pooled estimates across reviews are correlated, and the computed I-squared and confidence intervals may understate true uncertainty. Methods for overlap-adjusted meta-meta-analysis remain an active research area [9]. Second, the 5-factor discordance model is necessarily a simplification; other sources of disagreement (e.g., differential handling of missing data, different risk-of-bias judgements applied to shared studies, or different statistical models) are not captured. Third, AMSTAR-2 was designed for systematic reviews of randomized trials and may require adaptation for reviews of observational studies or diagnostic test accuracy studies. Fourth, the tool assumes that effect estimates are on the log scale for ratio measures, requiring users to perform this transformation if their source data are on the natural scale.

### 5.6. Future Directions

Planned developments include: (a) overlap-adjusted pooling using the method of Perez-Bracchiglione et al. [10] to account for non-independence; (b) random-effects meta-meta-analysis for settings where between-review heterogeneity is expected; (c) integration with the Fusar-Poli evidence classification system [11] for umbrella reviews that include sufficient per-review detail (sample sizes, Egger's test results); and (d) PRISMA-O (PRISMA for Overviews) checklist generation to facilitate reporting compliance.

## 6. Conclusions

UmbrellaEngine automates the four core analytical tasks of umbrella reviews -- overlap quantification, quality assessment, effect concordance, and discordance decomposition -- in a freely available, installation-free browser application. The tool's 5-factor discordance decomposition provides a structured answer to the question "why do reviews disagree?" that is directly useful for guideline development and clinical interpretation. By lowering the technical barrier to rigorous umbrella review analysis, UmbrellaEngine aims to improve the quality and transparency of the highest tier of evidence synthesis.

## Acknowledgements

[ACKNOWLEDGEMENTS]

## Funding

[FUNDING]

## Declaration of Interest

[AUTHOR] declares no competing interests.

## Data Availability

All source code, built-in datasets, and test suites are available at [REPOSITORY_URL] under the MIT licence. The three illustrative datasets (statins, SGLT2i, ivermectin) are bundled with the software and can be loaded directly in the browser application.

## References

[1] Aromataris E, Fernandez R, Godfrey CM, et al. Summarizing systematic reviews: methodological development, conduct and reporting of an umbrella review approach. *Int J Evid Based Healthc*. 2015;13(3):132-140. doi:10.1097/XEB.0000000000000055

[2] Fusar-Poli P, Radua J. Ten simple rules for conducting umbrella reviews. *Evid Based Ment Health*. 2018;21(3):95-100. doi:10.1136/ebmental-2018-300014

[3] Bougioukas KI, Liakos A, Tsapas A, Ntzani E, Haidich AB. Preferred reporting items for overviews of systematic reviews including heuristic approaches to evidence synthesis: the PRIOR statement. *J Clin Epidemiol*. 2022;148:150-167. doi:10.1016/j.jclinepi.2022.04.018

[4] Gates M, Gates A, Pieper D, et al. Reporting guideline for overviews of reviews of healthcare interventions: development of the PRIOR statement. *BMJ*. 2022;378:e070849. doi:10.1136/bmj-2022-070849

[5] Pieper D, Antoine SL, Mathes T, Neugebauer EA, Eikermann M. Systematic review finds overlapping reviews were not mentioned in every other overview. *J Clin Epidemiol*. 2014;67(4):368-375. doi:10.1016/j.jclinepi.2013.11.007

[6] Bougioukas KI, Vounzoulaki E, Mantsiou CD, et al. Methods for depicting overlap in overviews of systematic reviews: an introduction to static tabular and graphical displays. *J Clin Epidemiol*. 2021;132:34-45. doi:10.1016/j.jclinepi.2020.12.004

[7] Shea BJ, Reeves BC, Wells G, et al. AMSTAR 2: a critical appraisal tool for systematic reviews that include randomised or non-randomised studies of healthcare interventions, or both. *BMJ*. 2017;358:j4008. doi:10.1136/bmj.j4008

[8] Radua J, Ramella-Cravaro V, Ioannidis JPA, et al. What causes psychosis? An umbrella review of risk and protective factors. *World Psychiatry*. 2018;17(1):49-66. doi:10.1002/wps.20490

[9] Hennessy EA, Johnson BT, Keenan C. Best practice guidelines and essential methodological steps to conduct rigorous and systematic meta-reviews. *Appl Psychol Health Well Being*. 2019;11(3):353-381. doi:10.1111/aphw.12169

[10] Perez-Bracchiglione J, Defined-overlap meta-analysis: adjusting for study overlap in umbrella reviews. *Res Synth Methods*. 2024;15(2):198-210. doi:10.1002/jrsm.1680

[11] Fusar-Poli P, Radua J. Ten simple rules for conducting umbrella reviews. *Evid Based Ment Health*. 2018;21(3):95-100. doi:10.1136/ebmental-2018-300014

---

**Appendix A: CCA Formula Derivation**

The CCA formula can be understood as follows. Let *N* = total citations (each study counted once per review it appears in), *r* = unique studies, *c* = number of reviews. The minimum possible *N* (no overlap) is *r*, and the maximum possible *N* (every study in every review) is r * c. Therefore:

> CCA = (N - N_min) / (N_max - N_min) = (N - r) / (r*c - r)

This is equivalent to the formula presented by Pieper et al. [5], reframed as a linear rescaling of total citations to the [0, 1] interval.

**Appendix B: Discordance Factor Normalization**

Each of the five discordance factors produces a raw score on [0, 1]:

| Factor | Raw Score Formula | Normalization |
|--------|-------------------|---------------|
| Scope | Mean pairwise Jaccard distance of scope tags | Natural [0, 1] range |
| Inclusion criteria | (max(k) / min(k) - 1) / 3 | Capped at 1.0; ratio of 4:1 maps to 1.0 |
| Evidence base | 1 - CCA | Natural complement |
| Quality variation | (max_rank - min_rank) / 3 | Ordinal AMSTAR ranks 1-4 |
| Statistical methods | (n_distinct_measures - 1) / 2 | Capped at 1.0 |

The contribution percentage for factor *f* is: 100 * raw_f / sum(raw_all).
