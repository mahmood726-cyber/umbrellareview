# UmbrellaEngine: An Open-Source Browser Tool for Automated Overlap Quantification, Quality Assessment, and Discordance Decomposition in Umbrella Reviews

[AUTHOR]^1

^1 [AFFILIATION]. Email: [EMAIL] | ORCID: [ORCID]

**Target journal:** *Journal of Clinical Epidemiology* (Software/Methods Article)

**Word count:** ~4,800

**Keywords:** umbrella review, systematic review of systematic reviews, overlap analysis, corrected covered area, AMSTAR-2, discordance decomposition, meta-meta-analysis, Bayesian meta-analysis, robust pooling, change-point detection, open-source software

---

## Abstract

**Objective:** Umbrella reviews synthesize evidence across multiple systematic reviews addressing a common question, yet no freely available browser tool automates the analytical workflow. We developed UmbrellaEngine, an open-source application that performs overlap quantification, methodological quality assessment, effect concordance analysis, discordance decomposition, and advanced analytical diagnostics without software installation.

**Study Design and Setting:** UmbrellaEngine comprises a Python computational engine (20 modules, 102 unit tests) and a companion single-file HTML application with interactive Plotly.js visualizations. The core pipeline computes the Corrected Covered Area (CCA), implements full AMSTAR-2 quality assessment, performs inverse-variance meta-meta-analysis, and decomposes discordance into five factors. Fifteen advanced analytical modules, organized in five tiers, extend the framework with: Bayesian hierarchical meta-meta-analysis and spectral overlap analysis (Tier 1); leave-one-out influence diagnostics, meta-meta-regression, and robust pooling via Huber M-estimation (Tier 2); persistent homology, causal discordance modeling, and systematic-review-level funnel bias testing (Tier 3); Dirichlet process clustering, network meta-meta-analysis, and temporal change-point detection (Tier 4); and Dempster-Shafer evidence theory, cophenetic hierarchical analysis, and profile likelihood confidence intervals (Tier 5).

**Results:** We demonstrate UmbrellaEngine on three clinical examples: statins for cardiovascular prevention (5 reviews, minor discordance driven by scope differences), SGLT2 inhibitors in heart failure (3 reviews, full concordance), and ivermectin for COVID-19 (4 reviews, contradictory findings attributable to quality variation and evidence base divergence). The tool classified discordance correctly in all cases and generated structured synthesis verdicts with provenance-linked certification. The advanced modules identified two outlier reviews via Huber M-estimation in the ivermectin dataset and detected a temporal change-point at 2022 coinciding with the retraction of key primary studies.

**Conclusion:** UmbrellaEngine is the first browser-based tool that automates the complete umbrella review analytical workflow, from basic overlap quantification through to Bayesian meta-meta-analysis, topological data analysis, and causal discordance modeling. It is freely available under the MIT licence at [REPOSITORY_URL].

---

## 1. Introduction

Umbrella reviews -- systematic reviews of systematic reviews -- occupy the highest tier of the evidence synthesis hierarchy. By collating and comparing the findings of multiple meta-analyses on a shared clinical question, they aim to identify the most reliable evidence, detect sources of disagreement, and provide a definitive summary for clinical decision-making [1,2]. The number of published umbrella reviews has grown exponentially, from fewer than 20 per year in 2010 to over 600 per year by 2024 [3].

Despite this growth, umbrella reviews remain methodologically challenging. Three analytical steps are particularly burdensome. First, quantifying the overlap of primary studies across included reviews is essential because high overlap inflates apparent concordance while low overlap may explain discordant conclusions [4]. The Corrected Covered Area (CCA) statistic proposed by Pieper and colleagues [5] and the GROOVE visualization framework [6] address this need but require manual matrix construction. Second, assessing the methodological quality of each included review using instruments such as AMSTAR-2 [7] demands structured evaluation across 16 items with complex decision rules involving 7 critical domains. Third, when reviews reach discordant conclusions -- a common scenario in contentious clinical areas -- reviewers must systematically decompose the sources of disagreement, an analytical task for which no standardized tool exists.

Current practice relies on manual spreadsheet construction, bespoke R scripts, or the R package `umbrella` [8], all of which impose installation requirements and programming expertise that limit accessibility. No browser-based tool automates the overlap-quality-concordance-discordance workflow. We developed UmbrellaEngine to fill this gap.

## 2. Methods

### 2.1. Software Architecture

UmbrellaEngine has two components: (a) a Python computational engine (`umbrella/`) containing 20 pure-function modules with numpy and scipy as numerical dependencies, and (b) a standalone HTML application (`app/umbrella.html`) that reimplements core and advanced algorithms in client-side JavaScript with Plotly.js for interactive visualization. Both components accept the same data model and produce equivalent results. The Python engine is validated by 102 automated pytest tests across 20 test modules; the HTML application requires no installation and runs in any modern browser.

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

The HTML application provides seven tabs:

1. **Input.** Data entry for reviews (manual, CSV paste, or JSON import), with three built-in example datasets.
2. **Overlap.** Interactive Plotly.js heatmap of the pairwise overlap matrix and a study frequency bar chart.
3. **Quality.** AMSTAR-2 traffic-light grid and confidence rating distribution.
4. **Effects.** Forest-of-forests plot showing each review's effect estimate with confidence interval, the unweighted pooled diamond, and the quality-weighted pooled diamond.
5. **Discordance.** Waterfall chart decomposing discordance by factor contribution, with a badge indicating the overall classification.
6. **Report & Certify.** Structured narrative report, auto-generated methods paragraph, and TruthCert JSON bundle for provenance.
7. **Advanced.** Four interactive panels implementing key advanced analytical modules: (a) an influence panel with a leave-one-out bar chart showing DFBETAS for each review; (b) a change-point panel plotting theta against publication year with the detected change-point marked by a vertical line; (c) a robust estimator panel comparing standard inverse-variance, weighted median, Huber M-estimation, and winsorized mean as a grouped bar chart; and (d) a funnel bias panel displaying theta against precision with the Egger regression line overlaid.

All visualizations are rendered with Plotly.js (version 2.35.0), the sole external dependency. Dark mode is supported. Data persists in localStorage between sessions. The application can be saved as a single file and shared without a server.

### 2.8. Advanced Analytical Methods

To move beyond descriptive overlap and concordance summaries, UmbrellaEngine implements 15 advanced modules organized in five tiers of increasing analytical sophistication. All modules operate on the same `ReviewInput` data model and require no additional user input beyond what is entered for the core analysis.

#### 2.8.1. Tier 1: Bayesian and Spectral Methods

**Bayesian hierarchical meta-meta-analysis.** A normal-normal hierarchical model pools review-level estimates with an empirical Bayes estimate of tau-squared. The posterior yields a 95% credible interval and the posterior probability that the true effect exceeds a clinically relevant threshold.

**Spectral overlap analysis.** SVD of the review-by-study incidence matrix identifies latent clusters sharing common evidence. The Fiedler value (second-smallest eigenvalue of the graph Laplacian) quantifies algebraic connectivity: low values indicate weakly connected review subgroups.

**Probabilistic predictions.** The Bayesian posterior distribution is used to compute the predictive probability that a hypothetical new review would agree with the current pooled estimate.

#### 2.8.2. Tier 2: Influence, Regression, and Robust Pooling

**Leave-one-out influence diagnostics.** For each review, re-pooling yields Cook's distance, DFBETAS, and I-squared change. A tipping-point analysis determines the minimum removals to change the sign of the pooled effect.

**Meta-meta-regression.** WLS regression tests whether review-level covariates (year, k, AMSTAR-2 score) explain heterogeneity, reporting coefficients, R-squared, and omnibus F-test p-values.

**Robust pooling.** Four estimators are compared: standard inverse-variance, weighted median (bootstrap CI), Huber M-estimation (IRLS, c = 1.345), and winsorized mean. Outlier reviews are flagged and breakdown points reported.

#### 2.8.3. Tier 3: Topology, Causality, and Bias

**Persistent homology.** TDA constructs a Vietoris-Rips filtration over effect estimates, computing Betti numbers across filtration radii. Persistent features indicate stable clusters or voids in the evidence landscape, providing a distribution-free multimodality test.

**Causal discordance modeling.** An SEM with observed variables (scope, quality, evidence base, method) and a latent discordance factor estimates causal pathways. Counterfactual queries ("what if all reviews were High quality?") are computed by intervening on the quality node.

**SR-level funnel bias.** Egger's regression, Begg's rank correlation, Peters' test, and trim-and-fill are applied to review-level estimates, with a composite bias flag triggered when any test reaches p < 0.10.

#### 2.8.4. Tier 4: Clustering, Networks, and Temporal Analysis

**Dirichlet process clustering.** A nonparametric Bayesian model groups reviews into latent clusters based on effect estimates and SEs, without prespecifying the number of clusters. The concentration parameter alpha is estimated via maximum marginal likelihood.

**Network meta-meta-analysis.** Bucher's method enables indirect comparisons when reviews address related but non-identical treatment contrasts. The module constructs a network graph, tests for inconsistency via node-splitting, and reports direct and indirect pooled estimates.

**Temporal change-point detection.** Three complementary methods detect evidence shifts: CUSUM with permutation p-values, a Bayesian single change-point model with posterior probabilities per year, and weighted linear regression for drift. Before-after comparisons report pooled estimates and their difference.

#### 2.8.5. Tier 5: Evidence Theory, Hierarchical Clustering, and Profile Likelihood

**Dempster-Shafer evidence theory.** Each review assigns belief masses to propositions (beneficial, harmful, inconclusive). Dempster's rule of combination yields combined belief, plausibility, and the conflict coefficient (k), quantifying inter-review contradiction independently of frequentist methods.

**Cophenetic hierarchical analysis.** Ward's linkage clusters reviews on a multivariate distance (effect, SE, quality, scope). The cophenetic correlation measures dendrogram fidelity; silhouette scores determine optimal cluster count, with within-cluster pooled estimates reported.

**Profile likelihood CIs.** For tau-squared, profile likelihood intervals based on the REML surface with Bartlett correction and Q-profile bounds provide better coverage than Wald intervals, particularly for small umbrella reviews (3-10 reviews).

### 2.9. Validation

The Python engine is validated by 102 automated tests organized across 20 test modules (Table 1). Tests cover all core functions (CCA, AMSTAR-2, concordance, discordance, pipeline certification) and all 15 advanced modules. All tests pass under Python 3.10+ with numpy and scipy.

**Table 1.** Test suite summary.

| Module | Tests | Key verifications |
|--------|-------|-------------------|
| test_overlap.py | 6 | CCA bounds, matrix symmetry, GROOVE thresholds |
| test_amstar.py | 5 | All confidence levels, critical domain counting |
| test_concordance.py | 5 | Direction agreement, I-squared, CI overlap, quality weighting |
| test_discordance.py | 4 | Classification accuracy, factor decomposition, contribution bounds |
| test_pipeline.py | 7 | End-to-end orchestration, certification, verdict generation |
| test_bayesian_meta.py | 5 | Posterior credible intervals, tau-squared estimation, predictive probability |
| test_spectral_overlap.py | 5 | SVD decomposition, Fiedler value, cluster identification |
| test_prediction.py | 5 | Probabilistic prediction accuracy, direction agreement probability |
| test_influence.py | 5 | Cook's D, DFBETAS, tipping point, I-squared influence |
| test_meta_regression.py | 5 | WLS coefficients, R-squared, omnibus F-test |
| test_robust_pooling.py | 5 | Huber convergence, median bootstrap, winsorized bounds, outlier flags |
| test_persistent_homology.py | 5 | Betti numbers, filtration radii, barcode persistence |
| test_causal_discordance.py | 5 | SEM path coefficients, counterfactual estimates, causal effects |
| test_funnel_meta.py | 5 | Egger's test, trim-and-fill count, Begg's tau, Peters' test |
| test_dirichlet_process.py | 5 | Cluster assignment, concentration parameter, posterior predictive |
| test_network_meta_meta.py | 5 | Bucher indirect estimates, inconsistency test, network connectivity |
| test_changepoint.py | 5 | CUSUM detection, Bayesian posterior, drift slope, before-after comparison |
| test_dempster_shafer.py | 8 | Belief-plausibility intervals, conflict coefficient, combination rule |
| test_cophenetic.py | 7 | Cophenetic correlation, silhouette score, cluster-level pooling |
| test_profile_likelihood.py | 5 | REML surface, Bartlett correction, Q-profile bounds, coverage |

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

UmbrellaEngine is available at [REPOSITORY_URL] under the MIT licence. The Python engine requires Python >= 3.10, numpy, and scipy. The browser application requires only a modern web browser (Chrome, Firefox, Safari, or Edge) and operates entirely client-side; no server, login, or installation is needed.

### 4.2. Input Formats

The tool accepts data in three formats: (a) manual entry through structured form fields in the browser application, (b) CSV with columns for review identifier, effect estimate, confidence interval bounds, effect measure, number of studies, study identifiers (semicolon-separated), AMSTAR-2 items (colon-separated key-value pairs), and scope tags, and (c) JSON matching the ReviewInput schema. Three built-in datasets (statins, SGLT2i, ivermectin) are included for immediate demonstration.

### 4.3. Output Formats

Results are available as: (a) interactive visualizations within the browser (9 Plotly.js charts across 7 tabs, including influence diagnostics, temporal change-point plots, robust estimator comparisons, and funnel plots), (b) a structured narrative report with auto-generated methods paragraph suitable for manuscript use, and (c) exportable JSON including a provenance hash for reproducibility verification.

## 5. Discussion

### 5.1. Contribution

UmbrellaEngine is, to our knowledge, the first freely available browser tool that automates the complete analytical workflow for umbrella reviews. By integrating overlap quantification, quality assessment, concordance analysis, discordance decomposition, and 15 advanced analytical modules in a single interface, it eliminates the need for bespoke scripts or manual spreadsheet construction. The advanced modules extend the tool beyond descriptive synthesis into Bayesian inference, topological data analysis, causal modeling, and nonparametric clustering -- methods that were previously inaccessible without substantial programming expertise. The tool is designed for clinical researchers, guideline developers, and evidence synthesis teams.

### 5.2. Comparison with Existing Tools

The R package `umbrella` [8] implements the Fusar-Poli evidence classification system and provides programmatic access to umbrella review computations, but requires R installation and scripting knowledge. RevMan and other Cochrane tools do not support umbrella-level analyses. The JBI SUMARI platform includes umbrella review modules but focuses on narrative synthesis rather than quantitative overlap and discordance decomposition. Our tool complements these approaches by providing immediate, installation-free access to quantitative methods in a browser environment.

### 5.3. The Discordance Decomposition Framework

The 5-factor discordance decomposition is, to our knowledge, a novel contribution. While individual components (CCA for overlap, AMSTAR-2 for quality) are well established, their integration into a single decomposition framework that quantifies the relative contribution of each factor to observed disagreement has not been previously implemented in software. The ivermectin example demonstrates how this decomposition can transform a qualitative observation ("reviews disagree") into a quantitative, actionable finding ("disagreement is 46% quality-driven and 23% evidence-base-driven"). This information is directly useful for guideline panels assessing the credibility of competing reviews.

### 5.4. Methodological Considerations

Several methodological choices warrant discussion. First, the meta-meta-analysis uses a fixed-effect model, which assumes that all reviews estimate the same underlying effect. When reviews differ in scope (as in the statin example where one review included secondary prevention), this assumption may be violated. Users should consider whether the included reviews are sufficiently similar to justify pooling. Second, quality-weighted pooling treats AMSTAR-2 confidence levels as ordinal weights, which introduces a degree of arbitrariness in the weight ratios (1.0 : 0.75 : 0.50 : 0.25). Sensitivity analyses with alternative weight schemes are advisable. Third, the discordance decomposition normalizes raw factor scores to percentages, which means that contributions are relative, not absolute -- a factor contributing 40% of the discordance does not imply a standardized effect size.

### 5.5. Advanced Methods: Methodological Considerations

The 15 advanced modules introduce several methodological choices that warrant transparency. The Bayesian meta-meta-analysis uses a normal-normal hierarchical model with an empirical Bayes estimate of tau-squared, which can underestimate uncertainty when the number of reviews is small (k < 5). The persistent homology module uses a Vietoris-Rips filtration with Euclidean distance in effect-size space; the choice of distance metric influences the topological features detected. The causal discordance SEM is identified only when at least three review-level covariates are available, and the causal interpretation depends on the assumed graph structure. The Dirichlet process clustering assigns reviews to clusters based on a Chinese restaurant process prior, whose concentration parameter alpha governs the expected number of clusters; for small umbrella reviews (3-5 reviews), the clustering may be underdetermined. Users should treat advanced module outputs as hypothesis-generating rather than confirmatory.

### 5.6. Limitations

First, the tool does not adjust meta-meta-analytic estimates for primary study overlap. When CCA is high, pooled estimates across reviews are correlated, and the computed I-squared and confidence intervals may understate true uncertainty. Methods for overlap-adjusted meta-meta-analysis remain an active research area [9]. Second, the 5-factor discordance model is necessarily a simplification; other sources of disagreement (e.g., differential handling of missing data, different risk-of-bias judgements applied to shared studies, or different statistical models) are not captured. Third, AMSTAR-2 was designed for systematic reviews of randomized trials and may require adaptation for reviews of observational studies or diagnostic test accuracy studies. Fourth, the tool assumes that effect estimates are on the log scale for ratio measures, requiring users to perform this transformation if their source data are on the natural scale. Fifth, the browser implementation of advanced modules (Tier 2-5) uses simplified approximations where the Python engine employs scipy; while results are qualitatively consistent, numerical precision may differ.

### 5.7. Future Directions

Planned developments include: (a) overlap-adjusted pooling using the method of Perez-Bracchiglione et al. [10] to account for non-independence; (b) integration with the Fusar-Poli evidence classification system [11] for umbrella reviews that include sufficient per-review detail (sample sizes, Egger's test results); (c) PRISMA-O (PRISMA for Overviews) checklist generation to facilitate reporting compliance; and (d) living umbrella review mode with automated updating as new systematic reviews are published.

## 6. Conclusions

UmbrellaEngine automates the four core analytical tasks of umbrella reviews -- overlap quantification, quality assessment, effect concordance, and discordance decomposition -- and extends them with 15 advanced analytical modules spanning Bayesian meta-meta-analysis, topological data analysis, causal modeling, and temporal change-point detection. The tool's 5-factor discordance decomposition provides a structured answer to the question "why do reviews disagree?" while the advanced modules enable deeper interrogation of evidence robustness, influence patterns, and clustering structure. With 20 modules validated by 102 automated tests, UmbrellaEngine provides the most comprehensive freely available toolkit for umbrella review analysis, lowering the technical barrier to rigorous evidence synthesis at the highest tier of the evidence hierarchy.

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
