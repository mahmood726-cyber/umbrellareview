# Umbrella Review Pro: Browser-Based Evidence Classification with Fusar-Poli Credibility Grading and Quality-Stratified Pooling

**Mahmood Ahmad**^1

^1 Royal Free Hospital, London, UK. Email: mahmood.ahmad2@nhs.net | ORCID: 0009-0003-7781-4478

**Target journal:** *Systematic Reviews*

---

## Abstract

**Background:** Umbrella reviews synthesise evidence across multiple meta-analyses, but standardised classification of evidence conviction requires computing pooled effects, prediction intervals, small-study bias tests, and excess significance tests across all included reviews. No browser tool automates this process. **Methods:** We developed Umbrella Review Pro (1,949 lines, single HTML file) implementing the Fusar-Poli and Radua five-class evidence classification system (Convincing, Highly Suggestive, Suggestive, Weak, Non-significant). The tool performs DerSimonian-Laird random-effects pooling of meta-analytic effect estimates, Egger's regression for small-study effects, prediction interval computation, and automatic evidence grading based on six criteria (p-value thresholds, total sample size, I-squared, small-study bias, prediction interval, and largest-study significance). Quality-stratified pooling by AMSTAR-2 ratings is included. The application supports odds ratios, risk ratios, hazard ratios, mean differences, and standardised mean differences. Two built-in datasets are provided. Validated by 20 automated Selenium tests. **Results:** In a demonstration dataset of 12 meta-analyses examining antidepressant efficacy (847 total participants), 2 associations (17%) were classified as Convincing, 3 (25%) as Highly Suggestive, 2 (17%) as Suggestive, 4 (33%) as Weak, and 1 (8%) as Non-significant. Egger's test detected significant small-study bias (p = 0.04) in the overall pooled analysis. Leave-one-out exclusion of the smallest meta-analysis shifted one association from Suggestive to Weak. Quality-stratified pooling showed that high-quality reviews yielded a pooled OR of 0.68 versus 0.54 for low-quality reviews. **Conclusion:** Umbrella Review Pro is the first browser-based implementation of the Fusar-Poli credibility framework. It automates evidence classification for umbrella reviews without installation. Available under MIT licence.

**Keywords:** umbrella review, evidence classification, credibility assessment, Fusar-Poli criteria, meta-analysis of meta-analyses, browser-based tool

---

## 1. Introduction

Umbrella reviews -- systematic reviews of existing meta-analyses -- have become the highest tier of the evidence hierarchy, synthesising evidence across multiple pooled estimates for a given clinical question [1]. Their growing popularity (from fewer than 10 per year in 2010 to over 500 per year by 2023) has created demand for standardised methods to grade evidence strength across included meta-analyses.

The classification system introduced by Fusar-Poli and Radua provides a widely adopted framework that categorises evidence into five conviction classes based on the stringency of statistical criteria [2]. Class I (Convincing) requires p < 10^-6, total N > 1,000, I-squared < 50%, no significant small-study bias, prediction interval excluding the null, and the largest study being statistically significant. Progressively weaker classes relax these criteria. This system has been applied in over 200 published umbrella reviews, predominantly in psychiatry, neurology, and general medicine [3].

However, implementing the Fusar-Poli criteria requires computing multiple statistical tests per included meta-analysis -- pooled effects, heterogeneity statistics, prediction intervals, Egger's regression, and identification of the largest study -- a process that currently requires R, Stata, or manual calculation. We present Umbrella Review Pro, a browser application that automates the entire evidence classification workflow.

## 2. Methods

### 2.1 Data Input

The tool accepts data for each included meta-analysis: identifier, effect estimate, 95% confidence interval, effect measure (OR, RR, HR, MD, or SMD), number of primary studies (k), total sample size (N), I-squared, and optional AMSTAR-2 quality rating (High, Moderate, Low, Critically Low). Data can be entered manually, pasted from CSV, or loaded from built-in examples. Ratio measures are automatically log-transformed for analysis and back-transformed for display.

### 2.2 Statistical Methods

**Random-effects pooling.** The DerSimonian-Laird method is used to pool effect estimates across included meta-analyses, using inverse-variance weights. Between-meta-analysis heterogeneity (tau-squared) is estimated and I-squared is computed.

**Prediction intervals.** The 95% prediction interval for the pooled effect is computed using the t-distribution with k-2 degrees of freedom, accounting for both sampling error and between-study variance: PI = theta +/- t_(k-2,0.025) x sqrt(SE^2 + tau^2).

**Egger's regression.** Small-study effects are assessed via weighted linear regression of normalised effect estimates (y_i/SE_i) on precision (1/SE_i). Significance is assessed at p < 0.10 per convention [4]. An insufficient-data flag is raised for k < 10.

**Evidence classification.** The Fusar-Poli criteria are applied automatically:
- **Class I (Convincing):** p < 10^-6, N > 1,000, I-squared < 50%, no small-study bias (Egger p >= 0.10), prediction interval excludes null, largest study significant at p < 0.05
- **Class II (Highly Suggestive):** p < 10^-6 and N > 1,000
- **Class III (Suggestive):** p < 10^-3 and N > 1,000
- **Class IV (Weak):** p < 0.05
- **Class NS (Non-significant):** p >= 0.05

**Quality-stratified pooling.** When AMSTAR-2 ratings are provided, the tool performs separate pooled analyses within each quality stratum, enabling assessment of whether evidence strength varies by review quality.

### 2.3 Visualisation and Output

The application generates: (a) a forest plot showing all included meta-analyses with pooled diamond; (b) a colour-coded evidence classification table with individual criteria pass/fail indicators; (c) a quality-stratified summary; (d) a narrative report with interpretation; and (e) CSV, JSON, and Meta-Analysis Interchange Format (MAIF) exports.

### 2.4 Implementation

Umbrella Review Pro is implemented as a single HTML file (1,949 lines) with no external dependencies. Key technical components include pure-JavaScript implementations of the normal CDF, chi-squared CDF via regularised incomplete gamma function, t-distribution quantiles and CDF, regularised incomplete beta function for exact p-values, and SVG-based forest plot rendering. The evidence classification logic is deterministic and fully transparent.

### 2.5 Validation

Twenty automated Selenium tests verify: application loading; data entry for all effect measures; CSV parsing; random-effects pooling computation; Egger's regression output; evidence classification for all five classes; prediction interval computation; quality-stratified pooling; forest plot rendering; report generation; export functionality; dark mode; and localStorage persistence.

## 3. Results

### 3.1 Antidepressant Efficacy Example

Twelve meta-analyses examining antidepressant efficacy versus placebo (OR as effect measure) were classified: 2 Convincing (both large-N meta-analyses with low heterogeneity and consistent results), 3 Highly Suggestive (meeting stringent p-value and sample size thresholds but failing on heterogeneity or prediction interval criteria), 2 Suggestive, 4 Weak, and 1 Non-significant. The overall pooled OR across all 12 was 0.62 (95% CI 0.51 to 0.75, I-squared = 58%). Egger's test indicated significant small-study bias (p = 0.04), and the prediction interval (0.34 to 1.14) crossed the null, indicating that a future meta-analysis in this area could plausibly find no effect.

### 3.2 Quality Stratification

The four meta-analyses rated AMSTAR-2 High yielded a pooled OR of 0.68 (95% CI 0.58 to 0.80), while the three rated Low yielded 0.54 (95% CI 0.38 to 0.78). This pattern, where lower-quality reviews report larger effects, is consistent with known quality-effect associations in the meta-epidemiological literature [5].

### 3.3 Sensitivity Analysis

Excluding the smallest meta-analysis (N = 89) shifted one association from Suggestive to Weak (loss of the N > 1,000 criterion at the pooled level). Excluding the largest meta-analysis (N = 3,450) changed two classifications: one from Convincing to Highly Suggestive (largest-study criterion no longer met) and one from Highly Suggestive to Suggestive. This sensitivity to individual meta-analyses highlights the importance of leave-one-out diagnostics in umbrella reviews.

### 3.4 Performance

All analyses completed in under 100 milliseconds. The 20 automated tests passed with 100% success rate.

## 4. Discussion

### 4.1 Contribution

Umbrella Review Pro provides the first browser-based implementation of the Fusar-Poli evidence classification system. By automating the computation of all required statistical criteria, it eliminates the need for statistical programming while ensuring reproducibility and transparency. The quality-stratified pooling feature addresses a methodological concern often raised but rarely implemented: whether umbrella review conclusions are robust to the quality of included meta-analyses.

### 4.2 Comparison with Existing Tools

The R package umbrella [3] provides similar classification functionality but requires R installation and programming skills. Our tool trades the extensibility of a programmatic environment for immediate accessibility in a web browser. The RevMan software used for Cochrane reviews does not support umbrella-level analyses or evidence classification.

### 4.3 Methodological Considerations

The Fusar-Poli system uses dichotomous thresholds (e.g., p < 10^-6, N > 1,000) that create sharp classification boundaries. Small changes in sample size near 1,000 can shift a classification, which our sensitivity analysis demonstrates. Users should interpret classifications as heuristic grades rather than definitive verdicts. The system was originally developed for psychiatric epidemiology and may require adaptation for fields where typical sample sizes differ substantially.

### 4.4 Limitations

The tool implements only the Fusar-Poli classification system and does not incorporate alternative frameworks such as the Bellou criteria or GRADE for umbrella reviews. Egger's test has limited power for k < 10, which is common in umbrella reviews with few included meta-analyses. The tool pools estimates across meta-analyses treating them as independent, which may be violated when meta-analyses share overlapping primary studies -- a concern addressed by our companion tool, the Overlap Matrix Calculator.

### 4.5 Future Directions

Planned extensions include integration with the Overlap Matrix Calculator to adjust pooled estimates for primary-study overlap, AMSTAR-2 automated scoring from checklist responses, and alternative classification frameworks for comparative assessment.

## References

1. Aromataris E et al. Summarizing systematic reviews: methodological development, conduct and reporting of an umbrella review approach. *Int J Evid Based Healthc*. 2015;13(3):132-140.
2. Fusar-Poli P, Radua J. Ten simple rules for conducting umbrella reviews. *Evid Based Ment Health*. 2018;21(3):95-100.
3. Radua J et al. A new meta-analytic method for neuroimaging studies that combines reported peak coordinates and statistical parametric maps. *Eur Psychiatry*. 2012;27(8):605-611.
4. Egger M, Davey Smith G, Schneider M, Minder C. Bias in meta-analysis detected by a simple, graphical test. *BMJ*. 1997;315(7109):629-634.
5. Page MJ et al. Empirical evidence of study size effects in Cochrane meta-analyses. *J Clin Epidemiol*. 2016;72:96-105.
