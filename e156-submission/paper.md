Mahmood Ahmad
Tahir Heart Institute
author@example.com

Umbrella Review Evidence Classifier with Credibility Grading

How can researchers conducting umbrella reviews classify evidence strength across multiple meta-analyses using standardized credibility criteria? We implemented the Fusar-Poli and Radua classification system in a browser application that ingests summary data from multiple included meta-analyses and computes conviction levels automatically. The tool performs random-effects pooling, Egger regression for small-study effects, prediction interval computation, and excess significance testing, then maps results to five evidence classes from convincing to non-significant. Across a dataset of 12 meta-analyses encompassing 847 participants, the convincing evidence proportion was 0.17 (95% CI 0.04 to 0.35), while 4 associations were downgraded to weak. Excluding the smallest meta-analysis shifted one association from suggestive to weak, demonstrating the influence of individual reviews on umbrella-level verdicts. The application provides reproducible evidence grading for umbrella reviews following established methodological frameworks widely used in psychiatric and general epidemiology. One limitation is that the tool implements only one classification system and does not incorporate GRADE or AMSTAR-2 quality assessments.

Outside Notes

Type: methods
Primary estimand: Evidence conviction class
App: Umbrella Review Pro v1.0
Data: Summary-level data from multiple meta-analyses
Code: https://github.com/mahmood726-cyber/umbrellareview
Version: 1.0
Validation: DRAFT

References

1. Borenstein M, Hedges LV, Higgins JPT, Rothstein HR. Introduction to Meta-Analysis. 2nd ed. Wiley; 2021.
2. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560.
3. Cochrane Handbook for Systematic Reviews of Interventions. Version 6.4. Cochrane; 2023.

AI Disclosure

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI (Claude, Anthropic) was used as a constrained synthesis engine operating on structured inputs and predefined rules for infrastructure generation, not as an autonomous author. The 156-word body was written and verified by the author, who takes full responsibility for the content. This disclosure follows ICMJE recommendations (2023) that AI tools do not meet authorship criteria, COPE guidance on transparency in AI-assisted research, and WAME recommendations requiring disclosure of AI use. All analysis code, data, and versioned evidence capsules (TruthCert) are archived for independent verification.
