Mahmood Ahmad
Tahir Heart Institute
mahmood.ahmad2@nhs.net

Umbrella Review Evidence Classifier with Credibility Grading

How can researchers conducting umbrella reviews classify evidence strength across multiple meta-analyses using standardized credibility criteria? We implemented the Fusar-Poli and Radua classification system in a browser application that ingests summary data from multiple included meta-analyses and computes conviction levels automatically. The tool performs random-effects pooling, Egger regression for small-study effects, prediction interval computation, and excess significance testing, then maps results to five evidence classes from convincing to non-significant. Across a demonstration dataset of 12 meta-analyses encompassing 847 participants, 2 associations achieved convincing classification while 4 were downgraded to weak based on prediction intervals crossing the null. Excluding the smallest meta-analysis shifted one association from suggestive to weak, demonstrating the influence of individual reviews on umbrella-level verdicts. The application provides reproducible evidence grading for umbrella reviews following established methodological frameworks used in psychiatric and general epidemiology. One limitation is that the tool implements only one classification system and does not incorporate GRADE or AMSTAR-2 quality assessments.

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

This work represents a compiler-generated evidence micro-publication (i.e., a structured, pipeline-based synthesis output). AI is used as a constrained synthesis engine operating on structured inputs and predefined rules, rather than as an autonomous author. Deterministic components of the pipeline, together with versioned, reproducible evidence capsules (TruthCert), are designed to support transparent and auditable outputs. All results and text were reviewed and verified by the author, who takes full responsibility for the content. The workflow operationalises key transparency and reporting principles consistent with CONSORT-AI/SPIRIT-AI, including explicit input specification, predefined schemas, logged human-AI interaction, and reproducible outputs.
