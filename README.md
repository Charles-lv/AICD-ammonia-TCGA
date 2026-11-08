AICD–ammonia TCGA-COADREAD analysis
Code and derived input tables for the exploratory bulk transcriptomic analysis of TCGA colon and rectal adenocarcinoma (COADREAD; n = 598) reported in Supplementary S1 of the accompanying review manuscript (under review; the full citation will be added on publication).

Main result
In a multivariable Cox model for overall survival adjusted for age, sex, pathological stage and mismatch repair status (complete cases, n = 436), the eight-gene composite (AICD-PC1) had no independent prognostic value (HR 0.89, 95% CI 0.74–1.06; Wald P = 0.19). An association seen without adjustment (Kaplan–Meier split at a data-driven optimal cut-point) did not persist after adjustment. The result is reported as negative.

Scope and limitations
This is an exploratory, correlative analysis of bulk tumour transcripts. It cannot test a claim about where ammonia acts or on which cells: bulk expression averages over all cell types present, transcript abundance is not metabolic flux, and ammonia has no direct transcriptional readout. FOXP3, one of the eight genes, is expressed both by suppressive FOXP3^hi regulatory T cells and by non-suppressive FOXP3^lo T cells, which carry opposite prognostic associations in colorectal cancer (Saito et al., Nat Med 2016;22:679–684); a bulk FOXP3 transcript pools both. The analysis has not been replicated in an independent cohort and is not a biomarker study.

Repository contents
.
├── AICD_TCGA_analysis.R                          # single self-contained analysis script
├── AICD_genes_tpm.csv                            # input: 8 genes (rows) x 598 patients (TPM)
├── TCGA-COADREAD_cli_surv_binomial_patient.csv   # input: patient clinical + survival
├── outputs/                                      # created on run (figures + result tables)
├── LICENSE
└── README.md
Data source
Transcriptomic (TPM) and clinical/survival data for TCGA-COADREAD were obtained from the UCSC Xena platform (https://xenabrowser.net) and the NCI Genomic Data Commons (https://portal.gdc.cancer.gov). TCGA data are open access; the derived tables contain gene-level expression for eight genes and de-identified clinical variables only.

AICD_genes_tpm.csv: TPM for FOXP3, GLUD1, GLUL, SMS, CD8A, FOXM1, HADH and ATP1A1.
TCGA-COADREAD_cli_surv_binomial_patient.csv: overall survival (os, os_time), age, gender, pathologic_stage (Stage I/II vs III/IV), M_Stage and mismatch_repair_proteins (pMMR/dMMR).
Requirements and how to run
R ≥ 4.3 with CRAN packages tidyverse, data.table, survival, survminer, ggpubr, pheatmap, RColorBrewer, broom and scales (missing packages are installed by the script).

Rscript AICD_TCGA_analysis.R
All figures and result tables are written to ./outputs/, together with sessionInfo.txt recording the package versions used. A fixed seed (set.seed(123)) is set.

Analysis modules
Module	Analysis	Reported in the manuscript
1	Expression landscape; Spearman correlation among the 8 genes	No (exploratory)
2	Clinical-subgroup Wilcoxon comparisons (BH-adjusted)	No (exploratory)
3	Per-gene Kaplan–Meier (optimal cut-point) and univariate continuous Cox	No (exploratory)
4	Composite scores (4-gene z-score; 8-gene PCA, PC1) and multivariable Cox	Yes, AICD-PC1 multivariable result only (Supplementary S1)
5	log2(CD8A/FOXP3) ratio across clinical subgroups	No (exploratory)
Methodological notes
Expression values are analysed as log2(TPM + 1).
AICD-PC1 is the first principal component of the eight z-scored genes, sign-aligned so that higher values track higher overall expression, and enters Cox models as a continuous variable.
Kaplan–Meier dichotomisation uses maximally selected rank statistics (survminer::surv_cutpoint, minprop = 0.25); such data-driven cut-points are optimistic.
The multivariable model is restricted to patients with complete age, stage and mismatch repair data (n = 436).
License
Released under the MIT License (see LICENSE).
