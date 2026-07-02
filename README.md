# AICD–ammonia TCGA-COADREAD analysis

Reproducible R code for the in-cohort transcriptional analysis (TCGA colon and rectal
adenocarcinoma, COADREAD; *n* = 598) reported in the review:

> Yang F, Lv C. *The Ammonia Metabolic Immune Checkpoint: Spatially Stratified Cell Death,
> Regulatory T Cell Empowerment, and Resistance to Immune Checkpoint Inhibition in Solid Tumours.*

This exploratory analysis is provided for full transparency and reproducibility. It is a
hypothesis-illustrating, correlative analysis embedded in a narrative review; it is **not** an
independently validated biomarker study (see *Scope and limitations* below).

---

## Repository contents

```
.
├── AICD_TCGA_analysis.R        # single self-contained analysis script
├── AICD_genes_tpm.csv          # input: 8 genes (rows) x 598 patients (TPM)
├── TCGA-COADREAD_cli_surv_binomial_patient.csv   # input: patient clinical + survival
├── outputs/                    # created on run (figures + result tables)
└── README.md
```

The two input `.csv` files are derived tables from open-access TCGA data (see *Data source*),
included here so the analysis can be reproduced with a single command.

---

## Data source

Transcriptomic (TPM) and clinical/survival data are derived from **TCGA-COADREAD**, obtained
from the UCSC Xena platform (<https://xenabrowser.net>) and the NCI Genomic Data Commons
(<https://portal.gdc.cancer.gov>). TCGA data are open-access; the derived tables in this
repository contain gene-level expression for eight genes and de-identified clinical variables
only.

- `AICD_genes_tpm.csv` — log-linear TPM for FOXP3, GLUD1, GLUL, SMS, CD8A, FOXM1, HADH, ATP1A1.
- `TCGA-COADREAD_cli_surv_binomial_patient.csv` — overall survival (`os`, `os_time`), age,
  gender, `pathologic_stage` (Stage I/II vs III/IV), `M_Stage` (M0/M1), and
  `mismatch_repair_proteins` (pMMR/dMMR).

---

## Requirements

- **R ≥ 4.3**
- CRAN packages: `tidyverse`, `data.table`, `survival`, `survminer`, `ggpubr`, `pheatmap`,
  `RColorBrewer`, `broom`, `scales`

The script installs any missing packages automatically. Exact versions used for the published
analysis are recorded in `outputs/sessionInfo.txt` after a run.

> Note: the plotting theme requests the **Arial** font. On systems without Arial (common on
> Linux), figures still render but a font warning may appear; this does not affect any numerical
> result. To silence it, replace `base_family = "Arial"` with `base_family = ""`.

---

## How to run

```bash
# from the repository root (with the two .csv files present)
Rscript AICD_TCGA_analysis.R
```

All figures (PDF + 300-dpi PNG) and result tables (CSV) are written to `./outputs/`.
A fixed seed (`set.seed(123)`) is set; the analysis is deterministic.

---

## Analysis modules and mapping to the manuscript

| Module | Analysis | Manuscript item |
|--------|----------|-----------------|
| 1 | Expression landscape; Spearman correlation among the 8 genes | Figure 8A |
| 2 | Clinical-subgroup differential expression (Wilcoxon; BH-adjusted) across MMR / M-stage / pathologic stage | Figure 8B; Table 6 (expression) |
| 3 | Per-gene Kaplan–Meier (optimal cut-point) + univariate continuous Cox | Table 6 (survival) |
| 4 | Composite scores (4-gene Treg-axis z-score; 8-gene PCA signature, PC1) + **multivariable Cox** | Figure 8D (PC1 KM); Table 7 |
| 5 | log2(CD8A/FOXP3) effector-to-Treg ratio across clinical subgroups | Figure 8C |

**Key reproduced values:** GLUD1–SMS Spearman *ρ* = 0.39 (*p* = 1.2×10⁻²²); AICD-PC1 explains
27.0% of variance; PC1 Kaplan–Meier log-rank *p* = 1.9×10⁻⁴; log2(CD8A/FOXP3) ratio reduced in
dMMR vs pMMR (*p* = 0.007). In the **multivariable Cox model (*n* = 436** after excluding cases
with missing MMR or stage), the AICD-PC1 signature is **not** independently prognostic
(HR = 0.89, 95% CI 0.74–1.06, *p* = 0.19); advanced stage (HR = 4.65) and pMMR status
(HR = 5.15) dominate.

---

## Methodological notes (for transparency)

- Expression values are analysed as log₂(TPM + 1).
- Multiple comparisons within each subgroup contrast are controlled by the Benjamini–Hochberg
  method.
- Kaplan–Meier dichotomisation uses maximally selected rank statistics
  (`survminer::surv_cutpoint`, `minprop = 0.25`). Such data-driven cut-points are optimistic and
  require external validation; univariate Cox models are therefore reported on the **continuous**
  variables in parallel.
- The multivariable model is restricted to patients with complete stage and mismatch-repair
  status (*n* = 436).

## Scope and limitations

This is an exploratory, correlative analysis. Co-expression does not establish pathway activity,
cell-type specificity, ammonia utilisation, or causality. Optimal cut-points are optimistic, and
the findings have **not** been replicated in an independent cohort. The composite AICD signature
is not independently prognostic after adjustment. Results should be read as hypothesis-generating
support for the conceptual framework developed in the review, not as a validated biomarker.

---

## Citation

If you use this code, please cite the manuscript above and this repository.

## License

Released under the MIT License (see `LICENSE`).
