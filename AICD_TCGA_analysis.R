# ============================================================================
# AICD pan-cancer review — TCGA-COADREAD bioinformatic validation (R version)
# Manuscript: Yang F, Lv C. "The Ammonia Metabolic Immune Checkpoint..."
# Author of analysis script: Charles (Shengjing Hospital, China Medical University)
# Date: 2026-05
# R version baseline: R v4.3+; sessionInfo() saved at the end.
# ============================================================================
# Inputs (working directory):
#   AICD_genes_tpm.csv                            — 8 genes (rows) x 598 patients (TPM)
#   TCGA-COADREAD_cli_surv_binomial_patient.csv   — 598 patients clinical+survival
#   TCGA-COADREAD_group.csv                       — full-sample-ID T/N labels (701 rows)
# Outputs: ./outputs/  (figures PDF/PNG + CSV tables)
# ============================================================================

# ---------- Packages ----------
need <- c("tidyverse", "data.table", "survival", "survminer", "ggpubr",
          "pheatmap", "RColorBrewer", "broom", "scales")
ins  <- need[!need %in% installed.packages()[, "Package"]]
if (length(ins)) install.packages(ins, dependencies = TRUE)
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(survival); library(survminer)
  library(ggpubr); library(pheatmap); library(RColorBrewer); library(broom); library(scales)
})

set.seed(123)
OUT <- "outputs"; dir.create(OUT, showWarnings = FALSE)
theme_publication <- function() {
  theme_classic(base_family = "Arial", base_size = 11) +
    theme(plot.title = element_text(face = "plain", size = 11),
          axis.title = element_text(size = 10),
          legend.position = "right", legend.title = element_blank())
}

# ---------- Load data ----------
expr <- fread("AICD_genes_tpm.csv") |> as.data.frame()
rownames(expr) <- expr[[1]]; expr[[1]] <- NULL
expr <- expr[, colnames(expr) != "" & !is.na(colnames(expr))]
clin <- fread("TCGA-COADREAD_cli_surv_binomial_patient.csv") |> as.data.frame()
rownames(clin) <- clin[[1]]; clin[[1]] <- NULL

# log2(TPM + 1)
expr_log <- log2(as.matrix(expr) + 1)
gene_order <- c("FOXP3", "GLUD1", "GLUL", "SMS", "CD8A", "FOXM1", "HADH", "ATP1A1")
expr_log <- expr_log[gene_order, , drop = FALSE]
X <- t(expr_log) |> as.data.frame()    # samples x genes (598 x 8)

# Align with clinical
common <- intersect(rownames(X), rownames(clin))
X    <- X[common, ]
clin <- clin[common, ]
df   <- cbind(X, clin)
df$os_event   <- as.integer(df$os == "Dead")
df$os_time_mo <- df$os_time / 30.44
fwrite(df, file.path(OUT, "merged_data.csv"), row.names = TRUE)
cat(sprintf("Final cohort: n = %d\n", nrow(df)))

# ---------- Cohort summary (Table 1) ----------
sumrow <- function(x, cont = TRUE) {
  x <- na.omit(x)
  if (cont) sprintf("%.1f (%.1f-%.1f)", median(x), quantile(x, .25), quantile(x, .75))
  else paste(sprintf("%s %d (%.1f%%)", names(table(x)), table(x), prop.table(table(x))*100),
             collapse = ", ")
}
tab1 <- tibble::tribble(
  ~Variable, ~Value,
  "N",                                            as.character(nrow(df)),
  "Age, median (IQR), y",                         sumrow(df$age),
  "Sex",                                          sumrow(df$gender, FALSE),
  "Pathologic stage",                             sumrow(df$pathologic_stage, FALSE),
  "M stage",                                      sumrow(df$M_Stage, FALSE),
  "MMR status",                                   sumrow(df$mismatch_repair_proteins, FALSE),
  "OS events (Dead/Total)",                       sprintf("%d/%d (%.1f%%)",
                                                          sum(df$os_event), nrow(df),
                                                          sum(df$os_event)/nrow(df)*100),
  "Follow-up months, median (IQR)",               sumrow(df$os_time_mo)
)
fwrite(tab1, file.path(OUT, "Table1_cohort.csv"))
print(tab1)

# ============================================================================
# MODULE 1 — Expression landscape & Spearman correlation
# ============================================================================
# 1A: violin of 8 genes
g1A <- X |> rownames_to_column("Sample") |> pivot_longer(-Sample, names_to = "Gene", values_to = "log2TPM") |>
  mutate(Gene = factor(Gene, levels = gene_order)) |>
  ggplot(aes(Gene, log2TPM, fill = Gene)) +
  geom_violin(trim = FALSE, scale = "width", color = "black", linewidth = .3) +
  geom_boxplot(width = .12, outlier.size = .3, fill = "white") +
  scale_fill_brewer(palette = "Set2", guide = "none") +
  labs(x = NULL, y = expression(log[2]*"(TPM+1)"),
       title = paste0("AICD core genes in TCGA-COADREAD (n = ", nrow(df), ")")) +
  theme_publication()
ggsave(file.path(OUT, "Fig1A_expression_violin.pdf"), g1A, width = 7, height = 4)
ggsave(file.path(OUT, "Fig1A_expression_violin.png"), g1A, width = 7, height = 4, dpi = 300)

# 1B: Spearman correlation matrix
corr <- cor(X, method = "spearman")
pmat <- matrix(1, ncol(X), ncol(X), dimnames = list(colnames(X), colnames(X)))
for (i in seq_len(ncol(X))) for (j in seq_len(ncol(X))) if (i != j)
  pmat[i, j] <- suppressWarnings(cor.test(X[[i]], X[[j]], method = "spearman")$p.value)
fwrite(as.data.frame(corr) |> rownames_to_column("Gene"),
       file.path(OUT, "M1_spearman_rho.csv"))
fwrite(as.data.frame(pmat) |> rownames_to_column("Gene"),
       file.path(OUT, "M1_spearman_pval.csv"))

pdf(file.path(OUT, "Fig1B_correlation_heatmap.pdf"), width = 6, height = 5)
pheatmap(corr, color = colorRampPalette(rev(brewer.pal(7, "RdBu")))(100),
         breaks = seq(-1, 1, length.out = 101),
         display_numbers = round(corr, 2), number_color = "black",
         cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 10,
         main = "Spearman correlation among AICD core genes")
dev.off()

# ============================================================================
# MODULE 2 — Clinical-subgroup Wilcoxon test (Tumor internal DEG surrogate)
# ============================================================================
sub_test <- function(gene, col, g1, g2) {
  a <- df[[gene]][df[[col]] == g1]; a <- a[!is.na(a)]
  b <- df[[gene]][df[[col]] == g2]; b <- b[!is.na(b)]
  if (length(a) < 5 || length(b) < 5) return(NULL)
  tt <- wilcox.test(a, b, exact = FALSE)
  tibble(gene = gene, group = col, g1 = g1, n1 = length(a),
         g2 = g2, n2 = length(b),
         median_g1 = median(a), median_g2 = median(b),
         log2FC_medians = median(a) - median(b),
         W = unname(tt$statistic), p = tt$p.value)
}
combs <- list(c("mismatch_repair_proteins", "dMMR",          "pMMR"),
              c("M_Stage",                  "M1",            "M0"),
              c("pathologic_stage",         "Stage III/IV",  "Stage I/II"))
res2 <- map_dfr(gene_order, \(g) map_dfr(combs, \(co) sub_test(g, co[1], co[2], co[3])))
res2 <- res2 |> group_by(group) |> mutate(p_adj = p.adjust(p, method = "BH")) |> ungroup() |>
        mutate(sig = case_when(p_adj < .001 ~ "***", p_adj < .01 ~ "**",
                               p_adj < .05  ~ "*",   TRUE ~ "ns"))
fwrite(res2, file.path(OUT, "M2_subgroup_wilcoxon.csv"))

# Forest plot
g2 <- res2 |>
  mutate(comparison = factor(paste(g1, "vs", g2),
                             levels = c("dMMR vs pMMR", "M1 vs M0", "Stage III/IV vs Stage I/II")),
         gene = factor(gene, levels = gene_order)) |>
  ggplot(aes(log2FC_medians, gene, fill = p_adj < .05)) +
  geom_col(color = "black", linewidth = .3) +
  geom_vline(xintercept = 0, linewidth = .4) +
  geom_text(aes(label = sig, x = log2FC_medians + sign(log2FC_medians) * 0.02,
                hjust = ifelse(log2FC_medians >= 0, 0, 1)), size = 3.2) +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#d7263d"), guide = "none") +
  facet_wrap(~ comparison, scales = "free_x") +
  labs(x = expression("Median "*log[2]*"TPM difference"), y = NULL) +
  theme_publication()
ggsave(file.path(OUT, "Fig2_subgroup_forest.pdf"), g2, width = 11, height = 4)
ggsave(file.path(OUT, "Fig2_subgroup_forest.png"), g2, width = 11, height = 4, dpi = 300)

# ============================================================================
# MODULE 3 — Per-gene KM with optimal cutpoint + univariate Cox (continuous)
# ============================================================================
surv_df <- df |> select(all_of(gene_order), os_time, os_event) |> drop_na()
S <- Surv(surv_df$os_time, surv_df$os_event)

cut_results <- map_dfr(gene_order, function(g) {
  cp <- surv_cutpoint(surv_df, time = "os_time", event = "os_event",
                      variables = g, minprop = 0.25)
  cut_val <- cp$cutpoint$cutpoint
  hi <- surv_df[[g]] > cut_val
  lr <- survdiff(Surv(os_time, os_event) ~ hi, data = surv_df)
  km_p <- 1 - pchisq(lr$chisq, df = 1)
  cox_uni <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ ", g)), data = surv_df)
  sm <- summary(cox_uni)$coefficients[1, ]
  ci <- summary(cox_uni)$conf.int[1, ]
  tibble(gene = g, cutpoint = cut_val,
         n_high = sum(hi), n_low = sum(!hi),
         KM_logrank_p = km_p,
         HR_per_unit  = ci["exp(coef)"],
         HR_low       = ci["lower .95"],
         HR_high      = ci["upper .95"],
         Cox_p        = sm["Pr(>|z|)"])
})
fwrite(cut_results, file.path(OUT, "M3_single_gene_survival.csv"))
print(cut_results, n = 8)

# KM panel
km_list <- map(gene_order, function(g) {
  cut_val <- cut_results$cutpoint[cut_results$gene == g]
  dd <- surv_df |> mutate(group = ifelse(.data[[g]] > cut_val, "High", "Low"),
                          os_time_mo = os_time / 30.44)
  fit <- survfit(Surv(os_time_mo, os_event) ~ group, data = dd)
  ggsurvplot(fit, data = dd, palette = c("#d7263d", "#1f6feb"),
             conf.int = TRUE, pval = TRUE, pval.method = FALSE,
             risk.table = FALSE, legend.title = "", legend.labs = c("High", "Low"),
             xlab = "Months", ylab = "Overall survival", xlim = c(0, 120),
             ggtheme = theme_publication() + theme(plot.title = element_text(face = "italic")),
             title = g)$plot
})
g3 <- ggarrange(plotlist = km_list, ncol = 4, nrow = 2)
ggsave(file.path(OUT, "Fig3_KM_8genes.pdf"), g3, width = 14, height = 7)
ggsave(file.path(OUT, "Fig3_KM_8genes.png"), g3, width = 14, height = 7, dpi = 300)

# ============================================================================
# MODULE 4 — AICD Treg-axis composite score + multivariate Cox
# ============================================================================
axis_genes <- c("FOXP3", "GLUD1", "GLUL", "SMS")
z          <- scale(X[, axis_genes])
df$TregAxis_score        <- rowMeans(z)

# PCA-based 8-gene AICD signature
z8                       <- scale(X)
pc                       <- prcomp(z8, center = FALSE, scale. = FALSE)
sign_flip                <- sign(mean(pc$rotation[, 1]))
df$AICD_PC1              <- pc$x[, 1] * sign_flip
var_pc1                  <- (pc$sdev[1]^2) / sum(pc$sdev^2)
cat(sprintf("PC1 explains %.1f%% variance\n", var_pc1 * 100))

# CD8A/FOXP3 ratio (linear-TPM based, log2-transformed for symmetry)
df$CD8A_TPM              <- as.numeric(expr["CD8A", rownames(df)])
df$FOXP3_TPM             <- as.numeric(expr["FOXP3", rownames(df)])
df$log2_CD8A_FOXP3_ratio <- log2(df$CD8A_TPM + 1) - log2(df$FOXP3_TPM + 1)

# Optimal cutpoint for each composite
cp_axis  <- surv_cutpoint(df, "os_time", "os_event", "TregAxis_score",        minprop = .25)$cutpoint$cutpoint
cp_ratio <- surv_cutpoint(df, "os_time", "os_event", "log2_CD8A_FOXP3_ratio", minprop = .25)$cutpoint$cutpoint
cp_pc1   <- surv_cutpoint(df, "os_time", "os_event", "AICD_PC1",              minprop = .25)$cutpoint$cutpoint

km_score <- function(var, cut, title) {
  dd  <- df |> mutate(grp = ifelse(.data[[var]] > cut, "High", "Low"),
                      mo = os_time / 30.44) |> drop_na(all_of(c(var, "os_time", "os_event")))
  fit <- survfit(Surv(mo, os_event) ~ grp, data = dd)
  ggsurvplot(fit, data = dd, palette = c("#d7263d", "#1f6feb"),
             conf.int = TRUE, pval = TRUE, pval.method = FALSE,
             risk.table = FALSE, legend.title = "",
             legend.labs = c(paste0("High (n=", sum(dd$grp == "High"), ")"),
                             paste0("Low (n=",  sum(dd$grp == "Low"),  ")")),
             xlab = "Months", ylab = "Overall survival", xlim = c(0, 120),
             ggtheme = theme_publication(), title = title)$plot
}
g4 <- ggarrange(
  km_score("TregAxis_score",        cp_axis,  "Treg-axis score (FOXP3+GLUD1+GLUL+SMS)"),
  km_score("log2_CD8A_FOXP3_ratio", cp_ratio, expression(log[2]*"(CD8A / FOXP3)")),
  km_score("AICD_PC1",              cp_pc1,   "AICD signature (PC1, 8 genes)"),
  ncol = 3
)
ggsave(file.path(OUT, "Fig4_KM_scores.pdf"), g4, width = 13, height = 4)
ggsave(file.path(OUT, "Fig4_KM_scores.png"), g4, width = 13, height = 4, dpi = 300)

# Multivariate Cox
cox_df <- df |> select(os_time, os_event, TregAxis_score, log2_CD8A_FOXP3_ratio,
                       AICD_PC1, age, gender, pathologic_stage, mismatch_repair_proteins) |>
  mutate(gender_M       = as.integer(gender == "Male"),
         stage_advanced = as.integer(pathologic_stage == "Stage III/IV"),
         MMR_pMMR       = as.integer(mismatch_repair_proteins == "pMMR")) |>
  drop_na()
cat(sprintf("Multivariate Cox n = %d\n", nrow(cox_df)))

cox_multi <- map_dfr(c("TregAxis_score", "log2_CD8A_FOXP3_ratio", "AICD_PC1"), function(v) {
  fm <- as.formula(paste0("Surv(os_time, os_event) ~ ", v,
                          " + age + gender_M + stage_advanced + MMR_pMMR"))
  fit <- coxph(fm, data = cox_df)
  s   <- summary(fit)$conf.int
  p   <- summary(fit)$coefficients[, "Pr(>|z|)"]
  tibble(model    = v,
         variable = rownames(s),
         HR       = s[, "exp(coef)"],
         low      = s[, "lower .95"],
         high     = s[, "upper .95"],
         p        = p)
})
fwrite(cox_multi, file.path(OUT, "M4_multivariate_cox.csv"))
print(cox_multi)

# Univariate Cox of the 3 composite scores
uni_scores <- map_dfr(c("TregAxis_score", "log2_CD8A_FOXP3_ratio", "AICD_PC1"), function(v) {
  fit <- coxph(as.formula(paste0("Surv(os_time, os_event) ~ ", v)), data = df)
  s   <- summary(fit)$conf.int
  tibble(score = v, HR = s[, "exp(coef)"], low = s[, "lower .95"],
         high  = s[, "upper .95"], p = summary(fit)$coefficients[, "Pr(>|z|)"])
})
fwrite(uni_scores, file.path(OUT, "M4_univariate_cox_scores.csv"))
print(uni_scores)

# ============================================================================
# MODULE 5 — CD8A/FOXP3 ratio across clinical subgroups
# ============================================================================
plot_box <- function(col, g1, g2, title) {
  sub <- df |> filter(.data[[col]] %in% c(g1, g2)) |>
    mutate(grp = factor(.data[[col]], levels = c(g2, g1)))
  p   <- wilcox.test(sub$log2_CD8A_FOXP3_ratio[sub$grp == g1],
                     sub$log2_CD8A_FOXP3_ratio[sub$grp == g2],
                     exact = FALSE)$p.value
  ptxt <- if (p < 0.001) "p < 0.001" else paste0("p = ", formatC(p, format = "f", digits = 3))
  ggplot(sub, aes(grp, log2_CD8A_FOXP3_ratio, fill = grp)) +
    geom_boxplot(width = .5, outlier.size = .5, linewidth = .3) +
    geom_jitter(width = .15, color = "black", alpha = .35, size = .6) +
    scale_fill_manual(values = c("#1f6feb", "#d7263d"), guide = "none") +
    annotate("text", x = 1.5, y = max(sub$log2_CD8A_FOXP3_ratio, na.rm = TRUE) + .25,
             label = ptxt, size = 3.5) +
    labs(x = NULL, y = expression(log[2]*"(CD8A / FOXP3)"), title = title) +
    theme_publication()
}
g5 <- ggarrange(plot_box("mismatch_repair_proteins", "dMMR",         "pMMR",       "MMR status"),
                plot_box("M_Stage",                  "M1",           "M0",         "M stage"),
                plot_box("pathologic_stage",         "Stage III/IV", "Stage I/II", "Pathologic stage"),
                ncol = 3)
ggsave(file.path(OUT, "Fig5_ratio_subgroups.pdf"), g5, width = 11, height = 4)
ggsave(file.path(OUT, "Fig5_ratio_subgroups.png"), g5, width = 11, height = 4, dpi = 300)

m5 <- map_dfr(list(c("mismatch_repair_proteins", "dMMR",         "pMMR"),
                   c("M_Stage",                  "M1",           "M0"),
                   c("pathologic_stage",         "Stage III/IV", "Stage I/II")),
              function(co) {
                a <- df$log2_CD8A_FOXP3_ratio[df[[co[1]]] == co[2]]; a <- a[!is.na(a)]
                b <- df$log2_CD8A_FOXP3_ratio[df[[co[1]]] == co[3]]; b <- b[!is.na(b)]
                tibble(group = co[1], g1 = co[2], n1 = length(a),
                       g2    = co[3], n2 = length(b),
                       median_g1 = median(a), median_g2 = median(b),
                       p = wilcox.test(a, b, exact = FALSE)$p.value)
              })
fwrite(m5, file.path(OUT, "M5_ratio_vs_subgroups.csv"))
print(m5)

# ---------- Session info ----------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
cat("\n=== DONE ===\n")
cat("All outputs written to ./", OUT, "/\n", sep = "")
