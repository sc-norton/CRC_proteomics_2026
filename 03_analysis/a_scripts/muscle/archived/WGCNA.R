# =============================================================================
# WGCNA Analysis — Muscle & Plasma Proteomics
# =============================================================================
# Networks built on full dataset (SURV T1 + SURV T2 + CTL T1)
# Two trait analyses:
#   (1) SURVvCTL: all samples, traits = group + baseline body composition
#   (2) SURV change: SURV samples only, traits = delta scores
# =============================================================================

library(WGCNA)
library(ggplot2)
library(ggpubr)
library(ggsignif)
library(dplyr)

options(stringsAsFactors = FALSE)
enableWGCNAThreads()

setwd(rprojroot::find_rstudio_root_file())
getwd()
# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------
tissues     <- c("muscle", "plasma")

# Set after inspecting SFT plots — placeholders until then
softPowers  <- c(muscle = 3, plasma = 3)

# -----------------------------------------------------------------------------
# 1. Load imputed matrices
# -----------------------------------------------------------------------------
# Matrices saved as CSV with proteins as rows (rownames = UniProt accessions)
# and samples as columns. First column is the rownames column from write.csv().

expr_t1     <- list()
expr_t2     <- list()
expr_t1vctl <- list()
expr_full   <- list()

for (tissue in tissues) {
  
  path <- file.path("02_imputation", "c_data",
                    paste0(tissue, "_imputed_matrix.csv"))
  
  df  <- read.csv(path, row.names = 1, check.names = FALSE)
  # df is now proteins x samples with correct rownames
  
  # Transpose to samples x proteins for WGCNA
  mat <- t(as.matrix(df))
  
  # CTL samples have no T2 — t1vctl captures SURV_T1 + CTL_T1
  expr_t1[[tissue]]     <- mat[grepl("_T1", rownames(mat)), ]
  expr_t2[[tissue]]     <- mat[grepl("_T2", rownames(mat)), ]
  expr_t1vctl[[tissue]] <- mat[!grepl("_T2", rownames(mat)), ]
  expr_full[[tissue]]   <- mat
  
  cat(tissue, "— dimensions:\n")
  cat("  expr_t1:     ", nrow(expr_t1[[tissue]]),     "samples x",
      ncol(expr_t1[[tissue]]),     "proteins\n")
  cat("  expr_t2:     ", nrow(expr_t2[[tissue]]),     "samples x",
      ncol(expr_t2[[tissue]]),     "proteins\n")
  cat("  expr_t1vctl: ", nrow(expr_t1vctl[[tissue]]), "samples x",
      ncol(expr_t1vctl[[tissue]]), "proteins\n")
  cat("  expr_full:   ", nrow(expr_full[[tissue]]),   "samples x",   
      ncol(expr_full[[tissue]]),   "proteins\n\n")
}

# -----------------------------------------------------------------------------
# 2. Load metadata
# -----------------------------------------------------------------------------
meta        <- list()
meta_t1     <- list()
meta_t2     <- list()
meta_t1vctl <- list()
meta_full   <- list()     # ADD: all samples
meta_surv   <- list()     # ADD: all SURV samples (T1 + T2)

for (tissue in tissues) {
  
  meta_file <- ifelse(tissue == "muscle", "CRm_meta.csv", "CRp_meta.csv")
  
  meta[[tissue]] <- read.csv(
    file.path("00_input", tissue, meta_file),
    row.names = 1,
    check.names = FALSE
  )
  
  meta_t1[[tissue]]     <- meta[[tissue]][rownames(expr_t1[[tissue]]), ]
  meta_t2[[tissue]]     <- meta[[tissue]][rownames(expr_t2[[tissue]]), ]
  meta_t1vctl[[tissue]] <- meta[[tissue]][rownames(expr_t1vctl[[tissue]]), ]
  meta_full[[tissue]]   <- meta[[tissue]][rownames(expr_full[[tissue]]), ]    # ADD
  
  # SURV samples only (T1 + T2) — for timepoint analysis
  surv_rows <- rownames(expr_full[[tissue]])[
    grepl("_T1|_T2", rownames(expr_full[[tissue]])) &
      meta_full[[tissue]]$cancer == "SURV"
  ]
  meta_surv[[tissue]] <- meta[[tissue]][surv_rows, ]                          # ADD
}

# -----------------------------------------------------------------------------
# 3. Compute delta scores and attach to metadata
# -----------------------------------------------------------------------------

add_deltas <- function(m) {
  m$delta_weight_kg       <- m$post_weight_kg       - m$pre_weight_kg
  m$delta_BMI             <- m$post_BMI             - m$pre_BMI
  m$delta_LBM_kg          <- m$post_LBM_kg          - m$pre_LBM_kg
  m$delta_ALM_kg          <- m$post_ALM_kg          - m$pre_ALM_kg
  m$delta_ALM_index       <- m$post_ALM_index       - m$pre_ALM_index
  m$delta_fat_kg          <- m$post_fat_kg          - m$pre_fat_kg
  m$delta_fat_pct         <- m$post_fat_pct         - m$pre_fat_pct
  m$delta_sts_max_pwr     <- m$post_sts_max_pwr     - m$pre_sts_max_pwr
  m$delta_grip_lbs        <- m$post_grip_lbs        - m$pre_grip_lbs
  m$delta_chest_press_lbs <- m$post_chest_press_lbs - m$pre_chest_press_lbs
  m$delta_leg_ext_lbs     <- m$post_leg_ext_lbs     - m$pre_leg_ext_lbs
  m
}

meta_t1     <- lapply(meta_t1,     add_deltas)
meta_t2     <- lapply(meta_t2,     add_deltas)
meta_t1vctl <- lapply(meta_t1vctl, add_deltas)
meta_full   <- lapply(meta_full,   add_deltas)   # ADD
meta_surv   <- lapply(meta_surv,   add_deltas)   # ADD

# -----------------------------------------------------------------------------
# 4. Build trait matrices
# -----------------------------------------------------------------------------

# Traits for SURVvCTL analysis (all samples at baseline + CTL)
make_traits_vctl <- function(m) {
  data.frame(
    age           = m$age,
    pre_weight_kg = m$pre_weight_kg,
    pre_BMI       = m$pre_BMI,
    pre_ALM_kg    = m$pre_ALM_kg,
    pre_ALM_index = m$pre_ALM_index,
    pre_fat_pct   = m$pre_fat_pct,
    cancer_SURV   = as.integer(m$cancer == "SURV"),
    sex_Male      = as.integer(m$sex == "Male"),
    row.names     = rownames(m)
  )
}

# Traits for SURV-only change score analysis
make_traits_surv <- function(m) {
  data.frame(
    age                   = m$age,
    yrs_since_dx          = m$yrs_since_dx,
    yrs_since_chemo       = m$yrs_since_chemo,
    delta_weight_kg       = m$delta_weight_kg,
    delta_BMI             = m$delta_BMI,
    delta_LBM_kg          = m$delta_LBM_kg,
    delta_ALM_kg          = m$delta_ALM_kg,
    delta_ALM_index       = m$delta_ALM_index,
    delta_fat_pct         = m$delta_fat_pct,
    delta_sts_max_pwr     = m$delta_sts_max_pwr,
    delta_grip_lbs        = m$delta_grip_lbs,
    delta_chest_press_lbs = m$delta_chest_press_lbs,
    delta_leg_ext_lbs     = m$delta_leg_ext_lbs,
    supp_CRE              = as.integer(m$supp == "CRE"),
    sex_Male              = as.integer(m$sex == "Male"),
    row.names             = rownames(m)
  )
}

# ADD: Traits for timepoint analysis (all SURV samples, T1 + T2)
make_traits_timepoint <- function(m) {
  data.frame(
    timepoint = as.integer(grepl("_T2", rownames(m))),   # T1=0, T2=1
    supp_CRE  = as.integer(m$supp == "CRE"),
    sex_Male  = as.integer(m$sex == "Male"),
    row.names = rownames(m)
  )
}

traits_t1vctl <- lapply(meta_t1vctl, make_traits_vctl)
traits_t1     <- lapply(meta_t1,     make_traits_surv)
traits_surv   <- lapply(meta_surv,   make_traits_timepoint)   # ADD

# -----------------------------------------------------------------------------
# 5. Sample quality check
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  sampleTree <- hclust(dist(expr_full[[tissue]]), method = "average")
  plot(sampleTree,
       main = paste("Sample clustering —", tissue, "(FULL)"),
       sub  = "", xlab = "")
}

# -----------------------------------------------------------------------------
# 6. Soft thresholding
# -----------------------------------------------------------------------------

sft_results <- list()
powers      <- c(1:10, seq(12, 30, 2))

for (tissue in tissues) {
  cat("\n--- Soft thresholding:", tissue, "---\n")
  # Change expr_t1vctl → expr_full in pickSoftThreshold call
  sft_results[[tissue]] <- pickSoftThreshold(
    expr_full[[tissue]],       # CHANGED from expr_t1vctl
    powerVector = powers,
    networkType = "signed hybrid",
    corFnc      = "bicor",
    verbose     = 3
  )
  
  # And in the SFT plots, update the sample clustering call:
  sampleTree <- hclust(dist(expr_full[[tissue]]), method = "average")   # CHANGED
  plot(sampleTree,
       main = paste("Sample clustering —", tissue, "(all samples)"),    # CHANGED
       sub  = "", xlab = "")
  
  par(mfrow = c(1, 2))
  
  # Scale-free topology fit
  plot(sft_results[[tissue]]$fitIndices[, 1],
       -sign(sft_results[[tissue]]$fitIndices[, 3]) *
         sft_results[[tissue]]$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (R²)",
       type = "n",
       main = paste(tissue, "— Scale independence"))
  text(sft_results[[tissue]]$fitIndices[, 1],
       -sign(sft_results[[tissue]]$fitIndices[, 3]) *
         sft_results[[tissue]]$fitIndices[, 2],
       labels = powers, col = "red")
  abline(h = 0.80, col = "red")
  
  # Mean connectivity
  plot(sft_results[[tissue]]$fitIndices[, 1],
       sft_results[[tissue]]$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       type = "n",
       main = paste(tissue, "— Mean connectivity"))
  text(sft_results[[tissue]]$fitIndices[, 1],
       sft_results[[tissue]]$fitIndices[, 5],
       labels = powers, col = "red")
}

# Update softPowers in Config (section 0) after inspecting plots

# -----------------------------------------------------------------------------
# 7. Build networks
# -----------------------------------------------------------------------------

nets <- list()

for (tissue in tissues) {
  cat("\n--- Building network:", tissue, "---\n")
  
  nets[[tissue]] <- blockwiseModules(
    expr_full[[tissue]],       # CHANGED from expr_t1vctl
    power             = softPowers[[tissue]],
    networkType       = "signed hybrid",
    TOMType           = "signed",
    corType           = "bicor",
    minModuleSize     = 15,
    mergeCutHeight    = 0.25,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = file.path("02_imputation", "c_data",
                                  paste0(tissue, "_TOM")),
    verbose           = 3
  )
  
  cat("\n", tissue, "— module sizes:\n")
  print(table(nets[[tissue]]$colors))
  
  plotDendroAndColors(
    nets[[tissue]]$dendrograms[[1]],
    nets[[tissue]]$colors[nets[[tissue]]$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE,
    hang         = 0.03,
    addGuide     = TRUE,
    main         = paste(tissue, "— Module Dendrogram")
  )
}

# -----------------------------------------------------------------------------
# 8. Sanity checks — gene and sample coverage
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  surv_samples <- rownames(traits_t1[[tissue]])
  ME_samples   <- rownames(nets[[tissue]]$MEs)
  
  cat(tissue, "— SURV samples in MEs:", all(surv_samples %in% ME_samples), "\n")
  cat(tissue, "— missing SURV samples:", sum(!surv_samples %in% ME_samples), "\n")
  cat(tissue, "— total proteins in network:", length(nets[[tissue]]$colors), "\n")
  cat(tissue, "— proteins in expr_t1vctl:  ", ncol(expr_t1vctl[[tissue]]), "\n")
  cat(tissue, "— proteins in expr_full:     ", ncol(expr_full[[tissue]]), "\n\n")  # CHANGED
}

# -----------------------------------------------------------------------------
# 9. Module-trait heatmap helper
# -----------------------------------------------------------------------------

plot_module_trait_heatmap <- function(tissue, MEs, traits, n_samples,
                                      title, outfile) {
  moduleTraitCor  <- cor(MEs, traits, use = "p")
  moduleTraitPval <- corPvalueStudent(moduleTraitCor, n_samples)
  
  n_modules   <- ncol(MEs)
  bonf_cutoff <- 0.05 / n_modules
  
  textMatrix <- ifelse(
    moduleTraitPval < bonf_cutoff,
    paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPval, 1), ")**"),
    ifelse(
      moduleTraitPval < 0.05,
      paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPval, 1), ")*"),
      paste0(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPval, 1), ")")
    )
  )
  
  pdf(file = outfile, width = 10, height = 7)
  par(mar = c(8, 12, 3, 2))
  labeledHeatmap(
    Matrix      = moduleTraitCor,
    xLabels     = colnames(traits),
    yLabels     = colnames(MEs),
    colorLabels = FALSE,
    colors      = blueWhiteRed(50),
    textMatrix  = textMatrix,
    cex.text    = 0.6,
    cex.lab     = 0.8,
    main        = title
  )
  dev.off()
  
  invisible(list(cor = moduleTraitCor, pval = moduleTraitPval))
}

# Ensure output directory exists
dir.create(file.path("04_figures", "WGCNA"), recursive = TRUE,
           showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 10. Module-trait heatmaps: SURVvCTL (all samples)
# -----------------------------------------------------------------------------

mt_vctl <- list()

for (tissue in tissues) {
  MEs <- orderMEs(nets[[tissue]]$MEs)
  
  mt_vctl[[tissue]] <- plot_module_trait_heatmap(
    tissue    = tissue,
    MEs       = MEs,
    traits    = traits_t1vctl[[tissue]],
    n_samples = nrow(expr_full[[tissue]]),
    title     = paste(tissue, "— Module-Trait Relationships (SURV & CTL)"),
    outfile   = file.path("04_figures", "WGCNA",
                          paste0(tissue, "_module_trait_heatmap_SURVvCTL.pdf"))
  )
}

# -----------------------------------------------------------------------------
# 11. Module-trait heatmaps: SURV only (change scores)
# -----------------------------------------------------------------------------

mt_surv <- list()

for (tissue in tissues) {
  MEs_all      <- orderMEs(nets[[tissue]]$MEs)
  surv_samples <- rownames(traits_t1[[tissue]])
  MEs_surv     <- MEs_all[surv_samples, ]
  
  mt_surv[[tissue]] <- plot_module_trait_heatmap(
    tissue    = tissue,
    MEs       = MEs_surv,
    traits    = traits_t1[[tissue]],
    n_samples = nrow(MEs_surv),
    title     = paste(tissue, "— Module-Trait Relationships (SURV change scores)"),
    outfile   = file.path("04_figures", "WGCNA",
                          paste0(tissue, "_module_trait_heatmap_SURV_change.pdf"))
  )
}

# -----------------------------------------------------------------------------
# 11b. Module-trait heatmaps: SURV only (timepoint T1 vs T2)
# -----------------------------------------------------------------------------

mt_timepoint <- list()

for (tissue in tissues) {
  MEs_all      <- orderMEs(nets[[tissue]]$MEs)
  surv_samples <- rownames(traits_surv[[tissue]])
  MEs_surv     <- MEs_all[surv_samples, ]
  
  mt_timepoint[[tissue]] <- plot_module_trait_heatmap(
    tissue    = tissue,
    MEs       = MEs_surv,
    traits    = traits_surv[[tissue]],
    n_samples = nrow(MEs_surv),
    title     = paste(tissue, "— Module-Trait Relationships (SURV timepoint)"),
    outfile   = file.path("04_figures", "WGCNA",
                          paste0(tissue, "_module_trait_heatmap_SURV_timepoint.pdf"))
  )
}

# -----------------------------------------------------------------------------
# 11c. LMM-based ME-trait heatmap (mixed model approach)
# -----------------------------------------------------------------------------

library(lme4)
library(lmerTest)

# --- Helper: extract r-equivalent from LMM ---
# Uses semi-partial r² approach: sqrt(F * df_num / (F * df_num + df_den))
# signed by the fixed effect estimate direction
lmm_r_equiv <- function(ME, timepoint, group, subject) {
  
  df <- data.frame(
    ME        = ME,
    timepoint = timepoint,   # 0/1
    group     = group,       # "SURV"/"CTL"
    subject   = subject      # subject ID
  )
  
  # Remove rows with NA
  df <- df[complete.cases(df), ]
  
  # Fit LMM: ME ~ timepoint * group + (1 | subject)
  fit <- tryCatch(
    lmerTest::lmer(ME ~ timepoint * group + (1 | subject), data = df,
                   REML = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(c(r_timepoint = NA, r_group = NA, r_interaction = NA))
  
  coefs <- summary(fit)$coefficients
  
  # Extract t-values and convert to r-equivalent: r = t / sqrt(t² + df)
  # df approximated as n_obs - n_fixed_params
  n  <- nrow(df)
  p  <- length(fixef(fit))
  df_resid <- n - p
  
  t_to_r <- function(t) t / sqrt(t^2 + df_resid)
  
  r_group       <- if ("groupSURV" %in% rownames(coefs))
    t_to_r(coefs["groupSURV", "t value"]) else NA
  r_timepoint   <- if ("timepoint" %in% rownames(coefs))
    t_to_r(coefs["timepoint", "t value"]) else NA
  r_interaction <- if ("timepoint:groupSURV" %in% rownames(coefs))
    t_to_r(coefs["timepoint:groupSURV", "t value"]) else NA
  
  c(r_group       = r_group,
    r_timepoint   = r_timepoint,
    r_interaction = r_interaction)
}

lmm_pvals <- function(ME, timepoint, group, subject) {
  
  df <- data.frame(ME = ME, timepoint = timepoint,
                   group = group, subject = subject)
  df <- df[complete.cases(df), ]
  
  fit <- tryCatch(
    lmerTest::lmer(ME ~ timepoint * group + (1 | subject), data = df,
                   REML = FALSE),
    error = function(e) NULL
  )
  
  if (is.null(fit)) return(c(p_group = NA, p_timepoint = NA, p_interaction = NA))
  
  coefs <- summary(fit)$coefficients
  
  p_group       <- if ("groupSURV" %in% rownames(coefs))
    coefs["groupSURV", "Pr(>|t|)"] else NA
  p_timepoint   <- if ("timepoint" %in% rownames(coefs))
    coefs["timepoint", "Pr(>|t|)"] else NA
  p_interaction <- if ("timepoint:groupSURV" %in% rownames(coefs))
    coefs["timepoint:groupSURV", "Pr(>|t|)"] else NA
  
  c(p_group       = p_group,
    p_timepoint   = p_timepoint,
    p_interaction = p_interaction)
}

# --- Build combined heatmap matrix ---

build_combined_trait_matrix <- function(tissue, nets, expr_full,
                                        meta_full, meta_t1) {
  
  MEs      <- orderMEs(nets[[tissue]]$MEs)
  ME_names <- colnames(MEs)
  
  # Subject IDs — strip timepoint suffix to get subject ID
  # Assumes sample names like "CR006_T1", "CR006_T2"
  all_samples <- rownames(MEs)
  subject_ids <- sub("_T[12]$", "", all_samples)
  timepoint   <- as.integer(grepl("_T2", all_samples))
  group       <- meta_full[[tissue]]$cancer   # "SURV" or "CTL"
  
  # --- LMM r-equivalents (left columns) ---
  lmm_r_mat <- t(sapply(ME_names, function(me) {
    lmm_r_equiv(
      ME        = MEs[, me],
      timepoint = timepoint,
      group     = group,
      subject   = subject_ids
    )
  }))
  # rows = modules, cols = r_group, r_timepoint, r_interaction
  colnames(lmm_r_mat) <- c("Group\n(SURVvCTL)",
                           "Timepoint\n(Training)",
                           "Group×Time\n(Interaction)")
  
  lmm_p_mat <- t(sapply(ME_names, function(me) {
    lmm_pvals(
      ME        = MEs[, me],
      timepoint = timepoint,
      group     = group,
      subject   = subject_ids
    )
  }))
  colnames(lmm_p_mat) <- colnames(lmm_r_mat)
  
  # --- Pearson correlations (right columns) ---
  # Baseline phenotypes — SURV T1 only
  surv_t1_samples <- rownames(meta_t1[[tissue]])[
    meta_t1[[tissue]]$cancer == "SURV"
  ]
  MEs_surv_t1 <- MEs[surv_t1_samples, ]
  
  pearson_traits <- make_traits_surv(meta_t1[[tissue]][
    meta_t1[[tissue]]$cancer == "SURV", 
  ])
  
  pearson_r_mat <- cor(MEs_surv_t1, pearson_traits, use = "p")
  pearson_p_mat <- corPvalueStudent(pearson_r_mat, nrow(MEs_surv_t1))
  
  # --- Combine ---
  combined_r <- cbind(lmm_r_mat, pearson_r_mat)
  combined_p <- cbind(lmm_p_mat, pearson_p_mat)
  
  list(r = combined_r, p = combined_p)
}

# --- Run and plot ---

for (tissue in tissues) {
  
  cat("\n--- LMM + Pearson heatmap:", tissue, "---\n")
  
  res <- build_combined_trait_matrix(
    tissue    = tissue,
    nets      = nets,
    expr_full = expr_full,
    meta_full = meta_full,
    meta_t1   = meta_t1
  )
  
  n_modules   <- nrow(res$r)
  bonf_cutoff <- 0.05 / n_modules
  
  textMatrix <- ifelse(
    res$p < bonf_cutoff,
    paste0(signif(res$r, 2), "\n(", signif(res$p, 1), ")**"),
    ifelse(
      res$p < 0.05,
      paste0(signif(res$r, 2), "\n(", signif(res$p, 1), ")*"),
      paste0(signif(res$r, 2), "\n(", signif(res$p, 1), ")")
    )
  )
  
  # Add dividing line between LMM and Pearson columns via x-axis label spacing
  xlabs <- colnames(res$r)
  
  pdf(file.path("04_figures", "WGCNA",
                paste0(tissue, "_module_trait_heatmap_combined.pdf")),
      width = 14, height = 8)
  par(mar = c(10, 12, 3, 2))
  
  labeledHeatmap(
    Matrix      = res$r,
    xLabels     = xlabs,
    yLabels     = rownames(res$r),
    colorLabels = FALSE,
    colors      = blueWhiteRed(50),
    textMatrix  = textMatrix,
    cex.text    = 0.55,
    cex.lab     = 0.75,
    main        = paste(tissue,
                        "— Module-Trait Relationships (LMM + Pearson)")
  )
  
  # Draw vertical line separating LMM from Pearson columns
  abline(v = 3.5, lwd = 2, lty = 2, col = "black")
  
  dev.off()
}
# -----------------------------------------------------------------------------
# 12. Eigengene network heatmap
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  MEs <- orderMEs(nets[[tissue]]$MEs)
  
  pdf(file = file.path("04_figures", "WGCNA",
                       paste0(tissue, "_eigengene_network.pdf")),
      width = 12, height = 10)
  par(cex = 0.8)
  plotEigengeneNetworks(
    MEs,
    setLabels    = paste(tissue, "Eigengene Network"),
    marDendro    = c(0, 4, 2, 0),
    marHeatmap   = c(3, 4, 2, 2),
    cex.lab      = 0.8,
    xLabelsAngle = 90
  )
  dev.off()
}

# -----------------------------------------------------------------------------
# 13. ME vs trait scatter plots
# -----------------------------------------------------------------------------

# SURVvCTL — all samples
plot_ME_trait_vctl <- function(tissue, module, trait_name) {
  MEs   <- orderMEs(nets[[tissue]]$MEs)
  ME    <- MEs[, paste0("ME", module)]
  trait <- traits_t1vctl[[tissue]][[trait_name]]
  
  plot_df <- data.frame(
    ME    = ME,
    trait = trait,
    row.names = rownames(traits_t1vctl[[tissue]])
  )
  
  ggplot(plot_df, aes(x = ME, y = trait)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = module) +
    stat_cor(method = "pearson") +
    labs(x     = paste0("ME (", module, ")"),
         y     = trait_name,
         title = paste(tissue, "|", module, "vs", trait_name,
                       "(SURV & CTL)")) +
    theme_classic()
}

# SURV only — change scores
plot_ME_trait_surv <- function(tissue, module, trait_name) {
  MEs_all      <- orderMEs(nets[[tissue]]$MEs)
  surv_samples <- rownames(traits_t1[[tissue]])
  ME           <- MEs_all[surv_samples, paste0("ME", module)]
  trait        <- traits_t1[[tissue]][[trait_name]]
  
  plot_df <- data.frame(
    ME    = ME,
    trait = trait,
    row.names = surv_samples
  )
  
  ggplot(plot_df, aes(x = ME, y = trait)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", se = TRUE, color = module) +
    stat_cor(method = "pearson") +
    labs(x     = paste0("ME (", module, ")"),
         y     = trait_name,
         title = paste(tissue, "|", module, "vs", trait_name,
                       "(SURV change)")) +
    theme_classic()
}

# Example usage (update module colors after inspecting heatmaps):
plot_ME_trait_vctl("muscle", "brown", "cancer_SURV")
# plot_ME_trait_vctl("plasma", "blue",  "sex_Male")
plot_ME_trait_surv("muscle", "brown", "delta_ALM_kg")
# plot_ME_trait_surv("plasma", "blue",  "delta_ALM_kg")

# -----------------------------------------------------------------------------
# 14. ME violin plots by group
# -----------------------------------------------------------------------------

plot_ME_group <- function(tissue, module, group_var, meta_list) {
  MEs   <- orderMEs(nets[[tissue]]$MEs)
  ME    <- MEs[, paste0("ME", module)]
  group <- meta_list[[tissue]][[group_var]]
  
  plot_df <- data.frame(
    ME    = ME,
    group = factor(group),
    row.names = rownames(meta_list[[tissue]])
  )
  
  n_groups <- nlevels(plot_df$group)
  shades   <- colorRampPalette(c("white", module))(n_groups + 2)[2:(n_groups + 1)]
  names(shades) <- levels(plot_df$group)
  
  ggplot(plot_df, aes(x = group, y = ME, fill = group)) +
    geom_violin(width = 0.6, color = "black", trim = FALSE) +
    geom_boxplot(width = 0.03, color = "black", outlier.shape = NA) +
    geom_jitter(width = 0.08, size = 2, alpha = 0.6) +
    geom_signif(
      comparisons      = list(levels(plot_df$group)),
      test             = "t.test",
      map_signif_level = FALSE,
      tip_length       = 0.01,
      textsize         = 3.5,
      vjust            = -0.2,
      y_position       = max(plot_df$ME, na.rm = TRUE) * 1.1
    ) +
    scale_fill_manual(values = shades) +
    labs(x     = group_var,
         y     = paste0("ME (", module, ")"),
         title = paste(tissue, "|", module, "module by", group_var)) +
    theme_classic() +
    theme(legend.position = "none")
}

# Example usage (update module colors after inspecting heatmaps):
plot_ME_group("muscle", "brown", "cancer", meta_t1vctl)
# plot_ME_group("plasma", "blue",  "cancer", meta_t1vctl)

# -----------------------------------------------------------------------------
# 15. Hub protein extraction
# -----------------------------------------------------------------------------

get_module_info <- function(tissue, module_color) {
  
  module_genes <- names(nets[[tissue]]$colors)[
    nets[[tissue]]$colors == module_color
  ]
  
  MEs      <- orderMEs(nets[[tissue]]$MEs)
  kME      <- cor(expr_t1vctl[[tissue]], MEs, use = "p")
  kME_pval <- corPvalueStudent(kME, nrow(expr_t1vctl[[tissue]]))
  
  ME_col <- paste0("ME", module_color)
  
  data.frame(
    protein  = module_genes,
    kME      = kME[module_genes, ME_col],
    kME_pval = kME_pval[module_genes, ME_col],
    row.names = NULL
  ) |>
    dplyr::arrange(desc(kME))
}

# Example usage (update module colors after inspecting heatmaps):
brown_muscle <- get_module_info("muscle", "brown")
# blue_plasma  <- get_module_info("plasma", "blue")
# head(brown_muscle, 20)
# head(blue_plasma,  20)

# Write full lists — update module colors as needed
# out_dir_muscle <- file.path(base_dir, "05_analysis", "b_results",
#                             "muscle", "SURVvCTL", "WGCNA")
# out_dir_plasma <- file.path(base_dir, "05_analysis", "b_results",
#                             "plasma", "SURVvCTL", "WGCNA")
# dir.create(out_dir_muscle, recursive = TRUE, showWarnings = FALSE)
# dir.create(out_dir_plasma, recursive = TRUE, showWarnings = FALSE)
# write.csv(brown_muscle, file.path(out_dir_muscle,
#           "brown_module_membership_muscle.csv"), row.names = FALSE)
# write.csv(blue_plasma,  file.path(out_dir_plasma,
#           "blue_module_membership_plasma.csv"),  row.names = FALSE)

# -----------------------------------------------------------------------------
# 16. kME-logFC integration with limma contrast results
# -----------------------------------------------------------------------------
# Correlates per-protein module membership (kME) with per-protein logFC
# from each limma contrast. A strong correlation means the hub proteins of
# a module are also the most differentially abundant in that contrast.
#
# Interpretation:
#   Baseline_SURVvCTL r > 0 → hub proteins upregulated in SURV at baseline
#   Training_SURV     r < 0 → hub proteins decreased with training (reversal)
#   Both above together    → reversal module (cross-validates limma trajectories)
# -----------------------------------------------------------------------------

library(readxl)

# Adjust paths and column names to match your actual limma output files
limma_files <- list(
  muscle = file.path("03_analysis", "b_reports", "muscle_DA_results.xlsx"),
  plasma = file.path("03_analysis", "b_reports", "plasma_DA_results.xlsx")
)

contrasts_SURV <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")

get_kME_lfc_cor <- function(tissue, limma_file, contrasts) {
  
  MEs  <- orderMEs(nets[[tissue]]$MEs)
  kME  <- cor(expr_full[[tissue]], MEs, use = "p")
  
  # Load limma results — adjust sheet/column names as needed
  # Assumes columns: "contrast", "UniProt", "logFC"
  limma_dat <- read_xlsx(limma_file)
  
  results <- lapply(contrasts, function(ct) {
    
    lfc <- limma_dat |>
      dplyr::filter(contrast == ct) |>
      dplyr::select(UniProt, logFC) |>
      dplyr::distinct(UniProt, .keep_all = TRUE)
    
    lfc_vec <- setNames(lfc$logFC, lfc$UniProt)
    shared  <- intersect(rownames(kME), names(lfc_vec))
    
    if (length(shared) < 10) {
      warning("Fewer than 10 shared proteins for contrast: ", ct, " in ", tissue)
      return(NULL)
    }
    
    kME_sub <- kME[shared, , drop = FALSE]
    lfc_sub <- lfc_vec[shared]
    
    cor_vals <- apply(kME_sub, 2, function(me_vals) {
      cor(me_vals, lfc_sub, method = "pearson", use = "complete.obs")
    })
    
    data.frame(
      module   = gsub("^ME", "", names(cor_vals)),
      contrast = ct,
      r        = cor_vals,
      row.names = NULL
    )
  })
  
  dplyr::bind_rows(results)
}

kME_lfc <- list()

for (tissue in tissues) {
  kME_lfc[[tissue]] <- get_kME_lfc_cor(
    tissue     = tissue,
    limma_file = limma_files[[tissue]],
    contrasts  = contrasts_SURV
  )
  
  cat("\n--- kME-logFC correlations:", tissue, "---\n")
  print(
    kME_lfc[[tissue]] |>
      dplyr::arrange(contrast, desc(abs(r))) |>
      dplyr::group_by(contrast) |>
      dplyr::slice_head(n = 5)
  )
}

# Reversal modules: strong positive r for Baseline_SURVvCTL AND
#                   strong negative r for Training_SURV
reversal_modules <- lapply(tissues, function(tissue) {
  df <- kME_lfc[[tissue]]
  
  baseline <- df |> dplyr::filter(contrast == "Baseline_SURVvCTL") |>
    dplyr::select(module, r_baseline = r)
  training <- df |> dplyr::filter(contrast == "Training_SURV") |>
    dplyr::select(module, r_training = r)
  
  dplyr::inner_join(baseline, training, by = "module") |>
    dplyr::mutate(
      reversal_score = r_baseline - r_training,   # large = strong reversal
      reversal_flag  = r_baseline > 0.2 & r_training < -0.2
    ) |>
    dplyr::arrange(desc(reversal_score))
})
names(reversal_modules) <- tissues

cat("\n--- Reversal module candidates ---\n")
for (tissue in tissues) {
  cat("\n", tissue, ":\n")
  print(reversal_modules[[tissue]])
}

# Save kME-logFC results
for (tissue in tissues) {
  write.csv(
    kME_lfc[[tissue]],
    file.path("04_figures", "WGCNA",
              paste0(tissue, "_kME_logFC_correlations.csv")),
    row.names = FALSE
  )
  write.csv(
    reversal_modules[[tissue]],
    file.path("04_figures", "WGCNA",
              paste0(tissue, "_reversal_module_candidates.csv")),
    row.names = FALSE
  )
}