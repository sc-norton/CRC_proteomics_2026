# =============================================================================
# WGCNA Analysis — Muscle & Plasma Proteomics
# =============================================================================
# Adapted from Gai Gianetto et al. framework.
# Networks built on full dataset (SURV_T1 + SURV_T2 + CTL_T1).
# Three LMM contrasts (cancer_time factor, 3 levels):
#   (1) SURV_T1 vs CTL_T1  — baseline dysregulation
#   (2) SURV_T2 vs SURV_T1 — training effect within survivors
#   (3) SURV_T2 vs CTL_T1  — post-training dysregulation
# Pearson correlations stratified by group for baseline and change traits.
# =============================================================================

library(WGCNA)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()   # single-threaded for byte-reproducible TOM across re-runs
set.seed(26)

setwd(rprojroot::find_rstudio_root_file())
getwd()
# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------

tissues <- c("muscle", "plasma")

WGCNA_NETWORK_TYPE <- "signed hybrid"
WGCNA_TOM_TYPE     <- "signed"
WGCNA_COR_TYPE     <- "bicor"
WGCNA_MIN_MOD_SIZE <- 15L
WGCNA_MERGE_CUT_H  <- 0.25
WGCNA_R2_CUTOFF    <- 0.80      # relaxed from 0.85/0.90 for small-n proteomics

softPowers <- c(muscle = 3, plasma = 3)   # update after SFT inspection

DATA_DIR    <- "02_imputation/c_data"
PANEL_DIR   <- "04_figures/WGCNA"
REPORT_DIR  <- "04_figures/WGCNA/reports"

dir.create(PANEL_DIR,  recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

# Override cor() with WGCNA version during network computations
cor <- WGCNA::cor

# -----------------------------------------------------------------------------
# 1. Load imputed matrices
# -----------------------------------------------------------------------------

expr_full   <- list()   # all samples: SURV_T1 + SURV_T2 + CTL_T1
expr_t1     <- list()   # T1 only (SURV_T1 + CTL_T1)
expr_t2     <- list()   # T2 only (SURV_T2)
expr_t1vctl <- list()   # SURV_T1 + CTL_T1 (alias for t1 — kept for compatibility)

for (tissue in tissues) {

  path <- file.path(DATA_DIR, paste0(tissue, "_imputed_matrix.csv"))
  df   <- read.csv(path, row.names = 1, check.names = FALSE)
  mat  <- t(as.matrix(df))   # samples x proteins

  expr_full[[tissue]]   <- mat
  expr_t1[[tissue]]     <- mat[grepl("_T1", rownames(mat)), ]
  expr_t2[[tissue]]     <- mat[grepl("_T2", rownames(mat)), ]
  expr_t1vctl[[tissue]] <- mat[!grepl("_T2", rownames(mat)), ]

  cat(tissue, "— dimensions:\n")
  cat("  expr_full:   ", nrow(expr_full[[tissue]]),   "x", ncol(expr_full[[tissue]]),   "\n")
  cat("  expr_t1:     ", nrow(expr_t1[[tissue]]),     "x", ncol(expr_t1[[tissue]]),     "\n")
  cat("  expr_t2:     ", nrow(expr_t2[[tissue]]),     "x", ncol(expr_t2[[tissue]]),     "\n\n")

  # goodSamplesGenes check
  gsg <- goodSamplesGenes(expr_full[[tissue]], verbose = 0)
  if (!gsg$allOK) {
    expr_full[[tissue]] <- expr_full[[tissue]][gsg$goodSamples, gsg$goodGenes]
    message(sprintf("  %s: removed %d bad samples, %d bad proteins",
                    tissue,
                    sum(!gsg$goodSamples),
                    sum(!gsg$goodGenes)))
  }
}

# -----------------------------------------------------------------------------
# 2. Load metadata
# -----------------------------------------------------------------------------

meta        <- list()
meta_full   <- list()   # all samples
meta_t1     <- list()   # T1 only
meta_surv   <- list()   # SURV samples only (T1 + T2)

for (tissue in tissues) {

  meta_file <- ifelse(tissue == "muscle", "CRm_meta.csv", "CRp_meta.csv")

  meta[[tissue]] <- read.csv(
    file.path("00_input", tissue, meta_file),
    row.names    = 1,
    check.names  = FALSE
  )

  # Align to expr_full row order
  meta_full[[tissue]] <- meta[[tissue]][rownames(expr_full[[tissue]]), ]
  meta_t1[[tissue]]   <- meta[[tissue]][rownames(expr_t1[[tissue]]), ]

  # SURV subjects only (T1 + T2) — for timepoint and change analyses
  surv_rows <- rownames(expr_full[[tissue]])[
    meta_full[[tissue]]$cancer == "SURV"
  ]
  meta_surv[[tissue]] <- meta[[tissue]][surv_rows, ]

  # Derive cancer_time factor for LMM (3 levels)
  meta_full[[tissue]]$cancer_time <- factor(
    paste(meta_full[[tissue]]$cancer,
          ifelse(grepl("_T2", rownames(meta_full[[tissue]])), "T2", "T1"),
          sep = "_"),
    levels = c("CTL_T1", "SURV_T1", "SURV_T2")
  )

  # Timepoint within SURV
  meta_surv[[tissue]]$timepoint <- factor(
    ifelse(grepl("_T2", rownames(meta_surv[[tissue]])), "T2", "T1"),
    levels = c("T1", "T2")
  )

  cat(tissue, "— sample counts:\n")
  print(table(meta_full[[tissue]]$cancer_time))
  cat("\n")
}

# -----------------------------------------------------------------------------
# 3. Delta scores
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

meta_full <- lapply(meta_full, add_deltas)
meta_t1   <- lapply(meta_t1,   add_deltas)
meta_surv <- lapply(meta_surv, add_deltas)

# -----------------------------------------------------------------------------
# 4. Sample quality check
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  sampleTree <- hclust(dist(expr_full[[tissue]]), method = "average")
  png(file.path(REPORT_DIR, paste0(tissue, "_sample_tree.png")),
      width = 3000, height = 1200, res = 300)
  par(mar = c(0, 5, 3, 0))
  plot(sampleTree,
       main = paste(tissue, "— sample clustering (all samples)"),
       sub = "", xlab = "")
  dev.off()
}

# -----------------------------------------------------------------------------
# 5. Soft thresholding
# -----------------------------------------------------------------------------

sft_results <- list()
powers      <- c(1:10, seq(12, 30, 2))

for (tissue in tissues) {
  cat("\n--- Soft thresholding:", tissue, "---\n")

  sft_results[[tissue]] <- pickSoftThreshold(
    expr_full[[tissue]],
    powerVector = powers,
    networkType = WGCNA_NETWORK_TYPE,
    corFnc      = WGCNA_COR_TYPE,
    verbose     = 2
  )

  saveRDS(sft_results[[tissue]]$fitIndices,
          file.path(PANEL_DIR, paste0(tissue, "_sft_fitIndices.rds")))

  r2_vals   <- -sign(sft_results[[tissue]]$fitIndices$slope) *
                 sft_results[[tissue]]$fitIndices$SFT.R.sq
  power_idx <- which(r2_vals >= WGCNA_R2_CUTOFF)[1]
  auto_power <- if (!is.na(power_idx)) powers[power_idx] else NA

  message(sprintf("  %s: suggested power = %s (R² = %.3f)",
                  tissue,
                  ifelse(is.na(auto_power), "none found", auto_power),
                  ifelse(is.na(auto_power), NA, r2_vals[power_idx])))

  png(file.path(REPORT_DIR, paste0(tissue, "_SFT.png")),
      width = 3000, height = 1500, res = 300)
  par(mfrow = c(1, 2))

  plot(sft_results[[tissue]]$fitIndices[, 1], r2_vals,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (R²)",
       main = paste(tissue, "— Scale independence"), type = "n")
  text(sft_results[[tissue]]$fitIndices[, 1], r2_vals,
       labels = powers, cex = 0.9, col = "red")
  abline(h = WGCNA_R2_CUTOFF, col = "red", lty = 2)

  plot(sft_results[[tissue]]$fitIndices[, 1],
       sft_results[[tissue]]$fitIndices[, 5],
       xlab = "Soft Threshold (power)",
       ylab = "Mean Connectivity",
       main = paste(tissue, "— Mean connectivity"), type = "n")
  text(sft_results[[tissue]]$fitIndices[, 1],
       sft_results[[tissue]]$fitIndices[, 5],
       labels = powers, cex = 0.9, col = "red")
  dev.off()
}

# Update softPowers in Config after inspecting plots

# -----------------------------------------------------------------------------
# 6. Build networks
# -----------------------------------------------------------------------------

nets          <- list()
MEs           <- list()
module_colors <- list()

for (tissue in tissues) {
  cat("\n--- Building network:", tissue, "---\n")

  nets[[tissue]] <- blockwiseModules(
    expr_full[[tissue]],
    power             = softPowers[[tissue]],
    networkType       = WGCNA_NETWORK_TYPE,
    TOMType           = WGCNA_TOM_TYPE,
    corType           = WGCNA_COR_TYPE,
    minModuleSize     = WGCNA_MIN_MOD_SIZE,
    mergeCutHeight    = WGCNA_MERGE_CUT_H,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = file.path(DATA_DIR, paste0(tissue, "_TOM")),
    verbose           = 3
  )

  module_colors[[tissue]] <- nets[[tissue]]$colors
  n_modules <- length(unique(module_colors[[tissue]])) -
               ("grey" %in% module_colors[[tissue]])

  message(sprintf("  %s: %d modules (+ grey/unassigned)", tissue, n_modules))
  print(table(module_colors[[tissue]]))

  # Compute MEs on full dataset
  MEs[[tissue]] <- orderMEs(
    moduleEigengenes(expr_full[[tissue]],
                     colors = module_colors[[tissue]])$eigengenes
  )

  # Dendrogram
  png(file.path(REPORT_DIR, paste0(tissue, "_dendrogram.png")),
      width = 3000, height = 1500, res = 300)
  plotDendroAndColors(
    nets[[tissue]]$dendrograms[[1]],
    module_colors[[tissue]][nets[[tissue]]$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
    main = paste(tissue, "— Module Dendrogram")
  )
  dev.off()
}

# Save networks
for (tissue in tissues) {
  saveRDS(nets[[tissue]],          file.path(PANEL_DIR, paste0(tissue, "_network.rds")))
  saveRDS(MEs[[tissue]],           file.path(PANEL_DIR, paste0(tissue, "_MEs.rds")))
  saveRDS(module_colors[[tissue]], file.path(PANEL_DIR, paste0(tissue, "_module_colors.rds")))
  saveRDS(expr_full[[tissue]],     file.path(PANEL_DIR, paste0(tissue, "_datExpr.rds")))
}

# -----------------------------------------------------------------------------
# 7. kME (module membership)
# -----------------------------------------------------------------------------

kME_all <- list()

for (tissue in tissues) {
  kME_all[[tissue]] <- signedKME(expr_full[[tissue]], MEs[[tissue]])
  saveRDS(kME_all[[tissue]], file.path(PANEL_DIR, paste0(tissue, "_kME_all.rds")))
}

# -----------------------------------------------------------------------------
# 8. Hub protein extraction
# -----------------------------------------------------------------------------

hub_list <- list()

for (tissue in tissues) {
  hub_rows <- list()
  unique_mods <- setdiff(unique(module_colors[[tissue]]), "grey")

  for (mod in unique_mods) {
    mod_proteins <- names(module_colors[[tissue]])[module_colors[[tissue]] == mod]
    kme_col      <- paste0("kME", mod)
    if (!(kme_col %in% colnames(kME_all[[tissue]]))) next

    hub_rows[[mod]] <- tibble(
      protein  = rownames(kME_all[[tissue]]),
      kME      = kME_all[[tissue]][, kme_col]
    ) |>
      filter(protein %in% mod_proteins) |>
      arrange(desc(abs(kME))) |>
      mutate(module = mod, tissue = tissue)
  }

  hub_list[[tissue]] <- bind_rows(hub_rows)
  write_csv(hub_list[[tissue]],
            file.path(PANEL_DIR, paste0(tissue, "_hub_proteins.csv")))
}

# -----------------------------------------------------------------------------
# 9. LMM contrasts (cancer_time factor, 3 levels)
# -----------------------------------------------------------------------------
# Model: ME ~ cancer_time + (1 | pid)
# Reference: CTL_T1
# Contrasts via emmeans:
#   (1) SURV_T1 vs CTL_T1  — baseline dysregulation  (≈ Baseline_SURVvCTL)
#   (2) SURV_T2 vs SURV_T1 — training effect in SURV  (≈ Training_SURV)
#   (3) SURV_T2 vs CTL_T1  — post-training dysreg.    (≈ Training_SURVvCTL)

lmm_contrast_list <- list(
  Baseline_SURVvCTL  = c(-1,  1,  0),   # SURV_T1 - CTL_T1
  Training_SURV      = c( 0, -1,  1),   # SURV_T2 - SURV_T1
  Training_SURVvCTL  = c(-1,  0,  1)    # SURV_T2 - CTL_T1
)
# Level order: CTL_T1=1, SURV_T1=2, SURV_T2=3

lmm_results <- list()

for (tissue in tissues) {
  cat("\n--- LMM contrasts:", tissue, "---\n")

  all_mods  <- colnames(MEs[[tissue]])
  lmm_rows  <- list()
  lmm_data  <- as.data.frame(meta_full[[tissue]])
  lmm_data$sample_id   <- rownames(lmm_data)
  lmm_data$cancer_time <- factor(lmm_data$cancer_time,
                                  levels = c("CTL_T1", "SURV_T1", "SURV_T2"))
  # pid is subject ID — controls repeated measures random effect
  lmm_data$pid <- lmm_data$pid

  for (mod in all_mods) {
    lmm_data[["me_val"]] <- MEs[[tissue]][lmm_data$sample_id, mod]

    fit <- tryCatch(
      suppressWarnings(
        lmer(me_val ~ cancer_time + (1 | pid), data = lmm_data, REML = TRUE)
      ),
      error = function(e) {
        warning(sprintf("LMM failed for %s %s: %s", tissue, mod, e$message))
        NULL
      }
    )
    if (is.null(fit)) next

    singular <- isSingular(fit)
    if (singular) message(sprintf("  Singular fit: %s %s", tissue, mod))

    emm <- emmeans(fit, ~ cancer_time)

    for (cname in names(lmm_contrast_list)) {
      ctr <- contrast(emm, list(ctr = lmm_contrast_list[[cname]]))
      s   <- summary(ctr, ddf = "Kenward-Roger")

      tval   <- s$t.ratio
      df_k   <- s$df
      r_equiv <- sign(s$estimate) * sqrt(tval^2 / (tval^2 + df_k))

      lmm_rows <- c(lmm_rows, list(tibble(
        tissue   = tissue,
        module   = gsub("^ME", "", mod),
        contrast = cname,
        estimate = round(s$estimate, 5),
        SE       = round(s$SE, 5),
        df       = round(df_k, 2),
        t_ratio  = round(tval, 4),
        p_raw    = s$p.value,
        r_equiv  = round(r_equiv, 4),
        singular = singular
      )))
    }
  }

  lmm_results[[tissue]] <- bind_rows(lmm_rows)
  lmm_results[[tissue]]$p_bh <- p.adjust(lmm_results[[tissue]]$p_raw, method = "BH")

  message(sprintf("  %s: %d LMM tests (%d modules x %d contrasts)",
                  tissue, nrow(lmm_results[[tissue]]),
                  length(all_mods), length(lmm_contrast_list)))

  write_csv(lmm_results[[tissue]],
            file.path(PANEL_DIR, paste0(tissue, "_lmm_contrasts.csv")))
}

# -----------------------------------------------------------------------------
# 10. Pearson correlations — baseline traits (SURV T1 only)
# -----------------------------------------------------------------------------

baseline_trait_cols <- c(
  "age", "pre_weight_kg", "pre_BMI", "pre_ALM_kg",
  "pre_ALM_index", "pre_fat_pct"
)

bl_cor  <- list()
bl_pval <- list()

per_col_bh <- function(pmat) {
  bh_mat <- pmat
  for (j in seq_len(ncol(pmat))) bh_mat[, j] <- p.adjust(pmat[, j], method = "BH")
  bh_mat
}

for (tissue in tissues) {

  # SURV T1 only
  surv_t1_samples <- rownames(meta_t1[[tissue]])[
    meta_t1[[tissue]]$cancer == "SURV"
  ]
  MEs_surv_t1   <- MEs[[tissue]][surv_t1_samples, ]
  bl_traits_mat <- as.data.frame(meta_t1[[tissue]][surv_t1_samples, baseline_trait_cols])
  rownames(bl_traits_mat) <- surv_t1_samples

  bl_cor[[tissue]]  <- cor(MEs_surv_t1, bl_traits_mat, use = "pairwise.complete.obs")
  bl_pval_raw       <- matrix(NA_real_, nrow = nrow(bl_cor[[tissue]]),
                               ncol = ncol(bl_cor[[tissue]]),
                               dimnames = dimnames(bl_cor[[tissue]]))
  for (tr in colnames(bl_traits_mat)) {
    n_ok <- sum(!is.na(bl_traits_mat[, tr]))
    bl_pval_raw[, tr] <- corPvalueStudent(bl_cor[[tissue]][, tr], n_ok)
  }
  bl_pval[[tissue]] <- per_col_bh(bl_pval_raw)

  write_csv(as.data.frame(bl_cor[[tissue]])  |> rownames_to_column("module"),
            file.path(PANEL_DIR, paste0(tissue, "_baseline_trait_cor.csv")))
  write_csv(as.data.frame(bl_pval[[tissue]]) |> rownames_to_column("module"),
            file.path(PANEL_DIR, paste0(tissue, "_baseline_trait_pval_bh.csv")))
}

# -----------------------------------------------------------------------------
# 11. Pearson correlations — change traits (delta MEs vs delta phenotypes)
# -----------------------------------------------------------------------------
# delta ME = ME_T2 - ME_T1 per SURV subject
# Mirrors the reference script's delta_me approach

change_trait_cols <- c(
  "delta_ALM_kg", "delta_ALM_index", "delta_fat_pct",
  "delta_sts_max_pwr", "delta_grip_lbs",
  "delta_chest_press_lbs", "delta_leg_ext_lbs"
)

delta_ME  <- list()
ch_cor    <- list()
ch_pval   <- list()

for (tissue in tissues) {

  # Paired SURV subjects with both T1 and T2
  surv_meta <- meta_surv[[tissue]]
  t1_samples <- rownames(surv_meta)[grepl("_T1", rownames(surv_meta))]
  t2_samples <- rownames(surv_meta)[grepl("_T2", rownames(surv_meta))]

  # Strip timepoint to get subject IDs
  t1_pids <- sub("_T1$", "", t1_samples)
  t2_pids <- sub("_T2$", "", t2_samples)
  common_pids <- intersect(t1_pids, t2_pids)

  me_t1 <- MEs[[tissue]][paste0(common_pids, "_T1"), ]
  me_t2 <- MEs[[tissue]][paste0(common_pids, "_T2"), ]
  rownames(me_t1) <- common_pids
  rownames(me_t2) <- common_pids

  delta_ME[[tissue]] <- me_t2 - me_t1

  # Change traits from T1 metadata (delta scores already computed)
  surv_t1_meta <- meta_t1[[tissue]][
    meta_t1[[tissue]]$cancer == "SURV", change_trait_cols
  ]
  # Align to common_pids
  pid_to_t1_sample <- setNames(
    paste0(common_pids, "_T1"),
    common_pids
  )
  ch_traits_mat <- surv_t1_meta[pid_to_t1_sample[common_pids], , drop = FALSE]
  rownames(ch_traits_mat) <- common_pids

  ch_cor[[tissue]]  <- cor(delta_ME[[tissue]], ch_traits_mat,
                            use = "pairwise.complete.obs")
  ch_pval_raw <- matrix(NA_real_, nrow = nrow(ch_cor[[tissue]]),
                         ncol = ncol(ch_cor[[tissue]]),
                         dimnames = dimnames(ch_cor[[tissue]]))
  for (tr in colnames(ch_traits_mat)) {
    n_ok <- sum(!is.na(ch_traits_mat[, tr]))
    ch_pval_raw[, tr] <- corPvalueStudent(ch_cor[[tissue]][, tr], n_ok)
  }
  ch_pval[[tissue]] <- per_col_bh(ch_pval_raw)

  # Save
  saveRDS(delta_ME[[tissue]], file.path(PANEL_DIR, paste0(tissue, "_delta_ME.rds")))
  write_csv(as.data.frame(ch_cor[[tissue]])  |> rownames_to_column("module"),
            file.path(PANEL_DIR, paste0(tissue, "_change_trait_cor.csv")))
  write_csv(as.data.frame(ch_pval[[tissue]]) |> rownames_to_column("module"),
            file.path(PANEL_DIR, paste0(tissue, "_change_trait_pval_bh.csv")))

  message(sprintf("  %s: %d paired SURV subjects for delta ME",
                  tissue, length(common_pids)))
}

# Restore stats::cor
cor <- stats::cor

# -----------------------------------------------------------------------------
# 12. Combined heatmap — LMM r-equivalents + Pearson correlations
# -----------------------------------------------------------------------------

star_mat <- function(pmat) {
  ifelse(pmat < 0.001, "***",
  ifelse(pmat < 0.01,  "**",
  ifelse(pmat < 0.05,  "*", "")))
}

for (tissue in tissues) {

  # Pivot LMM r-equivalents to matrix (modules x contrasts)
  lmm_r_wide <- lmm_results[[tissue]] |>
    dplyr::select(module, contrast, r_equiv) |>
    pivot_wider(names_from = contrast, values_from = r_equiv) |>
    column_to_rownames("module")
  lmm_r_wide <- lmm_r_wide[gsub("^ME", "", colnames(MEs[[tissue]])), ]

  lmm_p_wide <- lmm_results[[tissue]] |>
    dplyr::select(module, contrast, p_bh) |>
    pivot_wider(names_from = contrast, values_from = p_bh) |>
    column_to_rownames("module")
  lmm_p_wide <- lmm_p_wide[gsub("^ME", "", colnames(MEs[[tissue]])), ]

  # Rename columns for display
  colnames(lmm_r_wide) <- c("SURV_T1\nvs CTL_T1",
                              "SURV_T2\nvs SURV_T1",
                              "SURV_T2\nvs CTL_T1")
  colnames(lmm_p_wide) <- colnames(lmm_r_wide)

  # Align Pearson matrices to same module order
  mod_order <- rownames(lmm_r_wide)
  bl_r <- bl_cor[[tissue]][paste0("ME", mod_order), , drop = FALSE]
  bl_p <- bl_pval[[tissue]][paste0("ME", mod_order), , drop = FALSE]
  ch_r <- ch_cor[[tissue]][paste0("ME", mod_order), , drop = FALSE]
  ch_p <- ch_pval[[tissue]][paste0("ME", mod_order), , drop = FALSE]
  rownames(bl_r) <- rownames(bl_p) <- mod_order
  rownames(ch_r) <- rownames(ch_p) <- mod_order

  # Combine: LMM | baseline Pearson | change Pearson
  combined_r <- as.matrix(cbind(lmm_r_wide, bl_r, ch_r))
  combined_p <- as.matrix(cbind(lmm_p_wide, bl_p, ch_p))

  text_mat <- paste(signif(combined_r, 2), star_mat(combined_p), sep = "\n")
  dim(text_mat) <- dim(combined_r)

  n_lmm  <- ncol(lmm_r_wide)
  n_bl   <- ncol(bl_r)
  n_ch   <- ncol(ch_r)

  png(file.path(PANEL_DIR, paste0(tissue, "_module_trait_heatmap_combined.png")),
      width = 4500, height = 3000, res = 300)
  par(mar = c(10, 12, 3, 3))

  labeledHeatmap(
    Matrix        = combined_r,
    xLabels       = colnames(combined_r),
    yLabels       = rownames(combined_r),
    colorLabels   = FALSE,
    colors        = blueWhiteRed(50),
    textMatrix    = text_mat,
    setStdMargins = FALSE,
    cex.text      = 0.5,
    cex.lab       = 0.65,
    zlim          = c(-1, 1),
    main          = paste(tissue,
                          "— Module-Trait (LMM r-equiv | Baseline Pearson | Change Pearson)\n* BH<0.05  ** BH<0.01  *** BH<0.001")
  )

  # Vertical lines separating the three sections
  abline(v = n_lmm + 0.5,         lwd = 2, lty = 2, col = "black")
  abline(v = n_lmm + n_bl + 0.5,  lwd = 2, lty = 2, col = "black")

  dev.off()
}

# -----------------------------------------------------------------------------
# 13. Reversal module candidates (from LMM results)
# -----------------------------------------------------------------------------
# Reversal: elevated at baseline (Baseline_SURVvCTL r > 0)
#           AND decreases with training (Training_SURV r < 0)
#           AND resolved post-training (Training_SURVvCTL near zero / non-sig)

reversal_modules <- list()

for (tissue in tissues) {

  df_wide <- lmm_results[[tissue]] |>
    dplyr::select(module, contrast, r_equiv, p_bh) |>
    pivot_wider(names_from  = contrast,
                values_from = c(r_equiv, p_bh))

  reversal_modules[[tissue]] <- df_wide |>
    mutate(
      reversal_score = r_equiv_Baseline_SURVvCTL - r_equiv_Training_SURV,
      reversal_flag  = r_equiv_Baseline_SURVvCTL >  0.2 &
                       r_equiv_Training_SURV      < -0.2
    ) |>
    arrange(desc(reversal_score))

  cat("\n--- Reversal module candidates:", tissue, "---\n")
  print(reversal_modules[[tissue]] |> filter(reversal_flag))

  write_csv(reversal_modules[[tissue]],
            file.path(PANEL_DIR, paste0(tissue, "_reversal_module_candidates.csv")))
}

# -----------------------------------------------------------------------------
# 14. kME-logFC integration with limma contrast results
# -----------------------------------------------------------------------------

library(readxl)

limma_files <- list(
  muscle = file.path("03_analysis", "b_reports", "muscle", "reformat", "all_data", "limma", "limma_results_all_muscle.xlsx"),
  plasma = file.path("03_analysis", "b_reports", "plasma", "reformat", "all_data", "limma", "limma_results_all_plasma.xlsx")
)

# Sheet indices and their corresponding contrast names
contrast_sheets <- list(
  Baseline_SURVvCTL = 11,
  Training_SURVvCTL = 12,
  Training_SURV     = 13
)

names(read_xlsx(limma_files[["muscle"]], sheet = 11))

kME_lfc <- list()

for (tissue in tissues) {
  
  if (!file.exists(limma_files[[tissue]])) {
    message("Limma file not found for ", tissue, " — skipping kME-logFC integration")
    next
  }
  
  kME_mat <- as.matrix(kME_all[[tissue]])
  rows    <- list()
  
  for (ct in names(contrast_sheets)) {
    
    sheet_dat <- read_xlsx(limma_files[[tissue]], sheet = contrast_sheets[[ct]])
    
    # Confirm expected columns exist — adjust names if needed
    # Expected: a UniProt accession column and a logFC column
    # Check with: names(read_xlsx(limma_files[["muscle"]], sheet = 11))
    lfc_df <- sheet_dat |>
      dplyr::select(uniprot, logFC) |>
      distinct(uniprot, .keep_all = TRUE) |>
      filter(!is.na(uniprot), !is.na(logFC))
    
    lfc_vec <- setNames(lfc_df$logFC, lfc_df$uniprot)
    shared  <- intersect(rownames(kME_mat), names(lfc_vec))
    
    if (length(shared) < 10) {
      warning(sprintf("< 10 shared proteins for %s %s", tissue, ct))
      next
    }
    
    kME_sub  <- kME_mat[shared, , drop = FALSE]
    lfc_sub  <- lfc_vec[shared]
    
    cor_vals <- apply(kME_sub, 2, function(me) {
      cor(me, lfc_sub, method = "pearson", use = "complete.obs")
    })
    
    rows[[ct]] <- tibble(
      tissue   = tissue,
      module   = gsub("^kME", "", names(cor_vals)),
      contrast = ct,
      r        = cor_vals
    )
  }
  
  kME_lfc[[tissue]] <- bind_rows(rows)
  
  cat("\n--- kME-logFC top correlations:", tissue, "---\n")
  print(
    kME_lfc[[tissue]] |>
      arrange(contrast, desc(abs(r))) |>
      group_by(contrast) |>
      slice_head(n = 5)
  )
  
  write_csv(kME_lfc[[tissue]],
            file.path(PANEL_DIR, paste0(tissue, "_kME_logFC_correlations.csv")))
}

# -----------------------------------------------------------------------------
# 15. SFT summary
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  fi <- sft_results[[tissue]]$fitIndices
  r2_vals <- -sign(fi$slope) * fi$SFT.R.sq

  sft_summary <- tibble(
    tissue            = tissue,
    selected_power    = softPowers[[tissue]],
    R_squared         = r2_vals[softPowers[[tissue]]],
    mean_connectivity = fi$mean.k.[softPowers[[tissue]]],
    n_proteins        = ncol(expr_full[[tissue]]),
    n_samples         = nrow(expr_full[[tissue]])
  )
  write_csv(sft_summary,
            file.path(PANEL_DIR, paste0(tissue, "_sft_summary.csv")))
}

message("\nDone.")
