# =============================================================================
# WGCNA Script 2 — ME-Trait Associations
# =============================================================================
# Loads the saved network from Script 1 and runs three targeted analyses.
#
# Design rationale:
#   CTL has only T1; SURV has T1 and T2. The three-level cancer_time factor
#   (CTL_T1 / SURV_T1 / SURV_T2) maps directly onto limma contrasts:
#     SURV_T1 vs CTL_T1  = Baseline_SURVvCTL
#     SURV_T2 vs SURV_T1 = Training_SURV
#     SURV_T2 vs CTL_T1  = Training_SURVvCTL
#
# Analysis 1 — LMM contrasts (cancer_time factor, emmeans)
#   Model: ME ~ cancer_time + (1 | pid)
#   Random effect accounts for repeated measures in SURV (T1 + T2 per subject)
#   CTL contributes as independent observations (single timepoint, no pairing)
#   r-equivalents extracted via Kenward-Roger df (lmerTest + emmeans)
#
# Analysis 2 — Pearson: baseline phenotypes (SURV_T1 MEs vs pre_ traits)
#   Only SURV_T1 samples used; CTL excluded (no training response expected)
#
# Analysis 3 — Pearson: functional outcomes (SURV_T1 MEs vs delta traits)
#   Delta ME approach: ME_T2 - ME_T1 per subject, correlated with delta phenotype
#
# Analysis 4 — kME-logFC integration with limma results
#
# Outputs saved to 04_figures/WGCNA/
# =============================================================================

library(WGCNA)
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(readxl)

options(stringsAsFactors = FALSE)
set.seed(26)

setwd(rprojroot::find_rstudio_root_file())

# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------

tissues <- c("muscle", "plasma")
OUT_DIR <- "04_figures/WGCNA"

# limma contrast sheets (adjust paths/sheet numbers to match your files)
limma_files <- list(
  muscle = file.path("03_analysis", "b_reports", "muscle", "reformat",
                     "all_data", "limma", "limma_results_all_muscle.xlsx"),
  plasma = file.path("03_analysis", "b_reports", "plasma", "reformat",
                     "all_data", "limma", "limma_results_all_plasma.xlsx")
)
contrast_sheets <- list(
  Baseline_SURVvCTL = 11,
  Training_SURVvCTL = 12,
  Training_SURV     = 13
)

# Baseline phenotype traits (SURV_T1 only)
baseline_trait_cols <- c(
  "age", "pre_weight_kg", "pre_BMI",
  "pre_ALM_kg", "pre_ALM_index", "pre_fat_pct"
)

# Functional outcome traits (delta scores — computed from T1 metadata)
delta_trait_cols <- c(
  "delta_ALM_kg", "delta_ALM_index", "delta_fat_pct",
  "delta_sts_max_pwr", "delta_grip_lbs",
  "delta_chest_press_lbs", "delta_leg_ext_lbs"
)

# LMM contrasts — numeric vectors over (CTL_T1, SURV_T1, SURV_T2)
lmm_contrast_list <- list(
  Baseline_SURVvCTL = c(-1,  1,  0),
  Training_SURV     = c( 0, -1,  1),
  Training_SURVvCTL = c(-1,  0,  1)
)

# -----------------------------------------------------------------------------
# 1. Load network outputs from Script 1
# -----------------------------------------------------------------------------

MEs           <- list()
module_colors <- list()
kME_all       <- list()
datExpr       <- list()
meta_full     <- list()

for (tissue in tissues) {
  MEs[[tissue]]           <- readRDS(file.path(OUT_DIR, paste0(tissue, "_MEs.rds")))
  module_colors[[tissue]] <- readRDS(file.path(OUT_DIR, paste0(tissue, "_module_colors.rds")))
  kME_all[[tissue]]       <- readRDS(file.path(OUT_DIR, paste0(tissue, "_kME.rds")))
  datExpr[[tissue]]       <- readRDS(file.path(OUT_DIR, paste0(tissue, "_datExpr.rds")))
}

meta_full <- readRDS(file.path(OUT_DIR, "meta_full.rds"))

# -----------------------------------------------------------------------------
# 2. Prepare metadata subsets
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

# SURV_T1 only (for baseline and delta trait correlations)
meta_surv_t1 <- lapply(tissues, function(tissue) {
  m <- meta_full[[tissue]]
  m[m$cancer == "SURV" & grepl("_T1", rownames(m)), ]
})
names(meta_surv_t1) <- tissues

# -----------------------------------------------------------------------------
# 3. Analysis 1 — LMM contrasts via emmeans (cancer_time factor)
# -----------------------------------------------------------------------------
# Model: ME ~ cancer_time + (1 | pid)
# Random effect on pid accounts for repeated SURV measures (T1 + T2)
# CTL contributes as independent single-timepoint observations
# Contrasts extracted via Kenward-Roger df for small-n validity

lmm_results <- list()

for (tissue in tissues) {
  cat(sprintf("\n--- LMM contrasts: %s ---\n", tissue))

  all_mods  <- colnames(MEs[[tissue]])
  lmm_data  <- as.data.frame(meta_full[[tissue]])
  lmm_data$sample_id   <- rownames(lmm_data)
  lmm_data$cancer_time <- factor(lmm_data$cancer_time,
                                  levels = c("CTL_T1", "SURV_T1", "SURV_T2"))

  lmm_rows <- list()

  for (mod in all_mods) {
    lmm_data[["me_val"]] <- MEs[[tissue]][lmm_data$sample_id, mod]

    fit <- tryCatch(
      suppressWarnings(
        lmer(me_val ~ cancer_time + (1 | pid),
             data = lmm_data, REML = TRUE)
      ),
      error = function(e) {
        warning(sprintf("LMM failed: %s %s: %s", tissue, mod, e$message))
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

      tval    <- s$t.ratio
      df_k    <- s$df
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
  lmm_results[[tissue]]$p_bh <- p.adjust(lmm_results[[tissue]]$p_raw,
                                           method = "BH")

  cat(sprintf("  %d LMM tests (%d modules x %d contrasts)\n",
              nrow(lmm_results[[tissue]]),
              length(all_mods),
              length(lmm_contrast_list)))

  write_csv(lmm_results[[tissue]],
            file.path(OUT_DIR, paste0(tissue, "_lmm_contrasts.csv")))
}

# -----------------------------------------------------------------------------
# 4. Analysis 2 — Pearson: baseline phenotypes (SURV_T1 only)
# -----------------------------------------------------------------------------

per_col_bh <- function(pmat) {
  bh_mat <- pmat
  for (j in seq_len(ncol(pmat)))
    bh_mat[, j] <- p.adjust(pmat[, j], method = "BH")
  bh_mat
}

bl_cor  <- list()
bl_pval <- list()

for (tissue in tissues) {

  surv_t1_samples <- rownames(meta_surv_t1[[tissue]])
  MEs_surv_t1     <- MEs[[tissue]][surv_t1_samples, ]

  traits_mat <- as.data.frame(
    meta_surv_t1[[tissue]][, intersect(baseline_trait_cols,
                                       names(meta_surv_t1[[tissue]]))]
  )
  rownames(traits_mat) <- surv_t1_samples

  bl_cor[[tissue]]  <- cor(MEs_surv_t1, traits_mat,
                            use = "pairwise.complete.obs")
  praw <- matrix(NA_real_, nrow = nrow(bl_cor[[tissue]]),
                 ncol = ncol(bl_cor[[tissue]]),
                 dimnames = dimnames(bl_cor[[tissue]]))
  for (tr in colnames(traits_mat)) {
    n_ok <- sum(!is.na(traits_mat[, tr]))
    praw[, tr] <- corPvalueStudent(bl_cor[[tissue]][, tr], n_ok)
  }
  bl_pval[[tissue]] <- per_col_bh(praw)

  write_csv(as.data.frame(bl_cor[[tissue]])  |> rownames_to_column("module"),
            file.path(OUT_DIR, paste0(tissue, "_baseline_trait_cor.csv")))
  write_csv(as.data.frame(bl_pval[[tissue]]) |> rownames_to_column("module"),
            file.path(OUT_DIR, paste0(tissue, "_baseline_trait_pval_bh.csv")))
}

# -----------------------------------------------------------------------------
# 5. Analysis 3 — Delta ME vs delta phenotypes (paired SURV subjects)
# -----------------------------------------------------------------------------
# Delta ME = ME_T2 - ME_T1 per SURV subject
# Correlate with delta phenotype traits from T1 metadata

delta_ME  <- list()
ch_cor    <- list()
ch_pval   <- list()

for (tissue in tissues) {

  all_samples  <- rownames(MEs[[tissue]])
  surv_t1_samp <- all_samples[
    meta_full[[tissue]]$cancer == "SURV" & grepl("_T1", all_samples)
  ]
  surv_t2_samp <- all_samples[
    meta_full[[tissue]]$cancer == "SURV" & grepl("_T2", all_samples)
  ]

  t1_pids     <- sub("_T1$", "", surv_t1_samp)
  t2_pids     <- sub("_T2$", "", surv_t2_samp)
  common_pids <- intersect(t1_pids, t2_pids)

  if (length(common_pids) == 0) {
    warning(sprintf("No paired SURV subjects found for %s — skipping delta ME", tissue))
    next
  }

  me_t1 <- MEs[[tissue]][paste0(common_pids, "_T1"), ]
  me_t2 <- MEs[[tissue]][paste0(common_pids, "_T2"), ]
  rownames(me_t1) <- common_pids
  rownames(me_t2) <- common_pids

  delta_ME[[tissue]] <- me_t2 - me_t1

  # Delta traits from SURV_T1 metadata
  t1_meta_surv <- meta_surv_t1[[tissue]]
  pid_to_t1    <- setNames(paste0(common_pids, "_T1"), common_pids)
  avail_delta  <- intersect(delta_trait_cols, names(t1_meta_surv))
  ch_traits    <- t1_meta_surv[pid_to_t1[common_pids], avail_delta, drop = FALSE]
  rownames(ch_traits) <- common_pids

  ch_cor[[tissue]]  <- cor(delta_ME[[tissue]], ch_traits,
                            use = "pairwise.complete.obs")
  praw_ch <- matrix(NA_real_, nrow = nrow(ch_cor[[tissue]]),
                    ncol = ncol(ch_cor[[tissue]]),
                    dimnames = dimnames(ch_cor[[tissue]]))
  for (tr in colnames(ch_traits)) {
    n_ok <- sum(!is.na(ch_traits[, tr]))
    praw_ch[, tr] <- corPvalueStudent(ch_cor[[tissue]][, tr], n_ok)
  }
  ch_pval[[tissue]] <- per_col_bh(praw_ch)

  saveRDS(delta_ME[[tissue]], file.path(OUT_DIR, paste0(tissue, "_delta_ME.rds")))
  write_csv(as.data.frame(ch_cor[[tissue]])  |> rownames_to_column("module"),
            file.path(OUT_DIR, paste0(tissue, "_delta_trait_cor.csv")))
  write_csv(as.data.frame(ch_pval[[tissue]]) |> rownames_to_column("module"),
            file.path(OUT_DIR, paste0(tissue, "_delta_trait_pval_bh.csv")))

  cat(sprintf("%s: %d paired SURV subjects for delta ME analysis\n",
              tissue, length(common_pids)))
}

# -----------------------------------------------------------------------------
# 6. Combined heatmap — LMM r-equivalents + Pearson (saved as PDF)
# -----------------------------------------------------------------------------

star_mat <- function(pmat) {
  ifelse(pmat < 0.001, "***",
  ifelse(pmat < 0.01,  "**",
  ifelse(pmat < 0.05,  "*", "")))
}

for (tissue in tissues) {

  # LMM: pivot to modules x contrasts matrix
  lmm_r_wide <- lmm_results[[tissue]] |>
    dplyr::select(module, contrast, r_equiv) |>
    pivot_wider(names_from = contrast, values_from = r_equiv) |>
    column_to_rownames("module")

  lmm_p_wide <- lmm_results[[tissue]] |>
    dplyr::select(module, contrast, p_bh) |>
    pivot_wider(names_from = contrast, values_from = p_bh) |>
    column_to_rownames("module")

  # Ensure consistent module order
  mod_order <- gsub("^ME", "", colnames(MEs[[tissue]]))
  lmm_r_wide <- lmm_r_wide[mod_order, names(lmm_contrast_list), drop = FALSE]
  lmm_p_wide <- lmm_p_wide[mod_order, names(lmm_contrast_list), drop = FALSE]

  colnames(lmm_r_wide) <- c("SURV_T1\nvs CTL_T1",
                              "SURV_T2\nvs SURV_T1",
                              "SURV_T2\nvs CTL_T1")
  colnames(lmm_p_wide) <- colnames(lmm_r_wide)

  # Pearson matrices aligned to same module order
  me_rownames <- paste0("ME", mod_order)
  bl_r <- bl_cor[[tissue]][me_rownames,  , drop = FALSE]
  bl_p <- bl_pval[[tissue]][me_rownames, , drop = FALSE]
  ch_r <- ch_cor[[tissue]][me_rownames,  , drop = FALSE]
  ch_p <- ch_pval[[tissue]][me_rownames, , drop = FALSE]
  rownames(bl_r) <- rownames(bl_p) <- mod_order
  rownames(ch_r) <- rownames(ch_p) <- mod_order

  combined_r <- as.matrix(cbind(lmm_r_wide, bl_r, ch_r))
  combined_p <- as.matrix(cbind(lmm_p_wide, bl_p, ch_p))

  text_mat <- paste(signif(combined_r, 2), star_mat(combined_p), sep = "\n")
  dim(text_mat) <- dim(combined_r)

  n_lmm <- ncol(lmm_r_wide)
  n_bl  <- ncol(bl_r)

  png(file.path(OUT_DIR, paste0(tissue, "_module_trait_heatmap_combined.png")),
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
      "— Module-Trait (LMM r-equiv | Baseline Pearson | Delta Pearson)\n* BH<0.05  ** BH<0.01  *** BH<0.001")
  )

  abline(v = n_lmm + 0.5,        lwd = 2, lty = 2, col = "black")
  abline(v = n_lmm + n_bl + 0.5, lwd = 2, lty = 2, col = "black")
  dev.off()
}

# -----------------------------------------------------------------------------
# 7. Analysis 4 — kME-logFC integration with limma results
# -----------------------------------------------------------------------------

kME_lfc <- list()

for (tissue in tissues) {

  if (!file.exists(limma_files[[tissue]])) {
    message(sprintf("Limma file not found for %s — skipping kME-logFC", tissue))
    next
  }

  kME_mat  <- as.matrix(kME_all[[tissue]])
  rows     <- list()

  for (ct in names(contrast_sheets)) {
    sheet_dat <- read_xlsx(limma_files[[tissue]],
                           sheet = contrast_sheets[[ct]])

    # Adjust column names to match your limma output
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

    cor_vals <- apply(kME_mat[shared, , drop = FALSE], 2, function(me) {
      cor(me, lfc_vec[shared], method = "pearson", use = "complete.obs")
    })

    rows[[ct]] <- tibble(
      tissue   = tissue,
      module   = gsub("^kME", "", names(cor_vals)),
      contrast = ct,
      r        = cor_vals
    )
  }

  kME_lfc[[tissue]] <- bind_rows(rows)
  write_csv(kME_lfc[[tissue]],
            file.path(OUT_DIR, paste0(tissue, "_kME_logFC_correlations.csv")))
}

# Reversal module candidates
reversal_modules <- lapply(tissues, function(tissue) {
  df <- kME_lfc[[tissue]]
  if (is.null(df) || nrow(df) == 0) return(NULL)

  baseline <- df |> filter(contrast == "Baseline_SURVvCTL") |>
    dplyr::select(module, r_baseline = r)
  training <- df |> filter(contrast == "Training_SURV") |>
    dplyr::select(module, r_training = r)

  inner_join(baseline, training, by = "module") |>
    mutate(
      reversal_score = r_baseline - r_training,
      reversal_flag  = r_baseline > 0.2 & r_training < -0.2
    ) |>
    arrange(desc(reversal_score))
})
names(reversal_modules) <- tissues

for (tissue in tissues) {
  if (!is.null(reversal_modules[[tissue]])) {
    write_csv(reversal_modules[[tissue]],
              file.path(OUT_DIR, paste0(tissue, "_reversal_modules.csv")))
    cat(sprintf("\n%s reversal module candidates:\n", tissue))
    print(reversal_modules[[tissue]] |> filter(reversal_flag))
  }
}

# Save all trait results for Script 3
saveRDS(list(
  lmm_results      = lmm_results,
  bl_cor           = bl_cor,
  bl_pval          = bl_pval,
  ch_cor           = ch_cor,
  ch_pval          = ch_pval,
  delta_ME         = delta_ME,
  kME_lfc          = kME_lfc,
  reversal_modules = reversal_modules
), file.path(OUT_DIR, "trait_results.rds"))

message("\nScript 2 complete. Proceed to wgcna_03_figures.R")
