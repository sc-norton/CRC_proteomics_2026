# =============================================================================
# WGCNA Script 1 — Network Construction
# =============================================================================
# Builds co-expression networks on the full dataset (SURV_T1 + SURV_T2 + CTL_T1)
# for muscle and plasma separately.
#
# Outputs (saved to 04_figures/WGCNA/):
#   - {tissue}_network.rds       — blockwiseModules output
#   - {tissue}_MEs.rds           — ordered module eigengenes (all samples)
#   - {tissue}_module_colors.rds — named vector of module assignments
#   - {tissue}_kME.rds           — signed kME matrix
#   - {tissue}_datExpr.rds       — expression matrix used for network
#   - {tissue}_sft_summary.csv   — soft threshold selection summary
#
# Run once; downstream scripts load these outputs.
# =============================================================================

library(WGCNA)
library(tidyverse)

options(stringsAsFactors = FALSE)
disableWGCNAThreads()   # single-threaded = byte-reproducible TOM
set.seed(26)

setwd(rprojroot::find_rstudio_root_file())
getwd()
# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------

tissues <- c("muscle", "plasma")

# Update after inspecting SFT plots (section 4)
softPowers <- c(muscle = 3, plasma = 3)

WGCNA_NETWORK_TYPE <- "signed hybrid"
WGCNA_TOM_TYPE     <- "signed"
WGCNA_COR_TYPE     <- "bicor"
WGCNA_MIN_MOD_SIZE <- 15L
WGCNA_MERGE_CUT_H  <- 0.25
WGCNA_R2_CUTOFF    <- 0.80

OUT_DIR <- "04_figures/WGCNA"
TOM_DIR <- "02_imputation/c_data"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Override cor() with WGCNA version during network computations
cor <- WGCNA::cor

# -----------------------------------------------------------------------------
# 1. Load imputed matrices
# -----------------------------------------------------------------------------

datExpr <- list()

for (tissue in tissues) {

  df  <- read.csv(
    file.path("02_imputation", "c_data",
              paste0(tissue, "_imputed_matrix.csv")),
    row.names    = 1,
    check.names  = FALSE
  )

  # Transpose: samples x proteins (WGCNA convention)
  mat <- t(as.matrix(df))

  # goodSamplesGenes check
  gsg <- goodSamplesGenes(mat, verbose = 0)
  if (!gsg$allOK) {
    mat <- mat[gsg$goodSamples, gsg$goodGenes]
    message(sprintf("  %s: removed %d bad samples, %d bad proteins after goodSamplesGenes",
                    tissue, sum(!gsg$goodSamples), sum(!gsg$goodGenes)))
  }

  datExpr[[tissue]] <- mat

  cat(sprintf("%s: %d samples x %d proteins\n",
              tissue, nrow(mat), ncol(mat)))
}

# -----------------------------------------------------------------------------
# 2. Load metadata and attach cancer_time factor
# -----------------------------------------------------------------------------

meta_full <- list()

for (tissue in tissues) {

  meta_file <- ifelse(tissue == "muscle", "CRm_meta.csv", "CRp_meta.csv")
  m <- read.csv(file.path("00_input", tissue, meta_file),
                row.names = 1, check.names = FALSE)

  # Align to datExpr row order
  m <- m[rownames(datExpr[[tissue]]), ]

  # cancer_time: three-level factor mapping directly to limma contrasts
  m$cancer_time <- factor(
    paste(m$cancer,
          ifelse(grepl("_T2", rownames(m)), "T2", "T1"),
          sep = "_"),
    levels = c("CTL_T1", "SURV_T1", "SURV_T2")
  )

  meta_full[[tissue]] <- m

  cat(sprintf("%s sample breakdown:\n", tissue))
  print(table(m$cancer_time))
  cat("\n")
}

# Save metadata for downstream scripts
saveRDS(meta_full, file.path(OUT_DIR, "meta_full.rds"))

# -----------------------------------------------------------------------------
# 3. Sample clustering (quality check)
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  sampleTree <- hclust(dist(datExpr[[tissue]]), method = "average")

  png(file.path(OUT_DIR, paste0(tissue, "_sample_tree.png")),
      width = 3000, height = 1200, res = 300)
  par(mar = c(0, 5, 3, 0))
  plot(sampleTree,
       main = paste(tissue, "— sample clustering (all samples)"),
       sub = "", xlab = "")
  dev.off()
}

# -----------------------------------------------------------------------------
# 4. Soft thresholding
# -----------------------------------------------------------------------------

sft_results <- list()
powers      <- c(1:10, seq(12, 30, 2))

for (tissue in tissues) {
  cat(sprintf("\n--- Soft thresholding: %s ---\n", tissue))

  sft_results[[tissue]] <- pickSoftThreshold(
    datExpr[[tissue]],
    powerVector = powers,
    networkType = WGCNA_NETWORK_TYPE,
    corFnc      = WGCNA_COR_TYPE,
    verbose     = 2
  )

  r2_vals   <- -sign(sft_results[[tissue]]$fitIndices$slope) *
                 sft_results[[tissue]]$fitIndices$SFT.R.sq
  auto_pow  <- powers[which(r2_vals >= WGCNA_R2_CUTOFF)[1]]

  message(sprintf("  %s: auto-suggested power = %s (R² = %.3f at selected power %d)",
                  tissue,
                  ifelse(is.na(auto_pow), "none found", auto_pow),
                  r2_vals[softPowers[[tissue]]],
                  softPowers[[tissue]]))

  # SFT plots
  png(file.path(OUT_DIR, paste0(tissue, "_SFT.png")),
      width = 3000, height = 1500, res = 300)
  par(mfrow = c(1, 2))

  plot(sft_results[[tissue]]$fitIndices[, 1], r2_vals,
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology R²",
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

  # Save SFT summary
  write_csv(
    tibble(
      tissue            = tissue,
      selected_power    = softPowers[[tissue]],
      R_squared         = r2_vals[softPowers[[tissue]]],
      mean_connectivity = sft_results[[tissue]]$fitIndices$mean.k.[softPowers[[tissue]]],
      n_proteins        = ncol(datExpr[[tissue]]),
      n_samples         = nrow(datExpr[[tissue]])
    ),
    file.path(OUT_DIR, paste0(tissue, "_sft_summary.csv"))
  )
}

# >>> STOP HERE and inspect SFT plots before proceeding <<<
# >>> Update softPowers in Config (section 0) then rerun from section 5 <<<

# -----------------------------------------------------------------------------
# 5. Build networks
# -----------------------------------------------------------------------------

nets          <- list()
MEs           <- list()
module_colors <- list()
kME_all       <- list()

for (tissue in tissues) {
  cat(sprintf("\n--- Building network: %s (power = %d) ---\n",
              tissue, softPowers[[tissue]]))

  nets[[tissue]] <- blockwiseModules(
    datExpr[[tissue]],
    power             = softPowers[[tissue]],
    networkType       = WGCNA_NETWORK_TYPE,
    TOMType           = WGCNA_TOM_TYPE,
    corType           = WGCNA_COR_TYPE,
    minModuleSize     = WGCNA_MIN_MOD_SIZE,
    mergeCutHeight    = WGCNA_MERGE_CUT_H,
    numericLabels     = FALSE,
    pamRespectsDendro = FALSE,
    saveTOMs          = TRUE,
    saveTOMFileBase   = file.path(TOM_DIR, paste0(tissue, "_TOM")),
    verbose           = 3
  )

  module_colors[[tissue]] <- nets[[tissue]]$colors
  n_mods <- length(unique(module_colors[[tissue]])) -
            ("grey" %in% module_colors[[tissue]])

  cat(sprintf("\n%s: %d modules detected (+ grey/unassigned)\n", tissue, n_mods))
  print(table(module_colors[[tissue]]))

  # Compute MEs on full dataset
  MEs[[tissue]] <- orderMEs(
    moduleEigengenes(datExpr[[tissue]],
                     colors = module_colors[[tissue]])$eigengenes
  )

  # Signed kME for all proteins
  kME_all[[tissue]] <- signedKME(datExpr[[tissue]], MEs[[tissue]])

  # Dendrogram
  png(file.path(OUT_DIR, paste0(tissue, "_dendrogram.png")),
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

# -----------------------------------------------------------------------------
# 6. Save all network outputs
# -----------------------------------------------------------------------------

for (tissue in tissues) {
  saveRDS(nets[[tissue]],          file.path(OUT_DIR, paste0(tissue, "_network.rds")))
  saveRDS(MEs[[tissue]],           file.path(OUT_DIR, paste0(tissue, "_MEs.rds")))
  saveRDS(module_colors[[tissue]], file.path(OUT_DIR, paste0(tissue, "_module_colors.rds")))
  saveRDS(kME_all[[tissue]],       file.path(OUT_DIR, paste0(tissue, "_kME.rds")))
  saveRDS(datExpr[[tissue]],       file.path(OUT_DIR, paste0(tissue, "_datExpr.rds")))
}

# Module assignment table (for joining to protein annotation)
for (tissue in tissues) {
  write_csv(
    tibble(
      protein       = names(module_colors[[tissue]]),
      module_color  = module_colors[[tissue]]
    ),
    file.path(OUT_DIR, paste0(tissue, "_module_assignments.csv"))
  )
}

# Hub proteins (top 10 per module by kME)
for (tissue in tissues) {
  hub_rows <- lapply(setdiff(unique(module_colors[[tissue]]), "grey"), function(mod) {
    kme_col <- paste0("kME", mod)
    if (!(kme_col %in% colnames(kME_all[[tissue]]))) return(NULL)
    tibble(
      protein = rownames(kME_all[[tissue]]),
      kME     = kME_all[[tissue]][, kme_col],
      module  = mod
    ) |>
      filter(protein %in% names(module_colors[[tissue]])[module_colors[[tissue]] == mod]) |>
      arrange(desc(kME)) |>
      slice_head(n = 10)
  })
  write_csv(bind_rows(hub_rows),
            file.path(OUT_DIR, paste0(tissue, "_hub_proteins.csv")))
}

# Restore stats::cor
cor <- stats::cor

message("\nScript 1 complete. Proceed to wgcna_02_traits.R")
