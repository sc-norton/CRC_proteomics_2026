# --- Fry Gene Set Analysis (standalone script) ---

# --- 0: Setup ---
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr,
               ggplot2, limma, writexl, msigdbr)

config <- list(
  report_dir = "03_analysis/b_reports/muscle",
  data_dir   = "03_analysis/c_data/muscle",
  SURV_RDS   = "01_normalization/c_data/muscle/SURV/stringent/BM15/02_DAList_normalized_SURV_stringent_BM15_muscle.RDS",
  imputed    = "02_imputation/c_data/mat_surv_imputed.csv"
)

fry_dir <- file.path(config$report_dir, "all_data", "limma", "fry")
dir.create(fry_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1: Load saved limma outputs ---
limma_out <- readRDS(
  file.path(config$report_dir, "SURVvCTL", "limma", "limma_out_SURV.RDS")
)
out_surv <- limma_out

# --- 2: Rebuild fit object ---
# We need the raw DAList and to refit — this is fast and avoids
# rerunning the full limma script. The DAList RDS is from normalization,
# not the limma script, so it's always available independently.

contrast_specs_surv <- c(
  "Baseline_SURVvCTL = SURV_T1 - CTL_T1",
  "Training_SURVvCTL = SURV_T2 - CTL_T1",
  "Training_SURV     = SURV_T2 - SURV_T1"
)

SURV <- readRDS(config$SURV_RDS)

mat_surv  <- as.matrix(SURV$data)
group_col <- "cancer_time"
block_col <- "pid"

dal_surv <- SURV
dal_surv <- add_design(dal_surv, paste0("~ 0 + ", group_col, " + (1 | ", block_col, ")"))
colnames(dal_surv$design$design_matrix) <-
  gsub(paste0("^", group_col), "", colnames(dal_surv$design$design_matrix))
dal_surv <- add_contrasts(dal_surv, contrasts_vector = contrast_specs_surv)

design_mat_surv   <- dal_surv$design$design_matrix
contrast_mat_surv <- dal_surv$design$contrast_matrix
block_var_surv    <- dal_surv$metadata[[block_col]]
aw_surv           <- arrayWeights(mat_surv, design_mat_surv)
dupcor_surv       <- duplicateCorrelation(mat_surv, design_mat_surv,
                                          block   = block_var_surv,
                                          weights = aw_surv)

fit_surv <- lmFit(mat_surv, design_mat_surv,
                  block       = block_var_surv,
                  correlation = dupcor_surv$consensus.correlation,
                  weights     = aw_surv) |>
  contrasts.fit(contrast_mat_surv) |>
  eBayes(robust = TRUE, trend = TRUE)

# --- 3: Map UniProt IDs to gene symbols ---
protein_ids <- rownames(mat_surv)

id_map <- out_surv$results_list[[1]] |>
  select(uniprot, gene) |>
  filter(!is.na(gene), gene != "") |>
  distinct()

gene_ids <- id_map$gene[match(protein_ids, id_map$uniprot)]

cat(sprintf("Proteins mapped to gene symbols: %d / %d\n",
            sum(!is.na(gene_ids)), length(protein_ids)))

# --- 4: Load MSigDB gene sets ---
msig_go <- msigdbr(species = "Homo sapiens",
                   collection    = "C5",
                   subcollection = "GO:BP")
msig_h  <- msigdbr(species = "Homo sapiens",
                   collection = "H")

gene_sets_go       <- split(msig_go$gene_symbol, msig_go$gs_name)
gene_sets_hallmark <- split(msig_h$gene_symbol,  msig_h$gs_name)

# --- 5: Build index lists ---
build_index <- function(gene_sets, ids, min_size = 10, max_size = 500) {
  idx <- lapply(gene_sets, function(gs) which(ids %in% gs))
  idx[sapply(idx, length) >= min_size & sapply(idx, length) <= max_size]
}

index_go       <- build_index(gene_sets_go,      gene_ids)
index_hallmark <- build_index(gene_sets_hallmark, gene_ids)

cat(sprintf("Gene sets retained — GO:BP: %d | Hallmark: %d\n",
            length(index_go), length(index_hallmark)))

# --- 6: Define contrasts ---
reversal_contrast <- c(Baseline_SURVvCTL = 1, Training_SURVvCTL = 0, Training_SURV = -1)

individual_contrasts <- list(
  disease  = c(Baseline_SURVvCTL = 1, Training_SURVvCTL = 0, Training_SURV = 0),
  exercise = c(Baseline_SURVvCTL = 0, Training_SURVvCTL = 0, Training_SURV = 1)
)

# --- 7: Fry runner function ---
run_fry <- function(fit, y, index, contrast, label, out_dir) {
  cat(sprintf("  Running fry: %s (%d gene sets)\n", label, length(index)))
  
  res <- fry(
    y        = y,
    index    = index,
    design   = fit$design,
    contrast = contrast
  )
  
  res$GeneSet <- rownames(res)
  res <- res[order(res$PValue), ]
  
  write.csv(res,
            file.path(out_dir, paste0("fry_", label, ".csv")),
            row.names = FALSE)
  
  cat(sprintf("    Significant (FDR<0.05): %d directional | %d mixed\n",
              sum(res$FDR < 0.05, na.rm = TRUE),
              sum(res$FDR.Mixed < 0.05, na.rm = TRUE)))
  
  invisible(res)
}

# --- 8: Run fry ---
fry_reversal_go <- run_fry(
  fit      = fit_surv,
  y        = mat_surv,
  index    = index_go,
  contrast = reversal_contrast,
  label    = "reversal_GOBP",
  out_dir  = fry_dir
)

fry_reversal_hallmark <- run_fry(
  fit      = fit_surv,
  y        = mat_surv,
  index    = index_hallmark,
  contrast = reversal_contrast,
  label    = "reversal_Hallmark",
  out_dir  = fry_dir
)

for (cname in names(individual_contrasts)) {
  run_fry(fit_surv, mat_surv, index_go,
          contrast = individual_contrasts[[cname]],
          label    = paste0(cname, "_GOBP"),
          out_dir  = fry_dir)
  run_fry(fit_surv, mat_surv, index_hallmark,
          contrast = individual_contrasts[[cname]],
          label    = paste0(cname, "_Hallmark"),
          out_dir  = fry_dir)
}

# --- Diagnose ---
cat(sprintf("NA values:  %d\n", sum(is.na(mat_surv))))
cat(sprintf("Inf values: %d\n", sum(is.infinite(mat_surv))))
cat(sprintf("NaN values: %d\n", sum(is.nan(mat_surv))))

# How many proteins affected?
na_rows  <- rowSums(is.na(mat_surv) | is.infinite(mat_surv) | is.nan(mat_surv))
cat(sprintf("Proteins with any missing/infinite values: %d / %d\n",
            sum(na_rows > 0), nrow(mat_surv)))

# =============================================================================
# Fry gene set analysis: cancer signature sets (PREvCTL) tested on POSTvPRE
# -----------------------------------------------------------------------------
# Approach:
#   1. Load saved limma outputs (no upstream re-running required)
#   2. Define gene sets from proteins significantly dysregulated in PREvCTL
#   3. Run fry() on those sets using the POSTvPRE contrast
#   4. Visualize and save results
# =============================================================================

pacman::p_load(limma, dplyr, tibble, ggplot2, patchwork, stringr)

# =============================================================================
# 0. USER-DEFINED PATHS AND SETTINGS — edit these
# =============================================================================

path_dalist  <- "01_normalization/c_data/SURV/stringent/BM15/02_DAList_normalized_SURV_stringent_BM15_muscle.RDS"
path_limma   <- "03_analysis/b_reports/SURVvCTL/limma/limma_out_SURV.RDS"
path_out_dir <- "03_analysis/b_reports/SURVvCTL/fry/cancer_signature"

padj_threshold <- 0.05   # significance cutoff for defining PREvCTL sets
fry_fdr_thresh <- 0.05   # FDR threshold for fry results

# Contrast names as they appear in your limma results
# Adjust these strings to match your actual contrast names
contrast_pre_ctl  <- "Baseline_SURVvCTL"   # PRE vs CTL contrast name
contrast_post_pre <- "Training_SURV"        # POST vs PRE contrast name
contrast_col_name <- "Training_SURV"        # column name in contrasts matrix

# Column names in your limma results table — adjust if different
col_accession <- "accession"    # UniProt accession or row identifier
col_gene      <- "gene_symbol"  # gene symbol column (may be NULL if not present)
col_contrast  <- "contrast"     # column identifying which contrast each row belongs to
col_padj      <- "adj.P.Val"    # adjusted p-value column
col_logfc     <- "logFC"        # log fold-change column

# =============================================================================
# 1. LOAD DATA
# =============================================================================

message("Loading DAList...")
da_list <- readRDS(path_dalist)

message("Loading limma results...")
lim_out <- readRDS(path_limma)

# Extract expression matrix
# Adjust slot name if your DAList uses a different structure
expr_mat <- da_list$data
message(sprintf("Expression matrix: %d proteins x %d samples", nrow(expr_mat), ncol(expr_mat)))

# Extract design and contrast matrix
design    <- lim_out$design
contrasts <- lim_out$contrasts

message("Contrasts available: ", paste(colnames(contrasts), collapse = ", "))

# =============================================================================
# 2. EXTRACT RESULTS TABLE
# =============================================================================

# lim_out$results is assumed to be a data frame with one row per protein per
# contrast, with columns: contrast, accession, gene_symbol, logFC, adj.P.Val
# Adjust the slot name below if your object uses a different name
# e.g. lim_out$result_table, lim_out$de_results, etc.

results_tbl <- lim_out$results

# Verify expected columns are present
required_cols <- c(col_accession, col_contrast, col_padj, col_logfc)
missing_cols  <- setdiff(required_cols, colnames(results_tbl))
if (length(missing_cols) > 0) {
  stop("Missing columns in results table: ", paste(missing_cols, collapse = ", "),
       "\nAvailable columns: ", paste(colnames(results_tbl), collapse = ", "))
}

message("Results table loaded: ", nrow(results_tbl), " rows")
message("Contrasts in results: ", paste(unique(results_tbl[[col_contrast]]), collapse = ", "))

# =============================================================================
# 3. DEFINE GENE SETS FROM PREvCTL
# =============================================================================

message(sprintf("\nDefining gene sets from %s (padj < %s)...", contrast_pre_ctl, padj_threshold))

pre_ctl <- results_tbl |>
  filter(.data[[col_contrast]] == contrast_pre_ctl,
         .data[[col_padj]]     <  padj_threshold)

sig_up   <- pre_ctl |> filter(.data[[col_logfc]] > 0) |> pull(.data[[col_accession]])
sig_down <- pre_ctl |> filter(.data[[col_logfc]] < 0) |> pull(.data[[col_accession]])

message(sprintf("  Proteins up in SURV at baseline   (Up set): n = %d", length(sig_up)))
message(sprintf("  Proteins down in SURV at baseline (Down set): n = %d", length(sig_down)))

if (length(sig_up) == 0 && length(sig_down) == 0) {
  stop("No significant proteins found in ", contrast_pre_ctl, " at padj < ", padj_threshold,
       ". Check contrast name and threshold.")
}

# Build named list for ids2indices
# Using All set for reference but Up/Down are the primary sets of interest
gs_custom <- list(
  CancerSig_Up   = sig_up,
  CancerSig_Down = sig_down,
  CancerSig_All  = c(sig_up, sig_down)
)

# Remove any empty sets
gs_custom <- gs_custom[sapply(gs_custom, length) > 0]
message("Gene sets to test: ", paste(names(gs_custom), collapse = ", "))

# =============================================================================
# 4. BUILD INDEX AGAINST EXPRESSION MATRIX ROWNAMES
# =============================================================================

message("\nBuilding gene set index...")

# If rownames of expr_mat are UniProt accessions, ids2indices matches directly.
# If rownames are something else (e.g. gene symbols), adjust identifiers= below.
gs_index <- limma::ids2indices(
  gene.sets   = gs_custom,
  identifiers = rownames(expr_mat)
)

# Report how many proteins from each set were matched
for (nm in names(gs_index)) {
  message(sprintf("  %s: %d proteins matched in expression matrix", nm, length(gs_index[[nm]])))
}

# Warn if any set has very few members
small_sets <- names(gs_index)[sapply(gs_index, length) < 5]
if (length(small_sets) > 0) {
  warning("These sets have < 5 matched proteins and fry results will be unreliable: ",
          paste(small_sets, collapse = ", "))
}

# =============================================================================
# 5. RUN FRY ON POSTvPRE CONTRAST
# =============================================================================

message(sprintf("\nRunning fry() on contrast: %s...", contrast_col_name))

# Extract the contrast vector for POSTvPRE
contrast_vec <- contrasts[, contrast_col_name]

fry_result <- fry(
  y        = expr_mat,
  index    = gs_index,
  design   = design,
  contrast = contrast_vec
)

# Add set name as explicit column and convert to tibble
fry_result <- fry_result |>
  rownames_to_column("gene_set") |>
  as_tibble() |>
  arrange(FDR)

message("\nFry results:")
print(fry_result)

# Interpretation helper column
fry_result <- fry_result |>
  mutate(
    interpretation = case_when(
      # Up set going Down = elevated proteins suppressed by training = reversal
      gene_set == "CancerSig_Up"   & Direction == "Down" & FDR < fry_fdr_thresh ~ "Reversal (suppressed)",
      # Down set going Up = depleted proteins restored by training = reversal
      gene_set == "CancerSig_Down" & Direction == "Up"   & FDR < fry_fdr_thresh ~ "Reversal (restored)",
      # Same direction as baseline dysregulation = exacerbation
      gene_set == "CancerSig_Up"   & Direction == "Up"   & FDR < fry_fdr_thresh ~ "Exacerbation",
      gene_set == "CancerSig_Down" & Direction == "Down" & FDR < fry_fdr_thresh ~ "Exacerbation",
      FDR >= fry_fdr_thresh ~ "Not significant",
      TRUE ~ "Significant (All set)"
    )
  )

message("\nInterpretation summary:")
print(fry_result |> select(gene_set, NGenes, Direction, PValue, FDR, interpretation))

# =============================================================================
# 6. BARCODE PLOTS — visualize within-set protein rank distributions
# =============================================================================

message("\nGenerating barcode plots...")

dir.create(path_out_dir, recursive = TRUE, showWarnings = FALSE)

# t-statistics for POSTvPRE contrast (used for barcode plot ranking)
# Adjust slot name if your limma output stores t-stats differently
t_stats <- lim_out$tstat[, contrast_col_name]

# Verify t-stats rownames match expr_mat rownames
if (!all(rownames(expr_mat) %in% names(t_stats))) {
  message("Note: some expression matrix rownames not found in t-stat vector — subsetting to common")
  common <- intersect(rownames(expr_mat), names(t_stats))
  t_stats <- t_stats[common]
}

pdf(file.path(path_out_dir, "fry_barcode_plots.pdf"),
    width = 10, height = 5, onefile = TRUE)

for (nm in names(gs_index)) {
  if (length(gs_index[[nm]]) < 2) next
  
  # Direction label for plot title
  this_dir <- fry_result |> filter(gene_set == nm) |> pull(Direction)
  this_fdr <- fry_result |> filter(gene_set == nm) |> pull(FDR) |> round(4)
  
  barcodeplot(
    statistics = t_stats,
    index      = gs_index[[nm]],
    main       = sprintf("%s  |  Direction: %s  |  FDR: %s", nm, this_dir, this_fdr),
    xlab       = sprintf("t-statistic (%s)", contrast_col_name),
    quantiles  = c(-1, 1),
    cex.main   = 0.9
  )
}

dev.off()
message("Barcode plots saved to: ", file.path(path_out_dir, "fry_barcode_plots.pdf"))

# =============================================================================
# 7. SUMMARY PLOT — dot plot of fry results
# =============================================================================

message("Generating summary dot plot...")

# Colour by interpretation
interp_colours <- c(
  "Reversal (suppressed)" = "#1D9E75",
  "Reversal (restored)"   = "#3B8BD4",
  "Exacerbation"          = "#D85A30",
  "Not significant"       = "#888780",
  "Significant (All set)" = "#7F77DD"
)

p_dot <- fry_result |>
  mutate(
    neg_log_fdr = -log10(FDR),
    label       = str_replace_all(gene_set, "_", " ")
  ) |>
  ggplot(aes(x = neg_log_fdr, y = reorder(label, neg_log_fdr),
             colour = interpretation, size = NGenes)) +
  geom_vline(xintercept = -log10(fry_fdr_thresh),
             linetype = "dashed", colour = "grey60", linewidth = 0.4) +
  geom_point(alpha = 0.9) +
  scale_colour_manual(values = interp_colours, name = "Interpretation") +
  scale_size_continuous(name = "N proteins", range = c(3, 8)) +
  labs(
    title    = "Fry gene set test: cancer signature sets",
    subtitle = sprintf("Sets defined from %s (padj < %s), tested on %s",
                       contrast_pre_ctl, padj_threshold, contrast_col_name),
    x        = "-log10(FDR)",
    y        = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title    = element_text(size = 12, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey40")
  )

ggsave(file.path(path_out_dir, "fry_dot_plot.pdf"),
       p_dot, width = 8, height = 4)
message("Dot plot saved to: ", file.path(path_out_dir, "fry_dot_plot.pdf"))

# =============================================================================
# 8. SAVE RESULTS TABLE
# =============================================================================

write.csv(fry_result,
          file.path(path_out_dir, "fry_results_cancer_signature.csv"),
          row.names = FALSE)
message("Results table saved to: ", file.path(path_out_dir, "fry_results_cancer_signature.csv"))

message("\nDone.")
