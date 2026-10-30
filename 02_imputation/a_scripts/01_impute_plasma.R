# =============================================================================
# Plasma Proteomics — Imputation via imp4p
# =============================================================================
# Input:  Normalized DAList RDS (proteins x samples, log2, NAs for missing)
# Output: (1) Imputed DAList RDS
#         (2) Plain imputed matrix RDS
#         (3) Plain imputed matrix CSV
#         (4) Report data RDS, consumed by 02_imputation_reports_plasma.R
#
# Method: impute.mi() — multiple imputation with MCAR/MNAR mixture model
#   - Mixture estimation: estim.mix()
#   - MNAR component:     impute.igcda() [Gaussian censored distribution]
#   - MCAR component:     impute.mle()   [EM algorithm; swap to "slsa" if slow]
#   - Iterations:         nb.iter (set to 3 for exploration, 10+ for manuscript)
#
# Reference: Giai Gianetto et al. (2020) doi:10.1101/2020.05.29.122770
# =============================================================================

library(imp4p)
library(dplyr)
# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

# NOTE: path updated to match normalization_plasma.R's current output location
# (SURV/ subfolder, "02_DAList_normalized_SURV_plasma.RDS") -- the old flat
# path "01_normalization/c_data/plasma/01_normalized_DAList_SURV.RDS" no
# longer exists now that plasma normalization writes into per-dataset
# subfolders like muscle does. Update this if your normalization output
# location differs.
PATH_DALIST  <- "01_normalization/c_data/plasma/02_normalized_data_plasma.RDS"
PATH_OUT_DIR <- "02_imputation/c_data/plasma"   # tissue subfolder, mirrors muscle layout

SEED       <- 26
NB_ITER    <- 3      # MI iterations: 3 for exploration, 10+ for final run
METHOD_MCAR <- "slsa" # "mle" (EM, default) or "slsa" (faster for large matrices)
METHOD_MNAR <- "igcda"
Q           <- 0.95  # quantile threshold for MNAR boundary
# -----------------------------------------------------------------------------
# 1. Load DAList and extract matrix
# -----------------------------------------------------------------------------

dal <- readRDS(PATH_DALIST)

# DAList stores normalized abundance in $data (proteins x samples)
mat_raw    <- dal$data
annotation <- dal$annotation

cat("Matrix dimensions (proteins x samples):", dim(mat_raw), "\n")
cat("Total missing values:", sum(is.na(mat_raw)), "\n")
cat("Percent missing:", round(100 * mean(is.na(mat_raw)), 2), "%\n")

# -----------------------------------------------------------------------------
# 2. Build condition vector from sample metadata
# -----------------------------------------------------------------------------

# Check available columns with: names(dal$samples)
# Adjust group/timepoint column names below to match your DAList $samples slot
meta <- dal$metadata

# Condition = group x timepoint (e.g. SURV_T1, CTL_T1, SURV_T2)
condition <- factor(paste(meta$cancer, meta$timepoint, sep = "_"))

cat("\nConditions detected:\n")
print(table(condition))

cat("\nMissing values per condition:\n")
print(tapply(colSums(is.na(mat_raw)), condition, sum))

# -----------------------------------------------------------------------------
# 3. Impute using impute.mi()
# -----------------------------------------------------------------------------
# impute.mi() runs the full MCAR/MNAR mixture pipeline internally:
#   1. Estimates mixture model per sample via estim.mix()
#   2. At each iteration, assigns each missing value MCAR/MNAR by Bernoulli
#      draws from estimated probabilities
#   3. MNAR values imputed via impute.igcda() (left-censored Gaussian)
#   4. MCAR values imputed via methodMCAR (MLE by default)
#   5. Final matrix is the average across nb.iter imputed matrices

set.seed(SEED)

mat_imputed <- impute.mi(
  tab          = mat_raw,
  conditions   = condition,
  nb.iter      = NB_ITER,
  methodMCAR   = METHOD_MCAR,
  methodMNAR   = METHOD_MNAR,
  q            = Q,
  selec        = nrow(mat_raw),   # prevents subsampling error when n < default 1000
  progress.bar = TRUE
)

# impute.mi() drops rownames and replaces colnames with integer indices
# Restore from mat_raw before any downstream use
rownames(mat_imputed) <- rownames(mat_raw)
colnames(mat_imputed) <- colnames(mat_raw)

cat("\nPost-imputation missing values:", sum(is.na(mat_imputed)), "\n")

# -----------------------------------------------------------------------------
# 4. Validation checks
# -----------------------------------------------------------------------------

# Dimensions and names preserved
stopifnot(identical(dim(mat_imputed), dim(mat_raw)))
stopifnot(identical(rownames(mat_imputed), rownames(mat_raw)))
stopifnot(identical(colnames(mat_imputed), colnames(mat_raw)))

# Observed values unchanged — imputation only fills missing positions
stopifnot(all(mat_imputed[!is.na(mat_raw)] == mat_raw[!is.na(mat_raw)]))

# Distribution check: imputed values should be left-shifted vs observed
# (confirms MNAR left-censoring component is operating correctly)
observed_vals <- mat_raw[!is.na(mat_raw)]
imputed_vals  <- mat_imputed[is.na(mat_raw)]

cat("\n--- Distribution shift (log2 scale) ---\n")
cat("  Observed mean:", round(mean(observed_vals), 3), "\n")
cat("  Imputed mean: ", round(mean(imputed_vals),  3), "\n")
cat("  Shift:        ", round(mean(observed_vals) - mean(imputed_vals), 3),
    "(should be positive)\n")
cat("  Observed sd:  ", round(sd(observed_vals),   3), "\n")
cat("  Imputed sd:   ", round(sd(imputed_vals),    3), "\n")

# -----------------------------------------------------------------------------
# 5. Missingness diagnostics
# -----------------------------------------------------------------------------

# Per-protein missingness classification
# MNAR-like: missing rate differs substantially across conditions
#            (suggests abundance-dependent dropout)
# MCAR-like: missing rate roughly equal across conditions
MNAR_DIFF_THRESH <- 0.2   # tunable; lower = more sensitive to condition differences

prot_miss_by_cond <- sapply(levels(condition), function(cond) {
  idx <- which(condition == cond)
  rowMeans(is.na(mat_raw[, idx, drop = FALSE]))
})

prot_class <- apply(prot_miss_by_cond, 1, function(x) {
  if (all(x == 0))                          return("complete")
  if (max(x) - min(x) >= MNAR_DIFF_THRESH) return("MNAR-like")
  return("MCAR-like")
})

cat("\n--- Per-protein missingness classification ---\n")
print(table(prot_class))
cat("  (MNAR threshold: condition missing-rate difference >=",
    MNAR_DIFF_THRESH, ")\n")

cat("\n--- Per-condition missingness summary ---\n")
miss_cond_summary <- data.frame(
  condition   = levels(condition),
  n_samples   = as.integer(table(condition)),
  n_missing   = as.integer(tapply(colSums(is.na(mat_raw)), condition, sum)),
  pct_missing = round(100 * tapply(colMeans(is.na(mat_raw)), condition, mean), 2)
)
print(miss_cond_summary, row.names = FALSE)

cat("\n--- Imputation provenance ---\n")
cat("  Package:        imp4p v", as.character(packageVersion("imp4p")), "\n")
cat("  Function:       impute.mi()\n")
cat("  MNAR method:    impute.igcda() [Gaussian censored distribution]\n")
cat("  MCAR method:   ", METHOD_MCAR, "\n")
cat("  MI iterations: ", NB_ITER, "\n")
cat("  q threshold:   ", Q, "\n")
cat("  Seed:          ", SEED, "\n")
cat("  Date:          ", as.character(Sys.Date()), "\n")

# -----------------------------------------------------------------------------
# 6. Inject imputed matrix back into DAList
# -----------------------------------------------------------------------------

dal_imputed      <- dal
dal_imputed$data <- mat_imputed

dal_imputed$imputation <- list(
  method      = "imp4p_impute.mi",
  methodMCAR  = METHOD_MCAR,
  methodMNAR  = METHOD_MNAR,
  nb_iter     = NB_ITER,
  q           = Q,
  seed        = SEED,
  date        = Sys.Date()
)

# -----------------------------------------------------------------------------
# 7. Save outputs
# -----------------------------------------------------------------------------

dir.create(PATH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# (1) Imputed DAList RDS
saveRDS(dal_imputed,
        file = file.path(PATH_OUT_DIR, "plasma_imputed_dalist.RDS"))

# (2) Plain imputed matrix RDS
saveRDS(mat_imputed,
        file = file.path(PATH_OUT_DIR, "plasma_imputed_matrix.RDS"))

# (3) Plain imputed matrix CSV (proteins as rows)
write.csv(mat_imputed,
          file = file.path(PATH_OUT_DIR, "plasma_imputed_matrix.csv"),
          row.names = TRUE)

# (4) Report data RDS -- bundles everything 02_imputation_reports_plasma.R
#     needs so plots/tables can be regenerated without rerunning impute.mi()
report_data <- list(
  config = list(
    path_dalist = PATH_DALIST,
    data_dir    = PATH_OUT_DIR,
    report_dir  = "02_imputation/b_reports/plasma"
  ),
  mat_raw           = mat_raw,
  mat_imputed       = mat_imputed,
  annotation        = annotation,
  meta              = meta,
  condition         = condition,
  prot_miss_by_cond = prot_miss_by_cond,
  prot_class        = prot_class,
  miss_cond_summary = miss_cond_summary,
  imputation_params = dal_imputed$imputation
)

saveRDS(report_data, file.path(PATH_OUT_DIR, "00_report_data_plasma.rds"))

cat("\nDone. Outputs written to:", PATH_OUT_DIR, "\n")
cat("  plasma_imputed_dalist.RDS\n")
cat("  plasma_imputed_matrix.RDS\n")
cat("  plasma_imputed_matrix.csv\n")
cat("  00_report_data_plasma.rds\n")