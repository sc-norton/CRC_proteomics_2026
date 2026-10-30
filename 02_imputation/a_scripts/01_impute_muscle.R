# =============================================================================
# Muscle Proteomics — Imputation via imp4p
# =============================================================================
# Input:  Normalized DAList RDS (proteins x samples, log2, NAs for missing)
# Output: (1) Imputed DAList RDS
#         (2) Plain imputed matrix RDS
#         (3) Plain imputed matrix CSV
#         (4) Report data RDS, consumed by 02_imputation_reports_muscle.R
#         (5) Convergence check diagnostic plot (if RUN_CONVERGENCE_CHECK)
#
# Method: impute.mi() — multiple imputation with MCAR/MNAR mixture model
#   - Mixture estimation: estim.mix()
#   - MNAR component:     impute.igcda() [Gaussian censored distribution]
#   - MCAR component:     impute.mle()   [EM algorithm; swap to "slsa" if slow]
#   - Iterations:         nb.iter (see convergence check, section 3, for
#                          how to pick this rather than guessing)
#
# Reference: Giai Gianetto et al. (2020) doi:10.1101/2020.05.29.122770
#
# NOTE ON nb.iter: impute.mi() averages across nb.iter stochastic imputation
# draws into a single point-estimate matrix (this is NOT classical
# Rubin's-rules multiple imputation, which would keep draws separate and
# pool downstream analyses) -- nb.iter = 3 is the package's OWN built-in
# default across its documentation, not a bare-minimum placeholder. There is
# no universally-cited "use N for a final run" recommendation from the
# package or its methods paper; the right way to choose nb.iter is an
# empirical Monte Carlo convergence check (section 3 below), not a round
# number picked from general MI heuristics.
# =============================================================================

library(imp4p)
library(dplyr)
# -----------------------------------------------------------------------------
# 0. Config
# -----------------------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

PATH_DALIST  <- "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS"
PATH_OUT_DIR <- "02_imputation/c_data/muscle"   # tissue subfolder, mirrors 01_normalization layout

SEED       <- 26
NB_ITER    <- 3      # Run the convergence check below FIRST, then set this
# based on those results rather than guessing.
METHOD_MCAR <- "slsa" # "mle" (EM, default) or "slsa" (faster for large matrices)
METHOD_MNAR <- "igcda"
Q           <- 0.95  # quantile threshold for MNAR boundary

# Convergence check config -- see section 3. Each value tested is a FULL
# impute.mi() run (the package's own docs flag even nb.iter=3 as "can be
# time consuming"), so this check's total runtime is roughly proportional
# to sum(CONVERGENCE_ITER_VALUES). Set RUN_CONVERGENCE_CHECK <- FALSE once
# you've picked and locked in a final NB_ITER, to skip re-running this on
# subsequent executions of this script.
RUN_CONVERGENCE_CHECK   <- FALSE
CONVERGENCE_ITER_VALUES <- c(3, 10, 20)

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

# Condition = group x timepoint (e.g. CRC_T1, CTL_T1, CRC_T2)
condition <- factor(paste(meta$cancer, meta$timepoint, sep = "_"))

cat("\nConditions detected:\n")
print(table(condition))

cat("\nMissing values per condition:\n")
print(tapply(colSums(is.na(mat_raw)), condition, sum))

# -----------------------------------------------------------------------------
# 3. Convergence check for nb.iter
# -----------------------------------------------------------------------------
# Runs impute.mi() at each value in CONVERGENCE_ITER_VALUES and compares the
# resulting imputed cells to each other. Observed (non-missing) cells are
# identical across every run by construction (impute.mi() never touches
# them), so only imputed cells are informative here.
#
# Same SEED is used before each call -- this isolates the effect of nb.iter
# itself on the AVERAGED result, rather than comparing runs that also differ
# by random starting point. High correlation and small mean_abs_diff between
# the two LARGEST tested values suggests the average has stabilized and
# going higher wouldn't meaningfully change results.
dir.create(PATH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)
convergence_mats <- list()

if (RUN_CONVERGENCE_CHECK) {
  cat("\n--- Convergence check: running impute.mi() at nb.iter =",
      paste(CONVERGENCE_ITER_VALUES, collapse = ", "), "---\n")
  
  for (n in CONVERGENCE_ITER_VALUES) {
    cat("\n  Running nb.iter =", n, "...\n")
    set.seed(SEED)
    m <- impute.mi(
      tab          = mat_raw,
      conditions   = condition,
      nb.iter      = n,
      methodMCAR   = METHOD_MCAR,
      methodMNAR   = METHOD_MNAR,
      q            = Q,
      progress.bar = TRUE
    )
    rownames(m) <- rownames(mat_raw)
    colnames(m) <- colnames(mat_raw)
    convergence_mats[[as.character(n)]] <- m
  }
  
  missing_idx <- is.na(mat_raw)
  
  convergence_summary <- data.frame(
    comparison    = character(0),
    cor           = numeric(0),
    mean_abs_diff = numeric(0)
  )
  
  if (length(CONVERGENCE_ITER_VALUES) >= 2) {
    for (i in 2:length(CONVERGENCE_ITER_VALUES)) {
      n1 <- as.character(CONVERGENCE_ITER_VALUES[i - 1])
      n2 <- as.character(CONVERGENCE_ITER_VALUES[i])
      v1 <- convergence_mats[[n1]][missing_idx]
      v2 <- convergence_mats[[n2]][missing_idx]
      convergence_summary <- rbind(convergence_summary, data.frame(
        comparison    = paste0(n1, " vs ", n2),
        cor           = cor(v1, v2),
        mean_abs_diff = mean(abs(v1 - v2))
      ))
    }
  }
  
  cat("\n--- Convergence summary (imputed cells only) ---\n")
  print(convergence_summary, row.names = FALSE)
  cat("\n  High correlation (~1) and small mean_abs_diff between the LAST two\n",
      "  comparisons suggest the averaged matrix has stabilized -- going higher\n",
      "  than the largest tested value likely wouldn't change results much.\n",
      "  If mean_abs_diff is still substantial at the highest values tested,\n",
      "  add a larger value to CONVERGENCE_ITER_VALUES and re-run this section.\n", sep = "")
  
  # Diagnostic scatter for the LAST (most relevant) comparison -- points on
  # the red y=x line indicate agreement between the two iteration counts.
  if (nrow(convergence_summary) > 0) {
    n_pair <- CONVERGENCE_ITER_VALUES[c(length(CONVERGENCE_ITER_VALUES) - 1, length(CONVERGENCE_ITER_VALUES))]
    v1 <- convergence_mats[[as.character(n_pair[1])]][missing_idx]
    v2 <- convergence_mats[[as.character(n_pair[2])]][missing_idx]
    
    pdf_path <- file.path(PATH_OUT_DIR, "convergence_check_muscle.pdf")
    pdf(pdf_path, width = 6, height = 6)
    plot(v1, v2,
         xlab = paste0("Imputed value, nb.iter = ", n_pair[1]),
         ylab = paste0("Imputed value, nb.iter = ", n_pair[2]),
         main = paste0("Convergence check: nb.iter ", n_pair[1], " vs ", n_pair[2]),
         pch = 16, cex = 0.4, col = adjustcolor("black", alpha.f = 0.3))
    abline(0, 1, col = "red", lty = 2)
    dev.off()
    cat("\n  Convergence diagnostic plot saved to:", pdf_path, "\n")
  }
  
  cat("\n  Set NB_ITER above based on the results before proceeding.\n")
} else {
  cat("\nRUN_CONVERGENCE_CHECK is FALSE -- skipping convergence check, proceeding directly to final imputation with NB_ITER =", NB_ITER, "\n")
}

# -----------------------------------------------------------------------------
# 4. Impute using impute.mi() [FINAL RUN]
# -----------------------------------------------------------------------------
# impute.mi() runs the full MCAR/MNAR mixture pipeline internally:
#   1. Estimates mixture model per sample via estim.mix()
#   2. At each iteration, assigns each missing value MCAR/MNAR by Bernoulli
#      draws from estimated probabilities
#   3. MNAR values imputed via impute.igcda() (left-censored Gaussian)
#   4. MCAR values imputed via methodMCAR (MLE by default)
#   5. Final matrix is the average across nb.iter imputed matrices
#
# If NB_ITER matches a value already computed during the convergence check
# above, that matrix is reused directly rather than re-running impute.mi()
# a second time with the same settings.

if (RUN_CONVERGENCE_CHECK && as.character(NB_ITER) %in% names(convergence_mats)) {
  cat("\nNB_ITER =", NB_ITER, "was already computed during the convergence check -- reusing that matrix.\n")
  mat_imputed <- convergence_mats[[as.character(NB_ITER)]]
} else {
  set.seed(SEED)
  mat_imputed <- impute.mi(
    tab          = mat_raw,
    conditions   = condition,
    nb.iter      = NB_ITER,
    methodMCAR   = METHOD_MCAR,
    methodMNAR   = METHOD_MNAR,
    q            = Q,
    progress.bar = TRUE
  )
  # impute.mi() drops rownames and replaces colnames with integer indices
  # Restore from mat_raw before any downstream use
  rownames(mat_imputed) <- rownames(mat_raw)
  colnames(mat_imputed) <- colnames(mat_raw)
}

cat("\nPost-imputation missing values:", sum(is.na(mat_imputed)), "\n")

# -----------------------------------------------------------------------------
# 5. Validation checks
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
# 6. Missingness diagnostics
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
# 7. Inject imputed matrix back into DAList
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
  date        = Sys.Date(),
  convergence_check_run    = RUN_CONVERGENCE_CHECK,
  convergence_iter_values  = if (RUN_CONVERGENCE_CHECK) CONVERGENCE_ITER_VALUES else NA
)

# -----------------------------------------------------------------------------
# 8. Save outputs
# -----------------------------------------------------------------------------

dir.create(PATH_OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# (1) Imputed DAList RDS
saveRDS(dal_imputed,
        file = file.path(PATH_OUT_DIR, "muscle_imputed_dalist.RDS"))

# (2) Plain imputed matrix RDS
saveRDS(mat_imputed,
        file = file.path(PATH_OUT_DIR, "muscle_imputed_matrix.RDS"))

# (3) Plain imputed matrix CSV (proteins as rows)
write.csv(mat_imputed,
          file = file.path(PATH_OUT_DIR, "muscle_imputed_matrix.csv"),
          row.names = TRUE)

# (4) Report data RDS -- bundles everything 02_imputation_reports_muscle.R
#     needs so plots/tables can be regenerated without rerunning impute.mi()
report_data <- list(
  config = list(
    path_dalist = PATH_DALIST,
    data_dir    = PATH_OUT_DIR,
    report_dir  = "02_imputation/b_reports/muscle"
  ),
  mat_raw              = mat_raw,
  mat_imputed          = mat_imputed,
  annotation           = annotation,
  meta                 = meta,
  condition            = condition,
  prot_miss_by_cond    = prot_miss_by_cond,
  prot_class           = prot_class,
  miss_cond_summary    = miss_cond_summary,
  imputation_params    = dal_imputed$imputation,
  convergence_summary  = if (RUN_CONVERGENCE_CHECK) convergence_summary else NULL
)

saveRDS(report_data, file.path(PATH_OUT_DIR, "00_report_data_muscle.rds"))

cat("\nDone. Outputs written to:", PATH_OUT_DIR, "\n")
cat("  muscle_imputed_dalist.RDS\n")
cat("  muscle_imputed_matrix.RDS\n")
cat("  muscle_imputed_matrix.csv\n")
cat("  00_report_data_muscle.rds\n")
if (RUN_CONVERGENCE_CHECK) cat("  convergence_check_muscle.pdf\n")