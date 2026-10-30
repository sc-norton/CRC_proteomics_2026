# === PROJECT: CRC Proteomics - Plasma =============================
# === Script: Normalization

# === PURPOSE: To load the raw data and metadata, assemble a DAList for use in
# proteoDA, filter the data for contaminants and missingness, and normalize
# for downstream analyses.

# === INPUTS: raw data, meta data

# === OUTPUTS: raw sample summary, normalization reports, flagged contaminants,
# pre & post QC reports, filtering effects, outlier diagnostics, final protein
# report, final sample report, normalized data (csv), normalized DA list,
# normalization method quality scores, and a consolidated 00_report_data.rds
# consumed by 02_normalization_reports_plasma.R

# --- 0: Setup --------------------------------------------
cat("\n>> 0 - Setup\n")
cat("\n Loading required packages...")

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, tibble)

cat("\n Packages successfully loaded.\n")

# Set working directory
setwd(rprojroot::find_rstudio_root_file())
getwd()

# Configuration
config <- list(
  # Data files - plasma
  raw_plasma  = "00_input/plasma/CRp_raw.csv",
  meta_plasma = "00_input/plasma/CRp_meta.csv",
  
  # Directories
  report_dir = "01_normalization/b_reports/plasma",
  data_dir   = "01_normalization/c_data/plasma",
  
  # Parameters
  miss_pct    = 50,
  outlier_k   = 3,        # NOTE: consensus_outlier below still uses n_flags >= 2
  # (plasma's original threshold) rather than muscle's
  # n_flags >= 3 -- kept as-is since this is a scientific
  # choice, not an architecture change.
  mahal_p     = 0.01,
  mad_k       = 3,
  norm_method = "cycloess"
)

dir.create(config$report_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir,   showWarnings = F, recursive = T)

build_output_path <- function(nm, base_dir) {
  out_dir <- file.path(base_dir, nm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  return(out_dir)
}

compute_missingness_stats <- function(d) {
  data <- d$data
  list(
    protein_missing = tibble(uniprot_id = rownames(data),
                             pct_missing = rowMeans(is.na(data)) * 100),
    sample_missing  = tibble(sample_id = colnames(data),
                             n_missing = colSums(is.na(data)))
  )
}

compute_pca <- function(data, log_transform = TRUE) {
  data_for_pca <- data
  for (j in seq_len(ncol(data_for_pca))) {
    nas <- is.na(data_for_pca[, j])
    if (any(nas)) data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
  if (log_transform) data_for_pca <- log2(data_for_pca + 1)
  prcomp(t(data_for_pca), center = TRUE, scale. = TRUE)
}

cat("\n Base directory: ", getwd(), "\n")

# --- 1: Load & Restructure Data ---------------------------------------
cat("\n>> 1 - Load & Restructure Data\n")

# Load Data
raw  <- read.csv(file.path(config$raw_plasma))
meta <- read.csv(file.path(config$meta_plasma))

# Rename annotation column names
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq", "sum_spectra")
colnames(raw)[1:6] <- annot_cols

# Make Uniprot_ID the row name
rownames(raw) <- raw$uniprot_id

# Split input data into protein intensity data and annotation data
intensity  <- raw[ , setdiff(names(raw), annot_cols)]
annotation <- raw[ , annot_cols]

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))
cat(sprintf("\n Meta data: %d samples\n", nrow(meta)))

# Rename metadata variables: SURV -> CRC (cancer survivor group label)
# -- mirrors the muscle normalization script, so downstream group labels
# (CTL / CRC_PRE / CRC_POST) and the reports script's recode_group_labels()
# stay consistent across tissues.
meta <- meta %>%
  mutate(
    cancer      = str_replace(cancer,      "^SURV$", "CRC"),
    cancer_time = str_replace(cancer_time, "^SURV",  "CRC")
  )

# Align row names of metadata with column names of data
cat(" Aligning data and metadata...\n")
rownames(meta) <- meta[ , 1]

meta <- meta[order(meta$sample_id), ]

cat(" Data aligned to metadata.\n")

# Validate alignment
cat(" Validating alignment...\n")

data_samples <- colnames(intensity)
meta_samples <- meta$sample_id

stopifnot(
  "Sample mismatch between data and metadata" =
    setequal(data_samples, meta_samples)
)

# Reorder intensity columns to match metadata row order
intensity <- intensity[, meta_samples]

cat("   Validation passed: all", length(meta_samples), "samples aligned.\n")

# Sample summary
cat("\n   Sample distribution:\n")
print(meta %>% count(cancer, supp, timepoint, cancer_time, supp_time))

sample_summary <- as.data.frame(meta %>%
                                  count(cancer, supp, timepoint, cancer_time, supp_time))

write.csv(sample_summary, file.path(config$report_dir, "00_raw_sample_summary_plasma.csv"))

# --- 2: Contaminant Identification (Keratins) ----------------------------
cat("\n>> 2 - Contaminant Identification\n")

# NOTE: plasma has no blood-contamination filtering step (unlike muscle) --
# plasma samples ARE blood, so the blood-vs-tissue marker logic used for
# muscle doesn't apply here. Keratin remains the only contaminant class.

# Identify keratins
keratins <- annotation %>%
  filter(grepl("keratin", description, ignore.case = TRUE)) %>%
  pull(uniprot_id)

# Join to annotation data
annotation <- annotation %>%
  mutate(is_keratin = uniprot_id %in% keratins)

n_flagged <- sum(annotation$is_keratin)
cat(sprintf("\n Proteins flagged as contaminants: %d\n", n_flagged))

cat("\n Flagged contaminants: \n")
print(annotation %>%
        filter(is_keratin) %>%
        dplyr::select(gene, uniprot_id, description))

# Export contaminants
flagged_contaminants <- annotation %>%
  filter(is_keratin) %>%
  mutate(flag_source = "Keratin")

write.csv(flagged_contaminants,
          file.path(config$report_dir, "01_flagged_contaminants_plasma.csv"),
          row.names = FALSE)
cat("\n Saved flagged contaminant CSV to:", config$report_dir, "\n")

# --- 3: Assemble DAList -----------------------------------------
cat("\n>> 3 - Assemble DAList \n")

# Check for duplicate uniprot_ids
cat(" Checking for duplicate uniprot ids...\n")
dup_ids <- annotation$uniprot_id[duplicated(annotation$uniprot_id)]

if (length(dup_ids) > 0) {
  cat(sprintf("   Found %d duplicate uniprot_ids — deduplicating by highest mean intensity.\n",
              length(dup_ids)))
  
  intensity_num    <- as.data.frame(lapply(intensity, as.numeric))
  annotation$row_mean <- rowMeans(intensity_num, na.rm = TRUE)
  
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  
  annotation    <- annotation[keep_idx, ]
  intensity     <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  
  cat(sprintf("   After deduplication: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids found.\n")
}

# Assemble DAList (all samples — CRC vs CTL, all supp arms pooled)
dal <- DAList(
  data       = intensity,
  metadata   = meta,
  annotation = annotation
)

cat(sprintf(" DAList: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# Single grouping column for missingness filtering / outlier detection /
# normalization diagnostics -- mirrors muscle's single-DAList approach.
group_col <- "cancer_time"

# --- 3b: Preserve raw (pre-filter) data for missingness/PCA reporting -----
cat("\n>> 3b - Snapshot raw data for reporting\n")
dal_raw <- dal   # DAList unfiltered, pre-outlier-removal

missingness_raw <- compute_missingness_stats(dal_raw)
pca_raw <- compute_pca(dal_raw$data, log_transform = TRUE)

# --- 4: Filter Data -------------------------------------
cat("\n>> 4 - Filter Data \n")

log_step <- function(log, step, n_before, n_after) {
  bind_rows(log, tibble(
    step      = step,
    n_before  = n_before,
    n_after   = n_after,
    n_removed = n_before - n_after
  ))
}

filter_log <- tibble(step = character(), n_before = integer(),
                     n_after = integer(), n_removed = integer())
n_start <- nrow(dal$data)
d <- dal

cat(" Step 1: Convert 0s to NA\n")
d <- zero_to_missing(d)

cat(" Step 2: Keratin contaminant removal\n")
n_before <- nrow(d$data)
d        <- filter_proteins_by_annotation(d, !is_keratin)
n_after  <- nrow(d$data)
filter_log <- log_step(filter_log, "Keratin contaminant filter", n_before, n_after)
cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

# --- 4a: Detection Proportion Assessment -------------------------------
cat("\n>> 4a - Detection proportion assessment\n")

# Assess protein detection within each biological group BEFORE applying
# the missingness filter. This allows direct comparison of:
#   1) >=50% detected in every group
#   2) >=5 samples detected in at least 1 group
#
# This is intentionally performed after the keratin filter but before the
# final missingness filter.

assess_detection_proportion <- function(d, grouping_column,
                                        min_reps = 5,
                                        min_prop = 0.50) {
  
  data_mat <- as.matrix(d$data)
  meta_df  <- d$metadata
  
  # Make sure metadata is aligned to data columns
  meta_df <- meta_df[match(colnames(data_mat), meta_df$sample_id), ]
  
  stopifnot(all(meta_df$sample_id == colnames(data_mat)))
  
  groups <- unique(meta_df[[grouping_column]])
  groups <- groups[!is.na(groups)]
  
  # Detection = non-missing value
  detected <- !is.na(data_mat)
  
  # Build one row per protein
  detection_df <- tibble(
    uniprot_id = rownames(data_mat)
  )
  
  # Calculate detection count and proportion for every group
  for (grp in groups) {
    
    sample_idx <- meta_df[[grouping_column]] == grp
    
    n_group <- sum(sample_idx)
    
    detected_n <- rowSums(detected[, sample_idx, drop = FALSE])
    
    detected_prop <- detected_n / n_group
    
    # Make safe column names for group labels
    grp_safe <- make.names(as.character(grp))
    
    detection_df[[paste0(grp_safe, "_n")]] <- detected_n
    detection_df[[paste0(grp_safe, "_total")]] <- n_group
    detection_df[[paste0(grp_safe, "_prop")]] <- detected_prop
    detection_df[[paste0(grp_safe, "_pct")]] <- detected_prop * 100
  }
  
  # Identify proportion columns
  prop_cols <- grep("_prop$", names(detection_df), value = TRUE)
  
  count_cols <- grep("_n$", names(detection_df), value = TRUE)
  
  # Overall detection
  detection_df <- detection_df %>%
    mutate(
      overall_n = rowSums(as.data.frame(detected)),
      overall_total = ncol(data_mat),
      overall_prop = overall_n / overall_total,
      overall_pct = overall_prop * 100,
      
      min_group_prop = apply(
        as.matrix(dplyr::select(., all_of(prop_cols))),
        1,
        min,
        na.rm = TRUE
      ),
      
      max_group_prop = apply(
        as.matrix(dplyr::select(., all_of(prop_cols))),
        1,
        max,
        na.rm = TRUE
      ),
      
      detection_range = max_group_prop - min_group_prop,
      
      # Filter 1:
      # Protein must be detected in >=50% of EVERY group
      pass_prop_50_each =
        min_group_prop >= min_prop,
      
      # Filter 2:
      # Protein must be detected in >=5 samples in AT LEAST ONE group
      pass_group_5_any =
        apply(
          as.matrix(dplyr::select(., all_of(count_cols))),
          1,
          function(x) any(x >= min_reps)
        ),
      
      # Does the protein's inclusion depend on filtering strategy?
      filter_disagreement =
        pass_prop_50_each != pass_group_5_any
    )
  
  # Add a useful interpretation category
  detection_df <- detection_df %>%
    mutate(
      detection_category = case_when(
        
        !pass_group_5_any ~
          "Low detection",
        
        pass_prop_50_each & pass_group_5_any ~
          "Robust detection",
        
        !pass_prop_50_each & pass_group_5_any &
          detection_range >= 0.50 ~
          "Strong group asymmetry",
        
        !pass_prop_50_each & pass_group_5_any ~
          "Moderate group asymmetry",
        
        TRUE ~
          "Other"
      )
    )
  
  # Add annotation
  detection_df <- detection_df %>%
    left_join(
      d$annotation %>%
        dplyr::select(uniprot_id, protein, gene, description),
      by = "uniprot_id"
    ) %>%
    relocate(uniprot_id, protein, gene, description)
  
  detection_df
}


# Run detection assessment BEFORE the missingness filter
detection_proportion <- assess_detection_proportion(
  d = d,
  grouping_column = group_col,
  min_reps = 5,
  min_prop = 0.50
)

cat("\n Detection proportion summary:\n")

cat(sprintf(
  "   Proteins assessed: %d\n",
  nrow(detection_proportion)
))

cat(sprintf(
  "   Pass >=50%% in EVERY group: %d\n",
  sum(detection_proportion$pass_prop_50_each, na.rm = TRUE)
))

cat(sprintf(
  "   Pass >=5 in AT LEAST 1 group: %d\n",
  sum(detection_proportion$pass_group_5_any, na.rm = TRUE)
))

cat(sprintf(
  "   Filter disagreement: %d proteins\n",
  sum(detection_proportion$filter_disagreement, na.rm = TRUE)
))

cat("\n Detection categories:\n")
print(
  detection_proportion %>%
    count(detection_category) %>%
    arrange(desc(n))
)

# Export full protein-level detection assessment
write.csv(
  detection_proportion,
  file.path(
    config$report_dir,
    "03_detection_proportion_assessment_plasma.csv"
  ),
  row.names = FALSE
)

# Export only proteins for which the two filtering approaches disagree
detection_filter_sensitive <- detection_proportion %>%
  filter(filter_disagreement) %>%
  arrange(desc(detection_range))

write.csv(
  detection_filter_sensitive,
  file.path(
    config$report_dir,
    "04_filter_sensitive_proteins_plasma.csv"
  ),
  row.names = FALSE
)

cat(sprintf(
  "\n   Exported %d filter-sensitive proteins.\n",
  nrow(detection_filter_sensitive)
))

cat(" Step 3: Missingness filter\n")
n_before <- nrow(d$data)
d        <- filter_proteins_by_proportion(d, min_prop = 0.5, grouping_column = group_col)
n_after  <- nrow(d$data)
filter_log <- log_step(filter_log, "Missingness filter", n_before, n_after)
cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

filter_log <- bind_rows(
  tibble(step = "Raw input", n_before = NA_integer_,
         n_after = n_start, n_removed = NA_integer_),
  filter_log
) %>% mutate(pct_retained = round(n_after / n_start * 100, 1))

dal <- d
cat("\n   Filtering summary:\n"); print(filter_log)

# --- 5: Outlier Detection & Removal -----------------------------------------
cat("\n>> 5 - Outlier Detection & Removal\n")

run_outlier_detection <- function(d, nm) {
  cat(sprintf("\n [%s] Running outlier detection...\n", nm))
  
  # --- Missingness ---
  pct_missing <- colMeans(is.na(d$data)) * 100
  
  paired_subjects <- d$metadata %>%
    dplyr::filter(cancer == "CRC") %>%
    dplyr::count(pid) %>%
    dplyr::filter(n == 2) %>%
    dplyr::pull(pid)
  
  delta_missing <- tibble(sample_id = character(), pid = character(),
                          pct_missing = numeric(), delta_missing = numeric())
  
  for (subj in paired_subjects) {
    rows  <- d$metadata %>% dplyr::filter(pid == subj)
    t1_id <- rows$sample_id[rows$timepoint == "T1"]
    t2_id <- rows$sample_id[rows$timepoint == "T2"]
    
    if (length(t1_id) == 1 && length(t2_id) == 1) {
      delta_missing <- dplyr::bind_rows(delta_missing, tibble(
        sample_id     = c(t1_id, t2_id),
        pid           = subj,
        pct_missing   = c(pct_missing[t1_id], pct_missing[t2_id]),
        delta_missing = pct_missing[t2_id] - pct_missing[t1_id]
      ))
    }
  }
  
  unpaired_ids <- setdiff(d$metadata$sample_id, delta_missing$sample_id)
  
  if (length(unpaired_ids) > 0) {
    delta_missing <- dplyr::bind_rows(delta_missing, tibble(
      sample_id     = unpaired_ids,
      pid           = d$metadata$pid[match(unpaired_ids, d$metadata$sample_id)],
      pct_missing   = pct_missing[unpaired_ids],
      delta_missing = NA_real_
    ))
  }
  
  # --- Thresholds ---
  miss_threshold <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
  delta_vals     <- delta_missing$delta_missing[!is.na(delta_missing$delta_missing)]
  
  if (length(delta_vals) > 2) {
    delta_threshold <- quantile(delta_vals, 0.75) + 1.5 * IQR(delta_vals)
    delta_lower     <- quantile(delta_vals, 0.25) - 1.5 * IQR(delta_vals)
  } else {
    delta_threshold <- Inf
    delta_lower     <- -Inf
  }
  
  delta_missing <- delta_missing %>%
    dplyr::mutate(
      miss_flag = pct_missing > miss_threshold |
        (!is.na(delta_missing) &
           (delta_missing > delta_threshold | delta_missing < delta_lower))
    )
  
  # --- PCA ---
  data_for_pca <- d$data
  
  for (j in seq_len(ncol(data_for_pca))) {
    nas <- is.na(data_for_pca[, j])
    if (any(nas)) data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
  }
  
  data_log2 <- log2(data_for_pca + 1)
  pca_res   <- prcomp(t(data_log2), center = TRUE, scale. = TRUE)
  
  pc_scores  <- pca_res$x[, 1:3]
  mahal_dist <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))
  
  pca_flags <- tibble(
    sample_id  = colnames(d$data),
    mahal_dist = mahal_dist,
    pca_flag   = mahal_dist > qchisq(0.99, df = 3)
  )
  
  # --- MAD ---
  sample_medians <- apply(log2(d$data + 1), 2, median, na.rm = TRUE)
  global_median  <- median(sample_medians)
  mad_val        <- mad(sample_medians)
  
  mad_flags <- tibble(
    sample_id     = names(sample_medians),
    sample_median = sample_medians,
    mad_deviation = abs(sample_medians - global_median),
    mad_flag      = abs(sample_medians - global_median) > 3 * mad_val
  )
  
  # --- Consensus ---
  # NOTE: plasma's original consensus rule (n_flags >= 2) is preserved here,
  # unlike muscle's stricter n_flags >= 3.
  outlier_diag <- delta_missing %>%
    dplyr::select(sample_id, pct_missing, delta_missing, miss_flag) %>%
    dplyr::left_join(pca_flags, by = "sample_id") %>%
    dplyr::left_join(mad_flags,  by = "sample_id") %>%
    dplyr::mutate(
      n_flags           = miss_flag + pca_flag + mad_flag,
      consensus_outlier = n_flags >= 2
    )
  
  n_outliers <- sum(outlier_diag$consensus_outlier)
  cat(sprintf("   Outliers detected: %d\n", n_outliers))
  
  return(list(
    outlier_diag = outlier_diag,
    pca_res      = pca_res,
    thresholds   = list(
      miss_threshold  = miss_threshold,
      delta_threshold = delta_threshold,
      delta_lower     = delta_lower
    )
  ))
}

run_outlier_removal <- function(d, outlier_diag, nm) {
  n_outliers <- sum(outlier_diag$consensus_outlier)
  
  if (n_outliers > 0) {
    outlier_ids <- outlier_diag %>% filter(consensus_outlier) %>% pull(sample_id)
    cat(sprintf("   [%s] Removing %d outlier sample(s): %s\n",
                nm, length(outlier_ids), paste(outlier_ids, collapse = ", ")))
    keep_idx   <- !colnames(d$data) %in% outlier_ids
    d$data     <- d$data[, keep_idx, drop = FALSE]
    d$metadata <- d$metadata[d$metadata$sample_id %in% colnames(d$data), ]
    cat(sprintf("   After outlier removal: %d samples remain\n", ncol(d$data)))
  } else {
    cat(sprintf("   [%s] No outlier samples removed.\n", nm))
  }
  
  d
}

# Run outlier detection
outlier_result <- run_outlier_detection(dal, "plasma")
outlier_diag   <- outlier_result$outlier_diag

# Save data before outlier removal
dal_pre_outlier <- dal   # snapshot before removal, for QC comparison

# Run outlier removal
dal <- run_outlier_removal(dal, outlier_diag, "plasma")

# --- 6: Normalize Data ----------------------------
cat("\n>> 6 - Normalize Data\n")

saveRDS(dal, file.path(config$data_dir, "00_DAList_pre_norm_plasma.rds"))

write_norm_report(dal, grouping_column = group_col,
                  output_dir = config$report_dir,
                  filename   = "02_norm_comparison_plasma.pdf", overwrite = TRUE)
write_qc_report(dal, color_column = group_col, label_column = "sample_id",
                output_dir = config$report_dir,
                filename   = "01_qc_report_pre_norm_plasma.pdf", overwrite = TRUE)

dal <- normalize_data(dal, norm_method = config$norm_method)

cat("\n Diagnostic: post-normalization data range (to check if already log-scale)\n")
rng <- range(dal$data, na.rm = TRUE)
cat(sprintf("  range: %.2f to %.2f\n", rng[1], rng[2]))

write_qc_report(dal, color_column = group_col, label_column = "sample_id",
                output_dir = config$report_dir,
                filename   = "03_qc_report_post_norm_plasma.pdf", overwrite = TRUE)

cat(sprintf("  Final: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

pca_postnorm <- compute_pca(dal$data, log_transform = FALSE) # TRUE if diagnostic shows linear scale

# --- 6b: Normalization Quality Scoring ---------------
cat("\n>> 6b - Normalization quality scoring\n")

compute_pcv <- function(mat, groups) {
  cvs <- unlist(lapply(unique(groups), function(grp) {
    sub <- mat[, groups == grp, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2 || mean(x) == 0) return(NA_real_)
      sd(x) / abs(mean(x))
    })
  }))
  median(cvs, na.rm = TRUE)
}

compute_pmad <- function(mat, groups) {
  mads <- unlist(lapply(unique(groups), function(grp) {
    sub <- mat[, groups == grp, drop = FALSE]
    apply(sub, 1, function(x) {
      x <- x[!is.na(x)]
      if (length(x) < 2) return(NA_real_)
      mad(x, constant = 1)
    })
  }))
  median(mads, na.rm = TRUE)
}

compute_cor <- function(mat, groups) {
  cors <- unlist(lapply(unique(groups), function(grp) {
    idx <- which(groups == grp)
    if (length(idx) < 2) return(NULL)
    cm <- cor(mat[, idx], use = "pairwise.complete.obs")
    cm[lower.tri(cm)]
  }))
  mean(cors, na.rm = TRUE)
}

NORMS <- c("log2", "median", "mean", "vsn", "quantile", "cycloess", "rlr", "gi")

cat("\n  Scoring normalisation methods...\n")

pre_dal <- readRDS(file.path(config$data_dir, "00_DAList_pre_norm_plasma.rds"))
grps    <- factor(pre_dal$metadata[[group_col]])

norm_scores <- purrr::map_dfr(NORMS, function(nm) {
  dal_n <- tryCatch(
    normalize_data(pre_dal, norm_method = nm),
    error = function(e) {
      cat(sprintf("    %s: FAILED (%s)\n", nm, e$message))
      NULL
    }
  )
  if (is.null(dal_n)) return(NULL)
  
  mat_n <- as.matrix(dal_n$data)
  
  pcv     <- compute_pcv(mat_n, grps)
  pmad    <- compute_pmad(mat_n, grps)
  cor_val <- compute_cor(mat_n, grps)
  
  cat(sprintf("    %s: PCV=%.4f  PMAD=%.4f  COR=%.4f\n", nm, pcv, pmad, cor_val))
  tibble(norm = nm,
         PCV = round(pcv, 4), PMAD = round(pmad, 4), COR = round(cor_val, 4))
})

norm_scores <- norm_scores |>
  mutate(
    PCV_rank        = rank(PCV),
    PMAD_rank       = rank(PMAD),
    COR_rank        = rank(-COR),
    norm_composite  = (PCV_rank + PMAD_rank + COR_rank) / 3
  ) |>
  arrange(norm_composite)

cat("\n  Normalisation rankings:\n")
print(norm_scores |> dplyr::select(norm, PCV, PMAD, COR, norm_composite))

write_csv(norm_scores,
          file.path(config$report_dir, "10_norm_quality_scores_plasma.csv"))

# --- 7: Export Normalized Data ----------------------
cat("\n>> 7 - Export Normalized Data\n")

export <- bind_cols(
  dal$annotation %>% dplyr::select(uniprot_id, protein, gene, description),
  dal$data
)

write.csv(export, file.path(config$data_dir, "01_normalized_data_plasma.csv"))
saveRDS(dal,       file.path(config$data_dir, "02_normalized_data_plasma.RDS"))

cat(sprintf("  Exported: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 8: Store Data for Reports --------------
report_data <- list(
  config                 = config,
  flagged_contaminants   = flagged_contaminants,
  outlier_diag           = outlier_diag,
  outlier_pca            = outlier_result$pca_res,
  outlier_thresholds     = outlier_result$thresholds,
  filter_log             = filter_log,
  dal_pre_outlier        = dal_pre_outlier,
  sample_summary         = sample_summary,
  dal                    = dal,
  dal_raw                = dal_raw,
  missingness_raw        = missingness_raw,
  pca_raw                = pca_raw,
  pca_postnorm           = pca_postnorm,
  norm_scores            = norm_scores
)

saveRDS(report_data, file.path(config$data_dir, "00_report_data.rds"))

writeLines(capture.output(sessionInfo()), file.path(config$data_dir, "sessionInfo.txt"))