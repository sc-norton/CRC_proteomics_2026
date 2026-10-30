# === PROJECT: CR Proteomics - Muscle =============================
# === Script: Normalization 

# === PURPOSE: To load the raw data and metadata, assemble a DAList for use in
# proteoDA, filter the data for contaminants and missingness, and normalize the
# for downstream analyses.

# === INPUTS: raw data, meta data, HPA muscle data, 

# === OUTPUTS: (for each dataset, CREvPLA and SURVvCTL)
# raw sample summary, normalization reports, flagged contaminants,
# pre & post QC reports, filtering effects, outlier diagnostics, final protein
# report, final sample report, normalized data (csv), normalized DA lists

# --- 0: Setup --------------------------------------------
cat("\n>> 0 - Setup\n")
cat("\n Loading required packages...")

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, org.Hs.eg.db, AnnotationDbi, tibble, ggpubr,
               Hmisc, ggnewscale)

cat("\n Packages successfully loaded.\n")

# Set working directory
setwd(rprojroot::find_rstudio_root_file())

# Configuration
config <- list(
  # Data files - muscle
  raw_muscle = "00_input/muscle/CRm_raw.csv",
  meta_muscle = "00_input/muscle/CRm_meta.csv",
  hpa_muscle = "00_input/muscle/HPA_muscle.csv",
  go_map_muscle = "00_input/muscle/GO_map.csv",
  
  # Data files - plasma
  raw_plasma = "00_input/plasma/CRp_raw.csv",
  meta_plasma = "00_input/plasma/CRp_meta.csv",
  
  # Directories
  report_dir = "01_normalization/b_reports",
  data_dir = "01_normalization/c_data",
  
  # Parameters
  miss_pct = 50,
  outlier_k = 3,
  mahal_p = 0.01,
  mad_k = 3,
  norm_method = "cycloess",
  B_M_thresholds = c(5, 10, 15, 20, Inf)
)

dir.create(config$report_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

build_output_path <- function(nm, base_dir) {
  
  dataset <- ifelse(grepl("^SURV", nm), "SURV", "CRE")
  fl      <- stringr::str_extract(nm, "none|moderate|stringent")
  bm      <- stringr::str_extract(nm, "BM[^_]+")
  
  out_dir <- file.path(base_dir, dataset, fl, bm)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(out_dir)
}

cat("\n Base directory: ", getwd(), "\n")

# --- 1: Load & Restructure Data ---------------------------------------
cat("\n>> 1 - Load & Restructure Data\n")

# Load Data
raw <- read.csv(file.path(config$raw_muscle))
meta <- read.csv(file.path(config$meta_muscle))
go_map <- read.csv(file.path(config$go_map_muscle))

# Rename annotation column names
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq", "sum_spectra")
colnames(raw)[1:6] <- annot_cols

# Make Uniprot_ID the row name
rownames(raw) <- raw$uniprot_id

# Split input data into protein intensity data and annotation data
intensity <- raw[ , setdiff(names(raw), annot_cols)] 
intensity.b <- raw[ , setdiff(names(raw), annot_cols)]
rownames(intensity.b) <- make.unique(raw$gene)
annotation <- raw[ ,annot_cols] 

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))
cat(sprintf("\n Meta data: %d samples\n", nrow(meta)))

# Join GO_map into annotation
annotation <- annotation %>%
  left_join(go_map, by = c("uniprot_id" = "From"))

# Align row names of metadata with column names of data
cat(" Aligning data and metadata...\n")
rownames(meta) <- meta[ ,1]

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
                                  count(cancer, 
                                        supp, 
                                        timepoint, 
                                        cancer_time, 
                                        supp_time))

write.csv(sample_summary, file.path(config$report_dir, "00_raw_sample_summary_muscle.csv"))

# --- 2a: Contaminant Assessment----------------------------------
cat("\n>> 2 - Contaminant Assessment\n")

# --- A: Human Protein Atlas - Skeletal Muscle
# Load Human Protein Atlas Muscle Protein data
HPAm <- read.csv(file.path(config$hpa_muscle))

# Select the Gene, Ensembl, and Uniprot columns
hpa_lookup <- HPAm %>%
  dplyr::select(Gene, Ensembl, Uniprot) %>% 
  distinct(Uniprot, .keep_all = T)

# Join to annotation data
annotation <- annotation %>%
  left_join(hpa_lookup, annotation, by = c("uniprot_id" = "Uniprot")) %>%
  mutate(in_hpa_muscle = uniprot_id %in% hpa_lookup$Uniprot)

n_flagged_hpa <- sum(!annotation$in_hpa_muscle)

cat(sprintf("\n Proteins flagged as contaminants (HPA): %d\n", n_flagged_hpa))

# --- 2c: Blood Protein Blacklist Filter ---
blood_blacklist_moderate <- c(
  # Albumin & transport
  "ALB", "TTR", "APOA1", "APOA2", "APOB", "APOC1", "APOC2", "APOC3", 
  "APOE", "APOH", "CLU",
  # Coagulation
  "FGA", "FGB", "FGG", "F2", "F5", "F7", "F10", "F13A1", "F13B",
  "SERPINC1", "SERPIND1", "PROS1", "PROC", "VWF",
  # Complement
  "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C4BPA", "C4BPB",
  "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  # Hemoglobin
  "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2",
  # Immunoglobulins (Fc regions / constant chains)
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM",
  "IGKC", "IGLC1", "IGLC2",
  # Acute phase / other plasma
  "CRP", "SAA1", "SAA2", "ORM1", "ORM2", "HP", "HPX", "TF",
  "A1BG", "A2M", "GC", "PLG", "KNG1", "ITIH1", "ITIH2", "ITIH3", "ITIH4",
  # Platelet markers 
  "PPBP", "PF4", "THBS1", "SELP"
)

library(dplyr)

# Dynamically pull true contaminants from annotation data
contaminant_immunoglobulins <- annotation %>% 
  filter(grepl("^IGH|^IGK|^IGL", gene)) %>% 
  pull(gene)

contaminant_keratins <- annotation %>% 
  filter(grepl("^KRT", gene)) %>% 
  pull(gene)

# Assemble the blood blacklist
blood_blacklist_stringent <- c(
  # Albumin & transport
  "ALB", "TTR", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3", 
  "APOC4", "APOD", "APOE", "APOH", "APOL1", "APOL2", "APOM", "CLU", "TF", "GC",
  # Coagulation
  "FGA", "FGB", "FGG", "F2", "F5", "F7", "F10", "F13A1", "F13B", "SERPINC1", 
  "SERPIND1", "PROS1", "PROC", "VWF", "PLG", "KNG1",
  # Complement System
  "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C4BPA", "C4BPB", "C5", "C6", 
  "C7", "C8A", "C8B", "C8G", "C9", "CFH", "CFB", "CFI", "SERPING1",
  # True Hemoglobins (explicitly named to avoid breaking growth factors)
  "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2", "HBM", "HBZ", "HBE1", "HBQ1",
  # Acute phase / Inter-alpha-trypsin inhibitors
  "CRP", "SAA1", "SAA2", "ORM1", "ORM2", "HP", "HPX", "A1BG", "A2M", "ITIH1", 
  "ITIH2", "ITIH3", "ITIH4", "SERPINA1",
  # Platelet / Neutrophil markers (Remove S100A8/A9 if studying active tissue injury)
  "PPBP", "PF4", "THBS1", "SELP", "LYZ", "DEFA3", "JCHAIN",
  # Dynamic vectors
  contaminant_immunoglobulins,
  contaminant_keratins
) %>% unique() # Ensures no double-counting


annotation <- annotation %>%
  mutate(in_moderate = gene %in% blood_blacklist_moderate,
         in_stringent = gene %in% blood_blacklist_stringent)

# Export contaminants
flagged_hpa <- annotation %>%
  filter(!in_hpa_muscle) %>%
  mutate(flag_source = "HPA")

flagged_blood_moderate <- annotation %>%
  filter(in_moderate) %>%
  mutate(flag_source = "Blood - moderate")

flagged_blood_stringent <- annotation %>%
  filter(in_stringent) %>%
  mutate(flag_source = "Blood - stringent")

# --- 2d: Define Filter Levels ------------------------------------
cat("\n>> 2d - Define Filter Levels\n")

# Level 1: No blood filtering
blood_filter_none <- c()

# Level 2: Moderate — most abundant/definite blood proteins only
blood_filter_moderate <- unique(blood_blacklist_moderate)

# Level 3: Stringent — full blacklist + HPA confident blood proteins
blood_filter_stringent <- unique(blood_blacklist_stringent)

filter_levels <- list(
  none       = blood_filter_none,
  moderate   = blood_filter_moderate,
  stringent  = blood_filter_stringent
)

annotation <- annotation %>%
  mutate(
    in_blood_none       = FALSE,
    in_blood_moderate   = gene %in% blood_filter_moderate,
    in_blood_stringent  = gene %in% blood_filter_stringent
  )

cat(sprintf("  None:       %d proteins to remove\n", length(blood_filter_none)))
cat(sprintf("  Moderate:   %d proteins to remove\n", length(blood_filter_moderate)))
cat(sprintf("  Stringent:  %d proteins to remove\n", length(blood_filter_stringent)))

# --- 2e: Contamination indexing ------------------------------------
cat("\n>> 2e - Contamination index\n")

blood_markers  <- c("ALB", "HBB", "HBA1")
muscle_markers <- c("MB", "CKM")
all_markers    <- c(blood_markers, muscle_markers)

marker_rows <- which(rownames(intensity.b) %in% all_markers)

if (length(marker_rows) > 0) {
  marker_intensity <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% all_markers], ], na.rm = TRUE)
  blood_intensity  <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% blood_markers], ], na.rm = TRUE)
  muscle_intensity <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% muscle_markers], ], na.rm = TRUE)
  total_intensity  <- colSums(intensity.b, na.rm = TRUE)
  blood_pct        <- (blood_intensity / total_intensity) * 100
  muscle_pct       <- (muscle_intensity / total_intensity) * 100
  B_M_ratio        <- (blood_pct / muscle_pct)
  
  
  marker_pct_df <- as.data.frame(t(intensity.b[marker_rows, ])) %>%
    mutate(across(everything(), ~ . / total_intensity * 100)) %>%
    rename_with(~ paste0(.x, "_pct")) %>%
    mutate(sample_id = colnames(intensity.b),
           blood_pct = blood_pct,
           muscle_pct = muscle_pct,
           B_M_ratio = B_M_ratio)
  
  cat(sprintf("  Blood %% of total intensity — min: %.2f%%, median: %.2f%%, max: %.2f%%\n",
              min(blood_pct), median(blood_pct), max(blood_pct)))
  
  cat(sprintf("  Muscle %% of total intensity — min: %.2f%%, median: %.2f%%, max: %.2f%%\n",
              min(muscle_pct), median(muscle_pct), max(muscle_pct)))
  
} else {
  stop("Markers not found in rownames of intensity.b — check gene names.")
}

# Join to metadata
meta <- meta %>%
  left_join(marker_pct_df, by = "sample_id")
rownames(meta) <- meta$sample_id

# --- 2f: Blood Score (correlation) ----------------------------
blood_score <- colMeans(intensity.b[rownames(intensity.b) %in% blood_markers, ], na.rm = TRUE)

cor_with_blood <- apply(intensity.b, 1, function(x) cor(x, blood_score, use = "pairwise"))
cor_df <- data.frame(
  gene = rownames(intensity.b),
  cor = cor_with_blood
)

top_genes <- cor_df %>%
  arrange(desc(abs(cor))) %>%
  dplyr::slice(1:50) %>%
  pull(gene)

# --- 2g: Top 5 Abundance proteins per sample ---------------------
top5_per_sample <- apply(intensity.b, 2, function(x) {
  names(sort(x, decreasing = TRUE, na.last = TRUE))[1:5]
})

top5_per_sample <- apply(intensity.b, 2, function(x) {
  idx <- order(x, decreasing = TRUE, na.last = TRUE)[1:5]
  tibble(gene = rownames(intensity.b)[idx],
         intensity = x[idx],
         rank = 1:5)
}) %>%
  bind_rows(.id = "sample_id")

# --- 3: Assemble DALists -----------------------------------------
cat("\n>> 3 - Assemble DAList \n")

# Check for duplicate uniprot_ids
cat(" Checking for duplicate uniprot ids...\n")
dup_ids <- annotation$uniprot_id[duplicated(annotation$uniprot_id)]

if (length(dup_ids) > 0) {
  cat(sprintf("   Found %d duplicate uniprot_ids — deduplicating by highest mean intensity.\n",
              length(dup_ids)))
  
  intensity_num <- as.data.frame(lapply(intensity, as.numeric))
  annotation$row_mean <- rowMeans(intensity_num, na.rm = TRUE)
  
  keep_idx <- annotation %>%
    mutate(row_idx = row_number()) %>%
    group_by(uniprot_id) %>%
    slice_max(row_mean, n = 1, with_ties = FALSE) %>%
    pull(row_idx)
  
  annotation <- annotation[keep_idx, ]
  intensity  <- intensity[keep_idx, ]
  annotation$row_mean <- NULL
  
  cat(sprintf("   After deduplication: %d unique proteins\n", nrow(annotation)))
} else {
  cat("   No duplicate uniprot_ids found.\n")
}

# Assemble DAList (all samples)
dal <- DAList(
  data       = intensity,
  metadata   = meta,
  annotation = annotation
)

# DAList (w/out CTL)
dal2 <- filter_samples(dal, dal$metadata$cancer != "CTL")

# Combine
# Create all combinations: dataset (SURV/CRE) x filter level (none/moderate/stringent)
dals <- list()

for (fl in names(filter_levels)) {
  for (bm in config$B_M_thresholds) {
    
    bm_label <- if (is.infinite(bm)) "none" else as.character(bm)
    nm       <- paste0(fl, "_BM", bm_label)
    
    anno_fl <- annotation %>%
      mutate(in_blood_blacklist = gene %in% filter_levels[[fl]])
    
    dal_fl  <- DAList(data = intensity, metadata = meta, annotation = anno_fl)
    dal2_fl <- filter_samples(dal_fl, dal_fl$metadata$cancer != "CTL")
    
    # store DAList + params together
    dals[[paste0("SURV_", nm)]] <- list(
      dal = dal_fl,
      blood_filter = fl,
      bm_thresh = bm
    )
    
    dals[[paste0("CRE_", nm)]] <- list(
      dal = dal2_fl,
      blood_filter = fl,
      bm_thresh = bm
    )
  }
}

cat(sprintf(" DAList assembled: %d proteins x %d samples \n", 
            nrow(dal$data),
            ncol(dal$data)))

cat(sprintf(" DAList2 assembled: %d proteins x %d samples \n", 
            nrow(dal2$data),
            ncol(dal2$data)))

# --- 4: Filter Data -------------------------------------
cat("\n>> 4 - Filter Data \n")

grouping_cols <- setNames(
  ifelse(grepl("^SURV", names(dals)), "cancer_time", "supp_time"),
  names(dals)
)

# Filter step logging function
log_step <- function(log, step, n_before, n_after) {
  bind_rows(log, tibble(
    step = step, 
    n_before = n_before, 
    n_after = n_after, 
    n_removed = n_before - n_after))
}

filter_by_bm_ratio <- function(d, bm_thresh, nm = NULL) {
  
  meta_df <- d$metadata
  
  # Build logical index safely
  keep_idx <- is.na(meta_df$B_M_ratio) | meta_df$B_M_ratio <= bm_thresh
  
  if (!is.logical(keep_idx) || length(keep_idx) != nrow(meta_df)) {
    stop(sprintf("[%s] Invalid keep_idx in BM filtering", nm))
  }
  
  keep_samples <- meta_df$sample_id[keep_idx]
  
  n_before <- ncol(d$data)
  
  # Apply filtering explicitly
  d$data     <- d$data[, colnames(d$data) %in% keep_samples, drop = FALSE]
  d$metadata <- d$metadata[d$metadata$sample_id %in% colnames(d$data), ]
  
  n_after <- ncol(d$data)
  
  cat(sprintf("   %d -> %d samples (%d removed; BM <= %s)\n",
              n_before, n_after, n_before - n_after, bm_thresh))
  
  return(d)
}

f.dals <- lapply(names(dals), function(nm) {
  entry <- dals[[nm]]
  d     <- entry$dal
  grp   <- grouping_cols[[nm]]
  
  fl        <- entry$blood_filter
  bm_thresh <- entry$bm_thresh
  
  if (is.null(bm_thresh) || is.na(bm_thresh)) {
    stop(sprintf("[%s] Invalid BM threshold", nm))
  }
  
  log     <- tibble(step = character(), n_before = integer(),    
                    n_after = integer(), n_removed = integer())
  n_start <- nrow(d$data)                                        
  
  # Step 1: Convert 0s to NA
  cat(sprintf(" [%s] Step 1: Convert 0s to NA\n", nm))
  d <- zero_to_missing(d)
  
  # Step 2: HPA muscle filter (all levels)
  cat(sprintf(" [%s] Step 2: HPA muscle contaminant removal\n", nm))
  n_before <- nrow(d$data)
  d <- filter_proteins_by_annotation(d, in_hpa_muscle)
  n_after  <- nrow(d$data)
  log <- log_step(log, "HPA muscle filter", n_before, n_after)
  cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  # Step 2b: Blood filter — conditional on filter level
  if (fl != "none") {
    cat(sprintf(" [%s] Step 2b: Blood filter (%s)\n", nm, fl))
    n_before <- nrow(d$data)
    d <- filter_proteins_by_annotation(d, !in_blood_blacklist)
    n_after  <- nrow(d$data)
    log <- log_step(log, paste0("Blood filter (", fl, ")"), n_before, n_after)
    cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  } else {
    cat(sprintf(" [%s] Step 2b: Blood filter skipped (none)\n", nm))
  }
  
  # Step 2c: Blood contamination filter (B:M ratio thresholding)
  cat(sprintf(" [%s] Step 2c: Blood contamination filter (B:M ratio)\n", nm))
  
  d <- filter_by_bm_ratio(d, bm_thresh, nm)
  
  # Step 3: Missingness filter
  cat(sprintf(" [%s] Step 3: Missingness filter\n", nm))
  n_before <- nrow(d$data)
  d <- filter_proteins_by_proportion(d, min_prop = 0.5, grouping_column = grp)
  n_after  <- nrow(d$data)
  log <- log_step(log, "Missingness filter", n_before, n_after)
  cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))
  
  log <- bind_rows(
    tibble(step = "Raw input", n_before = NA_integer_,
           n_after = n_start, n_removed = NA_integer_),
    log
  ) %>% mutate(pct_retained = round(n_after / n_start * 100, 1))
  
  list(dal = d, log = log)
})

names(f.dals) <- names(dals)
dals        <- lapply(f.dals, `[[`, "dal")
filter_logs <- lapply(f.dals, `[[`, "log")

# Print filter summaries
cat("\n   Filtering summary:\n")
lapply(names(filter_logs), function(nm) {
  cat(sprintf("\n  -- %s --\n", nm))
  print(filter_logs[[nm]])
})

# --- 5: Outlier Detection & Removal -----------------------------------------
cat("\n>> 6 - Outlier Detection & Removal\n")

run_outlier_detection <- function(d, nm) {
  cat(sprintf("\n [%s] Running outlier detection...\n", nm))
  
  # --- Missingness ---
  pct_missing <- colMeans(is.na(d$data)) * 100
  
  paired_subjects <- d$metadata %>%
    dplyr::filter(cancer == "SURV") %>%
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
  miss_threshold  <- quantile(pct_missing, 0.75) + 1.5 * IQR(pct_missing)
  delta_vals      <- delta_missing$delta_missing[!is.na(delta_missing$delta_missing)]
  
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
    if (any(nas)) {
      data_for_pca[nas, j] <- median(data_for_pca[, j], na.rm = TRUE)
    }
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
  outlier_diag <- delta_missing %>%
    dplyr::select(sample_id, pct_missing, delta_missing, miss_flag) %>%
    dplyr::left_join(pca_flags, by = "sample_id") %>%
    dplyr::left_join(mad_flags,  by = "sample_id") %>%
    dplyr::mutate(
      n_flags = miss_flag + pca_flag + mad_flag,
      consensus_outlier = n_flags >= 3
    )
  
  n_outliers <- sum(outlier_diag$consensus_outlier)
  cat(sprintf("   Outliers detected: %d\n", n_outliers))
  
  # Return everything needed for reporting
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

# Run Outlier Detection
outlier_diag <- lapply(names(dals),
                       function(nm) run_outlier_detection(dals[[nm]], nm)) |>
  setNames(names(dals))

# Run Outlier Removal

# Save data before outlier removal for quality control
dals_pre_outlier <- dals

dals <- lapply(names(dals), 
               function(nm) run_outlier_removal(dals[[nm]], outlier_diag[[nm]]$outlier_diag, nm)) |> 
  setNames(names(dals))

# --- 6: Normalize Data ----------------------------
cat("\n>> 10 - Normalize Data\n")

# Save pre-normalized DALists
lapply(names(dals), function(nm) {
  out_dir <- build_output_path(nm, config$data_dir)
  saveRDS(
    dals[[nm]],
    file.path(out_dir, sprintf("00_DAList_pre_norm_%s_muscle.rds", nm))
  )
})

# Pre-normalization reports
lapply(names(dals), function(nm) {
  report_dir <- build_output_path(nm, config$report_dir)
  write_norm_report(dals[[nm]], grouping_column = grouping_cols[nm],
                    output_dir = report_dir,
                    filename   = sprintf("02_norm_comparison_%s_muscle.pdf", nm),
                    overwrite  = TRUE)
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample_id",
                  output_dir = report_dir,
                  filename   = sprintf("01_qc_report_pre_norm_%s_muscle.pdf", nm),
                  overwrite  = TRUE)
})

# Normalize
dals <- lapply(dals, normalize_data, norm_method = config$norm_method)

# Post-norm reports
lapply(names(dals), function(nm) {
  report_dir <- build_output_path(nm, config$report_dir)
  write_qc_report(dals[[nm]], color_column = grouping_cols[nm], label_column = "sample_id",
                  output_dir = report_dir,
                  filename   = sprintf("03_qc_report_post_norm_%s_muscle.pdf", nm),
                  overwrite  = TRUE)
})

lapply(names(dals), function(nm) {
  cat(sprintf("  [%s] Final: %d proteins x %d samples\n", nm,
              nrow(dals[[nm]]$data), ncol(dals[[nm]]$data)))
})

# --- 7: Export Normalized Data ----------------------
cat("\n>> 13 - Export Normalized Data\n")

invisible(lapply(names(dals), function(nm) {
  d <- dals[[nm]]
  data_dir <- build_output_path(nm, config$data_dir)
  export <- bind_cols(
    d$annotation %>% dplyr::select(uniprot_id, protein, gene, description, Ensembl),
    d$data
  )
  
  write.csv(export, file.path(data_dir, sprintf("01_normalized_data_%s_muscle.csv", nm)))
  saveRDS(d,        file.path(data_dir, sprintf("02_DAList_normalized_%s_muscle.RDS", nm)))
  
  cat(sprintf("  [%s] Exported: %d proteins x %d samples\n", nm, nrow(d$data), ncol(d$data)))
}))

# --- 8: Store Data for Reports --------------
report_data <- list(
  config = config,
  flagged_hpa = flagged_hpa,
  flagged_blood_moderate = flagged_blood_moderate,
  flagged_blood_stringent = flagged_blood_stringent,
  marker_pct_df = marker_pct_df,
  outlier_diag = outlier_diag,
  filter_logs = filter_logs,
  dals_pre_outlier = dals_pre_outlier,
  sample_summary = sample_summary,
  dals = dals
)

saveRDS(report_data, file.path(config$data_dir, "00_report_data.rds"))

writeLines(capture.output(sessionInfo()), file.path(config$data_dir, "sessionInfo.txt"))
