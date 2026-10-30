# === PROJECT: CR Proteomics - Muscle =============================
# === Script: Normalization

# === PURPOSE: To load the raw data and metadata, assemble a DAList for use in
# proteoDA, filter the data for contaminants and missingness, and normalize
# for downstream analyses.

# === INPUTS: raw data, meta data, HPA muscle data,

# === OUTPUTS: raw sample summary, normalization reports, flagged contaminants,
# pre & post QC reports, filtering effects, outlier diagnostics, final protein
# report, final sample report, normalized data (csv), normalized DA list

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
getwd()

# Configuration
config <- list(
  # Data files - muscle
  raw_muscle    = "00_input/muscle/CRm_raw.csv",
  meta_muscle   = "00_input/muscle/CRm_meta.csv",
  hpa_muscle    = "00_input/muscle/HPA_muscle_annotations.csv",
  go_map_muscle = "00_input/muscle/GO_map.csv",
  
  # Directories
  report_dir = "01_normalization/b_reports/muscle",
  data_dir   = "01_normalization/c_data/muscle",
  
  # Parameters
  miss_pct    = 50,
  outlier_k   = 3,
  mahal_p     = 0.01,
  mad_k       = 3,
  norm_method = "cycloess",
  bm_thresh   = 15
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
raw    <- read.csv(file.path(config$raw_muscle))
meta   <- read.csv(file.path(config$meta_muscle))
go_map <- read.csv(file.path(config$go_map_muscle))

# Rename annotation column names
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq", "sum_spectra")
colnames(raw)[1:6] <- annot_cols

# Make Uniprot_ID the row name
rownames(raw) <- raw$uniprot_id

# Split input data into protein intensity data and annotation data
intensity   <- raw[ , setdiff(names(raw), annot_cols)]
intensity.b <- raw[ , setdiff(names(raw), annot_cols)]
rownames(intensity.b) <- make.unique(raw$gene)
annotation  <- raw[ , annot_cols]

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))
cat(sprintf("\n Meta data: %d samples\n", nrow(meta)))

# Join GO_map into annotation
annotation <- annotation %>%
  left_join(go_map, by = c("uniprot_id" = "From"))

# Rename metadata variables: SURV -> CRC (cancer survivor group label)
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

write.csv(sample_summary, file.path(config$report_dir, "00_raw_sample_summary_muscle.csv"))

# --- 2a: HPA-Based Contaminant Classification ------------------------------
cat("\n>> 2a - HPA-based contaminant classification\n")

HPAm <- read.csv(file.path(config$hpa_muscle))

hpa_lookup <- HPAm %>%
  dplyr::select(Gene, Ensembl, Uniprot) %>%
  distinct(Uniprot, .keep_all = TRUE)

annotation <- annotation %>%
  left_join(hpa_lookup, by = c("uniprot_id" = "Uniprot")) %>%
  mutate(in_hpa_muscle = uniprot_id %in% hpa_lookup$Uniprot)

cat(sprintf("  Proteins not found in HPA muscle reference: %d\n", sum(!annotation$in_hpa_muscle)))

# --- 2b: Curated Blood Blacklist (external knowledge tier) -----------------
cat("\n>> 2b - Curated blood blacklist\n")

# Dynamically pull true contaminants from annotation data
contam_ig <- annotation %>%
  filter(grepl("^IGH|^IGK|^IGL", gene)) %>%
  pull(gene)

contam_keratin <- annotation %>%
  filter(grepl("^KRT", gene)) %>%
  pull(gene)

contam_apo <- annotation %>%
  filter(grepl("^Apolipoprotein", description)) %>%
  pull(gene)

contam_hb <- annotation %>%
  filter(grepl("^Hemoglobin", description)) %>%
  pull(gene)

contam_ca <- c(
  # carbonic anhydrase 1 & 2
  "CA1", "CA2"
)

contam_etc <- c(
  # Hemopexin
  "HPX",
  # Platelet markers
  "JCHAIN"
)

contam_comp <- annotation %>%
  filter(grepl("^Complement", description)) %>%
  pull(gene)

contam_coag <- annotation %>%
  filter(grepl("^Coag", description)) %>%
  pull(gene)

contam_blood_top14 <- c(
  # 1. Albumin
  "ALB", 
  # 2. Immunoglobulin G (IgG)
  "IGHG1", "IGHG2", "IGHG3", "IGHG4", 
  # 3. Immunoglobulin A (IgA)
  "IGHA1", "IGHA2", 
  # 4. Immunoglobulin M (IgM)
  "IGHM", 
  # 5. Immunoglobulin D (IgD)
  "IGHD", 
  # 6. Immunoglobulin E (IgE)
  "IGHE", 
  # 7. Immunoglobulin kappa & lambda light chains
  "IGKC", "IGLC1", "IGLC2", "IGLC3", "IGLC6", "IGLC7", 
  # 8. Transferrin
  "TF", 
  # 9. Fibrinogen
  "FGA", "FGB", "FGG",
  #10. Alpha-2-macroglobulin
  "A2M",
  #11. Alpha-1-antitrypsin
  "SERPINA1",
  #12. Haptoglobin
  "HP",
  #13. Alpha-1-acid Glycoprotein
  "ORM1", "ORM2",
  #14. Apolipoprotein A-1
  "APOA1"
)

new_contam = c(contam_ig, contam_keratin, contam_blood_top14, contam_ca,
               contam_apo, contam_hb, contam_etc, contam_comp, contam_coag)

blood_blacklist <- c(
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
  # Platelet / Neutrophil markers
  "PPBP", "PF4", "THBS1", "SELP", "LYZ", "DEFA3", "JCHAIN",
  # Others
  "AHSG", "AHSP",
  # Dynamic vectors
  contam_ig
) %>% unique()

annotation <- annotation %>%
  mutate(in_blood_blacklist = gene %in% blood_blacklist)

annotation <- annotation %>%
  mutate(is_keratin = gene %in% contam_keratin)

# Export contaminants (diagnostic exports -- not filtering yet)
flagged_hpa <- annotation %>%
  filter(!in_hpa_muscle) %>%
  mutate(flag_source = "HPA")

flagged_keratin <- annotation %>%
  filter(is_keratin) %>%
  mutate(flag_scource = "Keratin")

flagged_blood <- annotation %>%
  filter(in_blood_blacklist) %>%
  mutate(flag_source = "Blood")

cat(sprintf("  Blood blacklist: %d proteins flagged\n", sum(annotation$in_blood_blacklist)))

# --- 2c: Per-Sample Contamination Index -------------------------------------
cat("\n>> 2c - Per-sample contamination index\n")

blood_markers  <- c("ALB", "HBB", "HBA1")
muscle_markers <- c("MB", "CKM")
all_markers    <- c(blood_markers, muscle_markers)

marker_rows <- which(rownames(intensity.b) %in% all_markers)

if (length(marker_rows) > 0) {
  marker_intensity <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% all_markers], ], na.rm = TRUE)
  blood_intensity  <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% blood_markers], ],  na.rm = TRUE)
  muscle_intensity <- colSums(intensity.b[marker_rows[rownames(intensity.b)[marker_rows] %in% muscle_markers], ], na.rm = TRUE)
  total_intensity  <- colSums(intensity.b, na.rm = TRUE)
  blood_pct        <- (blood_intensity / total_intensity) * 100
  muscle_pct       <- (muscle_intensity / total_intensity) * 100
  B_M_ratio        <- (blood_pct / muscle_pct)
  
  marker_pct_df <- as.data.frame(t(intensity.b[marker_rows, ])) %>%
    mutate(across(everything(), ~ . / total_intensity * 100)) %>%
    rename_with(~ paste0(.x, "_pct")) %>%
    mutate(sample_id  = colnames(intensity.b),
           blood_pct  = blood_pct,
           muscle_pct = muscle_pct,
           B_M_ratio  = B_M_ratio)
  
  marker_pct_df <- marker_pct_df %>%
    mutate(blood_log  = log1p(blood_pct),
           blood_z    = as.numeric(scale(blood_log)))
  
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

# --- 2d: Blood Correlation Index (data-driven tier) -------------------------
cat("\n>> 2d - Blood correlation index\n")

li <- log2(as.matrix(intensity.b) + 1)
li[!is.finite(li)] <- NA

anchor_rows  <- rownames(intensity.b) %in% blood_markers
blood_index  <- colMeans(li[anchor_rows, , drop = FALSE], na.rm = TRUE)
n_obs        <- rowSums(!is.na(li))

blood_cor <- suppressWarnings(as.vector(
  cor(t(li), blood_index, method = "spearman", use = "pairwise.complete.obs")
))
blood_cor[n_obs < 15] <- NA
names(blood_cor) <- rownames(intensity.b)

cor_df <- data.frame(gene = names(blood_cor), cor = blood_cor)

top_genes <- cor_df %>%
  arrange(desc(abs(cor))) %>%
  dplyr::slice(1:50) %>%
  pull(gene)

gene_by_uniprot <- setNames(annotation$gene, annotation$uniprot_id)
cor_by_gene     <- cor_df$cor[match(gene_by_uniprot[annotation$uniprot_id], cor_df$gene)]

annotation <- annotation %>%
  mutate(
    blood_cor            = cor_by_gene,
    is_blood_correlated  = coalesce(blood_cor > 0.45, FALSE)
  )

cat(sprintf("  Blood-correlated (data-driven, ungated): %d proteins flagged\n",
            sum(annotation$is_blood_correlated, na.rm = TRUE)))

# --- 2f: Top 5 Abundance proteins per sample --------------------------------
top5_per_sample <- apply(intensity.b, 2, function(x) {
  idx <- order(x, decreasing = TRUE, na.last = TRUE)[1:5]
  tibble(gene      = rownames(intensity.b)[idx],
         intensity = x[idx],
         rank      = 1:5)
}) %>%
  bind_rows(.id = "sample_id")

# --- 3: Assemble DAList -----------------------------------------
cat("\n>> 3 - Assemble DAList \n")

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

# Assemble DAList (all samples — SURV vs CTL)
dal <- DAList(
  data       = intensity,
  metadata   = meta,
  annotation = annotation
)

cat(sprintf(" DAList: %d proteins x %d samples\n", nrow(dal$data),  ncol(dal$data)))

check_collinearity <- function(dal, group_col) {
  df <- data.frame(blood_z = dal$metadata$blood_z,
                   group   = dal$metadata[[group_col]])
  fit <- aov(blood_z ~ group, data = df)
  cat(sprintf("\n blood_z ~ %s\n", group_col))
  print(summary(fit))
}

check_collinearity(dal, "cancer_time")

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

filter_by_bm_ratio <- function(d, bm_thresh, nm = NULL) {
  meta_df  <- d$metadata
  keep_idx <- is.na(meta_df$B_M_ratio) | meta_df$B_M_ratio <= bm_thresh
  
  if (!is.logical(keep_idx) || length(keep_idx) != nrow(meta_df)) {
    stop(sprintf("[%s] Invalid keep_idx in BM filtering", nm))
  }
  
  keep_samples <- meta_df$sample_id[keep_idx]
  n_before     <- ncol(d$data)
  
  d$data     <- d$data[, colnames(d$data) %in% keep_samples, drop = FALSE]
  d$metadata <- d$metadata[d$metadata$sample_id %in% colnames(d$data), ]
  
  n_after <- ncol(d$data)
  cat(sprintf("   %d -> %d samples (%d removed; B:M ratio <= %s)\n",
              n_before, n_after, n_before - n_after, bm_thresh))
  d
}

group_col <- "cancer_time"

filter_log <- tibble(step = character(), n_before = integer(),
                     n_after = integer(), n_removed = integer())
n_start <- nrow(dal$data)
d <- dal

cat(" Step 1: Convert 0s to NA\n")
d <- zero_to_missing(d)

cat(" Step 2: HPA muscle contaminant removal\n")
n_before <- nrow(d$data)
d        <- filter_proteins_by_annotation(d, in_hpa_muscle)
n_after  <- nrow(d$data)
filter_log <- log_step(filter_log, "HPA muscle filter", n_before, n_after)
cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

cat(" Step 3: Keratin filter\n")
n_before <- nrow(d$data)
d        <- filter_proteins_by_annotation(d, !is_keratin)
n_after  <- nrow(d$data)
filter_log <- log_step(filter_log, "Keratin filter", n_before, n_after)
cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

cat(" Step 4: Blood blacklist filter\n")
n_before <- nrow(d$data)
d        <- filter_proteins_by_annotation(d, !in_blood_blacklist)
n_after  <- nrow(d$data)
filter_log <- log_step(filter_log, "Blood blacklist filter", n_before, n_after)
cat(sprintf("   %d -> %d (%d removed)\n", n_before, n_after, n_before - n_after))

cat(sprintf(" Step 5: B:M ratio filter (threshold = %s)\n", config$bm_thresh))
d <- filter_by_bm_ratio(d, config$bm_thresh, "muscle")

# --- 4a: Detection Proportion Assessment -------------------------------
cat("\n>> 4a - Detection proportion assessment\n")

assess_detection_proportion <- function(d, grouping_column,
                                        min_reps = 5,
                                        min_prop = 0.50) {
  
  data_mat <- as.matrix(d$data)
  meta_df  <- d$metadata
  
  meta_df <- meta_df[match(colnames(data_mat), meta_df$sample_id), ]
  
  stopifnot(all(meta_df$sample_id == colnames(data_mat)))
  
  groups <- unique(meta_df[[grouping_column]])
  groups <- groups[!is.na(groups)]
  
  detected <- !is.na(data_mat)
  
  detection_df <- tibble(
    uniprot_id = rownames(data_mat)
  )
  
  for (grp in groups) {
    
    sample_idx <- meta_df[[grouping_column]] == grp
    
    n_group <- sum(sample_idx)
    
    detected_n <- rowSums(detected[, sample_idx, drop = FALSE])
    
    detected_prop <- detected_n / n_group
    
    grp_safe <- make.names(as.character(grp))
    
    detection_df[[paste0(grp_safe, "_n")]] <- detected_n
    detection_df[[paste0(grp_safe, "_total")]] <- n_group
    detection_df[[paste0(grp_safe, "_prop")]] <- detected_prop
    detection_df[[paste0(grp_safe, "_pct")]] <- detected_prop * 100
  }
  
  prop_cols <- grep("_prop$", names(detection_df), value = TRUE)
  
  count_cols <- grep("_n$", names(detection_df), value = TRUE)
  
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
      
      pass_prop_50_each =
        min_group_prop >= min_prop,
      
      pass_group_5_any =
        apply(
          as.matrix(dplyr::select(., all_of(count_cols))),
          1,
          function(x) any(x >= min_reps)
        ),
      
      filter_disagreement =
        pass_prop_50_each != pass_group_5_any
    )
  
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
  
  detection_df <- detection_df %>%
    left_join(
      d$annotation %>%
        dplyr::select(uniprot_id, protein, gene, description),
      by = "uniprot_id"
    ) %>%
    relocate(uniprot_id, protein, gene, description)
  
  detection_df
}


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

write.csv(
  detection_proportion,
  file.path(
    config$report_dir,
    "03_detection_proportion_assessment_muscle.csv"
  ),
  row.names = FALSE
)

detection_filter_sensitive <- detection_proportion %>%
  filter(filter_disagreement) %>%
  arrange(desc(detection_range))

write.csv(
  detection_filter_sensitive,
  file.path(
    config$report_dir,
    "04_filter_sensitive_proteins_muscle.csv"
  ),
  row.names = FALSE
)

cat(sprintf(
  "\n   Exported %d filter-sensitive proteins.\n",
  nrow(detection_filter_sensitive)
))

cat(" Step 6: Missingness filter (50% in each group)\n")
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
  
  sample_medians <- apply(log2(d$data + 1), 2, median, na.rm = TRUE)
  global_median  <- median(sample_medians)
  mad_val        <- mad(sample_medians)
  
  mad_flags <- tibble(
    sample_id     = names(sample_medians),
    sample_median = sample_medians,
    mad_deviation = abs(sample_medians - global_median),
    mad_flag      = abs(sample_medians - global_median) > 3 * mad_val
  )
  
  outlier_diag <- delta_missing %>%
    dplyr::select(sample_id, pct_missing, delta_missing, miss_flag) %>%
    dplyr::left_join(pca_flags, by = "sample_id") %>%
    dplyr::left_join(mad_flags,  by = "sample_id") %>%
    dplyr::mutate(
      n_flags           = miss_flag + pca_flag + mad_flag,
      consensus_outlier = n_flags >= 3
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
outlier_result <- run_outlier_detection(dal, "muscle")
outlier_diag   <- outlier_result$outlier_diag

# Save data before outlier removal
dal_pre_outlier <- dal   # snapshot before removal, for QC comparison

# Run outlier removal
dal <- run_outlier_removal(dal, outlier_diag, "muscle")

# --- 6: Normalize Data ----------------------------
cat("\n>> 6 - Normalize Data\n")

saveRDS(dal, file.path(config$data_dir, "00_DAList_pre_norm_muscle.rds"))

write_norm_report(dal, grouping_column = group_col,
                  output_dir = config$report_dir,
                  filename   = "02_norm_comparison_muscle.pdf", overwrite = TRUE)
write_qc_report(dal, color_column = group_col, label_column = "sample_id",
                output_dir = config$report_dir,
                filename   = "01_qc_report_pre_norm_muscle.pdf", overwrite = TRUE)

dal <- normalize_data(dal, norm_method = config$norm_method)

cat("\n Diagnostic: post-normalization data range (to check if already log-scale)\n")
rng <- range(dal$data, na.rm = TRUE)
cat(sprintf("  range: %.2f to %.2f\n", rng[1], rng[2]))

write_qc_report(dal, color_column = group_col, label_column = "sample_id",
                output_dir = config$report_dir,
                filename   = "03_qc_report_post_norm_muscle.pdf", overwrite = TRUE)

cat(sprintf("  Final: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

pca_postnorm <- compute_pca(dal$data, log_transform = FALSE) # TRUE if diagnostic shows linear scale

# --- 7: Export Normalized Data ----------------------
cat("\n>> 7 - Export Normalized Data\n")

export <- bind_cols(
  dal$annotation %>% dplyr::select(uniprot_id, protein, gene, description, Ensembl),
  dal$data
)

write.csv(export, file.path(config$data_dir, "01_normalized_data_muscle.csv"))
saveRDS(dal,       file.path(config$data_dir, "02_normalized_data_muscle.RDS"))

cat(sprintf("  Exported: %d proteins x %d samples\n", nrow(dal$data), ncol(dal$data)))

# --- 8: Store Data for Reports --------------
report_data <- list(
  config                 = config,
  flagged_hpa            = flagged_hpa,
  flagged_keratin        = flagged_keratin,
  flagged_blood          = flagged_blood,
  marker_pct_df          = marker_pct_df,
  outlier_result         = outlier_result,
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
  pca_postnorm           = pca_postnorm        
)

saveRDS(report_data, file.path(config$data_dir, "00_report_data.rds"))

writeLines(capture.output(sessionInfo()), file.path(config$data_dir, "sessionInfo.txt"))