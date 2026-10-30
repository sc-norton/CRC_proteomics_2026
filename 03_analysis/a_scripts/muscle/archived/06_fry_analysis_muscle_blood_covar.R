## =============================================================================
## 06_fry_reversal_signature_muscle_SURV.R
##
## Empirical signature-set fry() reversal test, adapted from colleague's
## panel_C_fry.R structure to your SURV/CTL design and proteoDA pipeline.
##
## Departs from the curated-database version (Hallmark/Reactome/GO:BP tested
## independently per contrast) in favor of:
##   1. Define Up/Down signature sets from proteins significant in the SOURCE
##      contrast (Baseline_SURVvCTL) -- i.e. the cancer dysregulation signature
##   2. Run fry() ONCE, on the TEST contrast (Training_SURV), against those
##      sets -- a single hypothesis test with a pre-specified expected
##      direction, rather than two independent fry() calls compared post-hoc
##   3. Report a circularity diagnostic: cor(t_source, t_test) across all
##      proteins, quantifying how much the source/test contrasts already
##      share variance before fry() is even run
##   4. Identify "driving proteins" within each signature set whose individual
##      t-statistic in the test contrast moves in the expected direction, and
##      run leading-edge ORA on that subset
##
## NOTE ON ASSUMPTIONS -- verify against your actual objects before trusting
## output. Anything marked >>> VERIFY <<< is inferred from pipeline convention.
## =============================================================================
getwd()
setwd("/Users/samnorton/Documents/CRC_pro")
getwd()

pacman::p_load(
  limma, proteoDA, dplyr, tidyr, purrr, tibble, readr, readxl,
  writexl, stringr, clusterProfiler
)

set.seed(26)

## -----------------------------------------------------------------------
## 0. Config
## -----------------------------------------------------------------------

config <- list(
  paths = list(
    imp_dal_rds    = "02_imputation/c_data/muscle/muscle_imputed_dalist.RDS",  
    limma_xlsx = "03_analysis/c_data/muscle/blood_covar/limma_results_all_muscle.xlsx",  # sheets 11-13 = SURV contrasts
    out_dir    = "03_analysis/c_data/muscle/blood_covar/fry_reversal",
    report_dir = "03_analysis/b_reports/muscle/blood_covar/fry_reversal"
  ),
  universe_n = 1958,
  sig_col    = "Pi.Val",   # primary significance criterion per pipeline convention
  sig_thresh = 0.05,
  min_set_size = 5,
  source_contrast = "Baseline_SURVvCTL",   # defines the cancer dysregulation signature
  test_contrast   = "Training_SURV",       # where reversal is tested (POST vs PRE)
  contrast_labels = c(
    Baseline_SURVvCTL = "PRE vs CTL",
    Training_SURVvCTL = "POST vs CTL",
    Training_SURV     = "POST vs PRE"
  ),
  block_col = "pid"   # >>> VERIFY <<< sample-level pairing variable
)

dir.create(config$paths$out_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(config$paths$report_dir, recursive = TRUE, showWarnings = FALSE)

build_output_path <- function(dir, filename) file.path(dir, filename) 

## -----------------------------------------------------------------------
## 1. Load DAList, expression matrix, and per-contrast result tables
## -----------------------------------------------------------------------

dal <- readRDS(config$paths$imp_dal_rds)

expr_mat <- as.matrix(dal$data)   # >>> VERIFY <<< rows = uniprot_id, cols = samples


n_na <- sum(!is.finite(expr_mat))
if (n_na > 0) {
  stop(sprintf(
    "expr_mat has %d non-finite values. fry() requires a complete matrix -- point config$paths$imputed_rds at the missForest output.",
    n_na
  ))
}

## Per-contrast result tables carrying t-statistic, logFC, and Pi.Val, keyed
## on the same uniprot_id used as expr_mat rownames. Sheets 11-13 = SURV
## contrasts per your file convention.
res_source <- read_excel(config$paths$limma_xlsx, sheet = "Baseline_SURVvCTL")
res_test   <- read_excel(config$paths$limma_xlsx, sheet = "Training_SURV")      

stopifnot(all(c("uniprot_id", "gene", "logFC", "t", "Pi.Val") %in% names(res_source))) 

## -----------------------------------------------------------------------
## 2. Design, contrasts, duplicateCorrelation (double pass, per convention)
## -----------------------------------------------------------------------

meta <- dal$metadata  
stopifnot(identical(colnames(expr_mat), meta$sample_id))  # >>> VERIFY <<< actual sample-id column name

meta$cancer_time <- factor(meta$cancer_time,            # >>> VERIFY <<< matches your WGCNA three-level factor convention
                           levels = c("CTL_T1", "SURV_T1", "SURV_T2"))
design <- model.matrix(~ 0 + cancer_time, data = meta)
colnames(design) <- gsub("^cancer_time", "", colnames(design))

block_id <- meta[[config$block_col]]

## Double duplicateCorrelation pass for stability, consistent with the
## primary DA model -- do not collapse this to a single pass.
corfit <- duplicateCorrelation(expr_mat, design, block = block_id)
cor_consensus <- corfit$consensus.correlation
message(sprintf("Within-subject consensus correlation (2-pass): %.4f", cor_consensus))

cm <- makeContrasts(
  Baseline_SURVvCTL = SURV_T1 - CTL_T1,
  Training_SURVvCTL = SURV_T2 - CTL_T1,
  Training_SURV     = SURV_T2 - SURV_T1,
  levels = design
)

## -----------------------------------------------------------------------
## 3. Circularity diagnostic
##
## cor(t_source, t_test) across ALL proteins (not just the signature set)
## quantifies how much the two contrasts already share variance before any
## selection happens. Report this alongside the fry result -- it's the
## quantitative version of the "hypothesis-driven, not unbiased" caveat.
## -----------------------------------------------------------------------

t_join <- res_source %>%
  dplyr::select(uniprot_id, t_source = t) %>%
  inner_join(res_test %>% dplyr::select(uniprot_id, t_test = t), by = "uniprot_id") %>%
  filter(uniprot_id %in% rownames(expr_mat))

circ_r <- stats::cor(t_join$t_source, t_join$t_test, use = "complete.obs")
message(sprintf("Circularity: r(t_%s, t_%s) = %.3f",
                config$source_contrast, config$test_contrast, circ_r))

## -----------------------------------------------------------------------
## 4. Define signature sets from the SOURCE contrast
## -----------------------------------------------------------------------

sig_source <- res_source %>%
  filter(.data[[config$sig_col]] < config$sig_thresh,
         uniprot_id %in% rownames(expr_mat))

sets <- list(
  up_ids   = sig_source$uniprot_id[sig_source$logFC > 0],
  down_ids = sig_source$uniprot_id[sig_source$logFC < 0]
)
sets$up   <- match(sets$up_ids,   rownames(expr_mat))
sets$down <- match(sets$down_ids, rownames(expr_mat))

message(sprintf(
  "%s signature (%s < %.2f): up = %d, down = %d",
  config$source_contrast, config$sig_col, config$sig_thresh,
  length(sets$up), length(sets$down)
))

## -----------------------------------------------------------------------
## 5. Run fry() once on the TEST contrast, against pre-specified direction
##
## Reversal logic: a protein UP at baseline (cancer > control) is expected
## to move DOWN during training if exercise reverses dysregulation, and vice
## versa. This is a single hypothesis test per set, not a post-hoc comparison
## of two independently-estimated Direction calls.
## -----------------------------------------------------------------------

run_fry_set <- function(idx, set_name, expected_direction) {
  if (length(idx) < config$min_set_size) {
    return(tibble(
      set = set_name, n = length(idx), direction = NA_character_,
      PValue = NA_real_, PValue.Mixed = NA_real_,
      expected = expected_direction, consistent = NA
    ))
  }
  res <- fry(
    expr_mat, index = idx, design = design,
    contrast = cm[, config$test_contrast],
    block = block_id, correlation = cor_consensus
  )
  tibble(
    set = set_name, n = length(idx), direction = res$Direction[1],
    PValue = res$PValue[1], PValue.Mixed = res$PValue.Mixed[1],
    expected = expected_direction, consistent = res$Direction[1] == expected_direction
  )
}

## A protein set UP at baseline: full reversal in training looks like Down;
## a set DOWN at baseline: full reversal looks like Up.
fry_up <- run_fry_set(sets$up,   "signature_up",   expected_direction = "Down") %>%
  mutate(sig_thresh_used = config$sig_thresh, cor_within = cor_consensus, circularity_r = circ_r)
fry_down <- run_fry_set(sets$down, "signature_down", expected_direction = "Up") %>%
  mutate(sig_thresh_used = config$sig_thresh, cor_within = cor_consensus, circularity_r = circ_r)

fry_reversal <- bind_rows(fry_up, fry_down) %>%
  mutate(
    source_contrast = config$source_contrast,
    test_contrast   = config$test_contrast,
    interpretation = case_when(
      is.na(consistent)               ~ "Set too small to test",
      PValue < 0.05 & consistent       ~ "Reversal supported",
      PValue < 0.05 & !consistent      ~ "Persistent / exacerbated (unexpected direction)",
      TRUE                             ~ "Not significant"
    )
  ) %>%
  relocate(source_contrast, test_contrast, set, n, expected, direction, consistent, interpretation)

write_csv(fry_reversal, build_output_path(config$paths$report_dir, "fry_reversal_signature_muscle_SURV.csv"))
saveRDS(fry_reversal, build_output_path(config$paths$out_dir, "fry_reversal_signature_muscle_SURV.RDS"))

print(fry_reversal, Inf)
## -----------------------------------------------------------------------
## 6. Driving proteins -- individual members whose test-contrast t-statistic
##    moves in the expected (reversal) direction
## -----------------------------------------------------------------------

driving_up <- res_test %>%
  filter(uniprot_id %in% sets$up_ids, t < 0) %>%      # expected Down in test contrast
  transmute(uniprot_id, gene, set = "signature_up",
            t_source = res_source$t[match(uniprot_id, res_source$uniprot_id)],
            t_test = t, logFC_test = logFC, Pi.Val_test = Pi.Val)

driving_down <- res_test %>%
  filter(uniprot_id %in% sets$down_ids, t > 0) %>%    # expected Up in test contrast
  transmute(uniprot_id, gene, set = "signature_down",
            t_source = res_source$t[match(uniprot_id, res_source$uniprot_id)],
            t_test = t, logFC_test = logFC, Pi.Val_test = Pi.Val)

driving_df <- bind_rows(driving_up, driving_down)
write_csv(driving_df, build_output_path(config$paths$report_dir, "fry_driving_proteins_muscle_SURV.csv"))

message(sprintf(
  "Driving proteins: %d / %d (up set), %d / %d (down set)",
  nrow(driving_up),   length(sets$up_ids),
  nrow(driving_down), length(sets$down_ids)
))

## -----------------------------------------------------------------------
## 7. Leading-edge ORA on driving proteins
##
## Manual analogue to GSEA leading-edge analysis: which pathways do the
## proteins actually moving in the expected direction belong to. Universe
## fixed at 1958 proteins per pipeline convention.
## -----------------------------------------------------------------------

## >>> VERIFY <<< swap in your existing GO-slim / Reactome TERM2GENE table
## rather than rebuilding one here -- this assumes you already have a
## `pathway_term2gene` object saved somewhere in the ORA portion of the
## pipeline (mirrors `build_pathway_collection()` in the colleague's script).
if (exists("pathway_term2gene")) {
  
  all_genes <- res_test$gene[res_test$uniprot_id %in% rownames(expr_mat)]
  
  run_ora_leading <- function(gene_vec, label) {
    if (length(gene_vec) < 5) return(tibble())
    tryCatch({
      enricher(
        gene          = gene_vec,
        universe      = all_genes,
        TERM2GENE     = pathway_term2gene,
        pAdjustMethod = "BH",
        pvalueCutoff  = 1,
        qvalueCutoff  = 1
      )@result %>%
        as_tibble() %>%
        filter(p.adjust < 0.1) %>%
        arrange(p.adjust) %>%
        mutate(set = label) %>%
        slice_head(n = 10)
    }, error = function(e) tibble())
  }
  
  ora_up   <- run_ora_leading(driving_up$gene,   "signature_up_driving")
  ora_down <- run_ora_leading(driving_down$gene, "signature_down_driving")
  
  ora_leading <- bind_rows(ora_up, ora_down)
  write_csv(ora_leading, build_output_path(config$paths$report_dir, "fry_leading_edge_ora_muscle_SURV.csv"))
} else {
  message("`pathway_term2gene` not found in environment -- skipping leading-edge ORA. Source your ORA setup script first, or point this at your existing TERM2GENE object.")
}

## -----------------------------------------------------------------------
## 8. Summary
## -----------------------------------------------------------------------

cat("\n=== fry reversal summary ===\n")
print(fry_reversal %>% dplyr::select(set, n, expected, direction, PValue, consistent, interpretation))
cat(sprintf("\nCircularity r(t_%s, t_%s) = %.3f\n", config$source_contrast, config$test_contrast, circ_r))
cat("\nOutputs:\n")
cat(sprintf("  %s\n", build_output_path(config$paths$report_dir, "fry_reversal_signature_muscle_SURV.csv")))
cat(sprintf("  %s\n", build_output_path(config$paths$report_dir, "fry_driving_proteins_muscle_SURV.csv")))
cat(sprintf("  %s\n", build_output_path(config$paths$report_dir, "fry_leading_edge_ora_muscle_SURV.csv")))
