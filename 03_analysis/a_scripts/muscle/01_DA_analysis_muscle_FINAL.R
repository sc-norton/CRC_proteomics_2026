# Limma DA Analysis

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, limma, writexl, msigdbr, fgsea, rlang)

# Configuration
config <- list(
  report_dir = "03_analysis/b_reports/muscle/limma",
  data_dir   = "03_analysis/c_data/muscle/limma",
  norm_RDS   = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

dir_tree <- list(
  supp = list(
    report = file.path(config$report_dir, "supp")
  ),
  cancer = list(
    report = file.path(config$report_dir, "cancer")
  )
)

invisible(lapply(dir_tree, function(ds_dirs) {
  lapply(ds_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}))

# --- 1: Load Data ----------------------
CRC <- readRDS(config$norm_RDS)

CRC$metadata[c("sample_id", "blood_pct")]

# --- 2: Limma Analysis Function --------

run_limma_analysis <- function(dal,
                               contrast_specs,
                               group_col   = "supp_time",
                               block_col   = "pid",
                               covariates  = NULL,
                               label       = "analysis",
                               lfc_thresh  = 0,
                               pval_thresh = 0.05,
                               dirs,
                               table_cols  = c("uniprot_id", "gene", "protein"),
                               title_col   = "gene") {
  
  cat(sprintf("\n>> Running limma: %s\n", label))
  
  # (a) Design & contrasts
  cat("  Step 1: Design & contrasts\n")
  
  cov_term <- if (!is.null(covariates)) {
    paste0(" + ", paste(covariates, collapse = " + "))
  } else ""
  
  design_formula <- paste0("~ 0 + ", group_col, cov_term,
                           " + (1 | ", block_col, ")")
  cat("    Formula:", design_formula, "\n")
  
  dal <- add_design(dal, design_formula)
  
  # --- repair covariate column names dropped by add_design() ---
  cn <- colnames(dal$design$design_matrix)
  empty_idx <- which(cn == "" | is.na(cn))
  
  if (length(empty_idx) > 0) {
    if (is.null(covariates) || length(empty_idx) != length(covariates)) {
      stop("Mismatch: ", length(empty_idx), " unnamed design column(s) but ",
           length(covariates), " covariate(s) supplied — inspect design_matrix manually.")
    }
    cn[empty_idx] <- covariates
    colnames(dal$design$design_matrix) <- cn
  }
  # ---------------------------------------------------------------
  
  colnames(dal$design$design_matrix) <-
    gsub(paste0("^", group_col), "", colnames(dal$design$design_matrix))
  dal <- add_contrasts(dal, contrasts_vector = contrast_specs)
  
  # (b) Fit model via proteoDA
  cat("  Step 2: Fit model\n")
  dal <- fit_limma_model(dal)
  
  # (c) Extract results
  cat("  Step 3: Extract & write results\n")
  dal <- extract_DA_results(dal,
                            pval_thresh       = pval_thresh,
                            lfc_thresh        = lfc_thresh,
                            adj_method        = "BH",
                            extract_intercept = FALSE)
  
  write_limma_tables(dal, output_dir = dirs$report, overwrite = TRUE)
  file.rename(
    file.path(dirs$report, "results.xlsx"),
    file.path(dirs$report, paste0("results_", label, ".xlsx"))
  )
  write_limma_plots(dal, grouping_column = group_col, output_dir = dirs$report,
                    table_columns = table_cols, title_column = title_col,
                    overwrite = TRUE)
  
  # (d) Read results back & compute transformed pi-values (ϖ = p^|logFC|)
  cat("  Step 4: Compute pi-values\n")
  per_contrast_dir <- file.path(dirs$report, "per_contrast_results")
  results_list <- lapply(setNames(names(dal$results), names(dal$results)), function(cname) {
    read_csv(file.path(per_contrast_dir, paste0(cname, ".csv")),
             show_col_types = FALSE)
  })
  
  results_list <- lapply(results_list, function(df) {
    df$uniprot <- df$uniprot_id %||% rownames(df)
    df$Pi.Val  <- (df$P.Value)^abs(df$logFC)
    df[, c("uniprot", setdiff(names(df), "uniprot"))]
  })
  
  # (e) DA classification
  cat("  Step 5: DA classification\n")
  results_list <- classify_DA(results_list,
                              lfc_thresh  = lfc_thresh,
                              pval_thresh = pval_thresh)
  da_log <- summarise_DA(results_list)
  
  write.csv(da_log,
            file.path(dirs$report, paste0("DA_criteria_", label, ".csv")),
            row.names = FALSE)
  
  # (f) Save RDS
  cat("  Step 6: Save RDS\n")
  saveRDS(list(results_list = results_list, da_log = da_log),
          file.path(dirs$report, paste0("limma_out_", label, ".RDS")))
  
  cat(sprintf("  Done: %s\n", label))
  list(results_list = results_list,
       da_log       = da_log)
}


# --- 3: DA Helper Functions ------------

classify_DA <- function(results_list, lfc_thresh = 0.6, pval_thresh = 0.05) {
  lapply(results_list, function(df) {
    classify <- function(lfc, pval, lfc_t, pval_t = pval_thresh)
      case_when(lfc >  lfc_t & pval < pval_t ~ "UP",
                lfc < -lfc_t & pval < pval_t ~ "DOWN",
                TRUE                         ~ "NDE")
    
    df$DE.adj   <- classify(df$logFC, df$adj.P.Val, lfc_thresh)
    df$DE.adj.1 <- classify(df$logFC, df$adj.P.Val, lfc_thresh, pval_t = 0.1)
    df$DE.nom   <- classify(df$logFC, df$P.Value,   lfc_thresh)
    df$DE.pi    <- classify(df$logFC, df$Pi.Val,    lfc_t = 0)
    df
  })
}

summarise_DA <- function(results_list) {
  cols   <- c(adj_05 = "DE.adj", adj_1 = "DE.adj.1", nom = "DE.nom", pi = "DE.pi")
  labels <- c(adj_05 = "FDR < 0.05 (BH)", adj_1 = "FDR < 0.10 (BH)", 
              nom = "Nominal P-val", pi = "Pi-value")
  do.call(rbind, lapply(names(cols), function(key) {
    do.call(rbind, lapply(names(results_list), function(cname) {
      df <- results_list[[cname]]
      data.frame(comparison = cname, method = labels[[key]],
                 n_UP   = sum(df[[cols[[key]]]] == "UP"),
                 n_DOWN = sum(df[[cols[[key]]]] == "DOWN"),
                 stringsAsFactors = FALSE)
    }))
  }))
}


# --- 4: Build Contrasts ----------------

contrast_specs_supp <- c(
  "Baseline_supp       = CRE_T1 - PLA_T1",
  "Post_training_supp  = CRE_T2 - PLA_T2",
  "Training_CRE        = CRE_T2 - CRE_T1",
  "Training_PLA        = PLA_T2 - PLA_T1",
  "Interaction_supp    = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

# contrast_specs_sex <- c(
#   "Baseline_sex        = Female_T1 - Male_T1",
#   "Post-training_sex   = Female_T2 - Male_T2",
#   "Training_F          = Female_T2 - Female_T1",
#   "Training_M          = Male_T2 - Male_T1",
#   "Interaction_sex     = (Female_T2 - Female_T1) - (Male_T2 - Male_T1)"
# )

contrast_specs_cancer <- c(
  "Cancer   = CRC_T1 - CTL_T1",
  "Residual = CRC_T2 - CTL_T1",
  "Training = CRC_T2 - CRC_T1"
)

# --- 5: Run Analysis -------------------

out_supp_cv <- run_limma_analysis(CRC, contrast_specs_supp,
                              group_col = "supp_time",
                              covariates = "blood_pct",
                              label     = "supp_cv",
                              dirs      = dir_tree$supp)

out_cancer_cv <- run_limma_analysis(CRC, contrast_specs_cancer,
                               group_col  = "cancer_time",
                               covariates = "blood_pct",
                               label      = "cancer_cv",
                               dirs       = dir_tree$cancer)

out_supp <- run_limma_analysis(CRC, contrast_specs_supp,
                                  group_col = "supp_time",
                                  label     = "supp",
                                  dirs      = dir_tree$supp)

out_cancer <- run_limma_analysis(CRC, contrast_specs_cancer,
                                    group_col  = "cancer_time",
                                    label      = "cancer",
                                    dirs       = dir_tree$cancer)

# out_cre_sex <- run_limma_analysis(CRE, contrast_specs_cre_sex,
#                                   group_col = "sex_time",
#                                   covariates = "blood_pct",
#                                   label     = "SEX",
#                                   dirs      = dir_tree$FvM)

# --- 6: Combine & Export Results -------

# Matrix sheet: identifiers + every participant sample column, built once
# directly from CRC (not from a results_list df) so it isn't tied to
# whichever group's contrast happened to produce it.
sample_cols <- colnames(CRC$data)

normalized_matrix <- bind_cols(
  CRC$annotation %>% dplyr::select(uniprot_id, protein, gene, description),
  CRC$data
) %>%
  dplyr::rename(uniprot = uniprot_id)   # match the "uniprot" naming your stats sheets use

# Stats-only columns for every contrast sheet
stats_cols <- c("uniprot", "protein", "gene", "Ensembl", "description",
                "logFC", "CI.L", "CI.R", "average_intensity", "t",
                "B", "P.Value", "adj.P.Val", "sig.PVal", "sig.FDR",
                "Pi.Val", "DE.adj", "DE.adj.1", "DE.nom", "DE.pi")

results_combined_cv <- c(out_supp_cv$results_list, out_cancer_cv$results_list)
results_combined    <- c(out_supp$results_list,    out_cancer$results_list)

results_combined_cv <- lapply(results_combined_cv, function(df) {
  as.data.frame(df)[, stats_cols, drop = FALSE]
})
results_combined <- lapply(results_combined, function(df) {
  as.data.frame(df)[, stats_cols, drop = FALSE]
})

workbook_sheets    <- c(list(normalized_matrix = normalized_matrix), results_combined)
workbook_sheets_cv <- c(list(normalized_matrix = normalized_matrix), results_combined_cv)

da_combined    <- rbind(out_supp$da_log,    out_cancer$da_log)
da_combined_cv <- rbind(out_supp_cv$da_log, out_cancer_cv$da_log)

write_xlsx(workbook_sheets,    file.path(config$data_dir, "limma_results_all_muscle.xlsx"))
write_xlsx(workbook_sheets_cv, file.path(config$data_dir, "limma_results_all_muscle_cv.xlsx"))
write.csv(da_combined,    file.path(config$report_dir, "DA_methods_combined_muscle.csv"), row.names = FALSE)
write.csv(da_combined_cv, file.path(config$report_dir, "DA_methods_combined_muscle_cv.csv"), row.names = FALSE)

# STOP HERE #
# --- 7: Export UP/DOWN Protein Lists ---

cat("\n>> Exporting UP/DOWN protein lists per comparison\n")

export_dir <- file.path(config$report_dir, "all_data", "limma", "DA_protein_lists")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

export_DA_lists_xlsx <- function(results_combined, de_col, suffix, export_dir) {
  sheets <- lapply(names(results_combined), function(cname) {
    df  <- results_combined[[cname]]
    sig <- df[df[[de_col]] %in% c("UP", "DOWN"), ]
    sig$direction <- sig[[de_col]]
    sig <- sig[order(sig$direction, sig$logFC), ]
    cat(sprintf("  %s -> UP: %d, DOWN: %d (%s)\n",
                cname,
                sum(sig$direction == "UP"),
                sum(sig$direction == "DOWN"),
                suffix))
    sig
  })
  names(sheets) <- names(results_combined)
  write_xlsx(sheets,
             file.path(export_dir, paste0("DA_protein_lists_", suffix, ".xlsx")))
}

export_DA_lists_xlsx(results_combined, "DE.pi",    "pi",      export_dir)
export_DA_lists_xlsx(results_combined, "DE.nom",   "p.raw",   export_dir)
export_DA_lists_xlsx(results_combined, "DE.adj",   "FDR_05",  export_dir)
export_DA_lists_xlsx(results_combined, "DE.adj.1", "FDR_10",  export_dir)

cat("\n  Protein lists saved to:", export_dir, "\n")
