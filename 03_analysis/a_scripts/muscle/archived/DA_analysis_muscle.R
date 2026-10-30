# Limma DA Analysis

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, limma, writexl, msigdbr, fgsea, rlang)

# Configuration
config <- list(
  report_dir = "03_analysis/b_reports/muscle",
  data_dir   = "03_analysis/c_data/muscle",
  SURV_RDS   = "01_normalization/c_data/muscle/SURV/02_DAList_normalized_SURV_muscle.RDS",
  CRE_RDS    = "01_normalization/c_data/muscle/CRE/02_DAList_normalized_CRE_muscle.RDS"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

dir_tree <- list(
  CREvPLA  = list(root  = file.path(config$report_dir, "CREvPLA"),
                  limma = file.path(config$report_dir, "CREvPLA", "limma")),
  SURVvCTL = list(root  = file.path(config$report_dir, "SURVvCTL"),
                  limma = file.path(config$report_dir, "SURVvCTL", "limma")),
  FvM      = list(root  = file.path(config$report_dir, "FvM"),
                  limma = file.path(config$report_dir, "FvM", "limma")),
  all_data = list(root  = file.path(config$report_dir, "all_data"),
                  limma = file.path(config$report_dir, "all_data", "limma"))
)

invisible(lapply(dir_tree, function(ds_dirs) {
  lapply(ds_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}))


# --- 1: Load Data ----------------------
SURV <- readRDS(config$SURV_RDS)
CRE  <- readRDS(config$CRE_RDS)


# --- 2: Limma Analysis Function --------

run_limma_analysis <- function(dal,
                               contrast_specs,
                               group_col   = "supp_time",
                               block_col   = "pid",
                               label       = "analysis",
                               lfc_thresh  = 0.6,
                               pval_thresh = 0.05,
                               dirs,
                               table_cols  = c("uniprot_id", "gene", "protein"),
                               title_col   = "gene") {
  
  cat(sprintf("\n>> Running limma: %s\n", label))
  
  # (a) Design & contrasts
  cat("  Step 1: Design & contrasts\n")
  dal <- add_design(dal, paste0("~ 0 + ", group_col, " + (1 | ", block_col, ")"))
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
  
  write_limma_tables(dal, output_dir = dirs$limma, overwrite = TRUE)
  file.rename(
    file.path(dirs$limma, "results.xlsx"),
    file.path(dirs$limma, paste0("results_", label, ".xlsx"))
  )
  write_limma_plots(dal, grouping_column = group_col, output_dir = dirs$limma,
                    table_columns = table_cols, title_column = title_col,
                    overwrite = TRUE)
  
  # (d) Read results back & compute transformed pi-values (ϖ = p^|logFC|)
  cat("  Step 4: Compute pi-values\n")
  per_contrast_dir <- file.path(dirs$limma, "per_contrast_results")
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
            file.path(dirs$limma, paste0("DA_criteria_", label, ".csv")),
            row.names = FALSE)
  
  # (f) Save RDS
  cat("  Step 6: Save RDS\n")
  saveRDS(list(results_list = results_list, da_log = da_log),
          file.path(dirs$limma, paste0("limma_out_", label, ".RDS")))
  
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
  labels <- c(adj_05 = "FDR < 0.05 (BH)", adj_1 = "FDR < 0.1 (BH)", 
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

contrast_specs_cre <- c(
  "Baseline_CREvPLA    = CRE_T1 - PLA_T1",
  "Training_CREvPLA    = CRE_T2 - PLA_T2",
  "Training_CRE        = CRE_T2 - CRE_T1",
  "Training_PLA        = PLA_T2 - PLA_T1",
  "Interaction_supp    = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

contrast_specs_cre_sex <- c(
  "Baseline_FvM        = Female_T1 - Male_T1",
  "Training_FvM        = Female_T2 - Male_T2",
  "Training_F          = Female_T2 - Female_T1",
  "Training_M          = Male_T2 - Male_T1",
  "Interaction_sex     = (Female_T2 - Female_T1) - (Male_T2 - Male_T1)"
)

contrast_specs_surv <- c(
  "Baseline_SURVvCTL   = SURV_T1 - CTL_T1",
  "Training_SURVvCTL   = SURV_T2 - CTL_T1",
  "Training_SURV       = SURV_T2 - SURV_T1"
)


# --- 5: Run Analysis -------------------

out_cre <- run_limma_analysis(CRE, contrast_specs_cre,
                              group_col = "supp_time",
                              label     = "CRE",
                              dirs      = dir_tree$CREvPLA)

out_cre_sex <- run_limma_analysis(CRE, contrast_specs_cre_sex,
                                  group_col = "sex_time",
                                  label     = "SEX",
                                  dirs      = dir_tree$FvM)

out_surv <- run_limma_analysis(SURV, contrast_specs_surv,
                               group_col = "cancer_time",
                               label     = "SURV",
                               dirs      = dir_tree$SURVvCTL)


# --- 6: Combine & Export Results -------

results_combined <- c(out_cre$results_list,
                      out_cre_sex$results_list,
                      out_surv$results_list)

results_combined <- lapply(results_combined, function(x) {
  df <- as.data.frame(x)
  df[, !names(df) %in% c("gene_symbol"), drop = FALSE]
})

da_combined <- rbind(out_cre$da_log,
                     out_cre_sex$da_log,
                     out_surv$da_log)

write_xlsx(results_combined,
           file.path(config$report_dir, "all_data", "limma", "limma_results_all_muscle.xlsx"))
write.csv(da_combined,
          file.path(config$report_dir, "all_data", "limma", "DA_methods_combined.csv"),
          row.names = FALSE)
write_xlsx(results_combined,
           file.path(config$data_dir, "limma_results_all_muscle.xlsx"))
write.csv(da_combined,
          file.path(config$data_dir, "DA_methods_combined.csv"),
          row.names = FALSE)


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
