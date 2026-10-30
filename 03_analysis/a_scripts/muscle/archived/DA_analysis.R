# Limma DA Analysis

# --- 0: Setup --------------------------
# Set working directory
setwd(rprojroot::find_rstudio_root_file())

# Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, limma, writexl, msigdbr, fgsea)

# Configuration
config <- list(
  # Directories
  report_dir = "03_analysis/b_reports/muscle",
  data_dir = "03_analysis/c_data/muscle",
  # Data files
  SURV_RDS = "01_normalization/c_data/muscle/SURV/stringent/BM15/02_DAList_normalized_SURV_stringent_BM15_muscle.RDS",
  CRE_RDS = "01_normalization/c_data/muscle/CRE/stringent/BM15/02_DAList_normalized_CRE_stringent_BM15_muscle.RDS"
)

dir.create(config$report_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

dir_tree <- list(
  CREvPLA = list(
    root         = file.path(config$report_dir, "CREvPLA"),
    limma        = file.path(config$report_dir, "CREvPLA", "limma")
  ),
  SURVvCTL = list(
    root         = file.path(config$report_dir, "SURVvCTL"),
    limma        = file.path(config$report_dir, "SURVvCTL", "limma")
  ),
  FvM = list(
    root         = file.path(config$report_dir, "FvM"),
    limma        = file.path(config$report_dir, "FvM", "limma")
  ),
  all_data = list(
    root         = file.path(config$report_dir, "all_data"),
    limma        = file.path(config$report_dir, "all_data", "limma")
  ))

# Create all directories at once
invisible(lapply(dir_tree, function(ds_dirs) {
  lapply(ds_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}))

# --- 1 : Load Data ----------------------
SURV <- readRDS(file.path(config$SURV_RDS)) 
  
CRE <- readRDS(file.path(config$CRE_RDS))

# --- 2:  Limma analysis function -----------

run_limma_analysis <- function(dal,
                               contrast_specs,
                               group_col    = "supp_time",
                               block_col    = "pid",
                               label        = "analysis",
                               lfc_thresh   = 0,
                               pval_thresh  = 0.05,
                               dirs,
                               table_cols   = c("uniprot_id", "gene", "protein"),
                               title_col    = "gene") {
  
  cat(sprintf("\n>> Running limma: %s\n", label))
  
  # (a) Design & contrasts
  cat("  Step 1: Design & contrasts\n")
  dal <- add_design(dal, paste0("~ 0 + ", group_col, " + (1 | ", block_col, ")"))
  colnames(dal$design$design_matrix) <-
    gsub(paste0("^", group_col), "", colnames(dal$design$design_matrix))
  dal <- add_contrasts(dal, contrasts_vector = contrast_specs)
  
  design_mat   <- dal$design$design_matrix
  contrast_mat <- dal$design$contrast_matrix
  block_var    <- dal$metadata[[block_col]]
  mat          <- as.matrix(dal$data)
  
  # (b) Array weights
  cat("  Step 2: Array weights\n")
  aw        <- arrayWeights(mat, design_mat)
  names(aw) <- colnames(mat)
  cat(sprintf("    Weights: %.2f - %.2f (median %.2f)\n",
              min(aw), max(aw), median(aw)))
  write_csv(tibble(sample_id = names(aw), array_weight = round(aw, 4)),
            file.path(dirs$limma, paste0("array_weights_", label, ".csv")))
  
  # (c) Fit model
  cat("  Step 3: Fit model\n")
  dupcor <- duplicateCorrelation(mat, design_mat, block = block_var, weights = aw)
  cat(sprintf("    Within-subject correlation: %.3f\n", dupcor$consensus.correlation))
  
  fit <- lmFit(mat, design_mat,
               block       = block_var,
               correlation = dupcor$consensus.correlation,
               weights     = aw) |>
    contrasts.fit(contrast_mat) |>
    eBayes(robust = TRUE, trend = TRUE)
  
  dal$eBayes_fit             <- fit
  dal$eBayes_fit$correlation <- dupcor$consensus.correlation
  
  # (d) Extract results
  cat("  Step 4: Extract & write results\n")
  dal <- extract_DA_results(dal, pval_thresh = pval_thresh,
                            lfc_thresh = 0, adj_method = "BH")
  write_limma_tables(dal, output_dir = dirs$limma, overwrite = TRUE)
  file.rename(
    file.path(dirs$limma, "results.xlsx"),
    file.path(dirs$limma, paste0("results_", label, ".xlsx"))
  )
  write_limma_plots(dal, grouping_column = group_col, output_dir = dirs$limma,
                    table_columns = table_cols, title_column = title_col,
                    overwrite = TRUE)
  
  # (e) Read results back & compute pi-values
  per_contrast_dir <- file.path(dirs$limma, "per_contrast_results")
  results_list <- lapply(setNames(names(dal$results), names(dal$results)), function(cname) {
    read_csv(file.path(per_contrast_dir, paste0(cname, ".csv")),
             show_col_types = FALSE)
  })
  
  results_list <- lapply(results_list, function(df) {
    df$uniprot    <- df$uniprot_id %||% rownames(df)
    df$Pi.Val     <- (df$P.Value)^abs(df$logFC)
    df[, c("uniprot", setdiff(names(df), "uniprot"))]
  })
  
  # (f) DA classification
  results_list <- classify_DA(results_list, lfc_thresh = lfc_thresh,
                              pval_thresh = pval_thresh)
  da_log       <- summarise_DA(results_list, lfc_thresh = lfc_thresh,
                               pval_thresh = pval_thresh)
  
  write.csv(da_log,
            file.path(dirs$limma, paste0("DA_criteria_", label, ".csv")),
            row.names = FALSE)
  
  # (g) Save results to RDS for later use
  cat("  Step 5: Saving results to RDS\n")
  saveRDS(list(results_list = results_list, da_log = da_log),
          file.path(dirs$limma, paste0("limma_out_", label, ".RDS")))
  
  cat(sprintf("  Done: %s\n", label))
  list(results_list = results_list,
       da_log       = da_log)
}


# --- 3:  DA helper functions -----------

classify_DA <- function(results_list, lfc_thresh = 0.6, pval_thresh = 0.05) {
  lapply(results_list, function(df) {
    classify <- function(lfc, pval)
      case_when(lfc >  lfc_thresh & pval < pval_thresh ~ "UP",
                lfc < -lfc_thresh & pval < pval_thresh ~ "DOWN",
                TRUE                                   ~ "NDE")
    classify_pi <- function(lfc, pval)
      case_when(lfc >  0 & pval < pval_thresh ~ "UP",
                lfc <  0 & pval < pval_thresh ~ "DOWN",
                TRUE                          ~ "NDE")
    
    df$DE.adj <- classify(df$logFC, df$adj.P.Val)
    df$DE.nom <- classify(df$logFC, df$P.Value)
    df$DE.pi  <- classify_pi(df$logFC, df$Pi.Val)
    df
  })
}


summarise_DA <- function(results_list, lfc_thresh = 0.6, pval_thresh = 0.05) {
  cols    <- c(adj = "DE.adj", nom = "DE.nom", pi = "DE.pi")
  labels  <- c(adj = "BH-adj P-val", nom = "Nominal P-val", pi = "Pi-value")
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

# --- 4:  Build Contrasts --------
contrast_specs_cre <- c(
  "Baseline_CREvPLA = CRE_T1 - PLA_T1",
  "Training_CREvPLA = CRE_T2 - PLA_T2",
  "Training_CRE     = CRE_T2 - CRE_T1",
  "Training_PLA     = PLA_T2 - PLA_T1",
  "Interaction_supp      = (CRE_T2 - CRE_T1) - (PLA_T2 - PLA_T1)"
)

contrast_specs_cre_sex <- c(
  "Baseline_FvM = Female_T1 - Male_T1",
  "Training_FvM = Female_T2 - Male_T2",
  "Training_F   = Female_T2 - Female_T1",
  "Training_M   = Male_T2 - Male_T1",
  "Interaction_sex  = (Female_T2 - Female_T1) - (Male_T2 - Male_T1)"
)

contrast_specs_surv <- c(
  "Baseline_SURVvCTL = SURV_T1 - CTL_T1",
  "Training_SURVvCTL = SURV_T2 - CTL_T1",
  "Training_SURV     = SURV_T2 - SURV_T1"
)

# --- 5: Run Analysis ---------------
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

# Combine results
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

# Export results
write_xlsx(results_combined,
           file.path(config$report_dir, "all_data", "limma", "limma_results_all_muscle.xlsx"))
write.csv(da_combined,
          file.path(config$report_dir, "all_data", "limma", "DA_methods_combined.csv"), row.names = FALSE)
write_xlsx(results_combined,
           file.path(config$data_dir, "limma_results_all_muscle.xlsx"))
write.csv(da_combined,
          file.path(config$data_dir, "DA_methods_combined.csv"), row.names = FALSE)

# Export UP and DOWN lists per comparison
cat("\n>> Exporting UP/DOWN protein lists per comparison\n")

export_dir <- file.path(config$report_dir, "all_data", "limma", "DA_protein_lists")
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

invisible(lapply(names(results_combined), function(cname) {
  df <- results_combined[[cname]]
  
  up_df   <- df[df$DE.pi == "UP",   ]
  down_df <- df[df$DE.pi == "DOWN", ]
  
  write.csv(up_df,
            file.path(export_dir, paste0(cname, "_pi_UP.csv")),
            row.names = FALSE)
  write.csv(down_df,
            file.path(export_dir, paste0(cname, "_pi_DOWN.csv")),
            row.names = FALSE)
  
  cat(sprintf("  %s -> UP: %d, DOWN: %d (Pi-value) \n", 
              cname, nrow(up_df), nrow(down_df)))
}))

invisible(lapply(names(results_combined), function(cname) {
  df <- results_combined[[cname]]
  
  up_df   <- df[df$DE.nom == "UP",   ]
  down_df <- df[df$DE.nom == "DOWN", ]
  
  write.csv(up_df,
            file.path(export_dir, paste0(cname, "_p.raw_UP.csv")),
            row.names = FALSE)
  write.csv(down_df,
            file.path(export_dir, paste0(cname, "_p.raw_DOWN.csv")),
            row.names = FALSE)
  
  cat(sprintf("  %s -> UP: %d, DOWN: %d (Raw p-value) \n", 
              cname, nrow(up_df), nrow(down_df)))
}))

invisible(lapply(names(results_combined), function(cname) {
  df <- results_combined[[cname]]
  
  up_df   <- df[df$DE.adj == "UP",   ]
  down_df <- df[df$DE.adj == "DOWN", ]
  
  write.csv(up_df,
            file.path(export_dir, paste0(cname, "_FDR_UP.csv")),
            row.names = FALSE)
  write.csv(down_df,
            file.path(export_dir, paste0(cname, "_FDR_DOWN.csv")),
            row.names = FALSE)
  
  cat(sprintf("  %s -> UP: %d, DOWN: %d (BH-adj p-value) \n", 
              cname, nrow(up_df), nrow(down_df)))
}))

cat("\n  Protein lists saved to:", export_dir, "\n")

