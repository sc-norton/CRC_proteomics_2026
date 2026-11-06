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

# --- 1: Load Data ----------------------
CRC <- readRDS(config$norm_RDS)

CRC$metadata[c("sample_id", "blood_pct", "B_M_ratio")]

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

# Covariate models. Every analysis below is run once per entry.
#   - Names become the label suffix and the report subfolder.
#     "none" -> label "supp"; "cv" -> "supp_cv"; "cv2" -> "supp_cv2"
#     (keeps your existing labels/filenames for the first two models).
#   - The first entry is the REFERENCE model the sensitivity analysis
#     (section 6b) compares every other model against.
#   - To also test both covariates jointly, add e.g.
#     cv_both = c("blood_pct", "COVARIATE_2")

covariate_models <- list(
  none = NULL,
  cv   = "blood_pct",
  cv2  = "B_M_ratio"  
)

ref_model <- names(covariate_models)[1]

analyses <- list(
  supp   = list(contrasts = contrast_specs_supp,   group_col = "supp_time"),
  cancer = list(contrasts = contrast_specs_cancer, group_col = "cancer_time")
)

# Guard: every covariate must exist and have no NAs. A covariate with NAs
# would silently drop those samples from THAT model only, so differences
# between runs would partly reflect different samples, not the covariate.
all_covs <- unique(unlist(covariate_models))
missing_cols <- setdiff(all_covs, colnames(CRC$metadata))
if (length(missing_cols) > 0) {
  stop("Covariate(s) not in metadata: ", paste(missing_cols, collapse = ", "))
}
na_counts <- vapply(CRC$metadata[all_covs], function(x) sum(is.na(x)), integer(1))
if (any(na_counts > 0)) {
  stop("Covariate(s) with missing values: ",
       paste(sprintf("%s (%d NA)", names(na_counts)[na_counts > 0],
                     na_counts[na_counts > 0]), collapse = ", "),
       ". Impute or subset first so every model uses the same samples.")
}

# Quick look at how each covariate is distributed across groups. A
# covariate that differs strongly by group will absorb part of the group
# effect -- worth knowing before interpreting a drop in DAP counts.
for (an in names(analyses)) {
  gc <- analyses[[an]]$group_col
  cat(sprintf("\n>> Covariates by %s\n", gc))
  print(
    CRC$metadata %>%
      group_by(.data[[gc]]) %>%
      summarise(across(all_of(all_covs),
                       list(mean = ~ mean(.x), sd = ~ sd(.x))),
                .groups = "drop")
  )
}

# Run every analysis x covariate model. Each run writes to its own
# subfolder (e.g. limma/supp/cv/) so proteoDA's per_contrast_results and
# plots from one run don't overwrite another's.
limma_runs <- imap(analyses, function(a, an) {
  imap(covariate_models, function(covs, cm) {
    label <- if (cm == ref_model) an else paste0(an, "_", cm)
    dirs  <- list(report = file.path(config$report_dir, an, cm))
    dir.create(dirs$report, recursive = TRUE, showWarnings = FALSE)
    run_limma_analysis(CRC, a$contrasts,
                       group_col  = a$group_col,
                       covariates = covs,
                       label      = label,
                       dirs       = dirs)
  })
})

# --- 6: Combine & Export Results -------

# Matrix sheet: identifiers + every participant sample column, built once
# directly from CRC (not from a results_list df) so it isn't tied to
# whichever group's contrast happened to produce it.
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

# All contrasts (supp + cancer) per covariate model
results_by_model <- map(set_names(names(covariate_models)), function(cm) {
  res <- do.call(c, unname(map(limma_runs, ~ .x[[cm]]$results_list)))
  map(res, ~ as.data.frame(.x)[, stats_cols, drop = FALSE])
})
da_by_model <- map(set_names(names(covariate_models)), function(cm) {
  bind_rows(map(limma_runs, ~ .x[[cm]]$da_log))
})

# --- 6a: Protein-level blood association ------------------------
# Per protein: Spearman r of normalized log2 abundance vs. each blood
# measure across samples. "_within" centers protein and blood measure on
# their cancer_time group means first, so between-group differences
# (incl. the biopsy-method confound) can't create a correlation.
# Descriptive only: paired T1/T2 samples aren't independent, so no p-values.

blood_vars <- c(blood_pct = "blood_pct", bm_ratio = "B_M_ratio")
adj_group  <- "cancer_time"   # covers every sample, incl. CTL
min_obs    <- 10

expr <- as.matrix(CRC$data)
meta <- CRC$metadata
stopifnot(all(colnames(expr) == meta$sample_id))

grp <- factor(meta[[adj_group]])
center_by_group <- function(x, g) x - ave(x, g, FUN = function(v) mean(v, na.rm = TRUE))
expr_c <- t(apply(expr, 1, center_by_group, g = grp))
n_obs  <- rowSums(!is.na(expr))

blood_assoc <- tibble(uniprot = CRC$annotation$uniprot_id)
for (nm in names(blood_vars)) {
  b <- meta[[blood_vars[[nm]]]]
  r_raw    <- as.vector(cor(t(expr),   b, method = "spearman", use = "pairwise.complete.obs"))
  r_within <- as.vector(cor(t(expr_c), center_by_group(b, grp),
                            method = "spearman", use = "pairwise.complete.obs"))
  r_raw[n_obs < min_obs] <- NA; r_within[n_obs < min_obs] <- NA
  blood_assoc[[paste0(nm, "_r_raw")]]    <- r_raw
  blood_assoc[[paste0(nm, "_r_within")]] <- r_within
}

results_by_model <- map(results_by_model, function(res_list) {
  map(res_list, ~ left_join(.x, blood_assoc, by = "uniprot"))
})

# Contrast-level check: does a contrast's t track blood association?
# A clearly positive or negative r means that contrast's signal partly
# moves with blood content.
blood_contrast_check <- imap_dfr(results_by_model, function(res_list, cm) {
  imap_dfr(res_list, function(df, cname) {
    tibble(model = cm, comparison = cname,
           r_t_vs_blood_pct = cor(df$t, df$blood_pct_r_within, method = "spearman",
                                  use = "complete.obs"),
           r_t_vs_bm_ratio  = cor(df$t, df$bm_ratio_r_within,  method = "spearman",
                                  use = "complete.obs"))
  })
})
write.csv(blood_contrast_check,
          file.path(config$report_dir, "blood_contrast_check_muscle.csv"), row.names = FALSE)

# File suffix: none -> "" (same filenames as before), others -> "_<name>"
file_suffix <- function(cm) if (cm == ref_model) "" else paste0("_", cm)

for (cm in names(covariate_models)) {
  write_xlsx(c(list(normalized_matrix = normalized_matrix), results_by_model[[cm]]),
             file.path(config$data_dir,
                       paste0("limma_results_all_muscle", file_suffix(cm), ".xlsx")))
  write.csv(da_by_model[[cm]],
            file.path(config$report_dir,
                      paste0("DA_methods_combined_muscle", file_suffix(cm), ".csv")),
            row.names = FALSE)
}

# Section 7 exports protein lists from the reference (no-covariate) model
results_combined <- results_by_model[[ref_model]]

# --- 6b: Covariate Sensitivity Analysis -

# Compares each covariate model against the reference model on:
#   (1) DAP counts per contrast and DA criterion
#   (2) Overlap of DAP calls (direction-aware): shared / lost / gained /
#       direction flips, plus Jaccard index
#   (3) Effect-size concordance across ALL proteins (not just DAPs):
#       correlation and slope of logFC, correlation of moderated t
# Counts alone can hide a lot: two models can each call 40 DAPs with only
# half of them in common. (2) and (3) show whether the covariate changes
# WHICH proteins are called and HOW MUCH estimates shift.

sens_dir <- file.path(config$report_dir, "sensitivity")
dir.create(sens_dir, recursive = TRUE, showWarnings = FALSE)

de_methods <- c("FDR < 0.05 (BH)" = "DE.adj",
                "FDR < 0.10 (BH)" = "DE.adj.1",
                "Nominal P-val"   = "DE.nom",
                "Pi-value"        = "DE.pi")
model_levels      <- names(covariate_models)
comparison_levels <- names(results_by_model[[ref_model]])

# One row per protein x contrast x model
all_res_long <- imap_dfr(results_by_model, function(res_list, cm) {
  imap_dfr(res_list, function(df, cname) {
    df %>%
      dplyr::transmute(model = cm, comparison = cname, uniprot, gene,
                       logFC, t, P.Value, adj.P.Val,
                       DE.adj, DE.adj.1, DE.nom, DE.pi)
  })
})

# (1) Counts
da_counts <- imap_dfr(da_by_model, ~ dplyr::mutate(.x, model = .y)) %>%
  dplyr::mutate(n_DAP      = n_UP + n_DOWN,
                model      = factor(model, levels = model_levels),
                comparison = factor(comparison, levels = comparison_levels))

da_counts_wide <- da_counts %>%
  dplyr::select(comparison, method, model, n_DAP) %>%
  tidyr::pivot_wider(names_from = model, values_from = n_DAP) %>%
  dplyr::arrange(method, comparison)

# (2) Overlap vs. reference, direction-aware
ref_calls_all <- all_res_long %>% dplyr::filter(model == ref_model)

overlap <- imap_dfr(de_methods, function(col, mlab) {
  ref <- ref_calls_all %>%
    dplyr::select(comparison, uniprot, ref_call = dplyr::all_of(col))
  all_res_long %>%
    dplyr::filter(model != ref_model) %>%
    dplyr::select(model, comparison, uniprot, call = dplyr::all_of(col)) %>%
    dplyr::inner_join(ref, by = c("comparison", "uniprot")) %>%
    dplyr::group_by(model, comparison) %>%
    dplyr::summarise(
      n_ref     = sum(ref_call != "NDE"),
      n_model   = sum(call != "NDE"),
      n_shared  = sum(ref_call != "NDE" & call == ref_call),
      n_lost    = sum(ref_call != "NDE" & call == "NDE"),   # DAP in ref only
      n_gained  = sum(ref_call == "NDE" & call != "NDE"),   # DAP in this model only
      n_flipped = sum(ref_call != "NDE" & call != "NDE" & call != ref_call),
      .groups   = "drop"
    ) %>%
    dplyr::mutate(
      method  = mlab,
      # NA when neither model calls anything for this contrast
      jaccard = dplyr::if_else(n_ref + n_model - n_shared > 0,
                               n_shared / (n_ref + n_model - n_shared),
                               NA_real_)
    )
}) %>%
  dplyr::mutate(model      = factor(model, levels = model_levels),
                comparison = factor(comparison, levels = comparison_levels)) %>%
  dplyr::select(method, comparison, model, dplyr::everything()) %>%
  dplyr::arrange(method, comparison, model)

# Protein-level detail: every protein that is a DAP in at least one model,
# with its call in each model side by side (FDR < 0.05 and nominal shown;
# swap in other DE.* columns if needed)
dap_detail <- all_res_long %>%
  dplyr::select(comparison, uniprot, gene, model, logFC, adj.P.Val, DE.adj, DE.nom) %>%
  dplyr::group_by(comparison, uniprot) %>%
  dplyr::filter(any(DE.adj != "NDE") | any(DE.nom != "NDE")) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(names_from  = model,
                     values_from = c(logFC, adj.P.Val, DE.adj, DE.nom)) %>%
  dplyr::arrange(factor(comparison, levels = comparison_levels), gene)

# (3) Effect-size concordance vs. reference, all proteins
ref_est <- ref_calls_all %>%
  dplyr::select(comparison, uniprot, ref_logFC = logFC, ref_t = t)

concordance <- all_res_long %>%
  dplyr::filter(model != ref_model) %>%
  dplyr::inner_join(ref_est, by = c("comparison", "uniprot")) %>%
  dplyr::group_by(model, comparison) %>%
  dplyr::summarise(
    n_proteins  = dplyr::n(),
    r_logFC     = stats::cor(logFC, ref_logFC, use = "complete.obs"),
    slope_logFC = unname(stats::coef(stats::lm(logFC ~ ref_logFC))[2]),
    r_t         = stats::cor(t, ref_t, use = "complete.obs"),
    median_abs_logFC_shift = stats::median(abs(logFC - ref_logFC), na.rm = TRUE),
    .groups     = "drop"
  ) %>%
  dplyr::mutate(model      = factor(model, levels = model_levels),
                comparison = factor(comparison, levels = comparison_levels)) %>%
  dplyr::arrange(comparison, model)

write_xlsx(
  list(counts_long  = da_counts,
       counts_wide  = da_counts_wide,
       overlap      = overlap,
       concordance  = concordance,
       dap_detail   = dap_detail),
  file.path(sens_dir, "covariate_sensitivity_muscle.xlsx")
)

# Plot: DAP counts per contrast, by model, one panel per DA criterion
p_counts <- ggplot(da_counts, aes(x = comparison, y = n_DAP, fill = model)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = n_DAP), position = position_dodge(width = 0.8),
            vjust = -0.3, size = 2.5) +
  facet_wrap(~ method, scales = "free_y") +
  labs(x = NULL, y = "Number of DAPs (UP + DOWN)", fill = "Covariate model",
       title = "DAP detection by covariate model — muscle") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(sens_dir, "dap_counts_by_model_muscle.pdf"), p_counts,
       width = 12, height = 8, device = cairo_pdf)

# Plot: logFC in each covariate model vs. reference, all proteins.
# Points on the dashed identity line = covariate didn't change the estimate.
scatter_df <- all_res_long %>%
  dplyr::filter(model != ref_model) %>%
  dplyr::inner_join(ref_est, by = c("comparison", "uniprot")) %>%
  dplyr::mutate(comparison = factor(comparison, levels = comparison_levels))

conc_labels <- concordance %>%
  dplyr::mutate(label = sprintf("r = %.2f\nslope = %.2f", r_logFC, slope_logFC))

p_scatter <- ggplot(scatter_df, aes(ref_logFC, logFC)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 0.5, alpha = 0.4) +
  geom_text(data = conc_labels, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.2, size = 2.6, inherit.aes = FALSE) +
  facet_grid(model ~ comparison) +
  labs(x = paste0("logFC (", ref_model, ")"), y = "logFC (covariate model)",
       title = "Effect-size concordance vs. no-covariate model — muscle") +
  theme_bw()

ggsave(file.path(sens_dir, "logFC_concordance_muscle.pdf"), p_scatter,
       width = 3 * length(comparison_levels),
       height = 3 * (length(model_levels) - 1) + 1,
       device = cairo_pdf, limitsize = FALSE)

cat("\n>> Covariate sensitivity summary (FDR < 0.05):\n")
print(overlap %>% dplyr::filter(method == "FDR < 0.05 (BH)"))
print(concordance)

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