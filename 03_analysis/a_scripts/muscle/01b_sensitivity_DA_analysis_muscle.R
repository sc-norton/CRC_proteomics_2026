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
  ),
  sens = list(
    report = file.path(config$report_dir, "sensitivity")
  )
)

invisible(lapply(dir_tree, function(ds_dirs) {
  lapply(ds_dirs, dir.create, recursive = TRUE, showWarnings = FALSE)
}))

# --- 1: Load Data ----------------------
CRC <- readRDS(config$norm_RDS)

CRC$metadata[c("sample_id", "blood_pct")]

# Arm-stratified grouping for the sensitivity analysis:
# CRE_T1, CRE_T2, PLA_T1, PLA_T2, CTL_T1
CRC$metadata$arm_time <- ifelse(as.character(CRC$metadata$cancer_time) == "CTL_T1",
                                "CTL_T1", as.character(CRC$metadata$supp_time))
print(table(CRC$metadata$cancer_time, CRC$metadata$arm_time, useNA = "ifany"))
stopifnot(setequal(unique(CRC$metadata$arm_time),
                   c("CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2", "CTL_T1")))

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

# Arm-adjusted equivalents of contrast_specs_cancer (equal weight per arm).
# Names match contrast_specs_cancer so the sensitivity section can pair them.
contrast_specs_cancer_adj <- c(
  "Cancer   = (CRE_T1 + PLA_T1)/2 - CTL_T1",
  "Residual = (CRE_T2 + PLA_T2)/2 - CTL_T1",
  "Training = ((CRE_T2 - CRE_T1) + (PLA_T2 - PLA_T1))/2"
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

# Sensitivity: arm-stratified model (CRE/PLA kept separate)
out_sens_cv <- run_limma_analysis(CRC, contrast_specs_cancer_adj,
                                  group_col  = "arm_time",
                                  covariates = "blood_pct",
                                  label      = "cancer_adj_cv",
                                  dirs       = dir_tree$sens)

out_sens <- run_limma_analysis(CRC, contrast_specs_cancer_adj,
                               group_col = "arm_time",
                               label     = "cancer_adj",
                               dirs      = dir_tree$sens)

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

# --- 6b: Sensitivity Analysis ----------
# Pooled CRC model (cancer_time) vs. arm-stratified model (arm_time),
# plus pooling diagnostics from the supplement model.

# (a) Pre-specified criteria ------------------------------------------
# Fix these BEFORE looking at the results.
sens_cfg <- list(
  hit_cols     = c(FDR_05 = "DE.adj", FDR_10 = "DE.adj.1",
                   nominal = "DE.nom", pi = "DE.pi"),
  plot_hit_col = "DE.pi",                 # which hit set to highlight in figures
  criteria     = list(r_logFC   = 0.90,   # Pearson r of logFC, all proteins
                      r_t       = 0.90,   # Pearson r of moderated t
                      slope     = c(0.80, 1.25),  # TLS slope, adjusted on pooled
                      dir_agree = 0.95),  # same sign among pooled hits
  panel_mm     = 178                      # PANEL_MD
)

# (b) Data prep -------------------------------------------------------

# Total-least-squares slope (symmetric: neither model is treated as error-free)
tls_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  v  <- prcomp(cbind(x[ok], y[ok]))$rotation[, 1]
  unname(v[2] / v[1])
}

# Per-protein table pairing pooled and arm-adjusted estimates, all contrasts
build_sens_table <- function(pooled, adjusted) {
  keep <- c("uniprot", "logFC", "t", "P.Value", "adj.P.Val",
            "DE.adj", "DE.adj.1", "DE.nom", "DE.pi")
  purrr::imap_dfr(pooled, function(p, cname) {
    a <- adjusted[[cname]]
    dplyr::inner_join(
      as.data.frame(p) %>% dplyr::select(gene, dplyr::all_of(keep)),
      as.data.frame(a) %>% dplyr::select(dplyr::all_of(keep)),
      by = "uniprot", suffix = c("_pool", "_adj")
    ) %>%
      dplyr::mutate(contrast    = cname,
                    sign_flip   = sign(logFC_pool) != sign(logFC_adj),
                    delta_logFC = logFC_adj - logFC_pool,
                    .before = 1)
  })
}

# One row per contrast x hit-definition, with pass/fail against sens_cfg
summarise_sens <- function(tab, hit_cols, crit) {
  tab %>%
    dplyr::group_by(contrast) %>%
    dplyr::group_modify(function(d, key) {
      base <- tibble::tibble(
        n_proteins = nrow(d),
        r_logFC    = cor(d$logFC_pool, d$logFC_adj, use = "complete.obs"),
        rho_logFC  = cor(d$logFC_pool, d$logFC_adj, use = "complete.obs",
                         method = "spearman"),
        slope_tls  = tls_slope(d$logFC_pool, d$logFC_adj),
        r_t        = cor(d$t_pool, d$t_adj, use = "complete.obs")
      )
      hits <- purrr::map_dfr(names(hit_cols), function(m) {
        hp <- d[[paste0(hit_cols[[m]], "_pool")]]
        ha <- d[[paste0(hit_cols[[m]], "_adj")]]
        is_hit    <- hp != "NDE"
        same_sign <- sign(d$logFC_adj) == sign(d$logFC_pool)
        set_p <- paste(d$uniprot, hp)[hp != "NDE"]   # protein + direction
        set_a <- paste(d$uniprot, ha)[ha != "NDE"]
        n_union <- length(union(set_p, set_a))
        tibble::tibble(
          hit_def      = m,
          n_hits_pool  = sum(is_hit),
          n_hits_adj   = sum(ha != "NDE"),
          dir_agree    = if (any(is_hit)) mean(same_sign[is_hit]) else NA_real_,
          retained_nom = if (any(is_hit))
            mean(d$P.Value_adj[is_hit] < 0.05 & same_sign[is_hit]) else NA_real_,
          jaccard      = if (n_union > 0)
            length(intersect(set_p, set_a)) / n_union else NA_real_
        )
      })
      dplyr::bind_cols(hits, base[rep(1, nrow(hits)), ])
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      pass_r_logFC = r_logFC >= crit$r_logFC,
      pass_r_t     = r_t >= crit$r_t,
      pass_slope   = dplyr::between(slope_tls, crit$slope[1], crit$slope[2]),
      pass_dir     = is.na(dir_agree) | dir_agree >= crit$dir_agree,
      robust       = pass_r_logFC & pass_r_t & pass_slope & pass_dir
    )
}

# Pooling diagnostics from the supplement model:
#   null-rate checks for Baseline_supp / Interaction_supp, and
#   CRE vs PLA training-response concordance
pool_diagnostics <- function(supp_res, pooled_training, hit_col) {
  null_rate <- purrr::map_dfr(c("Baseline_supp", "Interaction_supp"), function(cn) {
    p <- supp_res[[cn]]$P.Value
    tibble::tibble(
      contrast  = cn,
      n         = length(p),
      prop_p05  = mean(p < 0.05),         # expect ~0.05 under no difference
      pi0       = limma::propTrueNull(p), # expect ~1
      n_FDR05   = sum(supp_res[[cn]]$adj.P.Val < 0.05),
      n_FDR10   = sum(supp_res[[cn]]$adj.P.Val < 0.10)
    )
  })
  
  arm <- dplyr::inner_join(
    as.data.frame(supp_res$Training_CRE) %>% dplyr::select(uniprot, logFC),
    as.data.frame(supp_res$Training_PLA) %>% dplyr::select(uniprot, logFC),
    by = "uniprot", suffix = c("_CRE", "_PLA")
  )
  pooled_hits <- pooled_training$uniprot[pooled_training[[hit_col]] != "NDE"]
  same    <- sign(arm$logFC_CRE) == sign(arm$logFC_PLA)
  in_hits <- arm$uniprot %in% pooled_hits
  
  concordance <- tibble::tibble(
    r_logFC        = cor(arm$logFC_CRE, arm$logFC_PLA, use = "complete.obs"),
    rho_logFC      = cor(arm$logFC_CRE, arm$logFC_PLA, use = "complete.obs",
                         method = "spearman"),
    dir_agree_all  = mean(same, na.rm = TRUE),
    hit_def        = hit_col,
    n_pooled_hits  = sum(in_hits),
    dir_agree_hits = if (any(in_hits)) mean(same[in_hits], na.rm = TRUE) else NA_real_
  )
  
  list(null_rate   = null_rate,
       concordance = concordance,
       arm_table   = arm %>% dplyr::mutate(pooled_training_hit = in_hits))
}

# (c) Plot builders ---------------------------------------------------

# Pooled (x) vs adjusted (y); identity line; pooled hits coloured
plot_sens_scatter <- function(tab, metrics, var, hit_col) {
  x <- paste0(var, "_pool"); y <- paste0(var, "_adj")
  d <- tab %>% dplyr::mutate(hit = .data[[paste0(hit_col, "_pool")]])
  lbl <- metrics %>%
    dplyr::distinct(contrast, r_logFC, r_t, slope_tls) %>%
    dplyr::mutate(label = if (var == "logFC")
      sprintf("r = %.2f\nslope = %.2f", r_logFC, slope_tls)
      else sprintf("r = %.2f", r_t))
  
  ggplot(d, aes(.data[[x]], .data[[y]])) +
    geom_blank(aes(x = .data[[y]], y = .data[[x]])) +   # square axes per panel
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                linewidth = 0.3, colour = "grey40") +
    geom_point(data = dplyr::filter(d, hit == "NDE"),
               colour = "grey75", size = 0.5, alpha = 0.6) +
    geom_point(data = dplyr::filter(d, hit != "NDE"),
               aes(colour = hit), size = 1) +
    geom_text(data = lbl, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 2.6, inherit.aes = FALSE) +
    scale_colour_manual(values = c(UP = "#B2182B", DOWN = "#2166AC"),
                        name = paste0("Pooled hit\n(", hit_col, ")")) +
    facet_wrap(~ contrast, nrow = 1, scales = "free") +
    labs(x = paste(var, "(pooled CRC)"),
         y = paste(var, "(arm-adjusted)")) +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), aspect.ratio = 1)
}

# P-value histograms for the pooling-validity contrasts
plot_null_hist <- function(supp_res) {
  d <- purrr::map_dfr(c("Baseline_supp", "Interaction_supp"), function(cn)
    tibble::tibble(contrast = cn, P.Value = supp_res[[cn]]$P.Value))
  ggplot(d, aes(P.Value)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey60", colour = "white") +
    geom_hline(aes(yintercept = n / 20),
               data = d %>% dplyr::count(contrast),
               linetype = 2, colour = "#B2182B") +       # uniform expectation
    facet_wrap(~ contrast, nrow = 1) +
    labs(x = "Nominal P-value", y = "Proteins") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank())
}

# CRE vs PLA training logFC
plot_arm_concordance <- function(arm_table, conc) {
  ggplot(arm_table, aes(logFC_CRE, logFC_PLA)) +
    geom_blank(aes(x = logFC_PLA, y = logFC_CRE)) +
    geom_hline(yintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_vline(xintercept = 0, linewidth = 0.2, colour = "grey80") +
    geom_abline(slope = 1, intercept = 0, linetype = 2,
                linewidth = 0.3, colour = "grey40") +
    geom_point(data = dplyr::filter(arm_table, !pooled_training_hit),
               colour = "grey75", size = 0.5, alpha = 0.6) +
    geom_point(data = dplyr::filter(arm_table, pooled_training_hit),
               colour = "black", size = 1) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, size = 2.6,
             label = sprintf("r = %.2f\nsign agree (hits) = %.0f%%",
                             conc$r_logFC, 100 * conc$dir_agree_hits)) +
    labs(x = "Training logFC (CRE)", y = "Training logFC (PLA)") +
    theme_bw(base_size = 9) +
    theme(panel.grid = element_blank(), aspect.ratio = 1)
}

# (d) Run for both covariate settings --------------------------------

sens_runs <- list(
  cv = list(pooled = out_cancer_cv$results_list,
            adj    = out_sens_cv$results_list,
            supp   = out_supp_cv$results_list),
  nocv = list(pooled = out_cancer$results_list,
              adj    = out_sens$results_list,
              supp   = out_supp$results_list)
)

sens_out <- purrr::imap(sens_runs, function(r, v) {
  cat(sprintf("\n>> Sensitivity analysis: %s\n", v))
  out_dir <- dir_tree$sens$report
  
  tab     <- build_sens_table(r$pooled, r$adj)
  metrics <- summarise_sens(tab, sens_cfg$hit_cols, sens_cfg$criteria)
  diag    <- pool_diagnostics(r$supp, r$pooled$Training, sens_cfg$plot_hit_col)
  
  # Console summary for the highlighted hit definition
  plot_hit_name <- names(sens_cfg$hit_cols)[sens_cfg$hit_cols == sens_cfg$plot_hit_col]
  print(metrics %>%
          dplyr::filter(hit_def == plot_hit_name) %>%
          dplyr::select(contrast, r_logFC, slope_tls, r_t,
                        n_hits_pool, dir_agree, retained_nom, robust))
  print(diag$null_rate)
  print(diag$concordance)
  
  # Figures (cairo_pdf, PANEL_MD width)
  w <- sens_cfg$panel_mm
  ggsave(file.path(out_dir, sprintf("01_sens_logFC_scatter_%s_muscle.pdf", v)),
         plot_sens_scatter(tab, metrics, "logFC", sens_cfg$plot_hit_col),
         device = cairo_pdf, width = w, height = 70, units = "mm")
  ggsave(file.path(out_dir, sprintf("02_sens_t_scatter_%s_muscle.pdf", v)),
         plot_sens_scatter(tab, metrics, "t", sens_cfg$plot_hit_col),
         device = cairo_pdf, width = w, height = 70, units = "mm")
  ggsave(file.path(out_dir, sprintf("03_pooling_pval_hist_%s_muscle.pdf", v)),
         plot_null_hist(r$supp),
         device = cairo_pdf, width = w * 0.75, height = 60, units = "mm")
  ggsave(file.path(out_dir, sprintf("04_arm_concordance_%s_muscle.pdf", v)),
         plot_arm_concordance(diag$arm_table, diag$concordance),
         device = cairo_pdf, width = w / 2, height = 75, units = "mm")
  
  # Companion xlsx
  per_contrast <- split(tab, tab$contrast)
  names(per_contrast) <- paste0("proteins_", names(per_contrast))
  write_xlsx(c(list(metrics          = metrics,
                    pooling_null     = diag$null_rate,
                    pooling_concord  = diag$concordance,
                    arm_training_lfc = diag$arm_table),
               per_contrast),
             file.path(out_dir, sprintf("sensitivity_%s_muscle.xlsx", v)))
  
  list(table = tab, metrics = metrics, diagnostics = diag)
})

# Report bundle (avoids recomputation in report scripts)
saveRDS(list(sens_out = sens_out, sens_cfg = sens_cfg),
        file.path(config$data_dir, "00_sensitivity_report_data_muscle.rds"))

cat("\n  Sensitivity outputs saved to:", dir_tree$sens$report, "\n")

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