# === PROJECT: CR Proteomics - Muscle/Plasma =================================
# === Script: CRE/PLA Pooling Diagnostic

# === PURPOSE: Assess whether CRE and PLA arms can be validly pooled into a
# single SURV group for the SURV vs CTL comparison. Runs three checks:
#   (1) Randomization balance at baseline (Baseline_CREvPLA)
#   (2) Global supplement x time interaction rate (Interaction_supp)
#   (3) Arm-specific concordance for proteins actually driving the pooled
#       Training_SURVvCTL / Training_SURV results (the check that matters
#       most, since it's about the specific findings being reported, not
#       just proteome-wide interaction prevalence)
#
# This script reports evidence -- it does not render an automated pooling
# decision. Whether the evidence below supports pooling is a judgment call
# that should be made (and stated in methods) before looking at the pooled
# SURV vs CTL results themselves.

# === INPUTS:  limma_results_all_<tissue>.xlsx (sheets 1-5, 11-13)
# === OUTPUTS: summary tables (xlsx) + diagnostic plots (pdf/png)

# --- 0: Setup ----------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               patchwork, scales, writexl)

# --- Configuration -------------------------------------------------------
config <- list(
  limma_muscle = "03_analysis/c_data/muscle/limma",
  # NOTE: your CRE-specific fig_02 script uses
  # "03_analysis/b_reports/plasma/all_data/limma" (no "reformat"), while your
  # SURV fig_02 script uses the path below (with "reformat"). Confirm which
  # is correct for your plasma limma output before running -- left as the
  # "reformat" version here since that's the convention used everywhere else.
  limma_plasma = "03_analysis/c_data/plasma/limma",
  fig_dir      = "04_figures/fig_s2_supp_check/b_figures",
  data_dir     = "04_figures/fig_s2_supp_check/c_data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette ---------------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

tissue_colors <- c(Muscle = colors$pink_dark, Plasma = colors$aqua_dark)

# This script intentionally overrides the shared theme_cr() with a bordered-
# panel + subtitle variant suited to diagnostic/supplementary check plots.
theme_cr <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(face = "bold"),
      axis.title      = element_text(face = "bold", size = base_size),
      axis.text       = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = base_size),
      plot.title      = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle   = element_text(size = base_size - 1, color = colors$gray_dark),
      panel.border    = element_rect(colour = colors$gray_dark, fill = NA, linewidth = 0.6)
    )
}

pi_thresh <- 0.05

# --- 1: Load Data ----------------------------------------------------------
load_limma_muscle <- function(dir, tissue_label) {
  fpath <- file.path(dir, sprintf("limma_results_all_%s_cv.xlsx", tolower(tissue_label)))
  list(
    Baseline_supp  = read_excel(fpath, sheet = 2),
    Post_training_supp  = read_excel(fpath, sheet = 3),
    Training_CRE      = read_excel(fpath, sheet = 4),
    Training_PLA      = read_excel(fpath, sheet = 5),
    Interaction_supp  = read_excel(fpath, sheet = 6),
    Cancer = read_excel(fpath, sheet = 7),
    Residual = read_excel(fpath, sheet = 8),
    Training = read_excel(fpath, sheet = 9)
  )
}

load_limma_plasma <- function(dir, tissue_label) {
  fpath <- file.path(dir, sprintf("limma_results_all_%s.xlsx", tolower(tissue_label)))
  list(
    Baseline_supp  = read_excel(fpath, sheet = 2),
    Post_training_supp  = read_excel(fpath, sheet = 3),
    Training_CRE      = read_excel(fpath, sheet = 4),
    Training_PLA      = read_excel(fpath, sheet = 5),
    Interaction_supp  = read_excel(fpath, sheet = 6),
    Cancer = read_excel(fpath, sheet = 7),
    Residual = read_excel(fpath, sheet = 8),
    Training = read_excel(fpath, sheet = 9)
  )
}

limma_results <- list(
  Muscle = load_limma_muscle(config$limma_muscle, "Muscle"),
  Plasma = load_limma_plasma(config$limma_plasma, "Plasma")
)

# --- 2: Check 1 & 2 -- Baseline balance + Interaction rate ------------------
# Under a true null (no baseline imbalance / no interaction), the raw p<0.05
# rate should sit near 5% by construction. Pi.Val and FDR counts are reported
# alongside for context but don't have as clean a "5% under null" reference,
# since they're already adjusted/composite metrics.

compute_sig_summary <- function(df, contrast_label, tissue_label) {
  n_total <- nrow(df)
  tibble(
    tissue          = tissue_label,
    contrast        = contrast_label,
    n_total         = n_total,
    n_sig_pi        = sum(df$Pi.Val    < pi_thresh, na.rm = TRUE),
    pct_sig_pi      = round(100 * n_sig_pi / n_total, 2),
    n_sig_fdr05     = sum(df$adj.P.Val < 0.05, na.rm = TRUE),
    pct_sig_fdr05   = round(100 * n_sig_fdr05 / n_total, 2),
    n_sig_fdr10     = sum(df$adj.P.Val < 0.10, na.rm = TRUE),
    pct_sig_fdr10   = round(100 * n_sig_fdr10 / n_total, 2),
    n_sig_p05_raw   = sum(df$P.Value   < 0.05, na.rm = TRUE),
    pct_sig_p05_raw = round(100 * n_sig_p05_raw / n_total, 2)
  )
}

check_summary <- bind_rows(lapply(names(limma_results), function(tis) {
  bind_rows(
    compute_sig_summary(limma_results[[tis]]$Baseline_CREvPLA, "Baseline_CREvPLA (randomization check)", tis),
    compute_sig_summary(limma_results[[tis]]$Interaction_supp, "Interaction_supp (global interaction rate)", tis)
  )
}))

cat("\n=== Check 1 & 2: Baseline balance and interaction rate ===\n")
print(check_summary)

make_check_summary_plot <- function(summary_tbl) {
  ggplot(summary_tbl, aes(x = contrast, y = pct_sig_p05_raw, fill = tissue)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_hline(yintercept = 5, linetype = "dashed", color = colors$gray_dark) +
    geom_text(aes(label = sprintf("%.1f%%", pct_sig_p05_raw)),
              position = position_dodge(width = 0.7),
              vjust = -0.4, size = 3, fontface = "bold") +
    scale_fill_manual(values = tissue_colors, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "% proteins with raw p < 0.05",
         title = "Randomization Balance & Interaction Rate",
         subtitle = "Dashed line = 5% expected under a true null (balanced arms / no interaction)") +
    theme_cr() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1))
}

p_check_summary <- make_check_summary_plot(check_summary)

# --- 3: Check 3 -- Arm-specific concordance for the proteins that matter ---
# Restricts to proteins significant in the POOLED contrasts you actually
# report (Training_SURVvCTL and/or Training_SURV), then checks whether
# Training_CRE and Training_PLA agree in sign/magnitude for those specific
# proteins. This is more decision-relevant than the genome-wide interaction
# rate in Check 2 -- it asks whether pooling distorts your actual findings,
# not just whether an interaction exists somewhere in the proteome.

build_concordance_data <- function(tissue_results) {
  pooled_sig <- tissue_results$Training %>% 
    filter(Pi.Val < pi_thresh) %>% 
    pull(uniprot)
  
  cre_df <- tissue_results$Training_CRE %>%
    dplyr::select(uniprot, gene, logFC_CRE = logFC, PiVal_CRE = Pi.Val)
  pla_df <- tissue_results$Training_PLA %>%
    dplyr::select(uniprot, logFC_PLA = logFC, PiVal_PLA = Pi.Val)
  
  cre_df %>%
    inner_join(pla_df, by = "uniprot") %>%
    filter(uniprot %in% pooled_sig) %>%
    mutate(concordant = sign(logFC_CRE) == sign(logFC_PLA))
}

concordance_data <- lapply(limma_results, build_concordance_data)

compute_concordance_stats <- function(df, tissue_label) {
  if (nrow(df) < 3) {
    cat(sprintf("  [%s] Fewer than 3 pooled-significant proteins -- skipping cor.test\n", tissue_label))
    return(tibble(
      tissue = tissue_label, n = nrow(df), r = NA_real_, p_value = NA_real_,
      pct_concordant = if (nrow(df) > 0) round(100 * mean(df$concordant), 1) else NA_real_
    ))
  }
  
  ct <- cor.test(df$logFC_CRE, df$logFC_PLA)
  tibble(
    tissue         = tissue_label,
    n              = nrow(df),
    r              = unname(ct$estimate),
    p_value        = ct$p.value,
    pct_concordant = round(100 * mean(df$concordant), 1)
  )
}

concordance_stats <- bind_rows(lapply(names(concordance_data), function(tis) {
  compute_concordance_stats(concordance_data[[tis]], tis)
}))

cat("\n=== Check 3: Arm-specific concordance (pooled-significant proteins only) ===\n")
print(concordance_stats)

discordant_proteins <- bind_rows(lapply(names(concordance_data), function(tis) {
  concordance_data[[tis]] %>%
    filter(!concordant) %>%
    mutate(tissue = tis, divergence = abs(logFC_CRE - logFC_PLA)) %>%
    arrange(desc(divergence))
}))

cat(sprintf("\nTotal discordant proteins across tissues: %d\n", nrow(discordant_proteins)))
if (nrow(discordant_proteins) > 0) {
  cat("(sign-flipped between Training_CRE and Training_PLA -- worth reviewing individually\n")
  cat(" before trusting the pooled estimate for these specific proteins)\n")
  print(discordant_proteins %>% dplyr::select(tissue, gene, logFC_CRE, logFC_PLA, divergence) %>% head(20))
}

# --- 4: Concordance scatter plots -------------------------------------------

make_concordance_scatter <- function(df, tissue_label) {
  if (nrow(df) == 0) {
    return(ggplot() +
             labs(title = sprintf("Arm-Specific Concordance \u2014 %s", tissue_label),
                  subtitle = "No proteins significant in pooled contrasts") +
             theme_cr())
  }
  
  stats_row <- concordance_stats %>% filter(tissue == tissue_label)
  r_val     <- stats_row$r
  p_val     <- stats_row$p_value
  pct_conc  <- stats_row$pct_concordant
  n_prot    <- stats_row$n
  
  stat_label <- if (is.na(r_val)) {
    sprintf("n = %d (too few for r)\n%.1f%% sign-concordant", n_prot, pct_conc)
  } else {
    p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
    sprintf("r = %.2f, %s\n%.1f%% sign-concordant\nn = %d", r_val, p_label, pct_conc, n_prot)
  }
  
  ggplot(df, aes(x = logFC_CRE, y = logFC_PLA, color = concordant)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey40") +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_manual(values = c(`TRUE` = colors$blue, `FALSE` = colors$red),
                       labels = c(`TRUE` = "Concordant", `FALSE` = "Discordant"),
                       name = NULL) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = stat_label, size = 3, lineheight = 0.95, fontface = "bold",
             color = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    labs(x = "logFC \u2014 Training_CRE", y = "logFC \u2014 Training_PLA",
         title = sprintf("Arm-Specific Concordance \u2014 %s", tissue_label),
         subtitle = "Proteins significant in pooled Training_SURV") +
    theme_cr()
}

p_concordance_muscle <- make_concordance_scatter(concordance_data$Muscle, "Muscle")
p_concordance_plasma <- make_concordance_scatter(concordance_data$Plasma, "Plasma")

# --- 5: Assemble & Save ------------------------------------------------------

combined <- p_check_summary / (p_concordance_muscle + p_concordance_plasma) +
  plot_layout(heights = c(0.8, 1))

ggsave(file.path(config$fig_dir, "pooling_diagnostic_combined.pdf"),
       combined, width = 12, height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "pooling_diagnostic_combined.png"),
       combined, width = 12, height = 10, dpi = 300, bg = "white")

ggsave(file.path(config$fig_dir, "pooling_check_summary.pdf"),
       p_check_summary, width = 8, height = 5, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "pooling_concordance_muscle.pdf"),
       p_concordance_muscle, width = 6, height = 5, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "pooling_concordance_plasma.pdf"),
       p_concordance_plasma, width = 6, height = 5, device = cairo_pdf)

writexl::write_xlsx(
  list(
    Check1_2_Summary       = check_summary,
    Check3_Concordance     = concordance_stats,
    Discordant_Proteins    = discordant_proteins,
    Concordance_Data_All   = bind_rows(
      concordance_data$Muscle %>% mutate(tissue = "Muscle"),
      concordance_data$Plasma %>% mutate(tissue = "Plasma")
    )
  ),
  path = file.path(config$data_dir, "pooling_diagnostic_results.xlsx")
)

cat("\n>> Pooling diagnostic complete. Outputs saved to:\n")
cat("  ", config$fig_dir, "\n")
cat("  ", config$data_dir, "\n")

cat("\n=== Interpretation guide (this script reports evidence, not a verdict) ===\n")
cat("1. Baseline_CREvPLA: raw p<0.05 rate near 5% supports randomization balance.\n")
cat("   Substantially elevated -> investigate before pooling even at baseline.\n")
cat("2. Interaction_supp: raw p<0.05 rate near 5% supports no broad supplement-\n")
cat("   driven divergence. Elevated -> some proteome-wide interaction signal,\n")
cat("   though this alone doesn't invalidate pooling for your specific hits.\n")
cat("3. Concordance (the check that matters most): high r, high %-concordant,\n")
cat("   few discordant proteins among your pooled-significant hits supports\n")
cat("   pooling for THOSE results specifically. Any discordant proteins listed\n")
cat("   above are candidates for individual review/footnoting rather than\n")
cat("   trusting the pooled estimate at face value.\n")

# --- Diagnostic: why did proteins drop out of the concordance join? --------
diagnose_concordance_dropout <- function(tissue_results, tissue_label) {
  training_surv_sig <- tissue_results$Training %>%
    filter(Pi.Val < pi_thresh) %>%
    dplyr::select(uniprot, gene, logFC_SURV = logFC, PiVal_SURV = Pi.Val)
  
  cre_uniprots <- tissue_results$Training_CRE$uniprot
  pla_uniprots <- tissue_results$Training_PLA$uniprot
  
  training_surv_sig %>%
    mutate(
      in_cre_sheet = uniprot %in% cre_uniprots,
      in_pla_sheet = uniprot %in% pla_uniprots,
      status = case_when(
        in_cre_sheet & in_pla_sheet   ~ "Present in both (should survive join)",
        in_cre_sheet & !in_pla_sheet  ~ "Missing from Training_PLA sheet",
        !in_cre_sheet & in_pla_sheet  ~ "Missing from Training_CRE sheet",
        TRUE                          ~ "Missing from both (filtered from CRE dataset entirely)"
      )
    ) %>%
    count(status) %>%
    mutate(tissue = tissue_label, .before = 1)
}

dropout_summary <- bind_rows(
  diagnose_concordance_dropout(limma_results$Muscle, "Muscle"),
  diagnose_concordance_dropout(limma_results$Plasma, "Plasma")
)

cat("\n=== Why proteins drop out of the concordance join ===\n")
print(dropout_summary)

