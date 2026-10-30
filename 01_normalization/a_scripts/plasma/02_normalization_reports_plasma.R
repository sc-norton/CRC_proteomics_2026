# === PROJECT: CR Proteomics - Plasma =============================
# === Script: Normalization Reports

# === PURPOSE: Generate QC reports, contaminant summaries, filter logs,
# outlier diagnostics, CV plots, and protein/sample assessments from
# the saved normalization report data.

# === INPUTS:  00_report_data.rds (from normalization script)
# === OUTPUTS: plots, CSVs, and PDFs in report_dir / data_dir

# --- 0: Setup ------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, tibble, ggpubr, Hmisc, ggnewscale, scales,
               vegan, writexl)

# Configuration
config <- list(
  report_data = "01_normalization/c_data/plasma/00_report_data.rds",
  report_dir  = "01_normalization/b_reports/plasma",
  data_dir    = "01_normalization/c_data/plasma"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# Color palette
colors <- list(
  red          = "#C0392B", red_light    = "#E6A19A", red_dark     = "#7B241C",
  blue         = "#2471A3", blue_light   = "#A9CCE3", blue_dark    = "#154360",
  green        = "#4E8B6F", green_light  = "#A9C9B8", green_dark   = "#2E5E4E",
  pink         = "#D98C9C", pink_light   = "#F2C6CF", pink_dark    = "#A65C6B",
  purple       = "#7B5EA7", purple_light = "#C5B4DD", purple_dark  = "#4A3A73",
  yellow       = "#D4A017", yellow_light = "#F1D77A", yellow_dark  = "#9A6F0A",
  orange       = "#C4693A", orange_light = "#E5A07B", orange_dark  = "#8A4321",
  aqua         = "#5DAFA3", aqua_light   = "#A9D6CF", aqua_dark    = "#3A7F76",
  brown        = "#8C5A3C", brown_light  = "#C7A189", brown_dark   = "#5A3A28",
  red_pale     = "#FADBD8", blue_pale    = "#D6EAF8",
  gray_light   = "#F5F5F5", gray         = "#9E9E9E", gray_dark    = "#4A4A4A",
  black        = "#000000", white        = "#FFFFFF"
)

# Shared border theme -- appended to every plot in this script so all
# panels get a consistent border regardless of which theme_minimal()/
# theme_void() base each individual plot uses.
panel_border_theme <- theme(
  panel.border = element_rect(colour = colors$gray_dark, fill = NA, linewidth = 0.6)
)

# --- 1: Load Data --------------------------------------------------------
report_data <- readRDS(config$report_data)

names(report_data)

flagged_contaminants <- report_data$flagged_contaminants
sample_summary        <- report_data$sample_summary
outlier_diag           <- report_data$outlier_diag        # single tibble
outlier_pca             <- report_data$outlier_pca          # pca_res from outlier detection step
outlier_thresholds       <- report_data$outlier_thresholds   # list(miss_threshold, delta_threshold, delta_lower)
filter_log             <- report_data$filter_log
dal                     <- report_data$dal
dal_pre_outlier         <- report_data$dal_pre_outlier
annotation              <- dal$annotation
meta                    <- dal$metadata
dal_raw                 <- report_data$dal_raw
missingness_raw         <- report_data$missingness_raw
pca_raw                 <- report_data$pca_raw
pca_postnorm             <- report_data$pca_postnorm

# Post-outlier-removal, pre-normalization snapshot -- saved as its own RDS
# in the normalization script (not carried in report_data)
dal_prenorm <- readRDS(file.path(config$data_dir, "00_DAList_pre_norm_plasma.rds"))

# Simple output path helper (kept for any ad hoc subfolder needs)
build_output_path <- function(nm, base_dir) {
  out_dir <- file.path(base_dir, nm)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  return(out_dir)
}

# --- 1b: Removed Sample Set (for bolding flagged labels) ------------------

# Any sample present in raw data but absent from the final filtered/
# normalized DAList. NOTE: unlike muscle, plasma has no B:M ratio filter,
# so the only source of sample removal is consensus outlier detection --
# removed_samples and outlier_samples will be identical for plasma.
removed_samples <- setdiff(dal_raw$metadata$sample_id, dal$metadata$sample_id)

outlier_samples <- outlier_diag %>%
  filter(consensus_outlier) %>%
  pull(sample_id)

# Bold flagged sample IDs using base R plotmath expressions.
# (ggtext::element_markdown() cannot render markdown at angle != 0, so we use
# expression()/bquote() labels instead, which render bold + rotated natively.)
bold_flagged_labels <- function(x, flagged) {
  lapply(x, function(lbl) {
    if (lbl %in% flagged) bquote(bold(.(lbl))) else lbl
  })
}

# Relabel raw group values for display (underlying metadata untouched).
# cancer_time values are CTL_T1 / CRC_T1 / CRC_T2 (normalization script
# renames SURV -> CRC before assembling the DAList, mirroring muscle).
recode_group_labels <- function(x) {
  dplyr::case_when(
    x == "CTL_T1" ~ "CTL",
    x == "CRC_T1" ~ "CRC_PRE",
    x == "CRC_T2" ~ "CRC_POST",
    TRUE          ~ x
  )
}

# --- 2: Flagged Contaminants ---------------------------------------------

# Plasma has a single contaminant class (keratin) rather than muscle's
# HPA + blood-blacklist tiers, so there's no "all sources combined" summary
# needed here -- flagged_contaminants already is the complete flagged set.
write.csv(flagged_contaminants,
          file.path(config$report_dir, "01_flagged_contaminants_plasma.csv"),
          row.names = FALSE)
cat("\n Saved flagged contaminant CSV to:", config$report_dir, "\n")

# --- 3: Filter Log & Filtered Protein List --------------------------------

write_csv(filter_log, file.path(config$report_dir, "02_filtering_effects_plasma.csv"))

kept_uniprots    <- rownames(dal$data)
removed_uniprots <- setdiff(dal_raw$annotation$uniprot_id, kept_uniprots)

# reason is assigned by replaying the normalization script's actual filter
# order (keratin -> missingness), since plasma has only these two filter
# stages (no HPA/blood-blacklist tiers, no B:M ratio filter).
filtered_proteins <- dal_raw$annotation %>%
  filter(uniprot_id %in% removed_uniprots) %>%
  dplyr::select(uniprot_id, gene, description, is_keratin) %>%
  mutate(
    reason = case_when(
      is_keratin ~ "Keratin filter",
      TRUE       ~ "Missingness filter"
    )
  )

write_csv(filtered_proteins,
          file.path(config$report_dir, "03_filtered_proteins_plasma.csv"))

cat(" Saved: 02_filtering_effects_plasma.csv\n")
cat(" Saved: 03_filtered_proteins_plasma.csv\n")

# --- 4: Filter Plots -------------------------------------------------------

make_filter_plot <- function(log, label) {
  log_plot <- log %>% filter(!is.na(n_after))
  ggplot(mutate(log_plot, step = factor(step, levels = step)),
         aes(x = step, y = n_after)) +
    geom_col(fill = colors$blue, width = 0.5) +
    geom_text(aes(label = n_after), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "Number of Proteins",
         title = "Proteins Retained Through Filtering") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1)) +
    panel_border_theme
}

p_filter <- make_filter_plot(filter_log, "plasma")
ggsave(file.path(config$report_dir, "04_filtering_plot_plasma.pdf"),
       p_filter, width = 8, height = 5)
cat(" Saved: 04_filtering_plot_plasma.pdf\n")

make_protein_missingness_plot <- function(df, label) {
  ggplot(df, aes(x = pct_missing)) +
    geom_histogram(binwidth = 5, boundary = 0, fill = colors$blue, color = "white") +
    labs(x = "% Missing", y = "Number of Proteins",
         title = "Per-Protein Missingness (Raw)") +
    theme_minimal(base_size = 13) +
    panel_border_theme
}

make_sample_missingness_plot <- function(df, label, flagged_samples = character(0)) {
  ggplot(df, aes(x = reorder(sample_id, n_missing), y = n_missing)) +
    geom_col(fill = colors$orange, width = 0.7) +
    scale_x_discrete(labels = function(x) bold_flagged_labels(x, flagged_samples)) +
    labs(x = "Sample", y = "Number of Missing Proteins",
         title = "Per-Sample Missingness (Raw, Pre-Normalization)") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
    panel_border_theme
}

make_pca_plot <- function(pca_res, meta, group_col, label, title_suffix, group_colors = NULL,
                          run_permanova = TRUE, permanova_permutations = 999) {
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  pc_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(meta %>% dplyr::select(sample_id, group = all_of(group_col)), by = "sample_id") %>%
    mutate(group = recode_group_labels(group))
  
  p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(aes(fill = group), type = "norm", level = 0.95,
                 geom = "polygon", alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_point(size = 3) +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = title_suffix, color = NULL) +
    theme_minimal(base_size = 11) +
    panel_border_theme
  
  if (!is.null(group_colors)) {
    p <- p +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors, guide = "none")
  }
  
  if (run_permanova) {
    perm_df <- as.data.frame(pca_res$x) %>%
      rownames_to_column("sample_id") %>%
      left_join(meta %>% dplyr::select(sample_id, group = all_of(group_col)), by = "sample_id") %>%
      mutate(group = recode_group_labels(group)) %>%
      filter(!is.na(group))
    
    permanova_label <- tryCatch({
      dist_mat <- vegan::vegdist(perm_df %>% dplyr::select(-sample_id, -group),
                                 method = "euclidean")
      
      set.seed(42)
      res <- vegan::adonis2(dist_mat ~ group, data = perm_df,
                            permutations = permanova_permutations)
      
      r2_val  <- round(res$R2[1], 3)
      p_val   <- res$`Pr(>F)`[1]
      p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
      
      sprintf("PERMANOVA\nR\u00b2 = %.3f\n%s", r2_val, p_label)
    }, error = function(e) {
      message(sprintf("PERMANOVA failed for %s (%s): %s", label, title_suffix, e$message))
      NULL
    })
    
    if (!is.null(permanova_label)) {
      p <- p +
        annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                 label = permanova_label, size = 3, lineheight = 0.9,
                 color = colors$gray_dark, fill = colors$white,
                 label.size = 0.2)
    }
  }
  
  p
}

# --- 5: Outlier Diagnostics -------------------------------------------------

write_csv(outlier_diag,
          file.path(config$data_dir, "06_outlier_diagnostics_plasma.csv"))

miss_threshold  <- outlier_thresholds$miss_threshold
delta_threshold <- outlier_thresholds$delta_threshold
delta_lower     <- outlier_thresholds$delta_lower

sample_medians <- outlier_diag$sample_median
global_median  <- median(sample_medians, na.rm = TRUE)
mad_val        <- mad(sample_medians, na.rm = TRUE)

# A: Paired missingness
p1 <- ggplot(outlier_diag %>% filter(!is.na(delta_missing)),
             aes(x = pct_missing, y = delta_missing)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  geom_hline(yintercept = c(delta_lower, delta_threshold),
             linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "% Missing", y = "Delta Missing (T2 − T1)",
       title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  panel_border_theme

# B: PCA Mahalanobis
var_explained <- round(summary(outlier_pca)$importance[2, 1:2] * 100, 1)
pc_df <- as.data.frame(outlier_pca$x[, 1:2]) %>%
  rownames_to_column("sample_id") %>%
  left_join(outlier_diag %>% dplyr::select(sample_id, consensus_outlier), by = "sample_id")

p2 <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "B: PCA Outliers (Mahalanobis)") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none") +
  panel_border_theme

# C: MAD intensity
p3 <- ggplot(outlier_diag, aes(x = reorder(sample_id, sample_median), y = sample_median)) +
  geom_point(aes(color = consensus_outlier), size = 2.5) +
  geom_hline(yintercept = global_median, color = "black") +
  geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "Sample", y = "Median log2 intensity",
       title = "C: MAD Intensity Outliers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "none") +
  panel_border_theme

ggsave(
  file.path(config$report_dir, "06_outlier_diagnostic_plots_plasma.pdf"),
  p1 + p2 + p3 + plot_layout(ncol = 1),
  width = 10, height = 12
)
cat(" Saved: 06_outlier_diagnostic_plots_plasma.pdf\n")

# --- 6: CV Assessment -------------------------------------------------------

build_cv_data <- function(data, annot, meta) {
  
  # group_map tags every sample with BOTH its cancer_time-derived group
  # (CTL / CRC_PRE / CRC_POST) and its supp_time-derived group
  # (CRE_PRE / CRE_POST / PLA_PRE / PLA_POST) in one pass.
  group_map <- bind_rows(
    meta %>% filter(cancer == "CTL")          %>% mutate(group = "CTL"),
    meta %>% filter(cancer_time == "CRC_T1")  %>% mutate(group = "CRC_PRE"),
    meta %>% filter(cancer_time == "CRC_T2")  %>% mutate(group = "CRC_POST"),
    meta %>% filter(supp_time   == "CRE_T1")  %>% mutate(group = "CRE_PRE"),
    meta %>% filter(supp_time   == "CRE_T2")  %>% mutate(group = "CRE_POST"),
    meta %>% filter(supp_time   == "PLA_T1")  %>% mutate(group = "PLA_PRE"),
    meta %>% filter(supp_time   == "PLA_T2")  %>% mutate(group = "PLA_POST")
  ) %>%
    dplyr::select(sample_id, group) %>%
    filter(!is.na(group))
  
  rn   <- rownames(data)
  data <- as.data.frame(base::apply(data, 2, as.numeric))
  rownames(data) <- rn
  
  data %>%
    rownames_to_column("uniprot_id") %>%
    left_join(annot %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
    filter(!is.na(gene)) %>%
    dplyr::select(gene, any_of(group_map$sample_id)) %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
    left_join(group_map, by = "sample_id", relationship = "many-to-many") %>%
    filter(!is.na(group), !is.na(intensity)) %>%
    group_by(gene, group) %>%
    summarise(
      mean_int = mean(intensity, na.rm = TRUE),
      sd_int   = sd(intensity,   na.rm = TRUE),
      n        = sum(!is.na(intensity)),
      cv       = (sd_int / mean_int) * 100,
      .groups  = "drop"
    ) %>%
    filter(mean_int > 0, n >= 3)
}

# One build per data stage -- each already carries both the cancer_time
# and supp_time group tags, so SURV-side and CRE/PLA-side panels below are
# just different filters on the same three data frames.
cv_df_raw     <- build_cv_data(dal_raw$data,     dal_raw$annotation,     dal_raw$metadata)
cv_df_prenorm <- build_cv_data(dal_prenorm$data, dal_prenorm$annotation, dal_prenorm$metadata)
cv_df_norm    <- build_cv_data(dal$data,         dal$annotation,         dal$metadata)

make_cv_scatter <- function(cv_wide, x_group, y_group,
                            n_label   = 20,
                            delta_col = NULL,
                            label     = NULL) {
  
  missing <- setdiff(c(x_group, y_group), colnames(cv_wide))
  if (length(missing) > 0) stop("Groups not found in cv_wide: ", paste(missing, collapse = ", "))
  
  if (is.null(delta_col)) {
    delta_col <- paste0("cv_delta_", x_group, "v", y_group)
    cv_wide   <- cv_wide %>%
      mutate(!!delta_col := .data[[y_group]] - .data[[x_group]])
  }
  
  df        <- cv_wide %>% filter(!is.na(.data[[x_group]]), !is.na(.data[[y_group]]))
  cor_test  <- cor.test(df[[x_group]], df[[y_group]])
  r_val     <- unname(cor_test$estimate)
  p_val     <- cor_test$p.value
  p_label   <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  delta_lim <- max(abs(df[[delta_col]]), na.rm = TRUE)
  df_label  <- df %>% slice_max(abs(.data[[delta_col]]), n = n_label)
  
  ggplot(df, aes(x = .data[[x_group]], y = .data[[y_group]],
                 colour = .data[[delta_col]])) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey40", linewidth = 0.6) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", colour = colors$black, fill = "grey80",
                alpha = 0.15, linewidth = 0.8) +
    geom_text_repel(
      data          = df_label,
      aes(label     = gene),
      colour        = "black",
      fontface      = "bold",
      size          = 2.8,
      max.overlaps  = Inf,
      box.padding   = 0.4,
      point.padding = 0.2,
      segment.color = "grey60",
      segment.size  = 0.3,
      show.legend   = FALSE,
      seed          = 42
    ) +
    scale_colour_gradientn(
      colours = c(colors$purple_dark, colors$purple, colors$orange,
                  colors$yellow, colors$orange, colors$purple, colors$purple_dark),
      values  = scales::rescale(c(-delta_lim, -delta_lim * 0.5, -delta_lim * 0.15,
                                  0,
                                  delta_lim * 0.15, delta_lim * 0.5, delta_lim)),
      limits  = c(-delta_lim, delta_lim),
      name    = sprintf("\u0394CV\n(%s \u2212 %s)", y_group, x_group)
    ) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = sprintf("r = %.2f\n%s", r_val, p_label),
             size = 3.2, lineheight = 0.95, fontface = "bold",
             colour = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    labs(
      x        = bquote(CV[.(x_group)]~"(%)"),
      y        = bquote(CV[.(y_group)]~"(%)"),
      title    = sprintf("Per-protein CV: %s vs %s", x_group, y_group),
      subtitle = sprintf("%s", label)
    ) +
    theme_minimal() +
    panel_border_theme
}

ALL_GROUPS <- c("CTL", "CRC_PRE", "CRC_POST",
                "CRE_PRE", "CRE_POST", "PLA_PRE", "PLA_POST")

build_cv_wide_with_deltas <- function(cv_df) {
  cv_df %>%
    filter(group %in% ALL_GROUPS) %>%
    dplyr::select(gene, group, cv) %>%
    pivot_wider(names_from = group, values_from = cv) %>%
    mutate(
      cv_delta_t1vt2   = CRC_POST - CRC_PRE,
      cv_delta_t1vctl  = CRC_PRE  - CTL,
      cv_delta_t2vctl  = CRC_POST - CTL,
      cv_delta_ct1vct2 = CRE_POST - CRE_PRE,
      cv_delta_pt1vpt2 = PLA_POST - PLA_PRE,
      cv_delta_ct1vpt1 = CRE_PRE  - PLA_PRE,
      cv_delta_ct2vpt2 = CRE_POST - PLA_PRE
    )
}

cv_norm <- build_cv_wide_with_deltas(cv_df_norm)
cv_raw  <- build_cv_wide_with_deltas(cv_df_raw)

# CV violin plot
GROUP_COLOURS <- c(
  CTL      = colors$green,
  CRC_PRE  = colors$purple,
  CRC_POST = colors$orange,
  CRE_PRE  = colors$purple_light,
  CRE_POST = colors$orange_light,
  PLA_PRE  = colors$purple_dark,
  PLA_POST = colors$orange_dark
)

GROUPS_SURV <- c("CTL", "CRC_PRE", "CRC_POST")
GROUPS_CRE  <- c("CRE_PRE", "CRE_POST", "PLA_PRE", "PLA_POST")

compute_cv_summary_stats <- function(df, groups, n_boot = 2000) {
  df <- df %>% filter(group %in% groups, !is.na(cv))
  
  pooled_skew <- with(df, mean((cv - mean(cv))^3) / sd(cv)^3)
  use_median  <- abs(pooled_skew) > 1
  
  cat(sprintf("  CV summary stat: pooled skewness = %.2f -> using %s \u00b1 95%% CI\n",
              pooled_skew, ifelse(use_median, "Median (bootstrap)", "Mean (normal-theory)")))
  
  set.seed(42)
  df %>%
    group_by(group) %>%
    group_modify(~ {
      x <- .x$cv
      n <- length(x)
      
      if (use_median) {
        boot_medians <- replicate(n_boot, median(sample(x, n, replace = TRUE)))
        tibble(
          n = n, stat_type = "Median",
          estimate = median(x),
          ci_lower = quantile(boot_medians, 0.025, names = FALSE),
          ci_upper = quantile(boot_medians, 0.975, names = FALSE),
          y_max = max(x)
        )
      } else {
        ci <- Hmisc::smean.cl.normal(x, conf.int = 0.95)
        tibble(
          n = n, stat_type = "Mean",
          estimate = ci[["Mean"]],
          ci_lower = ci[["Lower"]],
          ci_upper = ci[["Upper"]],
          y_max = max(x)
        )
      }
    }) %>%
    ungroup()
}

make_cv_violin <- function(cv_df, groups, label, n_boot = 2000) {
  df <- cv_df %>%
    filter(group %in% groups) %>%
    mutate(group = factor(group, levels = groups))
  
  stats_df <- compute_cv_summary_stats(df, groups, n_boot = n_boot) %>%
    mutate(group = factor(group, levels = groups))
  
  y_range  <- diff(range(df$cv, na.rm = TRUE))
  stats_df <- stats_df %>%
    mutate(
      label_y   = y_max + y_range * 0.08,
      label_txt = sprintf("%.1f [%.1f\u2013%.1f]",
                          estimate, ci_lower, ci_upper)
    )
  
  ggplot(df, aes(x = group, y = cv, fill = group, color = group)) +
    geom_violin(alpha = 0.5, linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8, linewidth = 0.5,
                 outlier.shape = NA, colour = colors$gray_dark) +
    geom_label(
      data        = stats_df,
      aes(x = group, y = label_y, label = label_txt),
      inherit.aes = FALSE,
      size        = 2.8, lineheight = 0.9, fontface = "bold",
      color       = colors$gray_dark, fill = colors$white, label.size = 0.2
    ) +
    scale_fill_manual(values   = GROUP_COLOURS[groups], guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS[groups], guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.2))) +
    labs(x = NULL, y = "CV (%)", title = "Per-protein CV by group") +
    theme_minimal() +
    panel_border_theme
}

cv_violin_specs <- list(
  list(data = cv_df_prenorm, groups = GROUPS_SURV, label = "pre-normalization",
       file = "07_cv_violin_SURV_prenorm_plasma.pdf"),
  list(data = cv_df_norm,    groups = GROUPS_SURV, label = "normalized",
       file = "07_cv_violin_SURV_normalized_plasma.pdf"),
  list(data = cv_df_prenorm, groups = GROUPS_CRE,  label = "pre-normalization",
       file = "07_cv_violin_CRE_prenorm_plasma.pdf"),
  list(data = cv_df_norm,    groups = GROUPS_CRE,  label = "normalized",
       file = "07_cv_violin_CRE_normalized_plasma.pdf")
)

walk(cv_violin_specs, function(s) {
  p <- make_cv_violin(s$data, s$groups, s$label)
  ggsave(file.path(config$report_dir, s$file), p, width = 6, height = 6, device = cairo_pdf)
  cat(sprintf(" Saved: %s\n", s$file))
})

# CV scatter plots
cv_scatter_specs <- list(
  list(data = cv_norm, x = "CRC_PRE", y = "CTL",      delta = "cv_delta_t1vctl", label = "normalized", n = 30),
  list(data = cv_norm, x = "CRC_PRE", y = "CRC_POST", delta = "cv_delta_t1vt2",  label = "normalized", n = 30),
  list(data = cv_raw,  x = "CRC_PRE", y = "CTL",      delta = "cv_delta_t1vctl", label = "raw",        n = 50),
  list(data = cv_raw,  x = "CRC_PRE", y = "CRC_POST", delta = "cv_delta_t1vt2",  label = "raw",        n = 50)
)

walk(cv_scatter_specs, function(s) {
  plot_name <- sprintf("08_cv_scatter_%s_v_%s_%s_plasma.pdf", s$x, s$y, s$label)
  ggsave(
    file.path(config$report_dir, plot_name),
    make_cv_scatter(s$data, s$x, s$y,
                    delta_col = s$delta,
                    label     = s$label,
                    n_label   = s$n),
    width = 7, height = 6, device = cairo_pdf
  )
  cat(sprintf(" Saved: %s\n", plot_name))
})

# --- 6b: Combined CV Report (scatters above violin) -----------------------
# Pairs each stage's two scatter plots (PRE vs CTL, PRE vs POST) side by
# side above that stage's violin plot.
# NOTE: the "raw" scatter data (dal_raw, pre-filtering) and the
# "pre-normalization" violin data (dal_prenorm, post-filtering/
# pre-normalization) are paired here as the two "before normalization"
# views, but come from two distinct upstream objects -- flagging in case
# an exact-object match matters more than the conceptual pairing.

make_cv_report <- function(cv_wide, cv_violin_df, groups, stage_label, n_label,
                           scatter_heights = c(1, 1)) {
  
  scatter_subtitle <- sub("^[A-Za-z0-9]+,\\s*", "", stage_label)
  
  p_scatter_ctl <- make_cv_scatter(cv_wide, "CRC_PRE", "CTL",
                                   delta_col = "cv_delta_t1vctl",
                                   label = scatter_subtitle, n_label = n_label)
  p_scatter_post <- make_cv_scatter(cv_wide, "CRC_PRE", "CRC_POST",
                                    delta_col = "cv_delta_t1vt2",
                                    label = scatter_subtitle, n_label = n_label)
  p_violin <- make_cv_violin(cv_violin_df, groups, stage_label)
  
  (p_scatter_ctl + p_scatter_post) / p_violin +
    plot_layout(heights = scatter_heights)
}

cv_report_normalized <- make_cv_report(cv_norm, cv_df_norm, GROUPS_SURV,
                                       stage_label = "SURV, normalized", n_label = 30)
cv_report_prenorm <- make_cv_report(cv_raw, cv_df_prenorm, GROUPS_SURV,
                                    stage_label = "SURV, pre-normalization", n_label = 30)

ggsave(file.path(config$report_dir, "14_cv_report_SURV_normalized_plasma.pdf"),
       cv_report_normalized, width = 12, height = 12, device = cairo_pdf)
ggsave(file.path(config$report_dir, "14_cv_report_SURV_prenorm_plasma.pdf"),
       cv_report_prenorm, width = 12, height = 12, device = cairo_pdf)

cat(" Saved: 14_cv_report_SURV_normalized_plasma.pdf\n")
cat(" Saved: 14_cv_report_SURV_prenorm_plasma.pdf\n")

# --- 7: Protein Assessment -----------------------------------------------

miss_pct_final <- rowMeans(is.na(dal$data))
protein_summary <- tibble(
  Total      = nrow(dal$data),
  Complete   = sum(miss_pct_final == 0),
  Incomplete = sum(miss_pct_final > 0)
)

write.csv(protein_summary,
          file.path(config$data_dir, "11_protein_report_plasma.csv"),
          row.names = FALSE)
cat(" Saved: 11_protein_report_plasma.csv\n")
print(protein_summary)

# --- 8: Sample Assessment -------------------------------------------------

compute_sample_report <- function(meta) {
  
  paired_subjects <- meta %>% count(pid) %>% filter(n == 2) %>% nrow()
  paired_cre      <- meta %>% filter(supp == "CRE") %>% count(pid) %>% filter(n == 2) %>% nrow()
  paired_pla      <- meta %>% filter(supp == "PLA") %>% count(pid) %>% filter(n == 2) %>% nrow()
  unpaired        <- meta %>% count(pid) %>% filter(n == 1) %>% nrow()
  
  t1_count <- sum(meta$timepoint == "T1", na.rm = TRUE)
  t2_count <- sum(meta$timepoint == "T2", na.rm = TRUE)
  
  tibble(
    total_samples = nrow(meta),
    n_CRC        = sum(meta$cancer == "CRC", na.rm = TRUE),
    n_CTL        = sum(meta$cancer == "CTL", na.rm = TRUE),
    n_CRE        = sum(meta$supp   == "CRE", na.rm = TRUE),
    n_PLA        = sum(meta$supp   == "PLA", na.rm = TRUE),
    n_male       = sum(meta$sex    == "Male",   na.rm = TRUE),
    n_female     = sum(meta$sex    == "Female", na.rm = TRUE),
    n_T1         = t1_count,
    n_T2         = t2_count,
    n_CRE_T1     = sum(meta$supp_time   == "CRE_T1",  na.rm = TRUE),
    n_CRE_T2     = sum(meta$supp_time   == "CRE_T2",  na.rm = TRUE),
    n_PLA_T1     = sum(meta$supp_time   == "PLA_T1",  na.rm = TRUE),
    n_PLA_T2     = sum(meta$supp_time   == "PLA_T2",  na.rm = TRUE),
    n_male_T1    = sum(meta$sex_time    == "Male_T1",   na.rm = TRUE),
    n_male_T2    = sum(meta$sex_time    == "Male_T2",   na.rm = TRUE),
    n_female_T1  = sum(meta$sex_time    == "Female_T1", na.rm = TRUE),
    n_female_T2  = sum(meta$sex_time    == "Female_T2", na.rm = TRUE),
    n_paired     = paired_subjects,
    n_paired_CRE = paired_cre,
    n_paired_PLA = paired_pla,
    n_unpaired   = unpaired,
    n_T1_only    = t1_count - paired_subjects,
    n_T2_only    = t2_count - paired_subjects
  )
}

sample_report <- compute_sample_report(dal$metadata)

write.csv(sample_report,
          file.path(config$report_dir, "12_sample_report_plasma.csv"),
          row.names = FALSE)
cat(" Saved: 12_sample_report_plasma.csv\n")
print(sample_report)

writeLines(capture.output(sessionInfo()),
           file.path(config$data_dir, "sessionInfo_reports.txt"))

cat("\n>> Reporting complete.\n")

# --- 9: Combined Normalization Summary -------------------------------------
cat("\n>> 9 - Combined normalization summary\n")

pca_group_colors <- c(
  CTL      = colors$green,
  CRC_PRE  = colors$purple,
  CRC_POST = colors$orange
)

p_protmiss   <- make_protein_missingness_plot(missingness_raw$protein_missing, "plasma")
p_samplemiss <- make_sample_missingness_plot(missingness_raw$sample_missing, "plasma",
                                             flagged_samples = removed_samples)
p_pca_raw    <- make_pca_plot(pca_raw, dal_raw$metadata, "cancer_time", "plasma",
                              "PCA — Raw", group_colors = pca_group_colors)
p_pca_norm   <- make_pca_plot(pca_postnorm, dal$metadata, "cancer_time", "plasma",
                              "PCA — Post-Normalization", group_colors = pca_group_colors)

# Individual saves
ggsave(file.path(config$report_dir, "09_protein_missingness_plasma.pdf"), p_protmiss, width = 7, height = 5)
ggsave(file.path(config$report_dir, "09_pca_raw_plasma.pdf"), p_pca_raw, width = 6, height = 5)
ggsave(file.path(config$report_dir, "09_pca_postnorm_plasma.pdf"), p_pca_norm, width = 6, height = 5)
ggsave(file.path(config$report_dir, "09_sample_missingness_plasma.pdf"), p_samplemiss, width = 10, height = 5)

# Combined single PDF
combined <- (p_filter + p_protmiss) /
  (p_pca_raw + p_pca_norm) /
  p_samplemiss +
  plot_layout(heights = c(1, 1, 1))

ggsave(file.path(config$report_dir, "10_normalization_summary_plasma.pdf"),
       combined, width = 14, height = 20)

cat(" Saved: 10_normalization_summary_plasma.pdf (+ individual panels)\n")

# --- 10: Combined Normalization Data Excel ---------------------------------
cat("\n>> 10 - Combined normalization_data.xlsx\n")

# Sheet: raw, pre-filtered, pre-normalized annotated data matrix
raw_matrix_export <- bind_cols(
  dal_raw$annotation %>%
    dplyr::select(uniprot_id, protein, gene, description),
  dal_raw$data
)

# Sheet: raw data summary (proteins + samples broken down by group)
raw_data_summary <- sample_summary %>%
  mutate(
    total_proteins = nrow(dal_raw$data),
    total_samples  = nrow(dal_raw$metadata)
  ) %>%
  relocate(total_proteins, total_samples, .before = 1)

# Sheet: final normalized annotated data matrix
normalized_matrix_export <- bind_cols(
  dal$annotation %>%
    dplyr::select(uniprot_id, protein, gene, description),
  dal$data
)

# Sheet: sample-level filtering summary. NOTE: unlike muscle, plasma has no
# B:M ratio filter -- the only sample-removing step is outlier detection --
# so this is a two-row log (Raw input -> Outlier removal), not three.
sample_filter_log <- tibble(
  step     = c("Raw input", "Outlier removal"),
  n_before = c(NA_integer_, nrow(dal_raw$metadata)),
  n_after  = c(nrow(dal_raw$metadata), nrow(dal$metadata))
) %>%
  mutate(
    n_removed    = n_before - n_after,
    pct_retained = round(n_after / nrow(dal_raw$metadata) * 100, 1)
  )

# Sheet: samples removed + reason. Plasma's only removal reason is
# consensus outlier (no B:M ratio filter to attribute removals to).
sample_removed_reasons <- outlier_diag %>%
  filter(sample_id %in% removed_samples) %>%
  dplyr::select(sample_id, pct_missing, delta_missing, mahal_dist, sample_median, n_flags) %>%
  mutate(
    reason = "Consensus outlier",
    detail = sprintf("n_flags = %d of 3 (missingness/PCA/MAD)", n_flags)
  ) %>%
  dplyr::select(sample_id, reason, detail, everything())

# README sheet -- describes every other sheet, in the order they appear in the workbook
readme <- tibble(
  Sheet = c(
    "Raw_Data_Matrix", "Raw_Data_Summary", "Normalized_Data_Matrix",
    "Protein_Filter_Log", "Filtered_Proteins", "Sample_Filter_Log",
    "Samples_Removed", "Outlier_Diagnostics", "Final_Sample_Summary",
    "Final_Protein_Summary"
  ),
  Description = c(
    "Raw, pre-filtered, pre-normalized annotated protein x sample intensity matrix (all proteins, all samples as loaded).",
    "Raw sample counts broken down by cancer / supp / timepoint / cancer_time / supp_time group, plus total raw protein and sample counts.",
    "Final filtered and normalized annotated protein x sample intensity matrix used for downstream analyses.",
    "Protein-level filtering log: number of proteins retained/removed at each filtering step (keratin, missingness).",
    "Proteins removed during filtering, with the reason (Keratin filter or Missingness filter).",
    "Sample-level filtering summary: number of samples retained/removed at each stage (raw input, outlier removal).",
    "Samples removed from the dataset (consensus outliers), with supporting diagnostic values.",
    "Full outlier diagnostics (missingness, PCA Mahalanobis distance, MAD intensity deviation, consensus flag) for every sample that reached the outlier-detection step.",
    "Final sample counts by group, sex, timepoint, and pairing status, for the filtered/normalized dataset.",
    "Final protein counts by completeness (fully complete vs. containing missing values) in the filtered/normalized dataset."
  )
)

normalization_data_workbook <- list(
  README                  = readme,
  Raw_Data_Matrix         = raw_matrix_export,
  Raw_Data_Summary        = raw_data_summary,
  Normalized_Data_Matrix  = normalized_matrix_export,
  Protein_Filter_Log      = filter_log,
  Filtered_Proteins       = filtered_proteins,
  Sample_Filter_Log       = sample_filter_log,
  Samples_Removed         = sample_removed_reasons,
  Outlier_Diagnostics     = outlier_diag,
  Final_Sample_Summary    = sample_report,
  Final_Protein_Summary   = protein_summary
)

write_xlsx(normalization_data_workbook,
           file.path(config$data_dir, "normalization_data_plasma.xlsx"))

cat(" Saved: normalization_data_plasma.xlsx (11 sheets)\n")