# === PROJECT: CR Proteomics - Plasma =============================
# === Script: Imputation Reports

# === PURPOSE: QC the imp4p multiple imputation -- missingness overview,
# per-protein MCAR/MNAR classification, observed-vs-imputed distribution
# diagnostics, and before/after PCA and per-protein CV comparisons.

# === INPUTS:  00_report_data_plasma.rds (from 01_imputation_plasma.R)
# === OUTPUTS: plots, CSVs, and XLSX companions in report_dir / data_dir

# --- 0: Setup ------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, tibble, stringr, ggplot2, ggrepel,
               patchwork, scales, vegan, writexl, Hmisc)

# Configuration
config <- list(
  report_data = "02_imputation/c_data/plasma/00_report_data_plasma.rds",
  report_dir  = "02_imputation/b_reports/plasma",
  data_dir    = "02_imputation/c_data/plasma"
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# Color palette -- matches 01/02_normalization_reports for cross-figure consistency
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

panel_border_theme <- theme(
  panel.border = element_rect(colour = colors$gray_dark, fill = NA, linewidth = 0.6)
)

CLASS_COLOURS <- c(
  "complete"  = colors$green,
  "MCAR-like" = colors$blue,
  "MNAR-like" = colors$orange
)

GROUP_COLOURS <- c(
  CTL      = colors$green,
  CRC_PRE  = colors$purple,
  CRC_POST = colors$orange
)

# Relabel raw condition values for display (underlying condition untouched)
recode_condition_labels <- function(x) {
  dplyr::case_when(
    x == "CTL_T1"  ~ "CTL",
    x == "CRC_T1" ~ "CRC_PRE",
    x == "CRC_T2" ~ "CRC_POST",
    TRUE           ~ as.character(x)
  )
}

# --- 1: Load Data ----------------------------------------------------------

report_data <- readRDS(config$report_data)

mat_raw           <- report_data$mat_raw
mat_imputed       <- report_data$mat_imputed
annotation        <- report_data$annotation
meta              <- report_data$meta
condition         <- report_data$condition
prot_miss_by_cond <- report_data$prot_miss_by_cond
prot_class        <- report_data$prot_class
miss_cond_summary <- report_data$miss_cond_summary
imputation_params <- report_data$imputation_params

GROUPS <- c("CTL", "CRC_PRE", "CRC_POST")

cat(sprintf(" Loaded: %d proteins x %d samples\n", nrow(mat_raw), ncol(mat_raw)))
cat(sprintf(" Imputation method: %s (MCAR: %s, MNAR: %s, iter: %d, seed: %d)\n",
            imputation_params$method, imputation_params$methodMCAR,
            imputation_params$methodMNAR, imputation_params$nb_iter,
            imputation_params$seed))

# --- 2: Missingness Overview -----------------------------------------------
cat("\n>> 2 - Missingness overview\n")

# Per-protein missingness (pre-imputation)
protein_missing <- tibble(
  uniprot_id  = rownames(mat_raw),
  pct_missing = rowMeans(is.na(mat_raw)) * 100
)

p_protein_miss <- ggplot(protein_missing, aes(x = pct_missing)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = colors$blue, color = "white") +
  labs(x = "% Missing", y = "Number of Proteins",
       title = "Per-Protein Missingness (Pre-Imputation)") +
  theme_minimal(base_size = 13) +
  panel_border_theme

ggsave(file.path(config$report_dir, "01_protein_missingness_plasma.pdf"),
       p_protein_miss, width = 7, height = 5)
cat(" Saved: 01_protein_missingness_plasma.pdf\n")

# Per-sample missingness, coloured by condition
sample_missing <- tibble(
  sample_id = colnames(mat_raw),
  n_missing = colSums(is.na(mat_raw)),
  pct_missing = colMeans(is.na(mat_raw)) * 100,
  group = recode_condition_labels(as.character(condition))
)

p_sample_miss <- ggplot(sample_missing,
                        aes(x = reorder(sample_id, pct_missing), y = pct_missing, fill = group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = GROUP_COLOURS, name = NULL) +
  labs(x = "Sample", y = "% Missing",
       title = "Per-Sample Missingness (Pre-Imputation)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "bottom") +
  panel_border_theme

ggsave(file.path(config$report_dir, "02_sample_missingness_plasma.pdf"),
       p_sample_miss, width = 10, height = 5)
cat(" Saved: 02_sample_missingness_plasma.pdf\n")

# Per-condition missingness summary
miss_cond_summary_export <- as_tibble(miss_cond_summary) %>%
  mutate(group = recode_condition_labels(condition))

write_csv(miss_cond_summary_export,
          file.path(config$data_dir, "03_condition_missingness_summary_plasma.csv"))

p_condition_miss <- ggplot(miss_cond_summary_export,
                           aes(x = factor(group, levels = GROUPS), y = pct_missing, fill = group)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d/%d samples", n_samples, n_missing)),
            vjust = -0.4, size = 3) +
  scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
  labs(x = NULL, y = "Mean % Missing per Sample",
       title = "Missingness by Condition (Pre-Imputation)") +
  theme_minimal(base_size = 13) +
  panel_border_theme

ggsave(file.path(config$report_dir, "03_condition_missingness_plasma.pdf"),
       p_condition_miss, width = 7, height = 5)
cat(" Saved: 03_condition_missingness_plasma.pdf\n")

# Per-protein MCAR/MNAR/complete classification
class_summary <- tibble(class = prot_class) %>%
  count(class) %>%
  mutate(percentage = round(n / sum(n) * 100, 1)) %>%
  arrange(match(class, c("complete", "MCAR-like", "MNAR-like")))

write_csv(class_summary,
          file.path(config$data_dir, "04_protein_missingness_classification_plasma.csv"))
write_xlsx(class_summary,
           file.path(config$data_dir, "04_protein_missingness_classification_plasma.xlsx"))

p_class <- ggplot(class_summary,
                  aes(x = 2, y = n, fill = factor(class, levels = c("complete", "MCAR-like", "MNAR-like")))) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(0.5, 2.5)) +
  scale_fill_manual(
    values = CLASS_COLOURS,
    labels = paste0(class_summary$class, " (", class_summary$percentage, "%)")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  panel_border_theme +
  labs(title = "Per-Protein Missingness Classification\n(MNAR threshold: \u0394 condition missing-rate \u2265 0.2)",
       fill = NULL)

ggsave(file.path(config$report_dir, "04_protein_missingness_classification_plasma.pdf"),
       p_class, width = 8, height = 6)
cat(" Saved: 04_protein_missingness_classification_plasma.pdf\n")

# --- 3: Distribution Diagnostics (Observed vs Imputed) ---------------------
cat("\n>> 3 - Distribution diagnostics\n")

sample_group_map <- tibble(
  sample_id = colnames(mat_raw),
  group     = recode_condition_labels(as.character(condition))
)

long_status <- tibble(
  uniprot_id = rep(rownames(mat_raw), times = ncol(mat_raw)),
  sample_id  = rep(colnames(mat_raw), each  = nrow(mat_raw)),
  observed   = as.vector(mat_raw),
  final      = as.vector(mat_imputed)
) %>%
  mutate(status = if_else(is.na(observed), "Imputed", "Observed"),
         value  = if_else(is.na(observed), final, observed)) %>%
  left_join(sample_group_map, by = "sample_id")

dist_stats_overall <- long_status %>%
  group_by(status) %>%
  summarise(n = n(), mean = mean(value), sd = sd(value), .groups = "drop")

shift_val <- dist_stats_overall$mean[dist_stats_overall$status == "Observed"] -
  dist_stats_overall$mean[dist_stats_overall$status == "Imputed"]

p_density_overall <- ggplot(long_status, aes(x = value, fill = status)) +
  geom_density(alpha = 0.5, colour = NA) +
  geom_vline(data = dist_stats_overall, aes(xintercept = mean, colour = status),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  scale_fill_manual(values = c(Observed = colors$blue, Imputed = colors$orange), name = NULL) +
  scale_colour_manual(values = c(Observed = colors$blue_dark, Imputed = colors$orange_dark)) +
  annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
           label = sprintf("Shift (obs \u2212 imp) = %.2f", shift_val),
           size = 3.2, fontface = "bold",
           colour = colors$gray_dark, fill = colors$white, label.size = 0.2) +
  labs(x = "log2 intensity", y = "Density",
       title = "Observed vs Imputed Value Distributions") +
  theme_minimal(base_size = 13) +
  panel_border_theme

ggsave(file.path(config$report_dir, "05_density_observed_v_imputed_plasma.pdf"),
       p_density_overall, width = 7, height = 5)
cat(" Saved: 05_density_observed_v_imputed_plasma.pdf\n")

dist_stats_by_group <- long_status %>%
  group_by(group, status) %>%
  summarise(n = n(), mean = mean(value), sd = sd(value), .groups = "drop")

write_csv(dist_stats_overall,  file.path(config$data_dir, "06_distribution_stats_overall_plasma.csv"))
write_csv(dist_stats_by_group, file.path(config$data_dir, "06_distribution_stats_by_condition_plasma.csv"))
write_xlsx(
  list(overall = dist_stats_overall, by_condition = dist_stats_by_group),
  file.path(config$data_dir, "06_distribution_stats_plasma.xlsx")
)

p_density_by_group <- ggplot(long_status, aes(x = value, fill = status)) +
  geom_density(alpha = 0.5, colour = NA) +
  facet_wrap(~ factor(group, levels = GROUPS)) +
  scale_fill_manual(values = c(Observed = colors$blue, Imputed = colors$orange), name = NULL) +
  labs(x = "log2 intensity", y = "Density",
       title = "Observed vs Imputed Value Distributions, by Condition") +
  theme_minimal(base_size = 12) +
  panel_border_theme

ggsave(file.path(config$report_dir, "06_density_observed_v_imputed_by_condition_plasma.pdf"),
       p_density_by_group, width = 10, height = 5)
cat(" Saved: 06_density_observed_v_imputed_by_condition_plasma.pdf\n")

# --- 4: PCA Before vs After Imputation --------------------------------------
cat("\n>> 4 - PCA before/after imputation\n")

compute_pca_matrix <- function(mat) {
  mat_fill <- mat
  for (j in seq_len(ncol(mat_fill))) {
    nas <- is.na(mat_fill[, j])
    if (any(nas)) mat_fill[nas, j] <- median(mat_fill[, j], na.rm = TRUE)
  }
  prcomp(t(mat_fill), center = TRUE, scale. = TRUE)
}

pca_before <- compute_pca_matrix(mat_raw)      # NAs median-filled, naive baseline
pca_after  <- compute_pca_matrix(mat_imputed)  # already complete, fill loop is a no-op

make_pca_plot <- function(pca_res, meta_df, group_vec, title_suffix,
                          run_permanova = TRUE, permanova_permutations = 999) {
  var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
  pc_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("sample_id") %>%
    left_join(tibble(sample_id = names(group_vec), group = unname(group_vec)), by = "sample_id")
  
  p <- ggplot(pc_df, aes(x = PC1, y = PC2, color = group)) +
    stat_ellipse(aes(fill = group), type = "norm", level = 0.95,
                 geom = "polygon", alpha = 0.12, color = NA, show.legend = FALSE) +
    geom_point(size = 3) +
    scale_color_manual(values = GROUP_COLOURS) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
         y = sprintf("PC2 (%.1f%%)", var_explained[2]),
         title = title_suffix, color = NULL) +
    theme_minimal(base_size = 11) +
    panel_border_theme
  
  if (run_permanova) {
    perm_df <- as.data.frame(pca_res$x) %>%
      rownames_to_column("sample_id") %>%
      left_join(tibble(sample_id = names(group_vec), group = unname(group_vec)), by = "sample_id")
    
    permanova_label <- tryCatch({
      dist_mat <- vegan::vegdist(perm_df %>% dplyr::select(-sample_id, -group), method = "euclidean")
      set.seed(42)
      res <- vegan::adonis2(dist_mat ~ group, data = perm_df, permutations = permanova_permutations)
      r2_val  <- round(res$R2[1], 3)
      p_val   <- res$`Pr(>F)`[1]
      p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
      sprintf("PERMANOVA\nR\u00b2 = %.3f\n%s", r2_val, p_label)
    }, error = function(e) {
      message(sprintf("PERMANOVA failed (%s): %s", title_suffix, e$message))
      NULL
    })
    
    if (!is.null(permanova_label)) {
      p <- p +
        annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
                 label = permanova_label, size = 3, lineheight = 0.9,
                 color = colors$gray_dark, fill = colors$white, label.size = 0.2)
    }
  }
  
  p
}

group_vec <- setNames(recode_condition_labels(as.character(condition)), colnames(mat_raw))

p_pca_before <- make_pca_plot(pca_before, meta, group_vec, "PCA \u2014 Pre-Imputation (median-filled)")
p_pca_after  <- make_pca_plot(pca_after,  meta, group_vec, "PCA \u2014 Post-Imputation (imp4p)")

ggsave(file.path(config$report_dir, "07_pca_pre_imputation_plasma.pdf"),  p_pca_before, width = 6, height = 5)
ggsave(file.path(config$report_dir, "07_pca_post_imputation_plasma.pdf"), p_pca_after,  width = 6, height = 5)

combined_pca <- p_pca_before + p_pca_after + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(config$report_dir, "07_pca_before_v_after_imputation_plasma.pdf"),
       combined_pca, width = 12, height = 6)
cat(" Saved: 07_pca_before_v_after_imputation_plasma.pdf (+ individual panels)\n")

# --- 5: CV Before vs After Imputation --------------------------------------
cat("\n>> 5 - CV before/after imputation\n")

build_cv_data <- function(mat, annot, group_vec) {
  rn      <- rownames(mat)
  mat_num <- as.data.frame(base::apply(mat, 2, as.numeric))
  rownames(mat_num) <- rn
  
  mat_num %>%
    rownames_to_column("uniprot_id") %>%
    left_join(annot %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
    filter(!is.na(gene)) %>%
    pivot_longer(-c(uniprot_id, gene), names_to = "sample_id", values_to = "intensity") %>%
    mutate(group = group_vec[sample_id]) %>%
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

build_cv_wide <- function(cv_df) {
  cv_df %>%
    dplyr::select(gene, group, cv) %>%
    pivot_wider(names_from = group, values_from = cv) %>%
    mutate(
      cv_delta_prevctl  = CRC_PRE  - CTL,
      cv_delta_prevpost = CRC_POST - CRC_PRE
    )
}

cv_df_before <- build_cv_data(mat_raw,     annotation, group_vec)
cv_df_after  <- build_cv_data(mat_imputed, annotation, group_vec)

cv_wide_before <- build_cv_wide(cv_df_before)
cv_wide_after  <- build_cv_wide(cv_df_after)

write_xlsx(
  list(cv_before_imputation = cv_wide_before, cv_after_imputation = cv_wide_after),
  file.path(config$data_dir, "08_cv_before_after_imputation_plasma.xlsx")
)

compute_cv_summary_stats <- function(df, n_boot = 2000) {
  pooled_skew <- with(df, mean((cv - mean(cv))^3) / sd(cv)^3)
  use_median  <- abs(pooled_skew) > 1
  
  set.seed(42)
  df %>%
    group_by(group) %>%
    group_modify(~ {
      x <- .x$cv
      n <- length(x)
      if (use_median) {
        boot_medians <- replicate(n_boot, median(sample(x, n, replace = TRUE)))
        tibble(estimate = median(x),
               ci_lower = quantile(boot_medians, 0.025, names = FALSE),
               ci_upper = quantile(boot_medians, 0.975, names = FALSE),
               y_max = max(x), stat_type = "Median")
      } else {
        ci <- Hmisc::smean.cl.normal(x, conf.int = 0.95)
        tibble(estimate = ci[["Mean"]], ci_lower = ci[["Lower"]], ci_upper = ci[["Upper"]],
               y_max = max(x), stat_type = "Mean")
      }
    }) %>%
    ungroup()
}

make_cv_violin <- function(cv_df, title_label) {
  df <- cv_df %>% mutate(group = factor(group, levels = GROUPS))
  stats_df <- compute_cv_summary_stats(df) %>%
    mutate(group = factor(group, levels = GROUPS))
  
  y_range  <- diff(range(df$cv, na.rm = TRUE))
  stats_df <- stats_df %>%
    mutate(label_y = y_max + y_range * 0.08,
           label_txt = sprintf("%.1f [%.1f\u2013%.1f]", estimate, ci_lower, ci_upper))
  
  ggplot(df, aes(x = group, y = cv, fill = group, color = group)) +
    geom_violin(alpha = 0.5, linewidth = 0.5, trim = FALSE) +
    geom_boxplot(width = 0.15, alpha = 0.8, linewidth = 0.5,
                 outlier.shape = NA, colour = colors$gray_dark) +
    geom_label(data = stats_df, aes(x = group, y = label_y, label = label_txt),
               inherit.aes = FALSE, size = 2.8, lineheight = 0.9, fontface = "bold",
               color = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    scale_fill_manual(values = GROUP_COLOURS, guide = "none") +
    scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.2))) +
    labs(x = NULL, y = "CV (%)", title = title_label) +
    theme_minimal() +
    panel_border_theme
}

p_cv_violin_before <- make_cv_violin(cv_df_before, "Per-protein CV \u2014 Pre-Imputation")
p_cv_violin_after  <- make_cv_violin(cv_df_after,  "Per-protein CV \u2014 Post-Imputation")

ggsave(file.path(config$report_dir, "08_cv_violin_pre_imputation_plasma.pdf"),
       p_cv_violin_before, width = 6, height = 6, device = cairo_pdf)
ggsave(file.path(config$report_dir, "08_cv_violin_post_imputation_plasma.pdf"),
       p_cv_violin_after, width = 6, height = 6, device = cairo_pdf)

combined_cv_violin <- p_cv_violin_before + p_cv_violin_after
ggsave(file.path(config$report_dir, "08_cv_violin_before_v_after_imputation_plasma.pdf"),
       combined_cv_violin, width = 12, height = 6, device = cairo_pdf)
cat(" Saved: 08_cv_violin_before_v_after_imputation_plasma.pdf (+ individual panels)\n")

make_cv_scatter <- function(cv_wide, x_group, y_group, delta_col, n_label = 30, label = NULL) {
  df        <- cv_wide %>% filter(!is.na(.data[[x_group]]), !is.na(.data[[y_group]]))
  cor_test  <- cor.test(df[[x_group]], df[[y_group]])
  r_val     <- unname(cor_test$estimate)
  p_val     <- cor_test$p.value
  p_label   <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  delta_lim <- max(abs(df[[delta_col]]), na.rm = TRUE)
  df_label  <- df %>% slice_max(abs(.data[[delta_col]]), n = n_label)
  
  ggplot(df, aes(x = .data[[x_group]], y = .data[[y_group]], colour = .data[[delta_col]])) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", colour = "grey40", linewidth = 0.6) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", colour = colors$black, fill = "grey80", alpha = 0.15, linewidth = 0.8) +
    geom_text_repel(data = df_label, aes(label = gene), colour = "black", fontface = "bold",
                    size = 2.8, max.overlaps = Inf, box.padding = 0.4, point.padding = 0.2,
                    segment.color = "grey60", segment.size = 0.3, show.legend = FALSE, seed = 42) +
    scale_colour_gradientn(
      colours = c(colors$purple_dark, colors$purple, colors$orange, colors$yellow,
                  colors$orange, colors$purple, colors$purple_dark),
      values  = scales::rescale(c(-delta_lim, -delta_lim * 0.5, -delta_lim * 0.15, 0,
                                  delta_lim * 0.15, delta_lim * 0.5, delta_lim)),
      limits  = c(-delta_lim, delta_lim),
      name    = sprintf("\u0394CV\n(%s \u2212 %s)", y_group, x_group)
    ) +
    annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.1,
             label = sprintf("r = %.2f\n%s", r_val, p_label),
             size = 3.2, lineheight = 0.95, fontface = "bold",
             colour = colors$gray_dark, fill = colors$white, label.size = 0.2) +
    labs(x = bquote(CV[.(x_group)]~"(%)"), y = bquote(CV[.(y_group)]~"(%)"),
         title = sprintf("Per-protein CV: %s vs %s", x_group, y_group),
         subtitle = label) +
    theme_minimal() +
    panel_border_theme
}

cv_scatter_specs <- list(
  list(data = cv_wide_before, x = "CRC_PRE", y = "CTL",      delta = "cv_delta_prevctl",  label = "pre-imputation",  n = 30),
  list(data = cv_wide_before, x = "CRC_PRE", y = "CRC_POST", delta = "cv_delta_prevpost", label = "pre-imputation",  n = 30),
  list(data = cv_wide_after,  x = "CRC_PRE", y = "CTL",      delta = "cv_delta_prevctl",  label = "post-imputation", n = 30),
  list(data = cv_wide_after,  x = "CRC_PRE", y = "CRC_POST", delta = "cv_delta_prevpost", label = "post-imputation", n = 30)
)

cv_scatter_plots <- map(cv_scatter_specs, function(s) {
  make_cv_scatter(s$data, s$x, s$y, delta_col = s$delta, label = s$label, n_label = s$n)
})

walk2(cv_scatter_specs, cv_scatter_plots, function(s, p) {
  plot_name <- sprintf("09_cv_scatter_%s_v_%s_%s_plasma.pdf", s$x, s$y, gsub("-", "_", s$label))
  ggsave(file.path(config$report_dir, plot_name), p, width = 7, height = 6, device = cairo_pdf)
  cat(sprintf(" Saved: %s\n", plot_name))
})

# --- 6: Composite Report (all plots, one PDF) -------------------------------
cat("\n>> 6 - Composite report\n")

# One multi-page PDF bundling every plot above, grouped into logical
# sections/pages so each panel stays legible (cramming everything onto a
# single page gets illegible fast with this many panels).

page_missingness <- (p_protein_miss + p_sample_miss) / (p_condition_miss + p_class) +
  plot_annotation(title = "Imputation QC \u2014 Missingness Overview",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

page_distribution <- (p_density_overall / p_density_by_group) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "Imputation QC \u2014 Observed vs Imputed Distributions",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

page_pca <- combined_pca +
  plot_annotation(title = "Imputation QC \u2014 PCA Before vs After",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

page_cv_violin <- combined_cv_violin +
  plot_annotation(title = "Imputation QC \u2014 Per-Protein CV Before vs After",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

page_cv_scatter <- wrap_plots(cv_scatter_plots, ncol = 2) +
  plot_annotation(title = "Imputation QC \u2014 Per-Protein CV, PRE vs CTL / PRE vs POST",
                  theme = theme(plot.title = element_text(face = "bold", size = 15)))

composite_path <- file.path(config$report_dir, "00_imputation_composite_report_plasma.pdf")

cairo_pdf(composite_path, width = 14, height = 10, onefile = TRUE)
print(page_missingness)
print(page_distribution)
print(page_pca)
print(page_cv_violin)
print(page_cv_scatter)
dev.off()

cat(sprintf(" Saved: %s (5-page composite: missingness, distributions, PCA, CV violin, CV scatter)\n",
            basename(composite_path)))

# --- 7: Wrap-up --------------------------------------------------------------

writeLines(capture.output(sessionInfo()),
           file.path(config$data_dir, "sessionInfo_imputation_reports.txt"))

cat("\n>> Imputation reporting complete.\n")