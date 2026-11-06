# ============================================================
# Blood contamination by group — Muscle
# CR Proteomics
#
# Compares blood_pct (blood-marker % of total intensity) and B_M_ratio
# (blood : muscle marker ratio) across:
#   - cancer_time : CTL_T1 / CRC_T1 / CRC_T2
#   - supp_time   : PLA_T1 / CRE_T1 / PLA_T2 / CRE_T2
#
# Model: value ~ group + (1 | pid), Kenward-Roger df, then the SAME named
# contrasts as the limma DA script (same names, same "to - from" sign),
# so results line up directly with the DA contrasts.
#
# Both measures are fit on the log scale by default (ratios/percentages
# are right-skewed and bounded at 0). Contrasts are then reported as
# fold-ratios: 1.30 = 30% higher in the "to" group.
#
# Multiplicity: contrasts are a small, pre-specified set (not all
# pairwise), so p-values are unadjusted -- set config$p_adjust to
# "holm" if you'd rather correct within each measure x grouping.
# ============================================================

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, dplyr, tidyr, purrr, tibble, ggplot2, patchwork,
               lmerTest, emmeans, pbkrtest, writexl)

source(file.path("04_figures", "a_scripts", "00_theme.R"))  # colors, theme_cr()

# --- 1: Config -------------------------
config <- list(
  norm_RDS    = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS",
  report_dir  = "03_analysis/b_reports/muscle/blood_by_group",
  subject_col = "pid",
  measures = c(
    blood_pct = "Blood markers (% of total intensity)",
    B_M_ratio = "Blood : muscle marker ratio"
  ),
  log_transform = c(blood_pct = TRUE, B_M_ratio = TRUE),
  p_adjust      = "none"
)
dir.create(config$report_dir, recursive = TRUE, showWarnings = FALSE)

# Contrast weights are named by group level, written as "to - from" to
# match the limma contrast_specs exactly.
groupings <- list(
  cancer = list(
    group_col = "cancer_time",
    levels    = c("CTL_T1", "CRC_T1", "CRC_T2"),
    contrasts = list(
      Cancer   = c(CTL_T1 = -1, CRC_T1 =  1, CRC_T2 = 0),  # CRC_T1 - CTL_T1
      Residual = c(CTL_T1 = -1, CRC_T1 =  0, CRC_T2 = 1),  # CRC_T2 - CTL_T1
      Training = c(CTL_T1 =  0, CRC_T1 = -1, CRC_T2 = 1)   # CRC_T2 - CRC_T1
    )
  ),
  supp = list(
    group_col = "supp_time",
    levels    = c("PLA_T1", "CRE_T1", "PLA_T2", "CRE_T2"),
    contrasts = list(
      Baseline_supp      = c(PLA_T1 = -1, CRE_T1 =  1, PLA_T2 =  0, CRE_T2 = 0),
      Post_training_supp = c(PLA_T1 =  0, CRE_T1 =  0, PLA_T2 = -1, CRE_T2 = 1),
      Training_CRE       = c(PLA_T1 =  0, CRE_T1 = -1, PLA_T2 =  0, CRE_T2 = 1),
      Training_PLA       = c(PLA_T1 = -1, CRE_T1 =  0, PLA_T2 =  1, CRE_T2 = 0),
      Interaction_supp   = c(PLA_T1 =  1, CRE_T1 = -1, PLA_T2 = -1, CRE_T2 = 1)
    )
  )
)

# --- 2: Load & pre-flight checks -------
dal  <- readRDS(config$norm_RDS)
meta <- as_tibble(dal$metadata)

needed  <- c(config$subject_col, names(config$measures),
             map_chr(groupings, "group_col"))
missing <- setdiff(needed, colnames(meta))
if (length(missing) > 0) {
  cat("Metadata columns available:\n"); print(colnames(meta))
  stop("Missing metadata column(s): ", paste(missing, collapse = ", "))
}

for (gn in names(groupings)) {
  g <- groupings[[gn]]
  cat(sprintf("\n>> %s levels observed:\n", g$group_col))
  print(table(meta[[g$group_col]], useNA = "ifany"))
  absent <- setdiff(g$levels, unique(meta[[g$group_col]]))
  if (length(absent) > 0) {
    stop(g$group_col, " is missing expected level(s): ", paste(absent, collapse = ", "))
  }
}

for (m in names(config$measures)) {
  if (config$log_transform[[m]] && any(meta[[m]] <= 0, na.rm = TRUE)) {
    stop(m, " has values <= 0 -- can't log-transform. Set log_transform[['", m, "']] = FALSE.")
  }
}

# --- 3: Data prep + model --------------

prep_measure_data <- function(meta, measure, g) {
  meta %>%
    dplyr::filter(.data[[g$group_col]] %in% g$levels, !is.na(.data[[measure]])) %>%
    dplyr::transmute(
      subject = .data[[config$subject_col]],
      group   = factor(.data[[g$group_col]], levels = g$levels),
      value   = .data[[measure]],
      y       = if (config$log_transform[[measure]]) log(value) else value
    )
}

fit_measure <- function(grouping, measure) {
  g      <- groupings[[grouping]]
  log_it <- config$log_transform[[measure]]
  d      <- prep_measure_data(meta, measure, g)
  
  fit <- lmerTest::lmer(y ~ group + (1 | subject), data = d)
  
  descriptives <- d %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      n          = dplyr::n(),
      n_subjects = dplyr::n_distinct(subject),
      mean       = mean(value),
      sd         = sd(value),
      median     = median(value),
      q1         = quantile(value, 0.25),
      q3         = quantile(value, 0.75),
      geo_mean   = if (log_it) exp(mean(log(value))) else NA_real_,
      .groups    = "drop"
    )
  
  anova_tab <- anova(fit, ddf = "Kenward-Roger") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")
  
  emm <- emmeans::emmeans(fit, ~ group, lmer.df = "kenward-roger")
  
  emm_df <- summary(emm) %>% as.data.frame()
  if (log_it) {
    emm_df <- emm_df %>%
      dplyr::mutate(response = exp(emmean),
                    response_lcl = exp(lower.CL), response_ucl = exp(upper.CL))
  }
  
  con_df <- emmeans::contrast(emm, method = purrr::map(g$contrasts, ~ .x[g$levels]),
                              adjust = config$p_adjust) %>%
    summary(infer = TRUE) %>%
    as.data.frame()
  if (log_it) {
    con_df <- con_df %>%
      dplyr::mutate(ratio = exp(estimate),
                    ratio_lcl = exp(lower.CL), ratio_ucl = exp(upper.CL))
  }
  
  list(grouping = grouping, measure = measure, log_it = log_it,
       data = d, descriptives = descriptives, anova = anova_tab,
       emmeans = emm_df, contrasts = con_df,
       singular = lme4::isSingular(fit))
}

runs    <- tidyr::expand_grid(grouping = names(groupings), measure = names(config$measures))
results <- purrr::pmap(runs, fit_measure)
names(results) <- paste(runs$grouping, runs$measure, sep = "_")

singular_fits <- names(results)[purrr::map_lgl(results, "singular")]
if (length(singular_fits) > 0) {
  warning("Singular fit (subject variance ~ 0) in: ", paste(singular_fits, collapse = ", "),
          ". Estimates are still valid but the random intercept isn't doing anything there.")
}

# --- 4: Plots --------------------------

contrast_caption <- function(res) {
  res$contrasts %>%
    dplyr::mutate(txt = if (res$log_it) {
      sprintf("%s: %.2f\u00d7 [%.2f, %.2f], p = %.3f",
              contrast, ratio, ratio_lcl, ratio_ucl, p.value)
    } else {
      sprintf("%s: %+.2f [%.2f, %.2f], p = %.3f",
              contrast, estimate, lower.CL, upper.CL, p.value)
    }) %>%
    dplyr::pull(txt) %>%
    paste(collapse = "\n")
}

make_blood_plot <- function(res) {
  y_lab <- config$measures[[res$measure]]
  if (res$log_it) y_lab <- paste0(y_lab, " (log scale)")
  
  p <- ggplot(res$data, aes(x = group, y = value)) +
    geom_boxplot(outlier.shape = NA, width = 0.55, fill = "grey90") +
    # paired lines connect each subject's T1 and T2; CTL subjects have one
    # sample so draw no line
    geom_line(aes(group = subject), color = "grey60", linewidth = 0.3, alpha = 0.7) +
    geom_point(size = 1.6, alpha = 0.8, color = colors$primary) +
    labs(x = NULL, y = y_lab,
         title = paste0(config$measures[[res$measure]], " — ", res$grouping),
         caption = contrast_caption(res)) +
    theme_cr() +
    theme(plot.caption = element_text(hjust = 0, size = 8))
  
  if (res$log_it) p <- p + scale_y_log10()
  p
}

plots <- purrr::map(results, make_blood_plot)

combined <- wrap_plots(plots, ncol = length(groupings), byrow = FALSE)
cairo_pdf(file.path(config$report_dir, "blood_by_group_muscle.pdf"),
          width = 6 * length(groupings), height = 5.5 * length(config$measures))
print(combined)
dev.off()

# --- 5: Export -------------------------

tag <- function(df, res) {
  dplyr::mutate(df, grouping = res$grouping, measure = res$measure, .before = 1)
}

write_xlsx(
  list(
    descriptives = purrr::map_dfr(results, ~ tag(.x$descriptives, .x)),
    anova        = purrr::map_dfr(results, ~ tag(.x$anova, .x)),
    emmeans      = purrr::map_dfr(results, ~ tag(.x$emmeans, .x)),
    contrasts    = purrr::map_dfr(results, ~ tag(.x$contrasts, .x))
  ),
  file.path(config$report_dir, "blood_by_group_muscle.xlsx")
)

cat("\n>> Contrasts:\n")
print(purrr::map_dfr(results, ~ tag(.x$contrasts, .x)) %>%
        dplyr::select(grouping, measure, contrast, dplyr::any_of(c("ratio", "estimate")),
                      p.value))
