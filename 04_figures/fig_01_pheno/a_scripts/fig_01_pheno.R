# Fig Phenotypic Overview -- DXA/phenotypic variables, CRE vs PLA (whole CRC cohort)
# Cohort = union of participants in the muscle AND plasma metadata files
# (muscle-only, plasma-only, and both), collapsed to one row per subject.
#
# Two analysis frameworks, each producing a PRE/POST figure and a change-score figure:
#
#   "ANOVA"  -- linear mixed model: value ~ group * time + (1|subject), Kenward-Roger df.
#               All available observations are used (subjects missing a timepoint still
#               contribute). Questions answered:
#                 Time (pooled POST - PRE)          -> change during RET
#                 Group x Time (difference in change) -> supplement effect
#               The group main effect (CRE vs PLA averaged over PRE and POST) is plotted
#               alongside Time and Group x Time. It averages a randomized baseline
#               difference with half the interaction, so it is not a supplement effect.
#
#   "ANCOVA" -- POST ~ group + PRE_c (PRE centered at the outcome's cohort mean), one row
#               per subject, paired subjects only. Questions answered:
#                 Adjusted CRE - PLA difference at POST -> supplement effect
#                 Adjusted pooled change               -> change during RET
#               The change-score model (POST - PRE) ~ group + PRE_c gives an identical
#               group difference; it is used only to put the adjusted means on the
#               change scale for the change figure.
#
# Percent change: descriptive only. Individual 100 * (POST - PRE) / PRE is added to the
# exported data and summarised (mean +/- SD) in the manuscript table for reporting.
# All inference is on absolute values.
#
# Outputs (04_figures/fig_01_pheno/b_figures):
#   fig_anova_prepost   -- model-estimated PRE/POST means by group (95% CI)
#   fig_anova_change    -- model-estimated POST - PRE change by group (95% CI)
#   fig_ancova_prepost  -- observed PRE means, baseline-adjusted POST means (95% CI)
#   fig_ancova_change   -- baseline-adjusted POST - PRE change by group (95% CI)
#   fig_anova_pooled_prepost -- mixed-model PRE/POST means pooled over CRE and PLA
#                               (95% CI), Time p-value bracket
# Diagnostics (04_figures/fig_01_pheno/d_diagnostics):
#   diag_qq_change, diag_qq_resid_lmm, diag_qq_resid_ancova -- normal Q-Q plots
# Plus an xlsx with the manuscript table, all model output, and normality tables.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               lme4, lmerTest, pbkrtest, emmeans, scales, writexl, ggtext)

# --- Configuration -----------------------------------------------------------
config <- list(
  data_files = c(
    muscle = "00_input/muscle/CRm_meta.csv",
    plasma = "00_input/plasma/CRp_meta.csv"    # <-- CHECK: adjust to your plasma metadata path
  ),
  fig_dir    = "04_figures/fig_01_pheno/b_figures",
  data_dir   = "04_figures/fig_01_pheno/c_data",
  diag_dir   = "04_figures/fig_01_pheno/d_diagnostics",
  
  subject_col   = "pid",
  group_col     = "supp",
  cancer_col    = "cancer",
  cancer_filter = "SURV",
  
  variables = c(
    weight_kg       = "Weight (kg)",
    BMI             = "BMI",
    LBM_kg          = "Lean Body Mass (kg)",
    ALM_kg          = "Appendicular Lean Mass (kg)",
    ALM_index       = "ALM Index",
    fat_kg          = "Fat Mass (kg)",
    fat_pct         = "Fat (%)",
    sts_max_pwr     = "Sit-to-Stand Max Power",
    grip_lbs        = "Grip Strength (lbs)",
    chest_press_lbs = "Chest Press (lbs)",
    leg_ext_lbs     = "Leg Extension (lbs)"
  ),
  
  prefix_pre  = "pre_",
  prefix_post = "post_"
)

# A calling script (e.g. fig_01_full_cohort.R) can source this file with a
# `config_override` list in the sourcing environment to swap data files and output
# folders without editing anything here.
if (exists("config_override", inherits = FALSE)) config <- modifyList(config, config_override)

dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$diag_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

group_colors <- c(CRE = colors$orange, PLA = colors$blue_dark)
pooled_color <- colors$gray_dark   # bars for the pooled (CRE + PLA) figure

# --- 1: Load & reshape data ---------------------------------------------------
expected_cols <- as.vector(outer(
  c(config$prefix_pre, config$prefix_post),
  names(config$variables),
  paste0
))
id_cols       <- c(config$subject_col, config$group_col, config$cancer_col)
required_cols <- c(id_cols, expected_cols)

read_meta <- function(path) {
  if (grepl("\\.xlsx?$", path)) read_excel(path) else read_csv(path, show_col_types = FALSE)
}

meta_list <- lapply(config$data_files, read_meta)

# --- Diagnostic check --------------------------------------------------------
cat("=== Raw data structure (per file) ===\n")
for (tis in names(meta_list)) {
  df <- meta_list[[tis]]
  cat("\n[", tis, "] ", config$data_files[[tis]], "\n", sep = "")
  cat("  Columns:", paste(colnames(df), collapse = ", "), "\n")
  cat("  N rows:", nrow(df), "\n")
  missing_req <- setdiff(required_cols, colnames(df))
  if (length(missing_req) > 0) {
    stop("[", tis, "] file is missing required columns: ",
         paste(missing_req, collapse = ", "),
         "\nCheck config$*_col / config$variables / prefixes. If the ", tis,
         " metadata doesn't carry the pre_/post_ phenotypic columns, point ",
         "config$data_files at a master phenotype file instead.")
  }
  cat("  N unique subjects:", n_distinct(df[[config$subject_col]]), "\n")
  cat("  Rows per subject (expect mostly 1-2):\n")
  print(table(table(df[[config$subject_col]])))
}

meta_all <- imap_dfr(meta_list, function(df, tis) {
  df %>%
    dplyr::select(all_of(required_cols)) %>%
    mutate(
      across(all_of(id_cols), as.character),
      across(all_of(expected_cols), as.numeric),
      tissue = tis
    )
})

conflict_tbl <- meta_all %>%
  group_by(across(all_of(config$subject_col))) %>%
  summarise(across(all_of(c(config$group_col, config$cancer_col, expected_cols)),
                   ~ n_distinct(.x, na.rm = TRUE)), .groups = "drop") %>%
  pivot_longer(-all_of(config$subject_col), names_to = "column", values_to = "n_distinct") %>%
  filter(n_distinct > 1)
cat("\nSubject x column conflicts across files/rows:", nrow(conflict_tbl), "\n")
if (nrow(conflict_tbl) > 0) {
  warning("Conflicting values found for some subjects -- first non-NA (muscle file priority) was kept. ",
          "Inspect `conflict_tbl`.")
  print(conflict_tbl, n = 50)
}

avail <- meta_all %>%
  group_by(across(all_of(config$subject_col))) %>%
  summarise(has_muscle = "muscle" %in% tissue,
            has_plasma = "plasma" %in% tissue, .groups = "drop") %>%
  mutate(tissue_avail = case_when(
    has_muscle & has_plasma ~ "both",
    has_muscle              ~ "muscle only",
    has_plasma              ~ "plasma only",
    TRUE                    ~ "n/a"
  ))

first_non_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x)) x[1] else x[NA_integer_]
}

raw_subj <- meta_all %>%
  group_by(across(all_of(config$subject_col))) %>%
  summarise(across(all_of(c(config$group_col, config$cancer_col, expected_cols)), first_non_na),
            .groups = "drop") %>%
  left_join(avail, by = config$subject_col)

cat("\nN subjects in union of files:", nrow(raw_subj), "\n")

pheno_subj <- raw_subj %>%
  filter(.data[[config$cancer_col]] == config$cancer_filter,
         .data[[config$group_col]] %in% c("CRE", "PLA"))

cat("\nGroup counts (", config$cancer_filter, ", CRE/PLA only):\n", sep = "")
print(table(pheno_subj[[config$group_col]]))
cat("\nTissue availability by group:\n")
print(table(pheno_subj$tissue_avail, pheno_subj[[config$group_col]]))

cohort_availability <- pheno_subj %>%
  count(group = .data[[config$group_col]], tissue_avail) %>%
  pivot_wider(names_from = tissue_avail, values_from = n, values_fill = 0) %>%
  mutate(total = rowSums(across(-group)))

cohort_roster <- pheno_subj %>%
  dplyr::select(subject = all_of(config$subject_col),
                group   = all_of(config$group_col),
                tissue_avail)

var_label_lookup <- tibble(
  variable       = names(config$variables),
  variable_label = factor(unname(config$variables), levels = unname(config$variables))
)

# Long format: subject x variable x timepoint (used by the LMM)
pheno_long <- pheno_subj %>%
  dplyr::select(subject = all_of(config$subject_col), group = all_of(config$group_col),
                tissue_avail, all_of(expected_cols)) %>%
  pivot_longer(cols = all_of(expected_cols), names_to = "colname", values_to = "value") %>%
  mutate(
    time = case_when(
      str_starts(colname, config$prefix_pre)  ~ "PRE",
      str_starts(colname, config$prefix_post) ~ "POST"
    ),
    variable = str_remove(colname, paste0("^", config$prefix_pre, "|^", config$prefix_post)),
    group    = factor(group, levels = c("PLA", "CRE")),
    time     = factor(time,  levels = c("PRE", "POST"))
  ) %>%
  dplyr::select(-colname) %>%
  left_join(var_label_lookup, by = "variable")

# Wide format: one row per subject x variable, paired subjects only (used by the ANCOVA).
# PRE_c is centered at each outcome's mean PRE, so adjusted means are evaluated at
# the cohort's mean baseline. PCT (individual % change) is descriptive only -- it is
# not used in any model.
ancova_data <- pheno_long %>%
  dplyr::select(subject, group, variable, variable_label, time, value) %>%
  pivot_wider(names_from = time, values_from = value) %>%
  filter(!is.na(PRE), !is.na(POST)) %>%
  group_by(variable) %>%
  mutate(PRE_c = PRE - mean(PRE),
         CHG   = POST - PRE,
         PCT   = if_else(PRE > 0, 100 * (POST - PRE) / PRE, NA_real_)) %>%
  ungroup()

cat("\n=== Long data check ===\n")
cat("Rows:", nrow(pheno_long), "\n")
print(pheno_long %>% group_by(variable_label) %>% summarise(n_na = sum(is.na(value)), .groups = "drop"))
cat("\n=== Paired subjects per outcome (ANCOVA n) ===\n")
print(ancova_data %>% count(variable_label, group) %>%
        pivot_wider(names_from = group, values_from = n, values_fill = 0))

# --- 1b: Normality assessment (pre-model) ---------------------------------------
# Raw values by group x time and change scores by group. These screen for skew,
# outliers, and data-entry problems before modeling. The models themselves assume
# normal RESIDUALS, which are checked after fitting (section 3b).
#
# With n ~5-10 per cell, Shapiro-Wilk has little power to detect non-normality and a
# non-significant p is not evidence of normality. Read the table together with the
# Q-Q plots (section 7b). Flag = "check" if SW p < 0.05 or |skew| > 1 (heuristics).
skewness_ <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 3) return(NA_real_)
  m <- mean(x); s <- sqrt(mean((x - m)^2))
  if (s == 0) return(NA_real_)
  mean((x - m)^3) / s^3
}

kurtosis_excess <- function(x) {
  x <- x[!is.na(x)]; n <- length(x)
  if (n < 4) return(NA_real_)
  m <- mean(x); s <- sqrt(mean((x - m)^2))
  if (s == 0) return(NA_real_)
  mean((x - m)^4) / s^4 - 3
}

sw_test <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 3 || length(unique(x)) < 3) return(c(W = NA_real_, p = NA_real_))
  t <- shapiro.test(x)
  c(W = unname(t$statistic), p = t$p.value)
}

# Tukey 1.5 x IQR rule (z-scores are uninformative at these n)
n_iqr_outliers <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 4) return(NA_integer_)
  q <- quantile(x, c(0.25, 0.75)); iqr <- diff(q)
  as.integer(sum(x < q[1] - 1.5 * iqr | x > q[2] + 1.5 * iqr))
}

normality_summary <- function(df, value_col, ...) {
  df %>%
    group_by(...) %>%
    summarise(
      n           = sum(!is.na(.data[[value_col]])),
      skew        = skewness_(.data[[value_col]]),
      kurt_excess = kurtosis_excess(.data[[value_col]]),
      sw_W        = sw_test(.data[[value_col]])[["W"]],
      sw_p        = sw_test(.data[[value_col]])[["p"]],
      n_outliers  = n_iqr_outliers(.data[[value_col]]),
      .groups = "drop"
    ) %>%
    mutate(flag = case_when(
      is.na(sw_p)                    ~ "too few",
      sw_p < 0.05 | abs(skew) > 1    ~ "check",
      TRUE                           ~ "ok"
    ))
}

norm_raw    <- normality_summary(pheno_long,  "value", variable, variable_label, group, time)
norm_change <- normality_summary(ancova_data, "CHG",   variable, variable_label, group)

cat("\n=== Normality: raw values (group x time) -- flagged ===\n")
print(norm_raw %>% filter(flag != "ok") %>% dplyr::select(-variable), n = Inf)
cat("\n=== Normality: change scores (by group) -- flagged ===\n")
print(norm_change %>% filter(flag != "ok") %>% dplyr::select(-variable), n = Inf)

# --- Shared helpers ------------------------------------------------------------
# Tidy an emmGrid (means or contrasts) into a data frame with a common
# `estimate` column and 95% CIs.
tidy_emm <- function(emm) {
  s <- as.data.frame(summary(emm, infer = c(TRUE, TRUE)))
  if ("emmean" %in% names(s)) s <- dplyr::rename(s, estimate = emmean)
  s
}

fmt_p <- function(p) {
  case_when(
    is.na(p)  ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE      ~ sprintf("%.3f", p)
  )
}

# Bold the line when raw p < 0.05 (pre-specified battery, uncorrected)
fmt_line <- function(name, p) {
  text <- paste0(name, ifelse(!is.na(p) & p < 0.001, " p < 0.001", paste0(" p = ", fmt_p(p))))
  ifelse(!is.na(p) & p < 0.05, paste0("**", text, "**"), text)
}

# --- 2: "ANOVA" -- mixed model per variable -----------------------------------
DF_METHOD <- "Kenward-Roger"
emm_df    <- tolower(DF_METHOD)
emmeans::emm_options(msg.interaction = FALSE)

fit_variable_model <- function(var_name) {
  df <- pheno_long %>% filter(variable == var_name, !is.na(value)) %>% droplevels()
  tryCatch(
    lmerTest::lmer(value ~ group * time + (1 | subject), data = df),
    error = function(e) e
  )
}

fitted_models <- setNames(lapply(names(config$variables), fit_variable_model),
                          names(config$variables))

# Omnibus F-tests: exported for the supplement only. In a 2x2 each term has 1 df,
# so time and group:time are the same tests as the plotted contrasts (F = t^2).
extract_anova_table <- function(m) {
  if (inherits(m, "error")) {
    return(tibble(term = c("group", "time", "group:time"),
                  F_value = NA_real_, df_num = NA_real_, df_den = NA_real_,
                  p.value = NA_real_, singular_fit = NA))
  }
  a <- anova(m, ddf = DF_METHOD)
  tibble(term = rownames(a), F_value = a$`F value`, df_num = a$NumDF,
         df_den = a$DenDF, p.value = a$`Pr(>F)`,
         singular_fit = lme4::isSingular(m, tol = 1e-4))
}

model_results <- imap_dfr(fitted_models, function(m, var_name) {
  extract_anova_table(m) %>% mutate(variable = var_name, .before = 1)
}) %>%
  left_join(var_label_lookup, by = "variable") %>%
  relocate(variable_label, .after = variable)

if (any(model_results$singular_fit, na.rm = TRUE)) {
  cat("\nNote: singular fit flagged for:",
      paste(unique(model_results$variable_label[which(model_results$singular_fit)]), collapse = ", "),
      "-- random-intercept variance estimated near zero.\n")
}
failed <- names(fitted_models)[vapply(fitted_models, inherits, logical(1), "error")]
if (length(failed)) cat("\nLMM failed for:", paste(failed, collapse = ", "), "\n")

# Contrasts from the same fitted model (positive = CRE > PLA, POST > PRE)
run_anova_effects <- function(m) {
  if (inherits(m, "error")) return(NULL)
  emm_gt <- emmeans::emmeans(m, ~ group * time, lmer.df = emm_df)
  emm_g  <- emmeans::emmeans(m, ~ group,        lmer.df = emm_df)
  emm_t  <- emmeans::emmeans(m, ~ time,         lmer.df = emm_df)
  
  list(
    # Cell means for the PRE/POST figure
    means = tidy_emm(emm_gt),
    # Within-group POST - PRE for the change figure
    change = tidy_emm(emmeans::contrast(emm_gt, "revpairwise", by = "group", adjust = "none")),
    # Pooled PRE and POST means (equal weight over CRE and PLA) for the pooled figure
    time_means = tidy_emm(emm_t),
    # Pooled POST - PRE (equal weight over CRE and PLA) = time effect
    pooled_change = tidy_emm(emmeans::contrast(emm_t, "revpairwise", adjust = "none")),
    # (POST - PRE in CRE) - (POST - PRE in PLA) = Group x Time
    diff_in_change = tidy_emm(emmeans::contrast(emm_gt, interaction = "revpairwise", adjust = "none")),
    # Exported, not plotted
    group_at_time = tidy_emm(emmeans::contrast(emm_gt, "revpairwise", by = "time", adjust = "none")),
    group_main    = tidy_emm(emmeans::contrast(emm_g,  "revpairwise", adjust = "none")),
    n_obs      = as.integer(stats::nobs(m)),
    n_subjects = as.integer(unname(lme4::ngrps(m)["subject"]))
  )
}

anova_fx <- lapply(fitted_models, run_anova_effects)

# Stack one component across variables
collect_fx <- function(fx_list, which_tbl) {
  imap_dfr(fx_list, function(x, var_name) {
    if (is.null(x) || is.null(x[[which_tbl]])) return(NULL)
    as_tibble(x[[which_tbl]]) %>% mutate(variable = var_name, .before = 1)
  }) %>%
    left_join(var_label_lookup, by = "variable") %>%
    relocate(variable_label, .after = variable)
}

anova_means          <- collect_fx(anova_fx, "means")
anova_time_means     <- collect_fx(anova_fx, "time_means")
anova_change         <- collect_fx(anova_fx, "change")
anova_pooled_change  <- collect_fx(anova_fx, "pooled_change")
anova_diff_in_change <- collect_fx(anova_fx, "diff_in_change")
anova_group_at_time  <- collect_fx(anova_fx, "group_at_time")
anova_group_main     <- collect_fx(anova_fx, "group_main")

# --- 3: ANCOVA per variable ----------------------------------------------------
run_ancova <- function(var_name) {
  df <- ancova_data %>% filter(variable == var_name)
  n_by <- table(df$group)
  if (any(n_by < 2)) {
    cat("ANCOVA skipped for", var_name, "-- fewer than 2 paired subjects in a group\n")
    return(NULL)
  }
  m_post <- lm(POST ~ group + PRE_c, data = df)
  m_chg  <- lm(CHG  ~ group + PRE_c, data = df)
  emm_post <- emmeans::emmeans(m_post, ~ group)   # evaluated at PRE_c = 0 (mean baseline)
  emm_chg  <- emmeans::emmeans(m_chg,  ~ group)
  
  list(
    post_means    = tidy_emm(emm_post),
    change        = tidy_emm(emm_chg),
    # Adjusted CRE - PLA (identical whether taken from the POST or the change model)
    group_diff    = tidy_emm(emmeans::contrast(emm_post, "revpairwise", adjust = "none")),
    # Adjusted change averaged with equal weight over groups, tested against 0
    pooled_change = tidy_emm(emmeans::emmeans(m_chg, ~ 1)),
    coefs = tibble(n_CRE = as.integer(n_by[["CRE"]]), n_PLA = as.integer(n_by[["PLA"]]),
                   b_PRE = unname(coef(m_post)[["PRE_c"]]),
                   r2    = summary(m_post)$r.squared),
    resid = tibble(subject = df$subject, group = df$group,
                   resid = unname(resid(m_post)))
  )
}

ancova_fx <- setNames(lapply(names(config$variables), run_ancova), names(config$variables))

ancova_post_means    <- collect_fx(ancova_fx, "post_means")
ancova_change        <- collect_fx(ancova_fx, "change")
ancova_group_diff    <- collect_fx(ancova_fx, "group_diff")
ancova_pooled_change <- collect_fx(ancova_fx, "pooled_change")
ancova_coefs         <- collect_fx(ancova_fx, "coefs")
ancova_resid         <- collect_fx(ancova_fx, "resid")

# --- 3b: Normality of model residuals ----------------------------------------------
# This is the assumption that actually matters for the LMM and ANCOVA inference.
# Residuals are pooled across groups within each outcome (one residual distribution
# is assumed per model). LMM residuals are conditional (observation-level).
lmm_resid <- imap_dfr(fitted_models, function(m, var_name) {
  if (inherits(m, "error")) return(NULL)
  mf <- model.frame(m)
  tibble(variable = var_name, group = mf$group, time = mf$time,
         resid = unname(resid(m)))
}) %>%
  left_join(var_label_lookup, by = "variable") %>%
  relocate(variable_label, .after = variable)

norm_resid <- bind_rows(
  normality_summary(lmm_resid,    "resid", variable, variable_label) %>% mutate(model = "LMM",    .before = 1),
  normality_summary(ancova_resid, "resid", variable, variable_label) %>% mutate(model = "ANCOVA", .before = 1)
)

cat("\n=== Normality: model residuals ===\n")
print(norm_resid %>% dplyr::select(-variable), n = Inf)

# PRE for the ANCOVA figure = pooled cohort baseline (mean PRE of the paired
# participants, t-based 95% CI), shown identically for both groups. This is the
# baseline value the adjusted POST means are evaluated at, so the PRE -> POST line
# for each group is exactly that group's adjusted change (matches the change figure).
# Group-specific observed PRE means would mix two reference points and misrepresent
# the adjusted change; they are in the Main_Table instead.
ancova_pre_means <- ancova_data %>%
  group_by(variable, variable_label) %>%
  summarise(n = n(), m = mean(PRE), s = sd(PRE), .groups = "drop") %>%
  mutate(estimate = m,
         lower.CL = m - qt(0.975, n - 1) * s / sqrt(n),
         upper.CL = m + qt(0.975, n - 1) * s / sqrt(n),
         time     = "PRE") %>%
  tidyr::crossing(group = c("PLA", "CRE")) %>%
  dplyr::select(variable, variable_label, group, time, estimate, lower.CL, upper.CL)

ancova_prepost <- bind_rows(
  ancova_pre_means,
  ancova_post_means %>% mutate(time = "POST") %>%
    dplyr::select(variable, variable_label, group, time, estimate, lower.CL, upper.CL)
) %>%
  mutate(group = factor(group, levels = c("PLA", "CRE")),
         time  = factor(time,  levels = c("PRE", "POST")))

cat("\n=== ANOVA: pooled change (time) ===\n");          print(anova_pooled_change, n = Inf)
cat("\n=== ANOVA: difference in change (G x T) ===\n");  print(anova_diff_in_change, n = Inf)
cat("\n=== ANCOVA: adjusted CRE - PLA ===\n");           print(ancova_group_diff, n = Inf)
cat("\n=== ANCOVA: pooled adjusted change ===\n");       print(ancova_pooled_change, n = Inf)

# --- 4: Annotations --------------------------------------------------------------
anova_p <- var_label_lookup %>%
  left_join(anova_group_main     %>% dplyr::select(variable, p_group = p.value), by = "variable") %>%
  left_join(anova_pooled_change  %>% dplyr::select(variable, p_time  = p.value), by = "variable") %>%
  left_join(anova_diff_in_change %>% dplyr::select(variable, p_int   = p.value), by = "variable")

ancova_p <- var_label_lookup %>%
  left_join(ancova_pooled_change %>% dplyr::select(variable, p_time  = p.value), by = "variable") %>%
  left_join(ancova_group_diff    %>% dplyr::select(variable, p_group = p.value), by = "variable")

annot_anova_prepost <- anova_p %>%
  mutate(label = paste0(fmt_line("Group", p_group), "<br>",
                        fmt_line("Time", p_time), "<br>",
                        fmt_line("Group\u00d7Time", p_int)))
annot_anova_change <- anova_p %>%
  mutate(label = paste0(fmt_line("Group (main)", p_group), "<br>",
                        fmt_line("Pooled change", p_time), "<br>",
                        fmt_line("CRE vs PLA", p_int)))
annot_ancova_prepost <- ancova_p %>%
  mutate(label = fmt_line("CRE vs PLA (adj. POST)", p_group))
annot_ancova_change <- ancova_p %>%
  mutate(label = paste0(fmt_line("Pooled adj. change", p_time), "<br>", fmt_line("CRE vs PLA", p_group)))

# --- 5: Plot functions -------------------------------------------------------------
facet_theme <- function() {
  theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      strip.background   = element_rect(fill = colors$gray_light, color = NA),
      panel.grid.major.y = element_line(color = colors$gray_light),
      panel.grid.minor   = element_blank()
    )
}

add_annotation <- function(annot) {
  geom_richtext(
    data = annot, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    size = 2.3, hjust = -0.05, vjust = 1.1, lineheight = 1.1,
    fill = NA, label.color = NA, label.padding = unit(0, "lines")
  )
}

# PRE/POST: bars = group estimate (95% CI), clustered by group with PRE and POST
# side by side (PRE lighter). Individual participants are jittered points joined
# PRE -> POST; each subject gets one fixed jitter offset so its line stays inside
# its own group's cluster rather than crossing the other group's bars.
prepost_pos <- tibble(
  group = c("PLA", "PLA", "CRE", "CRE"),
  time  = c("PRE", "POST", "PRE", "POST"),
  x     = c(1, 1.8, 2.8, 3.6)
)

make_prepost_plot <- function(est, raw, annot, title, subtitle,
                              bar_width = 0.7, jitter_width = 0.12) {
  est <- est %>%
    filter(!is.na(estimate)) %>%
    mutate(group = as.character(group), time = as.character(time)) %>%
    left_join(prepost_pos, by = c("group", "time")) %>%
    mutate(group = factor(group, levels = c("PLA", "CRE")),
           time  = factor(time,  levels = c("PRE", "POST")))
  
  set.seed(1)
  subj_jit <- raw %>% distinct(subject) %>%
    mutate(jit = runif(n(), -jitter_width, jitter_width))
  
  raw <- raw %>%
    mutate(group = as.character(group), time = as.character(time)) %>%
    left_join(prepost_pos, by = c("group", "time")) %>%
    left_join(subj_jit, by = "subject") %>%
    mutate(x = x + jit, group = factor(group, levels = c("PLA", "CRE")))
  
  ggplot(est, aes(x = x, y = estimate)) +
    geom_col(aes(fill = group, alpha = time), width = bar_width, color = colors$black) +
    geom_line(data = raw, aes(x = x, y = value, group = subject),
              inherit.aes = FALSE, color = colors$gray_dark, alpha = 0.35, linewidth = 0.3) +
    geom_point(data = raw, aes(x = x, y = value), inherit.aes = FALSE,
               shape = 21, size = 1.4, alpha = 0.85, fill = "white", color = colors$gray_dark) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.5) +
    add_annotation(annot) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_alpha_manual(values = c(PRE = 0.4, POST = 1), name = "Time") +
    scale_x_continuous(breaks = prepost_pos$x,
                       labels = paste0(prepost_pos$time, "\n", prepost_pos$group)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
    facet_wrap(~ variable_label, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    facet_theme() +
    theme(axis.text.x = element_text(size = 7, lineheight = 0.9))
}

# Change score: group estimate (95% CI) as bars from zero, observed individual
# POST - PRE changes jittered on top.
make_change_plot <- function(est, raw, annot, title, subtitle) {
  est <- est %>% filter(!is.na(estimate)) %>%
    mutate(group = factor(group, levels = c("PLA", "CRE")))
  
  ggplot(est, aes(x = group, y = estimate, fill = group)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = colors$gray_dark, linewidth = 0.4) +
    geom_col(width = 0.6, color = colors$black) +
    geom_point(data = raw, aes(x = group, y = CHG), inherit.aes = FALSE,
               position = position_jitter(width = 0.1, height = 0, seed = 1),
               shape = 21, size = 1.4, alpha = 0.85, fill = "white", color = colors$gray_dark) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.5) +
    add_annotation(annot) +
    scale_fill_manual(values = group_colors, name = "Group") +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.35))) +
    facet_wrap(~ variable_label, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Change (POST \u2212 PRE)") +
    facet_theme() +
    theme(legend.position = "none")
}

# Pooled PRE/POST: one bar per timepoint (model-estimated mean, 95% CI, equal weight
# over CRE and PLA), individual participants jittered and joined PRE -> POST, and a
# bracket labeled with the Time p-value (bold if p < 0.05).
fmt_bracket <- function(p) {
  text <- ifelse(!is.na(p) & p < 0.001, "p < 0.001", paste0("p = ", fmt_p(p)))
  ifelse(!is.na(p) & p < 0.05, paste0("**", text, "**"), text)
}

make_pooled_prepost_plot <- function(est, raw, p_tbl, title, subtitle,
                                     bar_width = 0.6, jitter_width = 0.12) {
  est <- est %>%
    filter(!is.na(estimate)) %>%
    mutate(time = factor(as.character(time), levels = c("PRE", "POST")),
           x    = if_else(time == "PRE", 1, 2))
  
  set.seed(1)
  subj_jit <- raw %>% distinct(subject) %>%
    mutate(jit = runif(n(), -jitter_width, jitter_width))
  
  raw <- raw %>%
    left_join(subj_jit, by = "subject") %>%
    mutate(x = if_else(time == "PRE", 1, 2) + jit)
  
  # Bracket sits just above the highest point or CI in each panel
  brk <- bind_rows(est %>% transmute(variable_label, y = upper.CL),
                   raw %>% transmute(variable_label, y = value)) %>%
    group_by(variable_label) %>%
    summarise(y_top = max(y, na.rm = TRUE), .groups = "drop") %>%
    left_join(p_tbl, by = "variable_label") %>%
    mutate(y_tick = y_top * 1.03,
           y_br   = y_top * 1.06,
           y_lab  = y_top * 1.08,
           label  = fmt_bracket(p_time))
  
  ggplot(est, aes(x = x, y = estimate)) +
    geom_col(aes(alpha = time), fill = pooled_color, width = bar_width, color = colors$black) +
    geom_line(data = raw, aes(x = x, y = value, group = subject),
              inherit.aes = FALSE, color = colors$gray_dark, alpha = 0.35, linewidth = 0.3) +
    geom_point(data = raw, aes(x = x, y = value), inherit.aes = FALSE,
               shape = 21, size = 1.4, alpha = 0.85, fill = "white", color = colors$gray_dark) +
    geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.5) +
    geom_segment(data = brk, aes(x = 1, xend = 2, y = y_br,   yend = y_br),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 1, xend = 1, y = y_tick, yend = y_br),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = brk, aes(x = 2, xend = 2, y = y_tick, yend = y_br),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_richtext(data = brk, aes(x = 1.5, y = y_lab, label = label), inherit.aes = FALSE,
                  size = 2.5, vjust = 0, fill = NA, label.color = NA,
                  label.padding = unit(0, "lines")) +
    scale_alpha_manual(values = c(PRE = 0.4, POST = 1), guide = "none") +
    scale_x_continuous(breaks = c(1, 2), labels = c("PRE", "POST")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    facet_wrap(~ variable_label, scales = "free_y") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    facet_theme()
}

# --- 6: Build the four figures -----------------------------------------------------
raw_long_obs <- pheno_long %>% filter(!is.na(value))
raw_long_paired <- ancova_data %>%
  dplyr::select(subject, group, variable_label, PRE, POST) %>%
  pivot_longer(c(PRE, POST), names_to = "time", values_to = "value") %>%
  mutate(time = factor(time, levels = c("PRE", "POST")))

p_anova_prepost <- make_prepost_plot(
  est   = anova_means %>% mutate(group = factor(group, levels = c("PLA", "CRE"))),
  raw   = raw_long_obs,
  annot = annot_anova_prepost,
  title    = "RET Outcomes: CRE vs PLA (mixed-model ANOVA)",
  subtitle = "Bars = model-estimated means \u00b1 95% CI; points/lines = individual participants (all available observations)"
)

p_anova_change <- make_change_plot(
  est   = anova_change,
  raw   = ancova_data,
  annot = annot_anova_change,
  title    = "Change with RET: CRE vs PLA (mixed-model ANOVA)",
  subtitle = "Model-estimated POST \u2212 PRE \u00b1 95% CI; points = observed change (paired participants)"
)

p_ancova_prepost <- make_prepost_plot(
  est   = ancova_prepost,
  raw   = raw_long_paired,
  annot = annot_ancova_prepost,
  title    = "RET Outcomes: CRE vs PLA (ANCOVA)",
  subtitle = "Bars: PRE = pooled baseline, POST = mean adjusted to that baseline (\u00b1 95% CI); points/lines = observed values, paired participants"
)

p_ancova_change <- make_change_plot(
  est   = ancova_change,
  raw   = ancova_data,
  annot = annot_ancova_change,
  title    = "Change with RET: CRE vs PLA (ANCOVA)",
  subtitle = "Baseline-adjusted POST \u2212 PRE \u00b1 95% CI; points = observed change (paired participants)"
)

p_anova_pooled_prepost <- make_pooled_prepost_plot(
  est   = anova_time_means,
  raw   = raw_long_obs,
  p_tbl = anova_p %>% dplyr::select(variable_label, p_time),
  title    = "RET Outcomes: Pooled CRE + PLA (mixed-model ANOVA)",
  subtitle = "Bars = model-estimated means \u00b1 95% CI, averaged over CRE and PLA; points/lines = individual participants; p = Time effect"
)

# --- 7: Save figures ------------------------------------------------------------
n_vars <- length(config$variables)
n_cols <- ceiling(sqrt(n_vars))
n_rows <- ceiling(n_vars / n_cols)
plot_width  <- max(8, 3.2 * n_cols)
plot_height <- max(5, 3.0 * n_rows)

save_fig <- function(p, name, dir = config$fig_dir) {
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
         width = plot_width, height = plot_height, device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), p,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
}

save_fig(p_anova_prepost,  "fig_anova_prepost")
save_fig(p_anova_change,   "fig_anova_change")
save_fig(p_ancova_prepost, "fig_ancova_prepost")
save_fig(p_ancova_change,  "fig_ancova_change")
save_fig(p_anova_pooled_prepost, "fig_anova_pooled_prepost")

# --- 7b: Normality diagnostic figures ---------------------------------------------
# Points close to the line = approximately normal. Look for systematic curvature
# (skew), S-shapes (heavy/light tails), or single points far off the line (outliers).
make_qq_plot <- function(df, sample_col, title, subtitle, by_group = TRUE) {
  p <- if (by_group) {
    ggplot(df, aes(sample = .data[[sample_col]], color = group)) +
      scale_color_manual(values = group_colors, name = "Group")
  } else {
    ggplot(df, aes(sample = .data[[sample_col]]))
  }
  p +
    stat_qq(size = 1.2, alpha = 0.8) +
    stat_qq_line(linewidth = 0.4) +
    facet_wrap(~ variable_label, scales = "free") +
    labs(title = title, subtitle = subtitle,
         x = "Theoretical quantiles", y = "Sample quantiles") +
    facet_theme()
}

save_fig(make_qq_plot(ancova_data, "CHG",
                      "Normality: change scores (POST \u2212 PRE)",
                      "Paired participants; one line per group"),
         "diag_qq_change", dir = config$diag_dir)
save_fig(make_qq_plot(lmm_resid, "resid",
                      "Normality: mixed-model residuals",
                      "Conditional residuals, pooled across groups and timepoints",
                      by_group = FALSE),
         "diag_qq_resid_lmm", dir = config$diag_dir)
save_fig(make_qq_plot(ancova_resid, "resid",
                      "Normality: ANCOVA residuals",
                      "POST ~ group + PRE; pooled across groups",
                      by_group = FALSE),
         "diag_qq_resid_ancova", dir = config$diag_dir)

# --- 8: Manuscript table ---------------------------------------------------------
pheno_summary <- pheno_long %>%
  group_by(variable, variable_label, group, time) %>%
  summarise(mean_val = mean(value, na.rm = TRUE), sd_val = sd(value, na.rm = TRUE),
            n = sum(!is.na(value)), .groups = "drop")

tbl_digits <- pheno_summary %>%
  group_by(variable) %>%
  summarise(d = ifelse(max(abs(mean_val), na.rm = TRUE) >= 20, 1L, 2L), .groups = "drop")

fx <- function(x, d) unname(mapply(function(xx, dd) formatC(xx, format = "f", digits = dd), x, d))
fx_ci <- function(est, lo, hi, d) paste0(fx(est, d), " [", fx(lo, d), ", ", fx(hi, d), "]")

ci_col <- function(tbl) {
  tbl %>% left_join(tbl_digits, by = "variable") %>%
    mutate(ci = fx_ci(estimate, lower.CL, upper.CL, d), p = fmt_p(p.value))
}

main_n <- pheno_long %>% filter(!is.na(value)) %>%
  group_by(variable, group) %>% summarise(n = n_distinct(subject), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = n, names_prefix = "n_", values_fill = 0)

main_npair <- ancova_data %>% count(variable, group) %>%
  pivot_wider(names_from = group, values_from = n, names_prefix = "npair_", values_fill = 0)

main_means <- pheno_summary %>% left_join(tbl_digits, by = "variable") %>%
  mutate(msd = paste0(fx(mean_val, d), " \u00b1 ", fx(sd_val, d)),
         col = paste(group, time, sep = "_")) %>%
  dplyr::select(variable, col, msd) %>%
  pivot_wider(names_from = col, values_from = msd)

t_base        <- ci_col(anova_group_at_time %>% filter(time == "PRE")) %>% transmute(variable, baseline_diff = ci)
t_chg         <- ci_col(anova_change) %>% dplyr::select(variable, group, ci) %>%
  pivot_wider(names_from = group, values_from = ci, names_prefix = "change_")
t_int         <- ci_col(anova_diff_in_change) %>% transmute(variable, diff_in_change = ci, diff_in_change_p = p)
t_pool        <- ci_col(anova_pooled_change)  %>% transmute(variable, pooled_change = ci, pooled_change_p = p)
t_ancova      <- ci_col(ancova_group_diff)    %>% transmute(variable, ancova_diff = ci, ancova_p = p)
t_ancova_pool <- ci_col(ancova_pooled_change) %>% transmute(variable, ancova_pooled = ci, ancova_pooled_p = p)

Main_Table <- var_label_lookup %>%
  left_join(main_n, by = "variable") %>% left_join(main_npair, by = "variable") %>%
  left_join(main_means, by = "variable") %>% left_join(t_base, by = "variable") %>%
  left_join(t_chg, by = "variable") %>% left_join(t_int, by = "variable") %>%
  left_join(t_pool, by = "variable") %>% left_join(t_ancova, by = "variable") %>%
  left_join(t_ancova_pool, by = "variable") %>%
  transmute(
    Outcome  = as.character(variable_label),
    n        = paste0(n_CRE, " / ", n_PLA),
    n_paired = paste0(npair_CRE, " / ", npair_PLA),
    CRE_PRE, CRE_POST, PLA_PRE, PLA_POST,
    baseline_diff, change_CRE, change_PLA,
    diff_in_change, diff_in_change_p,
    pooled_change, pooled_change_p,
    ancova_diff, ancova_p,
    ancova_pooled, ancova_pooled_p
  )

names(Main_Table) <- c(
  "Outcome", "n (CRE / PLA)", "n paired (CRE / PLA)",
  "CRE PRE (mean \u00b1 SD)", "CRE POST (mean \u00b1 SD)",
  "PLA PRE (mean \u00b1 SD)", "PLA POST (mean \u00b1 SD)",
  "Baseline diff, CRE \u2212 PLA [95% CI]",
  "Change (POST \u2212 PRE) CRE [95% CI]",
  "Change (POST \u2212 PRE) PLA [95% CI]",
  "Diff in change, CRE \u2212 PLA (Group\u00d7Time) [95% CI]", "Group\u00d7Time p",
  "Pooled change (POST \u2212 PRE) [95% CI]", "Pooled change p",
  "ANCOVA adj. diff, CRE \u2212 PLA [95% CI]", "ANCOVA p",
  "ANCOVA pooled adj. change (POST \u2212 PRE) [95% CI]", "ANCOVA pooled change p"
)

# Descriptive % change (mean +/- SD of individual % changes, paired participants),
# by group and pooled across all paired participants. Reporting only -- not used in
# any model.
pct_summary <- bind_rows(
  ancova_data %>% mutate(group = as.character(group)),
  ancova_data %>% mutate(group = "Pooled")
) %>%
  mutate(group = factor(group, levels = c("CRE", "PLA", "Pooled"))) %>%
  group_by(variable, variable_label, group) %>%
  summarise(n_paired   = sum(!is.na(PCT)),
            mean_pct   = mean(PCT,   na.rm = TRUE),
            sd_pct     = sd(PCT,     na.rm = TRUE),
            median_pct = median(PCT, na.rm = TRUE),
            .groups = "drop")

pct_cols <- c("% change CRE (mean \u00b1 SD)", "% change PLA (mean \u00b1 SD)",
              "% change Pooled (mean \u00b1 SD)")

t_pct <- pct_summary %>%
  transmute(Outcome = as.character(variable_label), group,
            msd = paste0(sprintf("%.1f", mean_pct), " \u00b1 ", sprintf("%.1f", sd_pct))) %>%
  pivot_wider(names_from = group, values_from = msd,
              names_glue = "% change {group} (mean \u00b1 SD)")

# Observed pooled (CRE + PLA) PRE and POST means, all available observations
pooled_obs <- pheno_long %>%
  group_by(variable, variable_label, time) %>%
  summarise(m = mean(value, na.rm = TRUE), s = sd(value, na.rm = TRUE), .groups = "drop") %>%
  left_join(tbl_digits, by = "variable") %>%
  transmute(Outcome = as.character(variable_label),
            col     = paste0("Pooled ", time, " (mean \u00b1 SD)"),
            msd     = paste0(fx(m, d), " \u00b1 ", fx(s, d))) %>%
  pivot_wider(names_from = col, values_from = msd)

pooled_obs_cols <- c("Pooled PRE (mean \u00b1 SD)", "Pooled POST (mean \u00b1 SD)")

Main_Table <- Main_Table %>%
  left_join(t_pct, by = "Outcome") %>%
  left_join(pooled_obs, by = "Outcome") %>%
  relocate(all_of(pct_cols), .after = "Change (POST \u2212 PRE) PLA [95% CI]") %>%
  relocate(all_of(pooled_obs_cols), .after = "PLA POST (mean \u00b1 SD)")

Main_Table_Notes <- tibble(Note = c(
  paste0("Mixed-model estimates: value ~ group * time + (1|subject); ", DF_METHOD,
         " df; all available observations."),
  "ANCOVA: POST ~ group + PRE (PRE centered at the cohort mean); paired participants only. The adjusted CRE \u2212 PLA difference is identical on the POST and change-score scales.",
  "Baseline diff is descriptive (no p). Pooled change = POST \u2212 PRE averaged with equal weight over CRE and PLA (mixed model: all observations; ANCOVA: paired participants, adjusted to the cohort mean baseline).",
  "% change = mean \u00b1 SD of individual 100 \u00d7 (POST \u2212 PRE) / PRE, paired participants only; Pooled = all paired participants combined (CRE + PLA). Descriptive; all inference is on absolute values.",
  "p-values are uncorrected for multiplicity (pre-specified outcome battery).",
  "No non-training comparison arm: pooled change describes change during RET, not an isolated causal effect of RET."
))

cat("\n=== Main table ===\n")
print(as.data.frame(Main_Table))

# --- 8b: ANCOVA main table (paired participants only) ---------------------------
# Every column is from the same paired sample the ANCOVA was fit to, so observed
# means, adjusted means, and tests all describe the same participants.
anc_obs <- ancova_data %>%
  pivot_longer(c(PRE, POST), names_to = "time", values_to = "value") %>%
  group_by(variable, group, time) %>%
  summarise(m = mean(value), s = sd(value), .groups = "drop") %>%
  left_join(tbl_digits, by = "variable") %>%
  mutate(msd = paste0(fx(m, d), " \u00b1 ", fx(s, d)),
         col = paste(group, time, sep = "_")) %>%
  dplyr::select(variable, col, msd) %>%
  pivot_wider(names_from = col, values_from = msd)

t_anc_post <- ci_col(ancova_post_means) %>% dplyr::select(variable, group, ci) %>%
  pivot_wider(names_from = group, values_from = ci, names_prefix = "adjpost_")
t_anc_chg  <- ci_col(ancova_change) %>% dplyr::select(variable, group, ci) %>%
  pivot_wider(names_from = group, values_from = ci, names_prefix = "adjchg_")
t_anc_info <- ancova_coefs %>%
  transmute(variable, b_pre = sprintf("%.2f", b_PRE), r2 = sprintf("%.3f", r2))

ANCOVA_Main_Table <- var_label_lookup %>%
  left_join(main_npair,    by = "variable") %>% left_join(anc_obs,    by = "variable") %>%
  left_join(t_anc_post,    by = "variable") %>% left_join(t_anc_chg,  by = "variable") %>%
  left_join(t_ancova,      by = "variable") %>% left_join(t_ancova_pool, by = "variable") %>%
  left_join(t_anc_info,    by = "variable") %>%
  transmute(
    Outcome  = as.character(variable_label),
    n_paired = paste0(npair_CRE, " / ", npair_PLA),
    CRE_PRE, CRE_POST, PLA_PRE, PLA_POST,
    adjpost_CRE, adjpost_PLA,
    adjchg_CRE, adjchg_PLA,
    ancova_diff, ancova_p,
    ancova_pooled, ancova_pooled_p,
    b_pre, r2
  )

names(ANCOVA_Main_Table) <- c(
  "Outcome", "n paired (CRE / PLA)",
  "CRE PRE (mean \u00b1 SD)", "CRE POST (mean \u00b1 SD)",
  "PLA PRE (mean \u00b1 SD)", "PLA POST (mean \u00b1 SD)",
  "Adj. POST CRE [95% CI]", "Adj. POST PLA [95% CI]",
  "Adj. change (POST \u2212 PRE) CRE [95% CI]", "Adj. change (POST \u2212 PRE) PLA [95% CI]",
  "Adj. diff, CRE \u2212 PLA [95% CI]", "Group p",
  "Pooled adj. change (POST \u2212 PRE) [95% CI]", "Pooled change p",
  "Baseline slope (b PRE)", "R\u00b2"
)

ANCOVA_Main_Table <- ANCOVA_Main_Table %>%
  left_join(t_pct, by = "Outcome") %>%
  relocate(all_of(pct_cols), .after = "Adj. change (POST \u2212 PRE) PLA [95% CI]")

ANCOVA_Table_Notes <- tibble(Note = c(
  "ANCOVA: POST ~ group + PRE (PRE centered at the cohort mean); paired participants only. All columns use the same paired sample.",
  "Adjusted POST and adjusted change are evaluated at the cohort mean baseline. Adj. diff (CRE \u2212 PLA) is identical on the POST and change scales and is the test of the supplement effect (Group p).",
  "Pooled adj. change = adjusted POST \u2212 PRE averaged with equal weight over CRE and PLA, tested against 0 (change during RET; no non-training control arm).",
  "% change = mean \u00b1 SD of individual 100 \u00d7 (POST \u2212 PRE) / PRE; Pooled = all paired participants combined. Descriptive only.",
  "Baseline slope: ~1 = participants shifted together; <1 = regression to the mean. R\u00b2 = variance in POST explained by the model.",
  "p-values are uncorrected for multiplicity (pre-specified outcome battery)."
))

cat("\n=== ANCOVA main table ===\n")
print(as.data.frame(ANCOVA_Main_Table))

# --- 9: Data export ---------------------------------------------------------------
drop_var <- function(x) dplyr::select(x, -variable)

write_xlsx(
  list(
    Main_Table             = Main_Table,
    Main_Table_Notes       = Main_Table_Notes,
    ANCOVA_Main_Table      = ANCOVA_Main_Table,
    ANCOVA_Table_Notes     = ANCOVA_Table_Notes,
    ANOVA_Means            = drop_var(anova_means),
    ANOVA_Pooled_Means     = drop_var(anova_time_means),
    ANOVA_Change_by_Group  = drop_var(anova_change),
    ANOVA_Pooled_Change    = drop_var(anova_pooled_change),
    ANOVA_Diff_in_Change   = drop_var(anova_diff_in_change),
    ANOVA_Group_at_Time    = drop_var(anova_group_at_time),
    ANOVA_Group_Main       = drop_var(anova_group_main),
    ANOVA_F_Tests          = drop_var(model_results),
    ANCOVA_PrePost_Means   = drop_var(ancova_prepost),
    ANCOVA_Change_by_Group = drop_var(ancova_change),
    ANCOVA_Group_Diff      = drop_var(ancova_group_diff),
    ANCOVA_Pooled_Change   = drop_var(ancova_pooled_change),
    ANCOVA_Model_Info      = drop_var(ancova_coefs),
    Pct_Change_Summary     = drop_var(pct_summary),
    Summary_Stats          = drop_var(pheno_summary),
    Cohort_Availability    = cohort_availability,
    Cohort_Roster          = cohort_roster,
    Long_Data              = drop_var(pheno_long),
    ANCOVA_Data            = drop_var(ancova_data),
    Normality_Raw          = drop_var(norm_raw),
    Normality_Change       = drop_var(norm_change),
    Normality_Residuals    = drop_var(norm_resid)
  ),
  path = file.path(config$data_dir, "phenotypic_overview_data.xlsx")
)

cat("\nSaved figures to:", config$fig_dir, "\n")
cat("Saved data to:", config$data_dir, "\n")

sessionInfo()
