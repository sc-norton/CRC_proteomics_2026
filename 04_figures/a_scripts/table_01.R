# =============================================================================
# table1_cohort_summary.R
# Table 1: cohort baseline characteristics, from raw (unfiltered) metadata.
#
# Source: CRm_meta.csv -- one row per subject per timepoint (T1/T2); baseline
# demographic/anthropometric values are duplicated identically across a
# subject's T1/T2 rows, so filtering to timepoint == "T1" gives exactly one
# row per subject (every subject, CTL and SURV, has a T1 row).
#
# Two between-group comparisons, each as its own p-value column:
#   1. SURV vs CTL   -- whole cohort (n = 27)
#   2. CRE vs PLA    -- within SURV only (n = 17; supp is NA for CTL)
#
# Continuous variables: mean +/- SD, Welch two-sample t-test.
# Categorical variables: n (%), Fisher's exact test.
# Variables not applicable to one group (e.g. cancer stage for CTL) are
# reported as "N/A" with no p-value, not silently dropped or errored on.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, purrr, tibble, writexl)

# --- 0: Config -------------------------------------------------------------
config <- list(
  meta_csv = "00_input/muscle/CRm_meta.csv",   # >>> VERIFY <<< actual path in your project
  out_dir  = "04_figures/data/table_01"
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1: Load raw metadata, reduce to one row per subject (baseline / T1)
# =============================================================================

meta_raw <- read_csv(config$meta_csv, show_col_types = FALSE)

baseline <- meta_raw %>% filter(timepoint == "T1")

cat(sprintf("Baseline cohort: %d subjects (%d SURV, %d CTL)\n",
            nrow(baseline),
            sum(baseline$cancer == "SURV"), sum(baseline$cancer == "CTL")))

# =============================================================================
# 2: Variable definitions
# =============================================================================
# type: "continuous" -> mean +/- SD + t-test
#       "categorical" -> n (%) + Fisher's exact test
# Units taken directly from column-name suffixes (_cm, _kg, _lbs) where
# explicit; left generic where the source column doesn't state a unit
# (e.g. sts_max_pwr) rather than guessing.

var_defs <- tribble(
  ~domain,          ~var,                   ~label,                          ~type,
  "Demographics",   "age",                  "Age (yrs)",                     "continuous",
  "Demographics",   "sex",                  "Sex",                           "categorical",
  "Demographics",   "stage",                "Cancer Stage",                  "categorical",
  "Demographics",   "yrs_since_dx",         "Years Since Diagnosis",         "continuous",
  "Demographics",   "yrs_since_chemo",      "Years Since Chemotherapy",      "continuous",
  "Anthropometry",  "pre_height_cm",        "Height (cm)",                   "continuous",
  "Anthropometry",  "pre_weight_kg",        "Weight (kg)",                   "continuous",
  "Anthropometry",  "pre_BMI",              "BMI (kg/m\u00b2)",              "continuous",
  "Body Composition","pre_LBM_kg",          "Lean Body Mass (kg)",           "continuous",
  "Body Composition","pre_ALM_kg",          "Appendicular Lean Mass (kg)",   "continuous",
  "Body Composition","pre_ALM_index",       "ALM Index (kg/m\u00b2)",        "continuous",
  "Body Composition","pre_fat_kg",          "Fat Mass (kg)",                 "continuous",
  "Body Composition","pre_fat_pct",         "Body Fat (%)",                  "continuous",
  "Function",       "pre_sts_max_pwr",      "Sit-to-Stand Max Power",        "continuous",
  "Function",       "pre_grip_lbs",         "Grip Strength (lbs)",           "continuous",
  "Function",       "pre_chest_press_lbs",  "Chest Press 1RM (lbs)",         "continuous",
  "Function",       "pre_leg_ext_lbs",      "Leg Extension 1RM (lbs)",       "continuous"
)

comparisons <- list(
  list(name = "SURV_vs_CTL", group_col = "cancer", levels = c("SURV", "CTL"),
       data = baseline),
  list(name = "CRE_vs_PLA", group_col = "supp", levels = c("CRE", "PLA"),
       data = filter(baseline, cancer == "SURV"))
)

# =============================================================================
# 3: Summary + test helpers -- both fail soft (NA/"N/A"), never error out
# =============================================================================

fmt_mean_sd <- function(x, digits = 1) {
  n <- sum(!is.na(x))
  if (n == 0) return("N/A")
  sprintf("%.*f \u00b1 %.*f", digits, mean(x, na.rm = TRUE), digits, sd(x, na.rm = TRUE))
}

fmt_p <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("< 0.001")
  sprintf("%.3f", p)
}

summarize_continuous_row <- function(df, var, label, domain, group_col, levels) {
  g1 <- df[[var]][df[[group_col]] == levels[1]]
  g2 <- df[[var]][df[[group_col]] == levels[2]]
  n1 <- sum(!is.na(g1)); n2 <- sum(!is.na(g2))
  
  p <- NA_real_
  if (n1 >= 2 && n2 >= 2) {
    p <- tryCatch(t.test(g1, g2)$p.value, error = function(e) NA_real_)
  }
  
  tibble(
    domain = domain, variable = label, level = NA_character_,
    group1_summary = fmt_mean_sd(g1), group1_n = n1,
    group2_summary = fmt_mean_sd(g2), group2_n = n2,
    p_value = p, p_display = fmt_p(p),
    note = if (n1 == 0 || n2 == 0) "Not available for one group" else NA_character_
  )
}

summarize_categorical_rows <- function(df, var, label, domain, group_col, levels) {
  g1 <- df[[var]][df[[group_col]] == levels[1]]
  g2 <- df[[var]][df[[group_col]] == levels[2]]
  n1_total <- sum(!is.na(g1)); n2_total <- sum(!is.na(g2))
  
  if (n1_total == 0 && n2_total == 0) {
    return(tibble(
      domain = domain, variable = label, level = NA_character_,
      group1_summary = "N/A", group1_n = 0, group2_summary = "N/A", group2_n = 0,
      p_value = NA_real_, p_display = "", note = "Not available for either group"
    ))
  }
  
  all_levels <- sort(unique(na.omit(c(g1, g2))))
  
  p <- NA_real_
  if (n1_total > 0 && n2_total > 0) {
    tab <- table(c(g1, g2), c(rep(levels[1], length(g1)), rep(levels[2], length(g2))),
                 useNA = "no")
    p <- tryCatch(fisher.test(tab)$p.value, error = function(e) NA_real_)
  }
  
  map_dfr(all_levels, function(lvl) {
    n1 <- sum(g1 == lvl, na.rm = TRUE); n2 <- sum(g2 == lvl, na.rm = TRUE)
    tibble(
      domain = domain, variable = label, level = as.character(lvl),
      group1_summary = if (n1_total > 0) sprintf("%d (%.0f%%)", n1, 100 * n1 / n1_total) else "N/A",
      group1_n = n1,
      group2_summary = if (n2_total > 0) sprintf("%d (%.0f%%)", n2, 100 * n2 / n2_total) else "N/A",
      group2_n = n2,
      p_value = if (lvl == all_levels[1]) p else NA_real_,   # p shown once per variable, on first level row
      p_display = if (lvl == all_levels[1]) fmt_p(p) else "",
      note = if (n1_total == 0 || n2_total == 0) "Not available for one group" else NA_character_
    )
  })
}

build_comparison_table <- function(var_defs, comp) {
  pmap_dfr(var_defs, function(domain, var, label, type) {
    if (type == "continuous") {
      summarize_continuous_row(comp$data, var, label, domain, comp$group_col, comp$levels)
    } else {
      summarize_categorical_rows(comp$data, var, label, domain, comp$group_col, comp$levels)
    }
  }) %>%
    rename_with(~ paste0(comp$levels[1], "_summary"), group1_summary) %>%
    rename_with(~ paste0(comp$levels[1], "_n"), group1_n) %>%
    rename_with(~ paste0(comp$levels[2], "_summary"), group2_summary) %>%
    rename_with(~ paste0(comp$levels[2], "_n"), group2_n) %>%
    rename_with(~ paste0("p_", comp$name), p_value) %>%
    rename_with(~ paste0("p_display_", comp$name), p_display) %>%
    rename_with(~ paste0("note_", comp$name), note)
}

# =============================================================================
# 4: Build both comparison tables, then join on domain/variable/level
# =============================================================================

tbl_survctl <- build_comparison_table(var_defs, comparisons[[1]])
tbl_crepla  <- build_comparison_table(var_defs, comparisons[[2]])

table1 <- full_join(tbl_survctl, tbl_crepla, by = c("domain", "variable", "level")) %>%
  select(
    domain, variable, level,
    SURV_summary, SURV_n, CTL_summary, CTL_n, p_display_SURV_vs_CTL, note_SURV_vs_CTL,
    CRE_summary, CRE_n, PLA_summary, PLA_n, p_display_CRE_vs_PLA, note_CRE_vs_PLA
  ) %>%
  rename(
    `p (SURV vs CTL)` = p_display_SURV_vs_CTL,
    `Note (SURV vs CTL)` = note_SURV_vs_CTL,
    `p (CRE vs PLA)` = p_display_CRE_vs_PLA,
    `Note (CRE vs PLA)` = note_CRE_vs_PLA
  )

cat("\n=== Table 1 preview ===\n")
print(table1, n = Inf)

# =============================================================================
# 5: Export
# =============================================================================

write_xlsx(
  list(
    "Table 1"          = table1,
    "Baseline raw data" = baseline
  ),
  path = file.path(config$out_dir, "table_01_cohort_summary.xlsx")
)

cat(sprintf("\nSaved to %s\n", file.path(config$out_dir, "table_01_cohort_summary.xlsx")))
