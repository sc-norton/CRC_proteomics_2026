# Fig Phenotypic Overview -- DXA/phenotypic variables, CRE vs PLA
# Mixed model: value ~ group * time + (1|subject), interaction term is the
# test of interest (does change over time differ by group?).
# Barplots: PRE vs POST bars per group, faceted by variable.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               patchwork, lme4, lmerTest, emmeans, scales, writexl, ggtext)

# --- Configuration -----------------------------------------------------------
# ADJUST THESE to match your actual file. Everything downstream depends on
# these names being correct -- the diagnostic block in section 1 will tell
# you immediately if something doesn't match rather than failing silently
# or silently dropping variables.
config <- list(
  data_file  = "00_input/muscle/CRm_meta.csv",   # CRm_meta.csv structure
  fig_dir    = "04_figures/fig_phenotypic_overview",
  data_dir   = "04_figures/data/fig_phenotypic_overview",

  subject_col   = "pid",       # unique per subject (sample_id has a _T1/_T2 suffix, pid doesn't)
  group_col     = "supp",      # "CRE" / "PLA"
  cancer_col    = "cancer",    # "SURV" / "CTL" -- the CRE/PLA sub-trial only exists within SURV
  cancer_filter = "SURV",

  # Named vector: names are the base variable name as it appears AFTER the
  # pre_/post_ prefix in the file, values are the display label used in
  # figures/tables. Add/remove rows here to change what's plotted -- height
  # is intentionally left out (doesn't meaningfully change PRE->POST).
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

  # Columns are named as prefix + variable, e.g. "pre_BMI" / "post_BMI".
  prefix_pre  = "pre_",
  prefix_post = "post_"
)

dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

group_colors <- c(CRE = colors$orange, PLA = colors$blue_dark)

adjp_thresh <- 0.05

# --- 1: Load & reshape data ---------------------------------------------------
raw <- if (grepl("\\.xlsx?$", config$data_file)) {
  read_excel(config$data_file)
} else {
  read_csv(config$data_file, show_col_types = FALSE)
}

# --- Diagnostic check --------------------------------------------------------
# Run this before trusting anything downstream. If subject_col/group_col/
# variable columns don't show up as expected, fix config above -- don't
# patch symptoms further down.
cat("=== Raw data structure ===\n")
cat("Columns:", paste(colnames(raw), collapse = ", "), "\n")
cat("N rows:", nrow(raw), "\n")

stopifnot(
  "group_col not found in data -- check config$group_col" = config$group_col %in% colnames(raw),
  "subject_col not found in data -- check config$subject_col" = config$subject_col %in% colnames(raw),
  "cancer_col not found in data -- check config$cancer_col" = config$cancer_col %in% colnames(raw)
)

# Each subject shows up once per timepoint row (T1/T2), but the pre_*/post_*
# columns are already identical across those duplicate rows -- they're not
# genuinely long-format data. Collapsing to one row per subject here avoids
# double-counting every subject in both the summary stats and the mixed
# model below (which would silently inflate n and distort SEs/p-values).
rows_per_subject <- raw %>% count(.data[[config$subject_col]])
cat("Rows per subject before collapsing (expect mostly 1-2):\n")
print(table(rows_per_subject$n))

raw_subj <- raw %>% distinct(.data[[config$subject_col]], .keep_all = TRUE)
cat("N subjects after collapsing duplicate timepoint rows:", nrow(raw_subj), "\n")

# Restrict to the SURV subjects who received a CRE/PLA supplement assignment
# -- CTL subjects have no supp value and no POST measurements, so they're
# out of scope for this comparison.
pheno_subj <- raw_subj %>%
  filter(.data[[config$cancer_col]] == config$cancer_filter,
         .data[[config$group_col]] %in% c("CRE", "PLA"))

cat("\nGroup counts (", config$cancer_filter, ", CRE/PLA only):\n", sep = "")
print(table(pheno_subj[[config$group_col]]))

expected_cols <- as.vector(outer(
  c(config$prefix_pre, config$prefix_post),
  names(config$variables),
  paste0
))
missing_cols <- setdiff(expected_cols, colnames(raw))
if (length(missing_cols) > 0) {
  stop("The following expected PRE/POST columns are missing from the data: ",
       paste(missing_cols, collapse = ", "),
       "\nCheck config$variables and config$prefix_pre/prefix_post.")
}
cat("All expected variable columns found:", paste(expected_cols, collapse = ", "), "\n")

# Reshape wide -> long. Each row becomes subject x variable x timepoint.
pheno_long <- pheno_subj %>%
  dplyr::select(subject = all_of(config$subject_col), group = all_of(config$group_col),
                all_of(expected_cols)) %>%
  pivot_longer(
    cols = all_of(expected_cols),
    names_to = "colname",
    values_to = "value"
  ) %>%
  mutate(
    time = case_when(
      str_starts(colname, config$prefix_pre)  ~ "PRE",
      str_starts(colname, config$prefix_post) ~ "POST"
    ),
    variable = str_remove(colname, paste0("^", config$prefix_pre, "|^", config$prefix_post)),
    variable_label = config$variables[variable],
    variable_label = factor(variable_label, levels = unname(config$variables)),
    group = factor(group, levels = c("PLA", "CRE")),
    time  = factor(time,  levels = c("PRE", "POST"))
  ) %>%
  dplyr::select(-colname)

cat("\n=== Long data check ===\n")
cat("Rows:", nrow(pheno_long), "\n")
cat("Missing values per variable:\n")
print(pheno_long %>% group_by(variable_label) %>% summarise(n_na = sum(is.na(value)), .groups = "drop"))
print(table(pheno_long$variable_label, pheno_long$time, pheno_long$group))

# --- 2: Mixed models per variable ---------------------------------------------
# value ~ group * time + (1|subject), fit once per variable and stored so
# the same fitted model feeds both the omnibus ANOVA table here and the
# post-hoc emmeans contrasts in section 2b -- keeps the two consistent
# rather than refitting separately.
fit_variable_model <- function(var_name) {
  df <- pheno_long %>% filter(variable == var_name) %>% droplevels()
  tryCatch(
    lmerTest::lmer(value ~ group * time + (1 | subject), data = df),
    error = function(e) e
  )
}

fitted_models <- setNames(
  lapply(names(config$variables), fit_variable_model),
  names(config$variables)
)

extract_anova_table <- function(m) {
  if (inherits(m, "error")) {
    return(tibble(term = c("group", "time", "group:time"),
                   F_value = NA_real_, df_num = NA_real_, df_den = NA_real_,
                   p.value = NA_real_, singular_fit = NA,
                   note = conditionMessage(m)))
  }
  is_singular <- lme4::isSingular(m, tol = 1e-4)
  a <- anova(m)
  tibble(
    term    = rownames(a),
    F_value = a$`F value`,
    df_num  = a$NumDF,
    df_den  = a$DenDF,
    p.value = a$`Pr(>F)`,
    singular_fit = is_singular
  )
}

var_label_lookup <- tibble(
  variable       = names(config$variables),
  variable_label = factor(unname(config$variables), levels = unname(config$variables))
)

model_results <- imap_dfr(fitted_models, function(m, var_name) {
  extract_anova_table(m) %>% mutate(variable = var_name, .before = 1)
}) %>%
  left_join(var_label_lookup, by = "variable") %>%
  relocate(variable_label, .after = variable) %>%
  # BH correction is applied within each term family separately (all
  # variables' group p-values together, all time p-values together, all
  # interaction p-values together) rather than pooling across term types,
  # since these are three distinct hypothesis families.
  group_by(term) %>%
  mutate(p.adjust = p.adjust(p.value, method = "BH")) %>%
  ungroup()

cat("\n=== ANOVA results (group, time, group:time) ===\n")
print(model_results)

if (any(model_results$singular_fit, na.rm = TRUE)) {
  cat("\nNote: singular fit flagged for:",
      paste(unique(model_results$variable_label[which(model_results$singular_fit)]), collapse = ", "),
      "-- random-intercept variance estimated near zero for these variables.\n")
}

# --- 2b: Post-hoc contrasts (emmeans) ------------------------------------------
# Three families of follow-up contrasts per variable, from the SAME fitted
# model as the omnibus ANOVA above (Satterthwaite df, consistent with
# lmerTest::anova()):
#   1. Simple effect of Time within each Group  -- paired PRE vs POST,
#      separately for CRE and PLA. Answers "did this group actually change."
#   2. Simple effect of Group within each Time   -- CRE vs PLA at baseline,
#      and CRE vs PLA at POST. Answers "were groups different to begin
#      with" / "different after."
#   3. Interaction contrast (difference-of-differences) -- (CRE POST-PRE) -
#      (PLA POST-PRE). Same hypothesis as the omnibus group:time F-test
#      above, but gives an estimate + CI rather than only a p-value.
# "revpairwise" with factor levels ordered PLA/CRE and PRE/POST gives an
# intuitive direction: positive = CRE > PLA, positive = POST > PRE.
# adjust = "none" matches the raw-p convention used for the main ANOVA
# annotations (planned, non-exploratory contrasts on a curated battery, not
# an open-ended screen). Change POSTHOC_ADJUST to "holm" or "tukey" to
# apply correction instead.
POSTHOC_ADJUST <- "none"

run_posthoc <- function(m) {
  if (inherits(m, "error")) {
    return(tibble(contrast_family = NA_character_, contrast = NA_character_,
                   estimate = NA_real_, SE = NA_real_, df = NA_real_,
                   t.ratio = NA_real_, p.value = NA_real_))
  }
  emm <- emmeans::emmeans(m, ~ group * time)

  simple_time <- emmeans::contrast(emm, "revpairwise", by = "group", adjust = POSTHOC_ADJUST) %>%
    as.data.frame() %>%
    mutate(contrast_family = "Simple effect of Time within Group")

  simple_group <- emmeans::contrast(emm, "revpairwise", by = "time", adjust = POSTHOC_ADJUST) %>%
    as.data.frame() %>%
    mutate(contrast_family = "Simple effect of Group within Time")

  interaction_contrast <- emmeans::contrast(emm, interaction = "revpairwise", adjust = POSTHOC_ADJUST) %>%
    as.data.frame() %>%
    mutate(contrast_family = "Interaction (difference-of-differences)")

  bind_rows(simple_time, simple_group, interaction_contrast) %>%
    relocate(contrast_family, .before = 1)
}

posthoc_results <- imap_dfr(fitted_models, function(m, var_name) {
  run_posthoc(m) %>% mutate(variable = var_name, .before = 1)
}) %>%
  left_join(var_label_lookup, by = "variable") %>%
  relocate(variable_label, .after = variable)

cat("\n=== Post-hoc contrasts (emmeans, raw p) ===\n")
print(posthoc_results, n = Inf)

# --- 3: Summary stats for plotting --------------------------------------------
pheno_summary <- pheno_long %>%
  group_by(variable, variable_label, group, time) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    n        = sum(!is.na(value)),
    se_val   = sd_val / sqrt(n),
    .groups = "drop"
  )

# Build the per-facet annotation text: one block per variable listing all
# three ANOVA terms. Uses raw p.value (not BH-adjusted) -- this is a
# curated, pre-specified phenotypic battery rather than a high-dimensional
# screen, so multiplicity correction across variables was judged
# unnecessary here (BH-adjusted values remain available in model_results
# and the ANOVA_Annotations/Model_Results sheets if needed later).
# Each line is bolded independently (via markdown, rendered by geom_richtext)
# only when that term's raw p < 0.05 -- non-significant terms stay plain.
fmt_p <- function(p) {
  case_when(
    is.na(p)  ~ "NA",
    p < 0.001 ~ "<0.001",
    TRUE      ~ sprintf("%.3f", p)
  )
}

fmt_line <- function(name, p) {
  text <- paste0(name, " p = ", fmt_p(p))
  ifelse(!is.na(p) & p < 0.05, paste0("**", text, "**"), text)
}

anova_annotations <- model_results %>%
  dplyr::select(variable, variable_label, term, p.value) %>%
  pivot_wider(names_from = term, values_from = p.value) %>%
  mutate(
    label = paste0(
      fmt_line("Group", group), "<br>",
      fmt_line("Time", time), "<br>",
      fmt_line("Group\u00d7Time", `group:time`)
    )
  )

# --- 4: Plot -------------------------------------------------------------
# PRE vs POST bars per group, faceted by variable. Dodged bars + error bars
# (mean +/- SE), individual subject points jittered behind at low alpha so
# the underlying spread is visible without competing with the bars.
p_overview <- ggplot(pheno_summary, aes(x = time, y = mean_val, fill = group)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65, color = colors$black) +
  geom_errorbar(
    aes(ymin = mean_val - se_val, ymax = mean_val + se_val, group = group),
    position = position_dodge(width = 0.75), width = 0.2, linewidth = 0.5
  ) +
  geom_point(
    data = pheno_long,
    aes(x = time, y = value, fill = group, group = group),
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    shape = 21, size = 1.4, alpha = 0.4, color = colors$gray_dark, inherit.aes = FALSE
  ) +
  geom_richtext(
    data = anova_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, size = 2.3,
    hjust = -0.05, vjust = 1.1, lineheight = 1.1,
    fill = NA, label.color = NA, label.padding = unit(0, "lines")
  ) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  facet_wrap(~ variable_label, scales = "free_y") +
  labs(title = "RET Outcome Variables: CRE vs PLA", x = NULL, y = NULL) +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(
    strip.background = element_rect(fill = colors$gray_light, color = NA),
    panel.grid.major.y = element_line(color = colors$gray_light),
    panel.grid.minor   = element_blank()
  )

# --- 5: Save main figure -------------------------------------------------------
n_vars <- length(config$variables)
n_cols <- ceiling(sqrt(n_vars))
n_rows <- ceiling(n_vars / n_cols)
plot_width  <- max(8, 3.2 * n_cols)
plot_height <- max(5, 3.0 * n_rows)

ggsave(file.path(config$fig_dir, "fig_phenotypic_overview.pdf"),
       p_overview, width = plot_width, height = plot_height, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_phenotypic_overview.png"),
       p_overview, width = plot_width, height = plot_height, dpi = 300, bg = "white")

# --- 6: Percent change figure ---------------------------------------------------
# (POST - PRE) / PRE * 100 per subject per variable. This collapses to a
# single value per subject -- no repeated-measures structure left -- so the
# group comparison here is a simple two-group Welch's t-test, not a mixed
# model. This figure is a descriptive companion to the main one: raw-unit
# bars can obscure relative magnitude of change when variables sit on very
# different scales (e.g. kg vs an index), % change puts everything on a
# common footing.
pheno_pctchange <- pheno_subj %>%
  dplyr::select(subject = all_of(config$subject_col), group = all_of(config$group_col),
                all_of(expected_cols)) %>%
  pivot_longer(cols = all_of(expected_cols), names_to = "colname", values_to = "value") %>%
  mutate(
    time = case_when(
      str_starts(colname, config$prefix_pre)  ~ "PRE",
      str_starts(colname, config$prefix_post) ~ "POST"
    ),
    variable = str_remove(colname, paste0("^", config$prefix_pre, "|^", config$prefix_post))
  ) %>%
  dplyr::select(-colname) %>%
  pivot_wider(names_from = time, values_from = value) %>%
  mutate(
    pct_change     = (POST - PRE) / PRE * 100,
    variable_label = config$variables[variable],
    variable_label = factor(variable_label, levels = unname(config$variables)),
    group          = factor(group, levels = c("PLA", "CRE"))
  ) %>%
  filter(!is.na(pct_change))

cat("\n=== Percent change data check ===\n")
cat("Rows:", nrow(pheno_pctchange), "\n")
print(pheno_pctchange %>% count(variable_label, group))

# Welch's t-test per variable. estimate = diff(tt$estimate) = CRE mean - PLA
# mean (factor levels are PLA then CRE), matching the "positive = CRE > PLA"
# direction convention used in the post-hoc contrasts above.
run_pctchange_test <- function(df) {
  df <- droplevels(df)
  if (n_distinct(df$group) < 2 || any(table(df$group) < 2)) {
    return(tibble(estimate = NA_real_, statistic = NA_real_, df = NA_real_, p.value = NA_real_))
  }
  tt <- tryCatch(t.test(pct_change ~ group, data = df), error = function(e) NULL)
  if (is.null(tt)) return(tibble(estimate = NA_real_, statistic = NA_real_, df = NA_real_, p.value = NA_real_))
  tibble(
    estimate  = unname(diff(tt$estimate)),
    statistic = unname(tt$statistic),
    df        = unname(tt$parameter),
    p.value   = tt$p.value
  )
}

pctchange_results <- pheno_pctchange %>%
  group_by(variable, variable_label) %>%
  group_modify(~ run_pctchange_test(.x)) %>%
  ungroup()

cat("\n=== Percent change group comparison (Welch t-test, CRE - PLA) ===\n")
print(pctchange_results)

pctchange_summary <- pheno_pctchange %>%
  group_by(variable, variable_label, group) %>%
  summarise(
    mean_pct = mean(pct_change, na.rm = TRUE),
    sd_pct   = sd(pct_change, na.rm = TRUE),
    n        = sum(!is.na(pct_change)),
    se_pct   = sd_pct / sqrt(n),
    .groups = "drop"
  )

# Reuses fmt_line() from section 3 -- single line per facet since there's
# only one contrast here (CRE vs PLA), bolded when raw p < 0.05.
pctchange_annotations <- pctchange_results %>%
  mutate(label = fmt_line("CRE vs PLA", p.value))

p_pctchange <- ggplot(pctchange_summary, aes(x = group, y = mean_pct, fill = group)) +
  geom_col(width = 0.6, color = colors$black) +
  geom_errorbar(aes(ymin = mean_pct - se_pct, ymax = mean_pct + se_pct),
                width = 0.15, linewidth = 0.5) +
  geom_point(
    data = pheno_pctchange,
    aes(x = group, y = pct_change, fill = group),
    position = position_jitter(width = 0.1), shape = 21, size = 1.4, alpha = 0.4,
    color = colors$gray_dark, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = colors$gray_dark, linewidth = 0.4) +
  geom_richtext(
    data = pctchange_annotations,
    aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, size = 2.3, hjust = -0.05, vjust = 1.1,
    fill = NA, label.color = NA, label.padding = unit(0, "lines")
  ) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.25))) +
  facet_wrap(~ variable_label, scales = "free_y") +
  labs(title = "Percent Change (PRE \u2192 POST): CRE vs PLA", x = NULL, y = "% Change") +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(
    strip.background   = element_rect(fill = colors$gray_light, color = NA),
    panel.grid.major.y = element_line(color = colors$gray_light),
    panel.grid.minor   = element_blank(),
    legend.position    = "none",
    axis.text.x        = element_text(size = 9)
  )

ggsave(file.path(config$fig_dir, "fig_phenotypic_pctchange.pdf"),
       p_pctchange, width = plot_width, height = plot_height, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_phenotypic_pctchange.png"),
       p_pctchange, width = plot_width, height = plot_height, dpi = 300, bg = "white")

# --- 7: Combined data export ----------------------------------------------------
write_xlsx(
  list(
    Model_Results     = model_results %>% dplyr::select(-variable),
    Posthoc_Contrasts = posthoc_results %>% dplyr::select(-variable),
    ANOVA_Annotations = anova_annotations %>% dplyr::select(-variable),
    Summary_Stats     = pheno_summary %>% dplyr::select(-variable),
    Long_Data         = pheno_long %>% dplyr::select(-variable),
    PctChange_Ttest   = pctchange_results %>% dplyr::select(-variable),
    PctChange_Summary = pctchange_summary %>% dplyr::select(-variable),
    PctChange_Data    = pheno_pctchange %>% dplyr::select(-variable)
  ),
  path = file.path(config$data_dir, "phenotypic_overview_data.xlsx")
)

cat("\nSaved figures to:", config$fig_dir, "\n")
cat("Saved data to:", config$data_dir, "\n")
