# Reformat CRC pheno data

# Fig 01 Pheno -> jamovi export
# Reads the Long_Data sheet written by fig_01_phenotypic_overview.R and reshapes it
# into jamovi-ready tables (one header row, one value per cell, blanks = missing).
# Run fig_01_phenotypic_overview.R first.
#
# Outputs (04_figures/fig_01_pheno/c_data/jamovi):
#   pheno_jamovi.xlsx with three sheets:
#     wide     -- one row per participant; <var>_pre, _post, _chg, _pct per outcome.
#                 Use for RM ANOVA, ANCOVA, paired t-tests, descriptives.
#     long     -- one row per participant x timepoint; one column per outcome.
#                 Use for mixed models (e.g., the GAMLj module).
#     codebook -- what every column is.
#   pheno_jamovi_wide.csv / pheno_jamovi_long.csv -- same data; jamovi opens these
#   directly with column names intact (File > Open).

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, writexl)

# --- Configuration -----------------------------------------------------------
config <- list(
  in_file = "04_figures/fig_01_pheno/c_data/phenotypic_overview_data.xlsx",
  out_dir = "04_figures/fig_01_pheno/c_data/jamovi",
  
  # Must match config$variables in fig_01_phenotypic_overview.R (names become
  # jamovi column names; labels are matched against the Long_Data sheet)
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
  )
)

dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- 1: Load ------------------------------------------------------------------
long_raw <- read_excel(config$in_file, sheet = "Long_Data")

cat("Long_Data columns:", paste(colnames(long_raw), collapse = ", "), "\n")
cat("Rows:", nrow(long_raw), "\n")

var_lookup <- tibble(variable       = names(config$variables),
                     variable_label = unname(config$variables))

unmatched <- setdiff(unique(long_raw$variable_label), var_lookup$variable_label)
if (length(unmatched)) {
  stop("Labels in Long_Data not found in config$variables: ",
       paste(unmatched, collapse = ", "),
       "\nUpdate config$variables to match fig_01_phenotypic_overview.R.")
}

dat <- long_raw %>%
  left_join(var_lookup, by = "variable_label") %>%
  mutate(subject  = as.character(subject),
         group    = factor(group, levels = c("PLA", "CRE")),
         time     = factor(time,  levels = c("PRE", "POST")),
         variable = factor(variable, levels = names(config$variables)))

# --- 2: Wide (one row per participant) -----------------------------------------
wide <- dat %>%
  mutate(time = tolower(as.character(time))) %>%
  pivot_wider(id_cols     = c(subject, group, tissue_avail),
              names_from  = c(variable, time),
              names_sep   = "_",
              values_from = value)

for (v in names(config$variables)) {
  pre  <- wide[[paste0(v, "_pre")]]
  post <- wide[[paste0(v, "_post")]]
  wide[[paste0(v, "_chg")]] <- post - pre
  wide[[paste0(v, "_pct")]] <- ifelse(!is.na(pre) & pre > 0, 100 * (post - pre) / pre, NA_real_)
}

# Columns grouped by outcome: pre, post, chg, pct
wide_cols <- c("subject", "group", "tissue_avail",
               as.vector(t(outer(names(config$variables),
                                 c("_pre", "_post", "_chg", "_pct"), paste0))))

wide <- wide %>%
  dplyr::select(all_of(wide_cols)) %>%
  arrange(group, subject)

# --- 3: Long (one row per participant x timepoint) ------------------------------
long_jamovi <- dat %>%
  dplyr::select(subject, group, tissue_avail, time, variable, value) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  dplyr::select(subject, group, tissue_avail, time, all_of(names(config$variables))) %>%
  arrange(group, subject, time)

# --- 4: Codebook --------------------------------------------------------------------
suffix_desc <- c(
  pre  = "PRE value",
  post = "POST value",
  chg  = "change, POST \u2212 PRE",
  pct  = "% change, 100 \u00d7 (POST \u2212 PRE) / PRE"
)

codebook <- bind_rows(
  tibble(
    sheet  = c("wide + long", "wide + long", "wide + long", "long"),
    column = c("subject", "group", "tissue_avail", "time"),
    description = c(
      "Participant ID (set to ID / nominal in jamovi)",
      "Supplement group: CRE or PLA",
      "Proteomics tissue availability: both, muscle only, plasma only",
      "Timepoint: PRE or POST"
    )
  ),
  expand_grid(variable = names(config$variables), suffix = names(suffix_desc)) %>%
    transmute(sheet  = "wide",
              column = paste0(variable, "_", suffix),
              description = paste0(unname(config$variables[variable]), " \u2014 ",
                                   unname(suffix_desc[suffix]))),
  var_lookup %>%
    transmute(sheet = "long", column = variable,
              description = paste0(variable_label, " (value at the row's timepoint)"))
)

# --- 5: Write -------------------------------------------------------------------------
as_chr_factors <- function(x) mutate(x, across(where(is.factor), as.character))

wide_out <- as_chr_factors(wide)
long_out <- as_chr_factors(long_jamovi)

write_xlsx(
  list(wide = wide_out, long = long_out, codebook = codebook),
  path = file.path(config$out_dir, "pheno_jamovi.xlsx")
)
write_csv(wide_out, file.path(config$out_dir, "pheno_jamovi_wide.csv"), na = "")
write_csv(long_out, file.path(config$out_dir, "pheno_jamovi_long.csv"), na = "")

cat("\nWide:", nrow(wide_out), "participants x", ncol(wide_out), "columns\n")
print(table(wide_out$group))
cat("Long:", nrow(long_out), "rows x", ncol(long_out), "columns\n")
cat("Saved to:", config$out_dir, "\n")
