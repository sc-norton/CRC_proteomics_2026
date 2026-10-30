# =============================================================================
# table_DE_proteins_by_tissue_contrast.R
# Excel export: differentially expressed proteins (Pi score < 0.05) for each
# tissue x contrast pairing, one sheet per pairing (6 total: 2 tissues x 3
# contrasts).
#
# Columns per sheet: uniprot_id, gene, log2FC, fold change, t, p-value,
# adjusted p-value, PI score.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, writexl)

# --- 0: Config -------------------------------------------------------------
config <- list(
  limma_muscle = "03_analysis/c_data/muscle/reformat/limma_results_all_muscle.xlsx",  # >>> VERIFY <<< path
  limma_plasma = "03_analysis/c_data/plasma/reformat/limma_results_all_plasma.xlsx",  # >>> VERIFY <<< path
  out_dir      = "04_figures/data/table_DE_proteins",
  sig_col      = "Pi.Val",
  sig_thresh   = 0.05
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

# Contrast sheet names as they appear in each limma workbook, paired with
# the display label used elsewhere in the pipeline (volcano/GSEA figures) --
# kept in this exact order (PRE v CTL, POST v CTL, POST v PRE) to match
# that established convention.
# >>> VERIFY <<< sheet identifiers -- some scripts in this pipeline have
# used sheet NAMES ("Baseline_SURVvCTL") and others sheet INDICES (11/12/13)
# for these same three contrasts; confirm which your actual workbook uses
# and adjust the `sheet = ` arguments in Section 1 accordingly if names
# don't match.
contrasts <- tibble::tribble(
  ~sheet_name,          ~label,
  "Baseline_SURVvCTL",  "PRE_v_CTL",
  "Training_SURVvCTL",  "POST_v_CTL",
  "Training_SURV",      "POST_v_PRE"
)

# =============================================================================
# 1: Load + filter -- one data frame per tissue x contrast pairing
# =============================================================================

# >>> VERIFY <<< uniprot ID column name -- this pipeline's own scripts have
# used BOTH "uniprot_id" and bare "uniprot" for what appears to be the same
# field across different files. Confirm which your actual limma sheets use
# and adjust `uniprot_col` below (or the two literal names inside
# build_de_table()) if it doesn't match.
uniprot_col <- "uniprot_id"

build_de_table <- function(xlsx_path, sheet_name) {
  df <- read_excel(xlsx_path, sheet = sheet_name)
  
  df %>%
    filter(.data[[config$sig_col]] < config$sig_thresh) %>%
    transmute(
      uniprot_id      = .data[[uniprot_col]],
      gene            = gene,
      log2FC          = logFC,
      `fold change`   = 2^logFC,
      t               = t,
      `p-value`       = P.Value,
      `adjusted p-value` = adj.P.Val,
      `PI score`      = .data[[config$sig_col]]
    ) %>%
    arrange(`PI score`)
}

tissue_files <- list(muscle = config$limma_muscle, plasma = config$limma_plasma)

de_tables <- list()
for (tissue in names(tissue_files)) {
  for (i in seq_len(nrow(contrasts))) {
    sheet_nm <- contrasts$sheet_name[i]
    label    <- contrasts$label[i]
    key      <- paste0(tissue, "_", label)
    
    de_tables[[key]] <- build_de_table(tissue_files[[tissue]], sheet_nm)
    
    cat(sprintf("%-25s: %d DE proteins (%s < %.2f)\n",
                key, nrow(de_tables[[key]]), config$sig_col, config$sig_thresh))
  }
}

# =============================================================================
# 2: Export -- 6 sheets, Excel sheet-name length capped at 31 chars
# =============================================================================
# muscle_POST_v_CTL etc. are all well under 31 chars as-is, so no truncation
# needed here -- but names() below still runs through a length guard in
# case contrast labels change later.

sheet_names <- names(de_tables)
if (any(nchar(sheet_names) > 31)) {
  stop("One or more sheet names exceed Excel's 31-character limit: ",
       paste(sheet_names[nchar(sheet_names) > 31], collapse = ", "))
}

write_xlsx(de_tables, path = file.path(config$out_dir, "DE_proteins_by_tissue_contrast.xlsx"))

cat(sprintf("\nSaved to %s\n", file.path(config$out_dir, "DE_proteins_by_tissue_contrast.xlsx")))