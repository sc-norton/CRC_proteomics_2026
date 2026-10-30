# Fig GSEA ORA Bar — top 10 pathways per contrast, ranked by significance
#
# Classic ORA-style horizontal bar chart: one panel per contrast, bars are
# the top 10 non-redundant pathways ranked by p.adjust (pooled across every
# database in gsea_pathway_data), bar length = -log10(p.adjust), bar fill =
# direction (Up/Down by NES sign). Loads the canonical consolidated GSEA
# data from 03_analysis/a_scripts/02_gsea_pathway_consolidation.R instead of
# recomputing the full-table/redundancy-filter logic that gsea_bubble_plot.R
# already centralized there.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               scales, writexl)

# --- Configuration -----------------------------------------------------------
config <- list(
  consolidated_rds = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  fig_dir          = "04_figures/fig_gsea_ora",
  data_dir         = "04_figures/fig_gsea_ora/data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

# Same Up/Down convention as DIR_COLORS in fig_02.R / NES_LOW/NES_HIGH in
# gsea_bubble_plot.R.
DIR_COLORS <- c(Up = colors$red, Down = colors$blue)

# --- 1: Load consolidated GSEA data ------------------------------------------
# One canonical object -- see 02_gsea_pathway_consolidation.R. Re-derive
# nothing here that already lives on gsea_pathway_data.
stopifnot(
  "Consolidated GSEA data not found -- run 03_analysis/a_scripts/02_gsea_pathway_consolidation.R first" =
    file.exists(config$consolidated_rds)
)
gsea_pathway_data <- readRDS(config$consolidated_rds)

gsea_supp_full         <- gsea_pathway_data$gsea_supp_full
representative_lookup  <- gsea_pathway_data$representative_lookup
adjp_thresh            <- gsea_pathway_data$params$adjp_thresh
comp_defs              <- gsea_pathway_data$params$comp_defs
comp_labels            <- gsea_pathway_data$params$comp_labels

cat("=== Consolidated GSEA data check ===\n")
cat("Rows:", nrow(gsea_supp_full), "\n")
cat("Tissues:", paste(levels(gsea_supp_full$tissue), collapse = ", "), "\n")
cat("Generated at:", format(gsea_pathway_data$params$generated_at), "\n")

TOP_N <- 10

# --- 2: Select top-N per tissue x contrast -----------------------------------
# Pooled across every database (all 7 in db_slots), restricted first to each
# database's redundancy-filter survivors (representative_lookup), then
# ranked by p.adjust across the pooled set -- same two-stage restrict-then-
# rank approach as select_pathways_per_contrast() in gsea_bubble_plot.R,
# just with TOP_N = 10 instead of 30 and no databases= subsetting (all
# databases are candidates here).
select_top_pathways <- function(tissue_name, contrast_name, top_n) {
  reps <- representative_lookup %>% filter(tissue == tissue_name)

  df <- gsea_supp_full %>%
    filter(tissue == tissue_name, comparison == contrast_name) %>%
    dplyr::group_by(database) %>%
    dplyr::group_modify(~ {
      db_name <- as.character(.y$database)
      rep_ids <- reps %>% filter(database == db_name) %>% pull(ID)
      if (length(rep_ids) == 0) return(.x)
      dplyr::filter(.x, ID %in% rep_ids)
    }) %>%
    dplyr::ungroup() %>%
    filter(p.adjust < adjp_thresh) %>%
    arrange(p.adjust)

  if (nrow(df) == 0) {
    warning("No significant pathways for ", tissue_name, " / ", contrast_name)
    return(df)
  }

  df %>% slice_head(n = top_n)
}

build_ora_selection <- function(tissue_name, contrasts) {
  purrr::map_dfr(contrasts, function(contrast_name) {
    select_top_pathways(tissue_name, contrast_name, TOP_N)
  })
}

# --- 3: Plot -------------------------------------------------------------
# One figure per tissue, one facet per contrast. pathway_key disambiguates
# same-label pathways across contrasts/databases (same trick as pathway_key
# in gsea_bubble_plot.R) and encodes rank order for a clean top-to-bottom
# bar chart within each free_y facet.
make_ora_plot <- function(tissue_name) {
  group_name <- names(comp_defs)[[1]]   # only SURV is wired up currently
  contrasts  <- comp_defs[[group_name]]

  sel <- build_ora_selection(tissue_name, contrasts)
  if (nrow(sel) == 0) stop("No data to plot for tissue = ", tissue_name)

  sel <- sel %>%
    mutate(
      comparison_label = factor(comp_labels[as.character(comparison)],
                                 levels = comp_labels[contrasts]),
      direction   = ifelse(NES > 0, "Up", "Down"),
      axis_label  = paste0(Description_clean, " (", database, ")")
    ) %>%
    arrange(comparison_label, p.adjust) %>%
    group_by(comparison_label) %>%
    mutate(pathway_key = paste(comparison_label, row_number(), sep = "::")) %>%
    ungroup()

  label_lookup <- setNames(sel$axis_label, sel$pathway_key)
  sel <- sel %>%
    mutate(pathway_key = factor(pathway_key, levels = rev(sel$pathway_key)))

  p <- ggplot(sel, aes(x = pathway_key, y = -log10(p.adjust), fill = direction)) +
    geom_col(width = 0.7) +
    coord_flip() +
    facet_wrap(~ comparison_label, scales = "free_y", ncol = length(contrasts)) +
    scale_x_discrete(labels = label_lookup) +
    scale_fill_manual(values = DIR_COLORS, name = "Direction") +
    labs(
      title = paste0("GSEA ORA — ", tissue_name, " — Top ", TOP_N,
                     " pathways per contrast"),
      x = NULL, y = expression(-log[10]~"(FDR)")
    ) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      strip.background = element_rect(fill = colors$gray_light, color = NA),
      # Default strip.clip = "inherit" clips strip text to the panel width,
      # which truncates contrast labels like "Training_SURVvCTL" when the
      # panel is narrow -- "off" lets them overflow visibly instead. The
      # rightmost facet has nowhere to overflow INTO but plot margin, so
      # pad plot.margin generously or its strip label clips against the
      # device edge instead of the panel edge.
      strip.clip        = "off",
      plot.margin       = margin(5.5, 40, 5.5, 40, unit = "pt"),
      panel.grid.major.x = element_line(color = colors$gray_light),
      panel.grid.minor.x = element_blank()
    )

  list(plot = p, data = sel)
}

# --- 4: Build, save, and export companion data --------------------------------
for (tissue_name in levels(gsea_supp_full$tissue)) {
  res <- make_ora_plot(tissue_name)

  cat("\n", tissue_name, ": ", nrow(res$data), " bars across ",
      dplyr::n_distinct(res$data$comparison_label), " contrasts\n", sep = "")

  n_contrasts <- dplyr::n_distinct(res$data$comparison_label)
  fig_width   <- 5 * n_contrasts
  fig_height  <- 6

  file_stub <- sprintf("gsea_ora_top%d_%s", TOP_N, tolower(tissue_name))

  ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
         res$plot, width = fig_width, height = fig_height, dpi = 300)
  ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
         res$plot, width = fig_width, height = fig_height, dpi = 300, bg = "white")

  print(res$plot)
}

ora_selection_data <- purrr::map_dfr(levels(gsea_supp_full$tissue), function(tissue_name) {
  make_ora_plot(tissue_name)$data %>%
    dplyr::select(tissue, comparison = comparison_label, database, ID,
                  Description = Description_clean, NES, direction,
                  pvalue, p.adjust, setSize)
})

writexl::write_xlsx(
  list(
    Muscle = filter(ora_selection_data, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma = filter(ora_selection_data, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$data_dir, sprintf("gsea_ora_top%d_percontrast_selection.xlsx", TOP_N))
)

cat("\nSaved GSEA ORA top-", TOP_N, " figure(s) to: ", config$fig_dir,
    "\nSaved companion selection data to: ", config$data_dir, "\n", sep = "")
