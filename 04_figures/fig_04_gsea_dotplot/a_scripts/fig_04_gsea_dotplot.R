# GSEA Dotplot Composite -- RECONSTRUCTED from a rendered PDF, not a
# recovered original script. Every choice below that isn't directly
# verifiable from the PDF (exact ranking metric for "top 20", y-axis sort
# order, dot size scale bounds, wrap width, whether the pooled redundancy
# filter was applied to this specific figure) is flagged inline as an
# assumption. Check these against your memory of the original before
# treating this as final.
#
# Structure read off the PDF:
#   - 3 panels: A) Muscle, top 20 pathways ranked by significance in the
#     "Cancer" (CRC vs CTL) comparison; B) Muscle, top 20 ranked by
#     significance in "Training" (CRC+RET vs CRC); C) Plasma, ALL
#     significant pathways (no top-N cap).
#   - Each panel: y = pathway description, x = comparison (Cancer /
#     Residual / Training, displayed as CRC vs CTL / CRC+RET vs CTL /
#     CRC+RET vs CRC), fill = NES, size = -log10(adj. p), black dot
#     outline only when significant in THAT comparison (thin/no outline
#     otherwise). Rows grouped into boxed strips by database, box height
#     proportional to pathway count, empty databases simply absent.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, scales,
               ggplot2, cowplot, forcats, writexl, ggh4x)

# --- Configuration -----------------------------------------------------------
# ADJUST fig_dir/data_dir to match your actual project structure -- these
# are placeholders, since the original config wasn't available to recover.
config <- list(
  # NOTE: "gsea_pathway_data.rds" -- NOT "_consolidated.rds". fig_02.R and
  # fig_02_supp.R both reference the wrong (nonexistent) "_consolidated"
  # filename; this uses what 02_gsea_pathway_consolidation.R actually saves.
  gsea_data = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  fig_dir   = "04_figures/fig_04_gsea_dotplot/b_figures",
  data_dir  = "04_figures/fig_04_gsea_dotplot/c_data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ----------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

# --- Shared constants --------------------------------------------------------
gsea_group  <- "SURV"  # matches the group column added by consolidation script
COMPARISONS <- c("Cancer", "Residual", "Training")
COMP_LABELS <- c(
  Cancer   = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training = "CRC+RET vs CRC"
)

DB_LEVELS <- c("Hallmark", "GO:BP", "GO:CC", "GO:MF", "KEGG", "Reactome")

TOP_N_A <- 20
TOP_N_B <- 20

# Assumption: dot-size legend bounds match the PDF's ticks (0, 5, 10, 15).
# oob = squish means any -log10(p) beyond 15 gets clamped to the largest
# dot rather than warning/dropping -- reasonable for a handful of very
# significant outliers, but check none of your real values blow well past
# 15 in a way that would be misleading to visually cap.
SIZE_LIMITS <- c(0, 15)
SIZE_BREAKS <- c(0, 5, 10, 15)
NES_LIMITS  <- c(-3, 3)

WRAP_WIDTH <- 32  # y-axis label wrap width, in characters -- tune by eye

# --- 1: Load Data ------------------------------------------------------------
if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run 02_gsea_pathway_consolidation.R first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

gsea_supp_full        <- gsea_pathway_data$gsea_supp_full        %>% filter(group == gsea_group)
representative_lookup <- gsea_pathway_data$representative_lookup %>% filter(group == gsea_group)

# --- 2: Pathway selection -----------------------------------------------------
# Restricts to the pooled-redundancy-filter survivors (representative_lookup)
# before ranking/selecting -- ASSUMPTION: applied uniformly to all 3 panels
# for consistency with how "significant pathway" is defined everywhere else
# in this pipeline (fig_02.R's GSEA_COUNT_MODE = "representative" default).
# If this figure was actually built on raw (unfiltered) significance,
# dropping the semi_join() below reverts to that.
select_pathways <- function(tissue_name, rank_contrast = NULL, top_n = Inf,
                            comparisons = COMPARISONS) {
  
  base <- gsea_supp_full %>%
    filter(tissue == tissue_name, comparison %in% comparisons,
           database != "GO Slim") %>%
    semi_join(representative_lookup, by = c("tissue", "database", "ID"))
  
  if (!is.null(rank_contrast)) {
    # Rank by significance (ascending p.adjust) WITHIN rank_contrast only,
    # among pathways significant in that comparison -- then pull every
    # comparison's row for those pathway IDs so all 3 dots can be drawn.
    keep <- base %>%
      filter(comparison == rank_contrast, significant) %>%
      arrange(p.adjust) %>%
      distinct(database, ID, Description_clean) %>%
      slice_head(n = top_n)
  } else {
    # "All significant pathways": kept if significant in ANY comparison.
    keep <- base %>%
      filter(significant) %>%
      distinct(database, ID, Description_clean)
  }
  
  base %>% semi_join(keep, by = c("database", "ID"))
}

pathways_A <- select_pathways("Muscle", rank_contrast = "Cancer",   top_n = TOP_N_A)
pathways_B <- select_pathways("Muscle", rank_contrast = "Training", top_n = TOP_N_B)
pathways_C <- select_pathways("Plasma", rank_contrast = NULL,       top_n = Inf)

# --- 3: Dot plot builder ------------------------------------------------------
make_gsea_dotplot <- function(df, title, subtitle, legend = FALSE) {
  
  df <- df %>%
    mutate(
      comparison  = factor(comparison, levels = COMPARISONS),
      database    = factor(database, levels = DB_LEVELS[DB_LEVELS %in% unique(database)]),
      neglog10_p  = -log10(pmax(p.adjust, 1e-300)),
      significant = factor(significant, levels = c(TRUE, FALSE))
    )
  
  # y-axis order: most significant (lowest p.adjust across whichever
  # comparisons this pathway appears in) at the top, within each database
  # facet -- ASSUMPTION: exact sort key not verifiable from the PDF alone.
  order_key <- df %>%
    group_by(database, Description_clean) %>%
    summarise(min_p = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    arrange(database, min_p)
  
  df <- df %>%
    mutate(Description_clean = factor(Description_clean,
                                      levels = rev(order_key$Description_clean)))
  
  p <- ggplot(df, aes(x = comparison, y = Description_clean)) +
    geom_point(aes(fill = NES, size = neglog10_p, colour = significant),
               shape = 21, stroke = 0.6) +
    scale_fill_gradient2(
      low = colors$blue, mid = "white", high = colors$red, midpoint = 0,
      limits = NES_LIMITS, oob = scales::squish, name = "NES"
    ) +
    scale_size_continuous(
      range = c(0.5, 7), limits = SIZE_LIMITS, breaks = SIZE_BREAKS,
      name = expression(-log[10]*"(adj. "*italic(p)*")")
    ) +
    scale_colour_manual(
      values = c(`TRUE` = "black", `FALSE` = "grey75"), guide = "none"
    ) +
    scale_x_discrete(labels = COMP_LABELS) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = WRAP_WIDTH)) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      strip.background   = ggh4x::element_part_rect(
        side = "l", colour = "black", linewidth = 1.1, fill = "white"
      ),
      strip.text.y       = element_text(angle = 0, face = "bold", size = 7,
                                        margin = margin(l = 4, r = 2)),
      strip.placement     = "outside",
      axis.text.x         = element_text(angle = 45, hjust = 1, face = "bold", size = 8),
      axis.text.y         = element_text(size = 7),
      panel.spacing        = unit(0.3, "lines"),
      panel.grid.major.x  = element_blank(),
      panel.grid.major.y  = element_line(colour = "grey85", linewidth = 0.3),
      panel.grid.minor    = element_blank(),
      panel.border        = element_blank(),
      plot.subtitle        = element_text(size = 8, face = "italic", colour = colors$gray_dark)
    )
  
  if (!legend) p <- p + guides(fill = "none", size = "none")
  
  p
}

p_A <- make_gsea_dotplot(pathways_A,
                         title    = "GSEA \u2014 Muscle",
                         subtitle = sprintf("Top %d pathways in CRC vs CTL", TOP_N_A))

p_B <- make_gsea_dotplot(pathways_B,
                         title    = "GSEA \u2014 Muscle",
                         subtitle = sprintf("Top %d pathways in CRC+RET vs CRC", TOP_N_B))

p_C <- make_gsea_dotplot(pathways_C,
                         title    = "GSEA \u2014 Plasma",
                         subtitle = "All significant pathways")

# --- 4: Shared legend ---------------------------------------------------------
# Single NES + size legend shared across all 3 panels (matches the PDF,
# which shows one legend block rather than one per panel).
get_legend_safe <- function(p) {
  lg <- tryCatch(cowplot::get_legend(p), error = function(e) NULL)
  if (!is.null(lg)) return(lg)
  gt <- ggplotGrob(p)
  leg_idx <- which(grepl("guide-box", gt$layout$name))
  if (length(leg_idx) == 0) {
    warning("get_legend_safe(): no guide-box grob found -- legend will be blank.")
    return(grid::nullGrob())
  }
  gt$grobs[[leg_idx[1]]]
}

legend_plot <- make_gsea_dotplot(pathways_C,
                                 title = NULL, subtitle = NULL, legend = TRUE) +
  theme(legend.position = "right",
        legend.box      = "vertical",
        legend.title    = element_text(size = 8, face = "bold"),
        legend.text     = element_text(size = 7))

shared_legend <- get_legend_safe(legend_plot)

# --- 5: Assemble & Save -------------------------------------------------------
# Relative panel widths eyeballed from the PDF (A and B roughly equal, C
# narrower given fewer rows/shorter database list) -- tune once rendered.
top_row <- plot_grid(
  p_A, p_B, p_C,
  labels     = c("A", "B", "C"),
  nrow       = 1,
  rel_widths = c(0.75, 0.75, 0.7)
)

fig_gsea_dotplot <- plot_grid(
  top_row, shared_legend,
  nrow        = 1,
  rel_widths  = c(1, 0.08)
)

ggsave(file.path(config$fig_dir, "fig_gsea_dotplot_composite.pdf"),
       fig_gsea_dotplot, width = 16, height = 11, units = "in", device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_gsea_dotplot_composite.png"),
       fig_gsea_dotplot, width = 16, height = 11, units = "in", dpi = 300)

ggsave(file.path(config$fig_dir, "fig_gsea_dotplot_A.pdf"), p_A, width = 6,   height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_gsea_dotplot_B.pdf"), p_B, width = 6,   height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_gsea_dotplot_C.pdf"), p_C, width = 4.5, height = 8,  device = cairo_pdf)

# --- 6: Save underlying data ---------------------------------------------------
write_xlsx(
  list(
    A_Muscle_CRCvCTL_top20     = pathways_A %>% dplyr::select(-tissue, -group),
    B_Muscle_CRCRETvCRC_top20 = pathways_B %>% dplyr::select(-tissue, -group),
    C_Plasma_all_significant   = pathways_C %>% dplyr::select(-tissue, -group)
  ),
  path = file.path(config$data_dir, "fig_gsea_dotplot_composite_data.xlsx")
)

message("Saved figure + data to: ", config$fig_dir)