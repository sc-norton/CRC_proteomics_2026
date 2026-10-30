# =============================================================================
# fig_05_stitch.R
# Master orchestrator for Figure 05: loads shared data ONCE, builds all four
# panels via their build_panel_*() functions, composes a 380mm-wide 2x2
# composite with one shared legend and uniform tags, and saves it.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, ggplot2, patchwork, cowplot)

source("04_figures/a_scripts/fig_05_style.R")

config <- list(
  limma_muscle_xlsx = "03_analysis/c_data/muscle/limma/limma_results_all_muscle_cv.xlsx",
  gsea_data         = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  gsea_muscle_dir   = "03_analysis/c_data/muscle/GSEA",
  fry_rds           = "03_analysis/c_data/muscle/fry_reversal/fry_reversal_signature_muscle.RDS",
  sankey_script     = "04_figures/a_scripts/fig_05_panel_D_heatmap_sankey.R",
  fig_dir           = "04_figures/fig_05",
  data_dir          = "04_figures/fig_05/data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1: Load shared data ONCE
# =============================================================================
# >>> VERIFY <<< sheet identifiers -- switched to name-based (matching the
# fry script's original convention) rather than the index-based
# sheet = 11 / sheet = 13 used in the prior fig_05b.R. Revert to indices
# below if your workbook's tabs aren't named exactly this way.

limma_cancer      <- read_excel(config$limma_muscle_xlsx, sheet = "Cancer")
limma_training <- read_excel(config$limma_muscle_xlsx, sheet = "Training")

if (!file.exists(config$gsea_data)) stop("Consolidated GSEA data not found at ", config$gsea_data)
gsea_pathway_data <- readRDS(config$gsea_data)

load_gsea_rds <- function(dir) {
  files <- list.files(dir, pattern = "^(GSEA|ORA)_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^(GSEA|ORA)_results_|\\.rds$", "", basename(files))
  setNames(lapply(files, readRDS), nms)
}
gsea_muscle <- load_gsea_rds(config$gsea_muscle_dir)

fry_reversal <- readRDS(config$fry_rds)

# =============================================================================
# 2: Composite geometry -- 380mm wide, A:B = 2:1, C:D = 1:2
# Explicit header bands reserve space for tag+title+subtitle ABOVE each
# content row, rather than drawing text on top of panel content.
# =============================================================================

COMP_W <- 380
A_W    <- round(COMP_W * 2/3); B_W <- COMP_W - A_W
C_W    <- round(COMP_W * 1/3); D_W <- COMP_W - C_W

HEADER_H <- 18   # tag + title + subtitle band, above each content row
LEGEND_H <- 14
ROW1_H   <- 150  # A/B content height
ROW2_H   <- 150  # C/D content height
COMP_H   <- HEADER_H * 2 + LEGEND_H + ROW1_H + ROW2_H   # 360

TAG_SZ <- round(scale_text(BASE_TAG, COMP_W))
TTL_SZ <- TAG_SZ
SUB_SZ <- round(TAG_SZ * 0.7)

# =============================================================================
# 3: Build panels
# =============================================================================

source("04_figures/a_scripts/fig_05_panel_A_protein.R")
res_A <- build_panel_A(limma_cancer, limma_training, gsea_muscle, 
                       panel_width_mm = A_W, ora_outside = list(rev_up = "NADH\nDehydrogenase Assembly"))

source("04_figures/a_scripts/fig_05_panel_B_pathway.R")
res_B <- build_panel_B(gsea_pathway_data, panel_width_mm = B_W)

source("04_figures/a_scripts/fig_05_panel_C_fry.R")
res_C <- build_panel_C(limma_cancer, limma_training, fry_reversal, panel_width_mm = C_W)

sankey_env <- new.env()
sankey_env$PANEL_SOURCE_MODE <- TRUE
sankey_env$PANEL_WIDTH_MM    <- D_W
sankey_env$scale_text        <- scale_text
source(config$sankey_script, local = sankey_env)

# Strip Panel D's internal title/subtitle here rather than editing the
# sankey script -- fig_ribbon_panel is a single ggplot (not a nested
# patchwork), so plain + works fine, no & needed.
p_D <- sankey_env$fig_ribbon_panel + labs(title = NULL, subtitle = NULL)

# Rebuild the same stat string the sankey script would have put in its own
# subtitle, from objects still present in sankey_env after sourcing.
sub_D <- sprintf(
  "%d / %d proteins (%.0f%%) mapped to %d pathway groups representing %d GO terms",
  dplyr::n_distinct(sankey_env$edges$gene), sankey_env$n_prot,
  100 * dplyr::n_distinct(sankey_env$edges$gene) / sankey_env$n_prot,
  sankey_env$n_terms, nrow(sankey_env$ora_df_final)
)

# =============================================================================
# 4: Shared A/B legend
# =============================================================================

shared_legend_grob <- build_shared_legend(pal_shared, panel_width_mm = (A_W + B_W) / 2)

# =============================================================================
# 5: Compose
# =============================================================================

layout <- paste(
  "############",   # header 1 -- tag/title/subtitle for A/B
  "AAAAAAAABBBB",   # A/B content
  "############",   # shared legend band
  "############",   # header 2 -- tag/title/subtitle for C/D
  "CCCCDDDDDDDD",   # C/D content
  sep = "\n"
)

fig_raw <- wrap_elements(full = res_A$plot_composite) +
  wrap_elements(full = res_B$plot_composite) +
  wrap_elements(full = res_C$plot_composite) +
  wrap_elements(full = p_D) +
  plot_layout(design = layout,
             widths  = rep(1, 12),
             heights = c(HEADER_H, ROW1_H, LEGEND_H, HEADER_H, ROW2_H))

# --- Cumulative mm-from-top for each band, converted to ggdraw's y-fraction
# (y = 1 at top, 0 at bottom) ------------------------------------------------
y_frac <- function(mm_from_top) 1 - mm_from_top / COMP_H

Y_TAG_AB <- y_frac(3)                      # tag/title cancer, header 1
Y_SUB_AB <- y_frac(3 + SUB_SZ * 0.6)       # subtitle, one line below
Y_LEG_BOTTOM <- y_frac(HEADER_H + ROW1_H + LEGEND_H)
LEG_H_FRAC   <- LEGEND_H / COMP_H
Y_TAG_CD <- y_frac(HEADER_H + ROW1_H + LEGEND_H + 3)
Y_SUB_CD <- y_frac(HEADER_H + ROW1_H + LEGEND_H + 3 + SUB_SZ * 0.6)

X_A <- 0.01;  X_B <- round(A_W / COMP_W, 3) + 0.02
X_C <- 0.01;  X_D <- round(C_W / COMP_W, 3) + 0.02
X_TTL <- 0.03

fig_05 <- ggdraw(fig_raw) +
  draw_plot(shared_legend_grob, x = 0.30, y = Y_LEG_BOTTOM, width = 0.4, height = LEG_H_FRAC) +
  draw_label("A", x = X_A, y = Y_TAG_AB, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_A$title,    x = X_A + X_TTL, y = Y_TAG_AB, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_A$subtitle, x = X_A + X_TTL, y = Y_SUB_AB, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("B", x = X_B, y = Y_TAG_AB, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_B$title,    x = X_B + X_TTL, y = Y_TAG_AB, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_B$subtitle, x = X_B + X_TTL, y = Y_SUB_AB, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("C", x = X_C, y = Y_TAG_CD, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_C$title,    x = X_C + X_TTL, y = Y_TAG_CD, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(res_C$subtitle, x = X_C + X_TTL, y = Y_SUB_CD, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
  draw_label("D", x = X_D, y = Y_TAG_CD, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("Protein-to-Pathway", x = X_D + X_TTL, y = Y_TAG_CD, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label(sub_D, x = X_D + X_TTL, y = Y_SUB_CD, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")

fig_05

save_scaled(fig_05, config$fig_dir, "fig_05_composite", width_mm = COMP_W, height_mm = COMP_H)

# Standalone panel saves (each with its own legend/title -- useful for review)
save_scaled(res_A$plot, config$fig_dir, "fig_05a_protein_reversal", width_mm = A_W, height_mm = ROW1_H)
save_scaled(res_B$plot, config$fig_dir, "fig_05b_pathway_reversal", width_mm = B_W, height_mm = ROW1_H)
save_scaled(res_C$plot, config$fig_dir, "fig_05c_fry_barcode",      width_mm = C_W, height_mm = ROW2_H)
save_scaled(p_D,        config$fig_dir, "fig_05d_heatmap_sankey",   width_mm = D_W, height_mm = ROW2_H)
