# =============================================================================
# fig_05_stitch.R
# Master orchestrator for Figure 05: loads shared data ONCE, builds all four
# panels via their build_panel_*() functions, composes a 380mm-wide 2x2
# composite with one shared legend and uniform tags, and saves it.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, ggplot2, patchwork, cowplot, writexl)

source("04_figures/fig_05_reversal/a_scripts/fig_05_style.R")

config <- list(
  limma_muscle_xlsx = "03_analysis/c_data/muscle/limma/limma_results_all_muscle_cv.xlsx",
  gsea_data         = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  gsea_muscle_dir   = "03_analysis/c_data/muscle/GSEA",
  fry_rds           = "03_analysis/c_data/muscle/fry_reversal/fry_reversal_signature_muscle.RDS",
  sankey_script     = "04_figures/fig_05_reversal/a_scripts/fig_05_panel_D_heatmap_sankey.R",
  fig_dir           = "04_figures/fig_05_reversal/b_figures",
  data_dir          = "04_figures/fig_05_reversal/c_data"
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

source("04_figures/fig_05_reversal/a_scripts/fig_05_panel_A_protein.R")
# ora_outside is no longer passed here -- ORA_OUTSIDE_LABELS at the top of
# fig_05_panel_A_protein.R is now the single place to edit which pathway
# labels render outside their bar. build_panel_A() defaults to that list.
res_A <- build_panel_A(limma_cancer, limma_training, gsea_muscle,
                       panel_width_mm = A_W)

source("04_figures/fig_05_reversal/a_scripts/fig_05_panel_B_pathway.R")
res_B <- build_panel_B(gsea_pathway_data, panel_width_mm = B_W)

source("04_figures/fig_05_reversal/a_scripts/fig_05_panel_C_fry.R")
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

# =============================================================================
# 6: Export the data behind every panel to a single Excel workbook
# One sheet per panel (A/B/C/D), plus a README sheet describing them. Panel A
# gets two sheets (protein-level + the ORA run on each of its 4 quadrants);
# Panel D's sheet is sourced from fig_05_panel_D_heatmap_sankey.R's own
# sheet1 (built in sankey_env, unconditionally -- see that script's section 15).
# =============================================================================

fry_up   <- filter(res_C$fry_summary, set == "signature_up")
fry_down <- filter(res_C$fry_summary, set == "signature_down")

panel_c_export <- res_C$t_rank |>
  mutate(
    fry_up_n              = fry_up$n,
    fry_up_PValue         = fry_up$PValue,
    fry_up_PValue_Mixed   = fry_up$PValue.Mixed,
    fry_up_consistent     = fry_up$consistent,
    fry_down_n            = fry_down$n,
    fry_down_PValue       = fry_down$PValue,
    fry_down_PValue_Mixed = fry_down$PValue.Mixed,
    fry_down_consistent   = fry_down$consistent,
    cor_within            = fry_reversal$cor_within[1]
  )

readme_df <- tibble::tibble(
  Sheet = c("Panel_A_Proteins", "Panel_A_ORA_by_Quadrant", "Panel_B_Pathways",
            "Panel_C_FryReversal", "Panel_D_Protein_Pathway"),
  Panel = c("A", "A", "B", "C", "D"),
  Description = c(
    paste("One row per protein quantified in both contrasts -- the Panel A protein-reversal scatter.",
          "logFC/Pi.Val/adj.P.Val for CRC vs CTL (_CRC) and CRC+RET vs CRC (_RET), plus the reversal",
          "significance category (sig: Both / CRC v CTL only / CRC+RET v CRC only / Neither)."),
    paste("GO Biological Process ORA results (jaccard-deduplicated, full set) run separately on the",
          "significant proteins in each of Panel A's four scatter quadrants -- Reversed (Up/Down RET)",
          "and Exacerbated (Up/Down). The flanking bar charts in Panel A plot only the top 6 terms",
          "per quadrant; this sheet has every term, tagged by the 'quadrant' column."),
    paste("GSEA pathway-level (NES) reversal data underlying Panel B: NES and p.adjust for",
          "CRC vs CTL (_PRE) and CRC+RET vs CRC (_RET) for each GO Slim / Hallmark pathway,",
          "restricted to pathways significant (p.adjust < 0.05) in at least one contrast."),
    paste("Proteins ranked by CRC+RET vs CRC t-statistic, with the running enrichment-score trace",
          "(es_up/es_down) used for the Panel C fry barcode plot. fry_up_*/fry_down_* columns carry",
          "the fry rotation-test result (p-value, mixed p-value, n, consistency) for the",
          "significant-UP and significant-DOWN signature sets; cor_within is the within-subject",
          "correlation used by the test."),
    paste("One row per significant protein underlying Panel D's heatmap + Sankey: logFC across all",
          "three contrasts, trajectory classification (Reversed/Attenuated/Exacerbated/Unclassified),",
          "significance pattern across contrasts, and which curated pathway groups / GO terms",
          "(from the Panel D ORA) it maps to.")
  )
)

write_xlsx(
  list(
    "README"                  = readme_df,
    "Panel_A_Proteins"        = res_A$df_rev,
    "Panel_A_ORA_by_Quadrant" = res_A$ora_quadrants,
    "Panel_B_Pathways"        = res_B$nes_df |> mutate(across(where(is.factor), as.character)),
    "Panel_C_FryReversal"     = panel_c_export,
    "Panel_D_Protein_Pathway" = sankey_env$sheet1 |> mutate(across(where(is.factor), as.character))
  ),
  path = file.path(config$data_dir, "fig_05_reversal_panel_data.xlsx")
)

message("Panel data workbook saved to: ",
        file.path(config$data_dir, "fig_05_reversal_panel_data.xlsx"))