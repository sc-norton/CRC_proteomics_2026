# Volcano + Heatmaps

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
# Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, ggvenn, patchwork, ggtext, cowplot, ggh4x, ggalluvial,
               scales, gtsummary, limma, Hmisc, nlme, ggpubr, tidytext, ggsankey,
               enrichplot)

# Configuration
config <- list(
  # Data files
  limma_muscle = "03_analysis/c_data/muscle/limma",
  limma_plasma = "03_analysis/c_data/plasma/limma",
  # Directories
  fig_dir = "04_figures/fig_03_volcano_heat/b_figures",
  data_dir = "04_figures/fig_03_volcano_heat/c_data"
)

dir.create(config$fig_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

# Color palette + theme
source("04_figures/a_scripts/00_theme.R")
source("04_figures/fig_03_volcano_heat/a_scripts/fig_03_style.R")
colors$tissue <- c(muscle = "#F2C6CF", plasma = "#A9D6CF")  # script-specific

make_strip <- function(fills, type = "x") {
  n <- length(fills)
  if (type == "x") {
    strip_themed(
      background_x = elem_list_rect(fill = fills, linewidth = rep(0, n)),
      text_x       = elem_list_text(face = rep("bold", n))
    )
  } else if (type == "xy") {
    strip_themed(
      background_x = elem_list_rect(fill = fills$x, linewidth = rep(0, length(fills$x))),
      background_y = elem_list_rect(fill = fills$y, linewidth = rep(0, length(fills$y))),
      text_x       = elem_list_text(face = rep("bold", length(fills$x))),
      text_y       = elem_list_text(face = rep("bold", length(fills$y)))
    )
  }
}

# Labels
# --- 0e: Comparison Labels ---------------------------------------
comp_labels <- c(
  "Cancer"     = "CRC vs CTL",
  "Residual"     = "CRC+RET vs CTL",
  "Training"         = "CRC+RET vs CRC",
  "Baseline_supp"      = "CRE_PRE vs PLA _RE",
  "Post_training_supp"      = "CRE_POST vs PLA_POST",
  "Training_CRE"          = "CRE_POST vs CRE_PRE",
  "Training_PLA"          = "PLA_POST vs PLA_PRE",
  "Interaction_supp"      = "Interaction\nCRE v PLA",
  "UP_both"               = "UP in CRC at PRE & POST",
  "DOWN_both"             = "DOWN in CRC at PRE & POST",
  "UP_PRE_only"           = "UP in CRC PRE only",
  "UP_POST_only"          = "UP in CRC POST only",
  "DOWN_PRE_only"         = "DOWN in CRC PRE only",
  "DOWN_POST_only"        = "DOWN in CRC POST only"
)

surv_labels <- c(
  "Cancer" = "CRC vs CTL",
  "Residual" = "CRC+RET vs CTL",
  "Training"     = "CRC+RET vs CRC"
)


# --- 1: Load Data --------------------------
limma_results <- list(
  muscle = list(
    Cancer = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle_cv.xlsx"), 
                                   sheet = "Cancer"),
    Residual = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle_cv.xlsx"),
                                   sheet = "Residual"),
    Training     = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle_cv.xlsx"),
                                   sheet = "Training")
  ),
  plasma = list(
    Cancer = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"), 
                                   sheet = "Cancer"),
    Residual = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"),
                                   sheet = "Residual"),
    Training    = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"),
                                   sheet = "Training")
  ))

# --- 2: Restructure Data -------------------
# --- 2a: Volcano data ---
df_vol <- bind_rows(
  lapply(c("Muscle", "Plasma"), function(tis) {
    res <- if (tis == "Muscle") limma_results$muscle else limma_results$plasma 
    bind_rows(lapply(names(res), function(nm) {
      res[[nm]] %>% mutate(comparison = nm, tissue = tis)
    }))
  })
) %>% mutate(tissue = factor(tissue, levels = c("Muscle", "Plasma")))

# --- 3: Volcano plot function ---------------
make_volcano_plot <- function(df_vol, selected_comparisons,
                              n_labels = 5, pval_thresh = 0.05, lfc_thresh = 0.6) {
  df <- df_vol %>%
    filter(comparison %in% selected_comparisons) %>%
    mutate(comparison = factor(comparison, levels = selected_comparisons))
  
  df_label <- df %>%
    mutate(nudge_x = ifelse(DE.pi == "UP", 2.5, -2.5), nudge_y = 0.5) %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, tissue, DE.pi) %>%
    slice_max(abs(logFC), n = n_labels) %>%
    ungroup()
  
  # Count UP and DOWN per facet
  df_counts <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, tissue, DE.pi) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      x     = ifelse(DE.pi == "UP",  Inf, -Inf),
      y     = Inf,
      hjust = ifelse(DE.pi == "UP",  1.2, -0.2),
      label = as.character(n),
      color = ifelse(DE.pi == "UP", colors$red, colors$blue)
    )
  
  # Background shading per tissue
  tissue_bg <- df %>%
    distinct(tissue) %>%
    mutate(fill_color = unname(colors$tissue[tissue]))
  
  strip <- make_strip(
    fills = list(x = rep(colors$light_gray, length(selected_comparisons)),
                 y = unname(colors$tissue)),
    type  = "xy"
  )
  
  ggplot(df, aes(logFC, -log10(P.Value), color = DE.pi)) +
    geom_rect(
      data        = tissue_bg,
      aes(fill    = fill_color),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      alpha       = 0.15,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_text_repel(data = df_label, aes(label = gene),
                    color = "black", fontface = "bold", size = 3,
                    show.legend = FALSE, max.overlaps = Inf,
                    box.padding = 0.4, point.padding = 0.3,
                    force = 2, force_pull = 0.01,
                    min.segment.length = 0.05, segment.color = "black",
                    arrow = arrow(length = unit(0.01, "npc")), seed = 42) +
    geom_text(
      data        = df_counts,
      aes(x = x, y = y, label = label, hjust = hjust),
      vjust       = 1.5,
      fontface    = "bold",
      size        = 3,
      color       = df_counts$color,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values   = c(UP = colors$red, DOWN = colors$blue),
                       na.value = colors$gray) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_cartesian(xlim = c(-4, 4)) +
    facet_grid2(tissue ~ comparison, strip = strip,
                labeller = labeller(comparison = comp_labels)) +
    labs(x = bquote(bold("log2(Fold Change)")), y = bquote(bold("-log10(p-value)")),
         color = "Direction") +
    theme_cr(base_size = 10)
}
# Build volcano
volcano_all <- make_volcano_plot(df_vol, c("Cancer", "Residual", "Training"))

ggsave(file.path(config$fig_dir, "panel_3A.pdf"), volcano_all, device = "pdf")
ggsave(file.path(config$fig_dir, "panel_3A.png"), volcano_all, device = "png")


# ============================================================
# --- 4: Heatmap data (top 20 up/down significant per contrast) ---
# ============================================================
# Reuses df_vol directly rather than reloading — same rows that feed
# the volcano panel feed these heatmaps, so selection stays consistent
# with what's plotted in 3A.
#
# For EACH of the 3 SURV contrasts x EACH tissue, selects the top 20
# UP and top 20 DOWN significant proteins from that contrast/tissue
# (DE.pi != "NDE", ranked by Pi.Val), then shows those same proteins
# across ALL 3 contrasts on the x-axis — so the heatmap sits next to
# the volcano panel for that same contrast/tissue and shows how the
# selected proteins behave everywhere, not just where they were picked.
#
# Two versions per contrast x tissue: logFC (raw magnitude) and
# row-wise z-score (normalizes each protein's own pattern across the
# 3 contrasts, independent of how large its logFC is).

pacman::p_load(writexl)

heat_config <- list(
  contrasts = c("Cancer", "Residual", "Training"),
  n_top     = 20
)

# Ensure Pi.Val is present (DE.pi implies it was derived from this, but
# compute defensively in case a sheet is missing the column)
if (!"Pi.Val" %in% names(df_vol)) {
  df_vol <- df_vol %>% mutate(Pi.Val = P.Value ^ abs(logFC))
}

# Find a stable protein ID column (beyond gene symbol) to disambiguate
# proteins that share a gene symbol (isoforms), and to join against
# dal$data rownames (UniProt accessions like "P02144").
df_vol <- df_vol %>%
  mutate(protein_key = uniprot)

# --- 4a: Group-level abundance (for z-score heatmaps) ----------------
# Z-score should be by GROUP (CTL, CRC_PRE, CRC_POST), not by contrast —
# this needs the normalized/imputed abundance matrix + sample metadata,
# not just the limma results, so it's a separate data source from df_vol.
#
# TODO: confirm these paths match your normalized DAList location (memory
# suggests something under 01_normalization/c_data/SURV/[stringency]/
# [BM_threshold]/) and that dalist$data / dalist$metadata, group_col, and
# group_levels below match your actual proteoDA DAList structure.

dalist_paths <- list(
  Muscle = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS",
  Plasma = "01_normalization/c_data/plasma/02_normalized_data_plasma.RDS"
)

group_col    <- "cancer_time"
group_levels <- c("CTL_T1", "CRC_T1", "CRC_T2")
group_labels <- c(CTL_T1 = "CTL", CRC_T1 = "CRC", CRC_T2 = "CRC+RET")

dalists <- list(
  Muscle = readRDS(dalist_paths$Muscle),
  Plasma = readRDS(dalist_paths$Plasma)
)

#' Select top n_top up + top n_top down significant proteins from one
#' contrast x tissue combo, ranked by Pi.Val (ascending = most significant)
select_top_proteins <- function(df_vol, selection_comparison, selection_tissue,
                                n_top = 20) {
  sel <- df_vol %>%
    filter(comparison == selection_comparison, tissue == selection_tissue,
           DE.pi != "NDE")
  
  up   <- sel %>% filter(DE.pi == "UP")   %>% arrange(Pi.Val) %>% slice_head(n = n_top)
  down <- sel %>% filter(DE.pi == "DOWN") %>% arrange(Pi.Val) %>% slice_head(n = n_top)
  
  if (nrow(up) < n_top) {
    message(sprintf("[%s / %s] Only %d significant UP proteins available (requested %d)",
                    selection_tissue, selection_comparison, nrow(up), n_top))
  }
  if (nrow(down) < n_top) {
    message(sprintf("[%s / %s] Only %d significant DOWN proteins available (requested %d)",
                    selection_tissue, selection_comparison, nrow(down), n_top))
  }
  
  bind_rows(
    up   %>% mutate(direction = factor("Up",   levels = c("Up", "Down"))),
    down %>% mutate(direction = factor("Down", levels = c("Up", "Down")))
  )
}

#' Select proteins significant in BOTH contrast_a (default:
#' Training_SURVvCTL, "CRC+RET vs CTL") and contrast_b (default:
#' Training_SURV, "CRC+RET vs CRC") for one tissue.
#'
#' Unlike select_top_proteins(), there's no top-N cap by default (top_n =
#' Inf) — being significant in two contrasts is already a strong filter,
#' matching the "all matching proteins, no arbitrary cap" convention this
#' pipeline already uses elsewhere (build_presence_sets(), da_only_list in
#' the GSEA script). Pass a finite top_n if the resulting set turns out too
#' large to be a useful heatmap.
#'
#' direction/logFC/Pi.Val are taken from contrast_b (Training_SURV) as the
#' canonical contrast for ranking and labeling, since that's the final
#' comparison of interest. A protein significant in both contrasts can
#' still move in OPPOSITE directions between them (e.g. up relative to CTL
#' but down relative to the CRC baseline) — concordant_direction reports
#' this rather than silently filtering it out, so you can see if any
#' discordant-direction proteins exist and decide what to do with them.
select_significant_in_both <- function(df_vol, selection_tissue,
                                       contrast_a = "Residual",
                                       contrast_b = "Training",
                                       top_n = Inf) {
  
  sig_a <- df_vol %>%
    filter(tissue == selection_tissue, comparison == contrast_a, DE.pi != "NDE") %>%
    dplyr::select(protein_key, gene, logFC_a = logFC, PiVal_a = Pi.Val, DEpi_a = DE.pi)
  
  sig_b <- df_vol %>%
    filter(tissue == selection_tissue, comparison == contrast_b, DE.pi != "NDE") %>%
    dplyr::select(protein_key, logFC_b = logFC, PiVal_b = Pi.Val, DEpi_b = DE.pi)
  
  both <- inner_join(sig_a, sig_b, by = "protein_key") %>%
    mutate(
      concordant_direction = sign(logFC_a) == sign(logFC_b),
      direction = factor(if_else(DEpi_b == "UP", "Up", "Down"), levels = c("Up", "Down")),
      logFC     = logFC_b,                # canonical logFC = contrast_b (Training_SURV)
      Pi.Val    = pmax(PiVal_a, PiVal_b)  # rank by the WORSE (larger) of the two p-values --
      # "significant in both" should mean genuinely both,
      # not just barely squeaking by on one
    )
  
  message(sprintf(
    "[%s] %d proteins significant in both %s and %s (%d concordant direction, %d discordant)",
    selection_tissue, nrow(both), contrast_a, contrast_b,
    sum(both$concordant_direction), sum(!both$concordant_direction)
  ))
  
  if (any(!both$concordant_direction)) {
    message(sprintf(
      "  NOTE: %d of those proteins move in OPPOSITE directions between %s and %s -- ",
      sum(!both$concordant_direction), contrast_a, contrast_b),
      "check the `concordant_direction` column in the exported selection sheet.")
  }
  
  if (is.finite(top_n)) {
    up   <- both %>% filter(direction == "Up")   %>% arrange(Pi.Val) %>% slice_head(n = top_n)
    down <- both %>% filter(direction == "Down") %>% arrange(Pi.Val) %>% slice_head(n = top_n)
    both <- bind_rows(up, down)
  }
  
  both
}

#' Build the display label for each selected protein — shared by both the
#' logFC and z-score builders so row order matches between the two heatmaps
build_gene_labels <- function(selected_proteins) {
  selected_proteins %>%
    distinct(protein_key, gene, direction, logFC) %>%
    mutate(gene_label = gene)
}

#' Top-to-bottom display order: strongest Up first, then strongest Down
build_gene_order <- function(sel_genes) {
  sel_genes %>%
    mutate(sort_key = if_else(direction == "Up", -logFC, logFC)) %>%
    arrange(direction, sort_key) %>%
    pull(gene_label)
}

#' Wide logFC matrix (long format for ggplot), across all 3 CONTRASTS,
#' for one tissue — pairs with the volcano panel for that contrast
build_logfc_matrix <- function(df_vol, selected_proteins, selection_tissue, contrasts) {
  
  sel_genes  <- build_gene_labels(selected_proteins)
  gene_order <- build_gene_order(sel_genes)
  
  df_vol %>%
    filter(tissue == selection_tissue, protein_key %in% selected_proteins$protein_key) %>%
    distinct(protein_key, comparison, logFC, Pi.Val) %>%
    left_join(sel_genes %>% dplyr::select(protein_key, gene_label), by = "protein_key") %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_order)),  # rev: first item plots at top
      comparison = factor(comparison, levels = contrasts),
      sig        = Pi.Val < 0.05
    )
}

#' Row-wise z-score matrix (long format for ggplot), across the 3 GROUPS
#' (CTL_T1, SURV_T1, SURV_T2) — computed from mean normalized abundance
#' per group, not from logFC. Same gene order/labels as the logFC heatmap
#' so the two line up when placed side by side.
build_group_zscore_matrix <- function(dalist, selected_proteins,
                                      group_col = "cancer_time",
                                      group_levels = c("CTL_T1", "CRC_T1", "CRC_T2")) {
  
  sel_genes  <- build_gene_labels(selected_proteins)
  gene_order <- build_gene_order(sel_genes)
  
  abund <- dalist$data      # protein x sample abundance matrix
  meta  <- dalist$metadata  # sample_id + group_col live here
  
  abund_sel <- abund[rownames(abund) %in% selected_proteins$protein_key, , drop = FALSE]
  
  group_means <- purrr::map_dfc(group_levels, function(g) {
    samp_ids <- meta$sample_id[meta[[group_col]] == g]
    samp_ids <- intersect(samp_ids, colnames(abund_sel))
    rowMeans(abund_sel[, samp_ids, drop = FALSE], na.rm = TRUE)
  })
  colnames(group_means) <- group_levels
  group_means$protein_key <- rownames(abund_sel)
  
  group_means %>%
    pivot_longer(-protein_key, names_to = "group", values_to = "abund_mean") %>%
    group_by(protein_key) %>%
    mutate(
      z = {
        v <- abund_mean
        s <- sd(v, na.rm = TRUE)
        if (is.na(s) || s == 0) rep(0, length(v)) else as.numeric(scale(v))
      }
    ) %>%
    ungroup() %>%
    left_join(sel_genes %>% dplyr::select(protein_key, gene_label), by = "protein_key") %>%
    mutate(
      gene_label = factor(gene_label, levels = rev(gene_order)),
      group      = factor(group, levels = group_levels)
    )
}

# --- 5: Heatmap plot functions ---------------------------------------

#' logFC heatmap, x-axis = contrast. Keeps the significance asterisk
#' since Pi.Val < 0.05 is a per-contrast property.
#' range_val: optional fixed color-scale range. Defaults to NULL, which
#' auto-computes from this df's own logFC values (original behavior,
#' unchanged for every existing caller). Pass an explicit range_val when
#' two heatmaps need an IDENTICAL color scale so their legends can be
#' collected into one (e.g. Muscle + Plasma side-by-side).
plot_logfc_heatmap <- function(df, title, range_val = NULL) {
  if (is.null(range_val)) range_val <- max(abs(df$logFC), na.rm = TRUE)
  ggplot(df, aes(x = comparison, y = gene_label, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      data = filter(df, sig),
      aes(label = "*"),
      color = "black", size = 3.5, vjust = 0.75
    ) +
    scale_fill_gradient2(
      low = colors$blue, mid = "white", high = colors$red,
      midpoint = 0, limits = c(-range_val, range_val), oob = scales::squish,
      name = "logFC"
    ) +
    scale_x_discrete(labels = comp_labels) +
    labs(x = NULL, y = NULL, title = title) +
    theme_cr(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      axis.text.y  = element_text(face = "italic"),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4)
    )
}

#' z-score heatmap, x-axis = group (CTL_T1, SURV_T1, SURV_T2). No
#' significance markers here — significance is a contrast-level property
#' and doesn't map cleanly onto a single group.
plot_group_zscore_heatmap <- function(df, title) {
  range_val <- max(abs(df$z), na.rm = TRUE)
  ggplot(df, aes(x = group, y = gene_label, fill = z)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_gradient2(
      low = colors$blue, mid = "white", high = colors$red,
      midpoint = 0, limits = c(-range_val, range_val), oob = scales::squish,
      name = "z-score"
    ) +
    scale_x_discrete(labels = group_labels) +
    labs(x = NULL, y = NULL, title = title) +
    theme_cr(base_size = 9) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      axis.text.y  = element_text(face = "italic"),
      panel.border = element_rect(color = "grey40", fill = NA, linewidth = 0.4)
    )
}

# ============================================================
# --- 5b: Cell sizing shared by all heatmaps in this section ---
# ============================================================
# Width now scales with number of columns instead of a fixed 3.5in,
# so cells stay narrow regardless of whether you're plotting 3 contrasts,
# 3 groups, or (with the strip) 1 membership column + N contrasts.
col_width_in  <- 0.4   # target width per data column
base_width_in <- 1.4   # room for y-axis labels + legend

# --- 6 (updated): Driver — build both heatmap versions for one tissue x contrast ---
# Only the ggsave width lines changed (fixed 3.5 -> dynamic), everything
# else is identical to before.

build_contrast_heatmap_pair <- function(df_vol, dalists, selection_comparison, selection_tissue,
                                        contrasts, n_top = 20,
                                        group_col = "cancer_time",
                                        group_levels = c("CTL_T1", "CRC_T1", "CRC_T2"),
                                        fig_dir = config$fig_dir, data_dir = config$data_dir) {
  
  selected <- select_top_proteins(df_vol, selection_comparison, selection_tissue, n_top)
  
  if (nrow(selected) == 0) {
    warning(sprintf("No significant proteins found for %s / %s — skipping.",
                    selection_tissue, selection_comparison))
    return(invisible(NULL))
  }
  
  logfc_df <- build_logfc_matrix(df_vol, selected, selection_tissue, contrasts)
  z_df     <- build_group_zscore_matrix(dalists[[selection_tissue]], selected,
                                        group_col = group_col, group_levels = group_levels)
  
  row_height_in  <- 0.22
  base_height_in <- 1.5
  n_rows         <- n_distinct(logfc_df$gene_label)
  heatmap_height <- base_height_in + n_rows * row_height_in
  
  width_logfc <- base_width_in + n_distinct(logfc_df$comparison) * col_width_in
  width_z     <- base_width_in + n_distinct(z_df$group) * col_width_in
  
  title_stub <- sprintf("%s \u2014 %s", selection_tissue, unname(comp_labels[selection_comparison]))
  
  p_logfc <- plot_logfc_heatmap(logfc_df, paste0(title_stub, " (logFC)"))
  p_z     <- plot_group_zscore_heatmap(z_df, paste0(title_stub, " (group z-score)"))
  
  base_name <- sprintf("heatmap_top%d_%s_%s", n_top, tolower(selection_tissue), selection_comparison)
  
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.pdf")), p_logfc, device = "pdf",
         width = width_logfc, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.png")), p_logfc, device = "png",
         width = width_logfc, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.pdf")), p_z, device = "pdf",
         width = width_z, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.png")), p_z, device = "png",
         width = width_z, height = heatmap_height)
  
  write_xlsx(
    list(logFC = logfc_df, zscore = z_df, selection = selected),
    path = file.path(data_dir, paste0(base_name, ".xlsx"))
  )
  
  list(logfc_plot = p_logfc, zscore_plot = p_z, selected = selected,
       logfc_data = logfc_df, zscore_data = z_df)
}

build_significant_both_heatmap_pair <- function(df_vol, dalists, selection_tissue,
                                                contrasts, top_n = Inf,
                                                contrast_a = "Residual",
                                                contrast_b = "Training",
                                                group_col = "cancer_time",
                                                group_levels = c("CTL_T1", "CRC_T1", "CRC_T2"),
                                                fig_dir = config$fig_dir, data_dir = config$data_dir) {
  
  selected <- select_significant_in_both(df_vol, selection_tissue, contrast_a, contrast_b, top_n)
  
  if (nrow(selected) == 0) {
    warning(sprintf("No proteins significant in both %s and %s for %s — skipping.",
                    contrast_a, contrast_b, selection_tissue))
    return(invisible(NULL))
  }
  
  logfc_df <- build_logfc_matrix(df_vol, selected, selection_tissue, contrasts)
  z_df     <- build_group_zscore_matrix(dalists[[selection_tissue]], selected,
                                        group_col = group_col, group_levels = group_levels)
  
  row_height_in  <- 0.22
  base_height_in <- 1.5
  n_rows         <- n_distinct(logfc_df$gene_label)
  heatmap_height <- base_height_in + n_rows * row_height_in
  
  width_logfc <- base_width_in + n_distinct(logfc_df$comparison) * col_width_in
  width_z     <- base_width_in + n_distinct(z_df$group) * col_width_in
  
  title_stub <- sprintf("%s \u2014 Sig. in both %s & %s", selection_tissue,
                        unname(comp_labels[contrast_a]), unname(comp_labels[contrast_b]))
  
  p_logfc <- plot_logfc_heatmap(logfc_df, paste0(title_stub, " (logFC)"))
  p_z     <- plot_group_zscore_heatmap(z_df, paste0(title_stub, " (group z-score)"))
  
  base_name <- sprintf("heatmap_sigboth_%s_%s_%s", tolower(selection_tissue), contrast_a, contrast_b)
  
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.pdf")), p_logfc, device = "pdf",
         width = width_logfc, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.png")), p_logfc, device = "png",
         width = width_logfc, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.pdf")), p_z, device = "pdf",
         width = width_z, height = heatmap_height)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.png")), p_z, device = "png",
         width = width_z, height = heatmap_height)
  
  write_xlsx(
    list(logFC = logfc_df, zscore = z_df, selection = selected),
    path = file.path(data_dir, paste0(base_name, ".xlsx"))
  )
  
  list(logfc_plot = p_logfc, zscore_plot = p_z, selected = selected,
       logfc_data = logfc_df, zscore_data = z_df)
}


# ============================================================
# --- 7c: Union heatmap — CRC vs CTL and/or CRC+RET vs CRC, ---
# --- with a membership strip instead of separate heatmaps  ---
# ============================================================
#' Union of top-N up/down proteins from contrast_a and top-N up/down from
#' contrast_b — NOT top-N of the combined pool. Each contrast contributes
#' its own top 20 up + top 20 down (by Pi.Val); the two sets are unioned,
#' so the result is capped at 4*n_top_per_contrast rows (fewer if there's
#' overlap or either contrast lacks n_top_per_contrast hits per direction).
#'
#' membership is computed against the FULL significant set for each
#' contrast (DE.pi != "NDE"), not just the top-N — so a protein that's
#' top-20 in contrast_b but only "significant, rank 40" in contrast_a
#' still correctly shows "Both" rather than being mislabeled contrast_b-only.
select_union_two_contrasts <- function(df_vol, selection_tissue,
                                       contrast_a = "Cancer",
                                       contrast_b = "Training",
                                       n_top_per_contrast = 20) {
  
  # Full significant sets, for membership tagging
  sig_a_full <- df_vol %>%
    filter(tissue == selection_tissue, comparison == contrast_a, DE.pi != "NDE") %>%
    dplyr::select(protein_key, gene, logFC_a = logFC, PiVal_a = Pi.Val, DEpi_a = DE.pi)
  
  sig_b_full <- df_vol %>%
    filter(tissue == selection_tissue, comparison == contrast_b, DE.pi != "NDE") %>%
    dplyr::select(protein_key, gene, logFC_b = logFC, PiVal_b = Pi.Val, DEpi_b = DE.pi)
  
  # Top-N up/down FROM EACH CONTRAST INDEPENDENTLY — this is what caps
  # the row count, not a post-union ranking
  top_a <- sig_a_full %>% arrange(PiVal_a) %>%
    group_by(DEpi_a) %>% slice_head(n = n_top_per_contrast) %>% ungroup()
  top_b <- sig_b_full %>% arrange(PiVal_b) %>%
    group_by(DEpi_b) %>% slice_head(n = n_top_per_contrast) %>% ungroup()
  
  selected_keys <- union(top_a$protein_key, top_b$protein_key)
  
  combined <- tibble(protein_key = selected_keys) %>%
    left_join(sig_a_full, by = "protein_key") %>%
    left_join(sig_b_full, by = "protein_key") %>%
    mutate(
      gene = dplyr::coalesce(gene.x, gene.y)
    ) %>%
    dplyr::select(-gene.x, -gene.y) %>%
    mutate(
      membership = dplyr::case_when(
        !is.na(DEpi_a) & !is.na(DEpi_b) ~ "Both",
        !is.na(DEpi_a)                  ~ unname(comp_labels[contrast_a]),
        !is.na(DEpi_b)                  ~ unname(comp_labels[contrast_b])
      ),
      membership = factor(membership, levels = c(unname(comp_labels[contrast_a]),
                                                 unname(comp_labels[contrast_b]),
                                                 "Both")),
      logFC     = if_else(!is.na(logFC_b), logFC_b, logFC_a),
      Pi.Val    = if_else(!is.na(PiVal_b), PiVal_b, PiVal_a),
      DE.pi     = if_else(!is.na(DEpi_b), DEpi_b, DEpi_a),
      direction = factor(if_else(DE.pi == "UP", "Up", "Down"), levels = c("Up", "Down"))
    )
  
  message(sprintf(
    "[%s] %d proteins in union of top-%d up/down from %s and top-%d up/down from %s (%d %s only, %d %s only, %d both)",
    selection_tissue, nrow(combined), n_top_per_contrast, contrast_a, n_top_per_contrast, contrast_b,
    sum(combined$membership == unname(comp_labels[contrast_a])), contrast_a,
    sum(combined$membership == unname(comp_labels[contrast_b])), contrast_b,
    sum(combined$membership == "Both")
  ))
  
  combined
}

#' Narrow single-column tile strip, one row per protein, fill = which
#' contrast(s) it's significant in. Shares gene_label factor levels with
#' the main heatmap so rows line up when patchworked side by side.
#'
#' No per-column title (the "Significant in" text above the column is
#' gone) — the legend guide itself is still named "Significant in" so the
#' color key is still labeled, just not repeated as a plot title. No
#' panel border. Column has no title so it needs no top margin either.
plot_membership_strip <- function(df, membership_colors) {
  ggplot(df, aes(x = 1, y = gene_label, fill = membership)) +
    geom_tile(color = "white", linewidth = 0.4) +
    scale_fill_manual(values = membership_colors, name = "Significant in", drop = FALSE) +
    scale_x_continuous(breaks = NULL) +
    labs(x = NULL, y = NULL, title = NULL) +
    theme_cr(base_size = 9) +
    theme(
      axis.text.x     = element_blank(),
      axis.ticks.x    = element_blank(),
      axis.text.y     = element_text(face = "italic", size = 6),
      panel.border    = element_blank(),
      legend.position = "bottom",
      legend.box      = "vertical"
    )
}

#' Sizing/legend fixes from before, PLUS:
#'   - strip_width_in shrunk further (1.6 -> 0.45in) now that it carries
#'     no title and no legend swatches of its own to make room for
#'   - panel.border stripped from both the logFC and z-score panels
#'     (matches the strip, which already has no border)
#'   - title no longer lives inside either column's ggtitle(); it's added
#'     once via patchwork's plot_annotation() AFTER combining, so it's
#'     centered over the full strip+heatmap width instead of just the
#'     right-hand column
build_union_heatmap_trio <- function(df_vol, dalists, selection_tissue,
                                     contrasts,
                                     contrast_a = "Cancer",
                                     contrast_b = "Training",
                                     n_top_per_contrast = 20,
                                     group_col = "cancer_time",
                                     group_levels = c("CTL_T1", "CRC_T1", "CRC_T2"),
                                     strip_width_in = 0.45,
                                     min_total_width_in = 4.5,
                                     fig_dir = config$fig_dir, data_dir = config$data_dir) {
  
  selected <- select_union_two_contrasts(df_vol, selection_tissue, contrast_a, contrast_b,
                                         n_top_per_contrast)
  
  if (nrow(selected) == 0) {
    warning(sprintf("No proteins in union for %s / %s / %s — skipping.",
                    selection_tissue, contrast_a, contrast_b))
    return(invisible(NULL))
  }
  
  gene_order <- build_gene_order(
    selected %>% mutate(gene_label = gene) %>% dplyr::select(gene_label, direction, logFC)
  )
  
  strip_df <- selected %>%
    mutate(gene_label = factor(gene, levels = rev(gene_order)))
  
  logfc_df <- build_logfc_matrix(df_vol, selected, selection_tissue, contrasts)
  z_df     <- build_group_zscore_matrix(dalists[[selection_tissue]], selected,
                                        group_col = group_col, group_levels = group_levels)
  
  membership_levels <- levels(selected$membership)
  membership_colors <- setNames(
    c("#66C2A5", "#FC8D62", "#8DA0CB")[seq_along(membership_levels)],
    membership_levels
  )
  
  row_height_in  <- 0.22
  base_height_in <- 1.5
  n_rows         <- n_distinct(logfc_df$gene_label)
  heatmap_height <- base_height_in + n_rows * row_height_in
  
  width_logfc <- base_width_in + n_distinct(logfc_df$comparison) * col_width_in
  width_z     <- base_width_in + n_distinct(z_df$group) * col_width_in
  
  # Narrow cells shouldn't starve the title/axis/legend of room — floor
  # the total figure width regardless of how few data columns there are
  total_width_logfc <- max(min_total_width_in, strip_width_in + width_logfc)
  total_width_z     <- max(min_total_width_in, strip_width_in + width_z)
  
  title_stub <- stringr::str_wrap(
    sprintf("%s — %s and/or %s", selection_tissue,
            unname(comp_labels[contrast_a]), unname(comp_labels[contrast_b])),
    width = 60
  )
  
  p_strip <- plot_membership_strip(strip_df, membership_colors)
  p_logfc <- plot_logfc_heatmap(logfc_df, title = NULL) +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom")
  p_z <- plot_group_zscore_heatmap(z_df, title = NULL) +
    theme(axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          panel.border = element_blank(),
          legend.position = "bottom")
  
  combo_logfc <- (p_strip + p_logfc +
                    plot_layout(widths = c(strip_width_in, width_logfc), guides = "collect")) +
    plot_annotation(title = paste0(title_stub, " (logFC)")) &
    theme(legend.position = "bottom", legend.box = "vertical",
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  combo_z <- (p_strip + p_z +
                plot_layout(widths = c(strip_width_in, width_z), guides = "collect")) +
    plot_annotation(title = paste0(title_stub, " (group z-score)")) &
    theme(legend.position = "bottom", legend.box = "vertical",
          plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
  
  base_name <- sprintf("heatmap_union_%s_%s_%s", tolower(selection_tissue), contrast_a, contrast_b)
  
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.pdf")), combo_logfc, device = "pdf",
         width = total_width_logfc, height = heatmap_height, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(base_name, "_logFC.png")), combo_logfc, device = "png",
         width = total_width_logfc, height = heatmap_height, dpi = 300, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.pdf")), combo_z, device = "pdf",
         width = total_width_z, height = heatmap_height, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(base_name, "_zscore.png")), combo_z, device = "png",
         width = total_width_z, height = heatmap_height, dpi = 300, limitsize = FALSE)
  
  write_xlsx(
    list(logFC = logfc_df, zscore = z_df, selection = selected),
    path = file.path(data_dir, paste0(base_name, ".xlsx"))
  )
  
  list(logfc_plot = combo_logfc, zscore_plot = combo_z, strip_plot = p_strip,
       selected = selected, logfc_data = logfc_df, zscore_data = z_df)
}

# --- Driver: run for both tissues -----------------------------------
union_results <- list()

for (tis in c("Muscle", "Plasma")) {
  message(sprintf("Building union heatmap trio: %s", tis))
  res <- build_union_heatmap_trio(
    df_vol, dalists, selection_tissue = tis,
    contrasts = heat_config$contrasts,
    contrast_a = "Cancer", contrast_b = "Training",
    group_col = group_col, group_levels = group_levels
  )
  union_results[[tis]] <- res
}

message("Done. Union heatmap trios (strip + logFC/z-score) written to: ",
        config$fig_dir, " / ", config$data_dir)


# ============================================================
# --- 8: Composite figure — volcano (A) + Muscle heatmap (B) + Plasma heatmap (C) ---
# ============================================================
# Panel A: the FULL volcano — both tissues, all 3 SURV contrasts (this
# is just volcano_all, already built in section 3).
# Panel B: Muscle union heatmap (strip + logFC) for the specified
# contrast pair (Baseline_SURVvCTL = "CRC vs CTL", Training_SURV =
# "CRC+RET vs CRC").
# Panel C: same, for Plasma. B and C sit SIDE BY SIDE (not stacked).
#
# LAYOUT: A on the left, (B | C) on the right (`volcano_plot | (panel_b | panel_c)`).
# Patchwork stretches the whole right-hand block to match A's full
# height automatically (single-row layout).
#
# SHARED LEGEND: Muscle and Plasma each have their own logFC range by
# default, so their color scales (and therefore their legends) would
# NOT be identical -- patchwork only collapses guides that match
# exactly. To get one shared legend instead of two redundant ones:
#   1. Compute a single shared_range from BOTH tissues' logFC values
#      and pass it into plot_logfc_heatmap(..., range_val = shared_range)
#      for both, so the color scales are identical.
#   2. The membership ("Significant in") scale is already identical
#      across tissues (same fixed factor levels/colors), so no extra
#      work is needed there.
#   3. guides = "collect" is set ONCE, at the (panel_b | panel_c) level
#      -- NOT inside build_composite_heatmap_data()/individual panels --
#      so patchwork gathers all four leaf legends (strip x2, heatmap x2)
#      and merges the two identical pairs into one strip legend + one
#      heatmap legend for the whole B+C block.
#
# TAGGING: each tissue's strip + heatmap is now built as ONE ggplot
# object (build_composite_heatmap_single_plot(), using
# ggnewscale::new_scale_fill() for the second fill scale and
# ggh4x::force_panelsizes() for the narrow-strip/wide-heatmap column
# widths) rather than two plots glued with patchwork. A single ggplot
# can only receive one tag, so the top-level plot_annotation(tag_levels
# = "A") now correctly assigns exactly 3 letters total: A (volcano),
# B (muscle), C (plasma) -- no tag_level = "keep" workaround needed.
#
# TEXT SIZE: axis/legend text in B and C is shrunk relative to the
# standalone heatmap versions so ~80 rows stays legible once compressed
# to match A's height.

#' Data-only prep for one tissue's composite heatmap panel (selection +
#' strip data + logFC matrix). Split out from plot-building so the
#' driver can compute a shared logFC color range across tissues BEFORE
#' any ggplot object is built.
build_composite_heatmap_data <- function(df_vol, selection_tissue, contrasts,
                                         contrast_a, contrast_b, n_top_per_contrast) {
  
  selected <- select_union_two_contrasts(df_vol, selection_tissue, contrast_a, contrast_b,
                                         n_top_per_contrast)
  
  if (nrow(selected) == 0) {
    warning(sprintf("No proteins in union for %s / %s / %s — skipping.",
                    selection_tissue, contrast_a, contrast_b))
    return(invisible(NULL))
  }
  
  gene_order <- build_gene_order(
    selected %>% mutate(gene_label = gene) %>% dplyr::select(gene_label, direction, logFC)
  )
  strip_df <- selected %>% mutate(gene_label = factor(gene, levels = rev(gene_order)))
  logfc_df <- build_logfc_matrix(df_vol, selected, selection_tissue, contrasts)
  
  membership_levels <- levels(selected$membership)
  membership_colors <- setNames(
    c("#66C2A5", "#FC8D62", "#8DA0CB")[seq_along(membership_levels)],
    membership_levels
  )
  
  list(selected = selected, strip_df = strip_df, logfc_df = logfc_df,
       membership_colors = membership_colors)
}

#' Build the strip + heatmap as a SINGLE ggplot per tissue (not two
#' patchwork-glued plots) -- this is what fixes the A/B/C/D/E over-
#' tagging: patchwork's automatic tagging descends into every leaf
#' plot of a composite, and tag_level = "keep" wasn't stopping that
#' descent, so each tissue's strip and heatmap were tagged separately.
#' A single ggplot object can only ever receive one tag.
#'
#' Two fill scales in one plot (categorical "Significant in" for the
#' strip column, continuous logFC gradient for the comparison columns)
#' via ggnewscale::new_scale_fill(). Column widths are now set with a
#' CONTINUOUS x-axis and explicit numeric positions/tile widths -- not
#' facets. (The earlier facet_grid2(space = "free_x") approach only
#' gave a fixed ~1:3 ratio tied to category counts, and the
#' ggh4x::force_panelsizes() attempt at exact control crashes in this
#' ggh4x version with "attempt to set an attribute on NULL".) Manual
#' positioning sidesteps both: strip_rel_width sets the strip's width
#' relative to a heatmap column (1 = same width; 0.25 = one quarter --
#' the current default), and strip_gap_rel inserts blank space between
#' the strip and the heatmap block so it reads as a separate labeling
#' strip rather than another heatmap column.
#'
#' shared_range: same logFC color scale passed to both tissues so their
#' legends can still collect into one downstream.
build_composite_heatmap_single_plot <- function(tissue_data, contrasts,
                                                heatmap_text_size, shared_range,
                                                strip_rel_width = 0.25,  # strip width relative to one heatmap column (1 = equal; was effectively ~1 before)
                                                strip_gap_rel = 0.3) {  # blank space between strip and heatmap block, same relative units
  
  n_contrasts <- length(contrasts)
  heat_rel_width <- 1  # each heatmap column's relative width unit
  
  # Layout, left to right: [strip][gap][heat_1][heat_2]...[heat_n],
  # laid out like a simple barplot -- cumulative left edges + half of
  # each segment's own width gives its center.
  seg_widths <- c(strip_rel_width, strip_gap_rel, rep(heat_rel_width, n_contrasts))
  left_edges <- cumsum(c(0, seg_widths[-length(seg_widths)]))
  centers    <- left_edges + seg_widths / 2
  
  strip_center <- centers[1]
  heat_centers <- centers[3:(2 + n_contrasts)]
  
  strip_part <- tissue_data$strip_df %>%
    dplyr::select(gene_label, membership) %>%
    mutate(x_pos = strip_center, tile_width = strip_rel_width)
  
  logfc_part <- tissue_data$logfc_df %>%
    mutate(x_pos = heat_centers[match(as.character(comparison), contrasts)],
           tile_width = heat_rel_width)
  
  x_breaks <- c(strip_center, heat_centers)
  x_labels <- c("", unname(comp_labels[contrasts]))
  
  ggplot(mapping = aes(y = gene_label)) +
    geom_tile(data = strip_part, aes(x = x_pos, width = tile_width, fill = membership),
              color = "white", linewidth = 0.4) +
    scale_fill_manual(values = tissue_data$membership_colors, name = "Significant in", drop = FALSE,
                      guide = guide_legend(keywidth = unit(7, "pt"), keyheight = unit(7, "pt"),
                                           direction = "vertical")) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = logfc_part, aes(x = x_pos, width = tile_width, fill = logFC),
              color = "white", linewidth = 0.4) +
    geom_text(data = filter(logfc_part, sig), aes(x = x_pos, label = "*"),
              color = "black", size = 3.5, vjust = 0.75) +
    scale_fill_gradient2(
      low = colors$blue, mid = "white", high = colors$red,
      midpoint = 0, limits = c(-shared_range, shared_range), oob = scales::squish,
      name = "logFC",
      guide = guide_colorbar(barwidth = unit(6, "pt"), barheight = unit(45, "pt"), direction = "vertical")
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(x = NULL, y = NULL) +
    theme_cr(base_size = 9) +
    theme(
      axis.text.y     = element_text(face = "italic", size = heatmap_text_size),
      axis.text.x     = element_text(size = heatmap_text_size, angle = 30, hjust = 1),
      axis.ticks.x    = element_blank(),
      panel.border    = element_blank(),
      panel.grid      = element_blank(),
      legend.text     = element_text(size = heatmap_text_size - 1),
      legend.title    = element_text(size = heatmap_text_size - 1)
    )
}

#' Driver: builds Panel A (full volcano, reused from volcano_all) +
#' Panel B (Muscle heatmap) + Panel C (Plasma heatmap) SIDE BY SIDE with
#' one shared legend, combines as A | (B | C) so B+C together
#' height-match A, and saves.
#'
#' EQUAL CELL THICKNESS: Muscle and Plasma almost never have the same
#' protein count, so simply stretching both to the same panel height
#' (the previous behavior) makes their tiles different thicknesses --
#' Plasma's fewer rows get stretched into visibly fatter cells. Instead:
#'   1. row_h is computed from MUSCLE alone: (composite_height_in -
#'      overhead_in) / rows_muscle -- i.e. the cell thickness Muscle
#'      needs to exactly fill the full composite height.
#'   2. Plasma is given a height of exactly `overhead_in + rows_plasma *
#'      row_h` -- the height IT needs at that SAME row_h -- via
#'      patchwork::area() with row-height weights, not by letting it
#'      stretch to fill the column.
#'   3. Whatever vertical space is left in Plasma's column after that
#'      goes to patchwork::guide_area(), which is where the collected
#'      shared legend renders instead of floating at the bottom of the
#'      whole figure.
#' This assumes Muscle has >= rows than Plasma; if not, a warning fires
#' and the layout falls back to the old equal-stretch behavior.
build_full_volcano_union_composite <- function(df_vol, dalists, volcano_plot,
                                               contrasts,
                                               contrast_a = "Cancer",
                                               contrast_b = "Training",
                                               n_top_per_contrast = 20,
                                               heatmap_text_size = 4,  # pt size for heatmap axis/legend text -- shrunk so ~80 rows stays readable at compressed height
                                               title_size = 11,
                                               panel_a_width = 1.4,   # relative width, volcano vs. the B+C block
                                               panel_bc_widths = c(1, 1),  # relative width split, Muscle vs. Plasma
                                               overhead_in = 0.9,     # non-tile vertical space per panel (x-axis text, spacing) -- tune if rows still look mismatched
                                               composite_width_in = 12,
                                               composite_height_in = 8,
                                               fig_dir = config$fig_dir, data_dir = config$data_dir) {
  
  muscle_data <- build_composite_heatmap_data(df_vol, "Muscle", contrasts, contrast_a, contrast_b,
                                              n_top_per_contrast)
  plasma_data <- build_composite_heatmap_data(df_vol, "Plasma", contrasts, contrast_a, contrast_b,
                                              n_top_per_contrast)
  
  if (is.null(muscle_data) || is.null(plasma_data)) {
    warning("Missing Muscle or Plasma data -- skipping full composite.")
    return(invisible(NULL))
  }
  
  # Shared color scale across both tissues -- required for the two
  # heatmap legends to collapse into one instead of showing twice.
  shared_range <- max(abs(c(muscle_data$logfc_df$logFC, plasma_data$logfc_df$logFC)), na.rm = TRUE)
  
  panel_b <- build_composite_heatmap_single_plot(muscle_data, contrasts, heatmap_text_size, shared_range)
  panel_c <- build_composite_heatmap_single_plot(plasma_data, contrasts, heatmap_text_size, shared_range)
  
  rows_muscle <- n_distinct(muscle_data$logfc_df$gene_label)
  rows_plasma <- n_distinct(plasma_data$logfc_df$gene_label)
  
  if (rows_plasma > rows_muscle) {
    warning("Plasma has more rows than Muscle -- equal cell thickness with a legend gap isn't possible ",
            "this direction; falling back to equal-stretch layout (both panels fill the full height).")
    
    bc_block <- (panel_b | panel_c) +
      plot_layout(widths = panel_bc_widths, guides = "collect")
  } else {
    # Cell thickness Muscle needs to exactly fill the full composite height
    row_h <- (composite_height_in - overhead_in) / rows_muscle
    # Height Plasma needs at that SAME row_h -- NOT stretched to fill the column
    plasma_row_height <- overhead_in + rows_plasma * row_h
    legend_row_height <- composite_height_in - plasma_row_height
    
    bc_design <- c(
      patchwork::area(t = 1, l = 1, b = 2, r = 1),  # Muscle: full column height
      patchwork::area(t = 1, l = 2, b = 1, r = 2),  # Plasma: top row only, sized to match Muscle's row_h
      patchwork::area(t = 2, l = 2, b = 2, r = 2)   # shared legend: leftover space under Plasma
    )
    
    bc_block <- (panel_b + panel_c + patchwork::guide_area()) +
      plot_layout(design = bc_design, widths = panel_bc_widths,
                  heights = c(plasma_row_height, legend_row_height), guides = "collect") &
      theme(legend.box = "horizontal",
            legend.direction = "vertical",
            legend.text  = element_text(size = heatmap_text_size - 1),
            legend.title = element_text(size = heatmap_text_size - 1))
  }
  
  # `|` places A and (B | C) in one patchwork row, so the right-hand
  # block is stretched to match A's full height automatically.
  composite <- (volcano_plot | bc_block) +
    plot_layout(widths = c(panel_a_width, 1)) +
    plot_annotation(
      title = sprintf("%s vs. %s", unname(comp_labels[contrast_a]), unname(comp_labels[contrast_b])),
      tag_levels = "A"
    ) &
    theme(plot.title = element_text(size = title_size, face = "bold", hjust = 0.5),
          plot.tag   = element_text(size = 11, face = "bold"))
  
  base_name <- sprintf("composite_volcano_full_union_%s_%s", contrast_a, contrast_b)
  
  ggsave(file.path(fig_dir, paste0(base_name, ".pdf")), composite, device = "pdf",
         width = composite_width_in, height = composite_height_in, limitsize = FALSE)
  ggsave(file.path(fig_dir, paste0(base_name, ".png")), composite, device = "png",
         width = composite_width_in, height = composite_height_in, dpi = 300, limitsize = FALSE)
  
  list(composite = composite, volcano = volcano_plot,
       muscle_panel = panel_b, plasma_panel = panel_c,
       muscle_selected = muscle_data$selected, plasma_selected = plasma_data$selected)
}

full_composite_result <- build_full_volcano_union_composite(
  df_vol, dalists, volcano_plot = volcano_all,
  contrasts = heat_config$contrasts,
  contrast_a = "Cancer", contrast_b = "Training"
)

message("Done. Full volcano + Muscle/Plasma union heatmap composite written to: ", config$fig_dir)


# =============================================================================
# --- 9: Harmonized composite -- matches fig_05's style (fig_03_style.R) ---
# =============================================================================
# Same content as Section 8 (Panel A = full volcano, Panel B = Muscle union
# heatmap, Panel C = Plasma union heatmap) but rebuilt to match the rest of
# the paper's figures:
#   - mm-native geometry (COMP_W / panel widths / heights), not inches
#   - wrap_elements(full = ...) per panel instead of raw patchwork "+"/"|"
#     composition -- this is what let fig_05 skip patchwork's automatic
#     tagging entirely; tags/titles/subtitles are drawn manually via
#     cowplot::draw_label() in an explicit header band above the content,
#     the same pattern fig_05_stitch.R uses for all four of its panels
#   - the shared Muscle+Plasma legend is extracted once via
#     build_heatmap_legend() (fig_03_style.R) and placed with
#     cowplot::draw_plot() at an explicit position, instead of
#     patchwork::guide_area() -- this also sidesteps the ggh4x crash and
#     the automatic-tagging issues Section 8 worked around manually
#   - fig_theme_scaled() + scale_text() drive all panel text sizing, so
#     axis text, titles, and legend text scale consistently with fig_05's
#     panels rather than using this figure's own ad hoc pt sizes
#
# Panel A keeps its OWN Direction (UP/DOWN) legend inside the panel --
# it's a different scale than the Muscle/Plasma "Significant in"/logFC
# legend and isn't a good candidate to merge with it, so only B and C get
# stripped of their own legends via strip_for_composite().

#' Volcano, sized for composite embedding via fig_theme_scaled()/scale_text()
#' instead of the fixed pt sizes make_volcano_plot() uses for the
#' standalone panel_3A.pdf. Otherwise identical logic to make_volcano_plot().
make_volcano_plot_scaled <- function(df_vol, selected_comparisons, panel_width_mm,
                                     n_labels = 5, strip_fill_colors = NULL) {
  df <- df_vol %>%
    filter(comparison %in% selected_comparisons) %>%
    mutate(comparison = factor(comparison, levels = selected_comparisons))
  
  df_label <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, tissue, DE.pi) %>%
    slice_max(abs(logFC), n = n_labels) %>%
    ungroup()
  
  df_counts <- df %>%
    filter(DE.pi != "NDE") %>%
    group_by(comparison, tissue, DE.pi) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(
      x     = ifelse(DE.pi == "UP",  Inf, -Inf),
      y     = Inf,
      hjust = ifelse(DE.pi == "UP",  1.2, -0.2),
      label = as.character(n),
      color = ifelse(DE.pi == "UP", colors$red, colors$blue)
    )
  
  tissue_bg <- df %>%
    distinct(tissue) %>%
    mutate(fill_color = unname(colors$tissue[tissue]))
  
  # strip_fill_colors, when supplied, keys by comparison name and lets
  # the volcano's facet labels match the heatmap's "Significant in"
  # strip colors -- falls back to plain light_gray if not given (e.g.
  # for the standalone panel_3A.pdf via make_volcano_plot(), unaffected).
  strip_x_fills <- if (!is.null(strip_fill_colors)) {
    unname(strip_fill_colors[selected_comparisons])
  } else {
    rep(colors$light_gray, length(selected_comparisons))
  }
  
  strip <- make_strip(
    fills = list(x = strip_x_fills,
                 y = unname(colors$tissue)),
    type  = "xy"
  )
  
  ggplot(df, aes(logFC, -log10(P.Value), color = DE.pi)) +
    geom_rect(
      data        = tissue_bg,
      aes(fill    = fill_color),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      alpha       = 0.15,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray) +
    geom_point(size = 1.2, alpha = 0.8) +
    geom_text_repel(data = df_label, aes(label = gene),
                    color = "black", fontface = "bold", size = scale_text(BASE_GENE_LABEL, panel_width_mm),
                    show.legend = FALSE, max.overlaps = Inf,
                    box.padding = 0.4, point.padding = 0.3,
                    force = 2, force_pull = 0.01,
                    min.segment.length = 0.05, segment.color = "black",
                    arrow = arrow(length = unit(0.01, "npc")), seed = 42) +
    geom_text(
      data        = df_counts,
      aes(x = x, y = y, label = label, hjust = hjust),
      vjust       = 1.5,
      fontface    = "bold",
      size        = scale_text(BASE_COUNT, panel_width_mm),
      color       = df_counts$color,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values   = c(UP = colors$red, DOWN = colors$blue),
                       na.value = colors$gray) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    coord_cartesian(xlim = c(-4, 4)) +
    facet_grid2(tissue ~ comparison, strip = strip,
                labeller = labeller(comparison = comp_labels)) +
    labs(x = bquote(bold("log2(Fold Change)")), y = bquote(bold("-log10(p-value)")),
         color = "Direction", title = NULL) +
    fig_theme_scaled(panel_width_mm)
}

#' Muscle/Plasma union heatmap, sized for composite embedding. Same layout
#' logic as the Section 8 version (manual continuous-x positioning for the
#' narrow strip / wide heatmap columns) but text sizing now comes from
#' fig_theme_scaled()/scale_text() and legend key/colorbar geometry from
#' BASE_LEGEND_KEY/BASE_COLORBAR_LEN/BASE_COLORBAR_WID, all scaled by
#' panel_width_mm to match the rest of the composite's typography.
build_composite_heatmap_scaled <- function(tissue_data, contrasts, panel_width_mm, shared_range,
                                           strip_rel_width = 0.25, strip_gap_rel = 0.3) {
  
  n_contrasts    <- length(contrasts)
  heat_rel_width <- 1
  
  seg_widths <- c(strip_rel_width, strip_gap_rel, rep(heat_rel_width, n_contrasts))
  left_edges <- cumsum(c(0, seg_widths[-length(seg_widths)]))
  centers    <- left_edges + seg_widths / 2
  
  strip_center <- centers[1]
  heat_centers <- centers[3:(2 + n_contrasts)]
  
  strip_part <- tissue_data$strip_df %>%
    dplyr::select(gene_label, membership) %>%
    mutate(x_pos = strip_center, tile_width = strip_rel_width)
  
  logfc_part <- tissue_data$logfc_df %>%
    mutate(x_pos = heat_centers[match(as.character(comparison), contrasts)],
           tile_width = heat_rel_width)
  
  x_breaks <- c(strip_center, heat_centers)
  x_labels <- c("", unname(comp_labels[contrasts]))
  
  key_pt      <- round(scale_text(BASE_LEGEND_KEY, panel_width_mm), 1)
  colorbar_ln <- round(scale_text(BASE_COLORBAR_LEN, panel_width_mm), 1)
  colorbar_wd <- round(scale_text(BASE_COLORBAR_WID, panel_width_mm), 1)
  
  ggplot(mapping = aes(y = gene_label)) +
    geom_tile(data = strip_part, aes(x = x_pos, width = tile_width, fill = membership),
              color = "white", linewidth = 0.4) +
    scale_fill_manual(values = tissue_data$membership_colors, name = "Significant in", drop = FALSE,
                      guide = guide_legend(keywidth = unit(key_pt, "pt"), keyheight = unit(key_pt, "pt"),
                                           direction = "vertical")) +
    ggnewscale::new_scale_fill() +
    geom_tile(data = logfc_part, aes(x = x_pos, width = tile_width, fill = logFC),
              color = "white", linewidth = 0.4) +
    geom_text(data = filter(logfc_part, sig), aes(x = x_pos, label = "*"),
              color = "black", size = scale_text(BASE_SIG_MARK, panel_width_mm), vjust = 0.75) +
    scale_fill_gradient2(
      low = colors$blue, mid = "white", high = colors$red,
      midpoint = 0, limits = c(-shared_range, shared_range), oob = scales::squish,
      name = "logFC",
      guide = guide_colorbar(barwidth = unit(colorbar_wd, "pt"), barheight = unit(colorbar_ln, "pt"),
                             direction = "vertical")
    ) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    labs(x = NULL, y = NULL, title = NULL) +
    fig_theme_scaled(panel_width_mm) +
    theme(
      axis.text.x  = element_text(angle = 30, hjust = 1),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    )
}

#' Driver: assembles Panel A (full volcano) + Panel B (Muscle) + Panel C
#' (Plasma) using the fig_05-style header-band + wrap_elements() +
#' draw_label() + extracted-legend pattern, and saves with save_scaled().
#'
#' Equal cell thickness (same logic as Section 8): row_h is set by Muscle's
#' row count filling the full content height; Plasma is drawn at the SAME
#' row_h (so its tiles match Muscle's thickness) and only occupies the
#' portion of the height it needs, leaving the remainder of its column for
#' the extracted shared legend.
build_fig03_composite <- function(df_vol, volcano_comparisons = heat_config$contrasts,
                                  contrast_a = "Cancer",
                                  contrast_b = "Training",
                                  n_top_per_contrast = 20,
                                  comp_w_mm = 380,          # matches fig_05's COMP_W -- keep the master figure width consistent across the paper
                                  header_h_mm = 18,         # matches fig_05's HEADER_H
                                  total_content_h_mm = 170, # matches fig_05's ROW1_H
                                  overhead_mm = 25,         # non-tile vertical space per heatmap panel (axis text, spacing) -- tune if rows still look mismatched
                                  a_width_frac = 0.583,     # volcano's share of comp_w_mm (was panel_a_width = 1.4 of 2.4 total in Section 8)
                                  fig_dir = config$fig_dir, data_dir = config$data_dir) {
  
  # --- Column geometry ---
  a_w_mm  <- round(comp_w_mm * a_width_frac)
  bc_w_mm <- comp_w_mm - a_w_mm
  b_w_mm  <- round(bc_w_mm / 2)   # Muscle
  c_w_mm  <- bc_w_mm - b_w_mm     # Plasma
  
  # --- Data prep (unchanged from Section 8) ---
  muscle_data <- build_composite_heatmap_data(df_vol, "Muscle", volcano_comparisons, contrast_a, contrast_b,
                                              n_top_per_contrast)
  plasma_data <- build_composite_heatmap_data(df_vol, "Plasma", volcano_comparisons, contrast_a, contrast_b,
                                              n_top_per_contrast)
  
  if (is.null(muscle_data) || is.null(plasma_data)) {
    warning("Missing Muscle or Plasma data -- skipping harmonized composite.")
    return(invisible(NULL))
  }
  
  shared_range <- max(abs(c(muscle_data$logfc_df$logFC, plasma_data$logfc_df$logFC)), na.rm = TRUE)
  
  rows_muscle <- n_distinct(muscle_data$logfc_df$gene_label)
  rows_plasma <- n_distinct(plasma_data$logfc_df$gene_label)
  
  if (rows_plasma > rows_muscle) {
    warning("Plasma has more rows than Muscle -- equal cell thickness with a legend gap isn't possible ",
            "this direction; content heights will fall back to matching total_content_h_mm for both.")
    plasma_content_h_mm <- total_content_h_mm
    legend_h_mm <- 0
  } else {
    row_h_mm <- (total_content_h_mm - overhead_mm) / rows_muscle
    plasma_content_h_mm <- overhead_mm + rows_plasma * row_h_mm
    legend_h_mm <- total_content_h_mm - plasma_content_h_mm
  }
  
  comp_h_mm <- header_h_mm + total_content_h_mm
  
  # --- Volcano facet strip colors -- match the heatmap's "Significant
  # in" strip exactly (pulled from muscle_data$membership_colors, not
  # re-hardcoded) so the same hue means the same contrast in both panel
  # types. Any contrast NOT part of the union heatmap (e.g. the unused
  # "CRC+RET vs CTL") gets colors$yellow instead.
  strip_fill_colors <- setNames(rep(colors$yellow, length(volcano_comparisons)), volcano_comparisons)
  strip_fill_colors[contrast_a] <- unname(muscle_data$membership_colors[unname(comp_labels[contrast_a])])
  strip_fill_colors[contrast_b] <- unname(muscle_data$membership_colors[unname(comp_labels[contrast_b])])
  
  # --- Build panels ---
  p_A <- make_volcano_plot_scaled(df_vol, volcano_comparisons, panel_width_mm = a_w_mm,
                                  strip_fill_colors = strip_fill_colors)
  p_B <- build_composite_heatmap_scaled(muscle_data, volcano_comparisons, panel_width_mm = b_w_mm, shared_range = shared_range)
  p_C <- build_composite_heatmap_scaled(plasma_data, volcano_comparisons, panel_width_mm = c_w_mm, shared_range = shared_range)
  
  # Extract the shared legend from a DEDICATED reference plot built at
  # a_w_mm (the volcano's width), not from p_C itself -- p_C is sized at
  # c_w_mm (Plasma's narrow column), and since guide_legend()/
  # guide_colorbar()'s key/colorbar dimensions are baked in at plot-build
  # time (via BASE_LEGEND_KEY/BASE_COLORBAR_LEN/WID scaled by
  # panel_width_mm inside build_composite_heatmap_scaled()), overriding
  # just the legend TEXT theme afterward wouldn't be enough to match the
  # volcano's legend size -- the key/colorbar geometry itself needs to be
  # built at that width too. Muscle and Plasma share an identical scale
  # (shared_range, same membership colors), so this reference plot's
  # actual data choice (Plasma here) doesn't affect the extracted legend.
  legend_ref_plot <- build_composite_heatmap_scaled(plasma_data, volcano_comparisons,
                                                    panel_width_mm = a_w_mm, shared_range = shared_range)
  shared_legend_grob <- build_heatmap_legend(legend_ref_plot, panel_width_mm = a_w_mm)
  
  # A keeps its own Direction legend (different scale, not merged);
  # B and C lose theirs -- the extracted legend covers both.
  p_A_composite <- p_A + labs(title = NULL)
  p_B_composite <- strip_for_composite(p_B)
  p_C_composite <- strip_for_composite(p_C)
  
  design <- c(
    patchwork::area(t = 2, l = 1, b = 3, r = 1),  # A: full content height (rows 2-3), col 1
    patchwork::area(t = 2, l = 2, b = 3, r = 2),  # B (Muscle): full content height, col 2
    patchwork::area(t = 2, l = 3, b = 2, r = 3)   # C (Plasma): row 2 only, col 3 -- row 3 left blank for the legend
  )
  
  fig_raw <- wrap_elements(full = p_A_composite) +
    wrap_elements(full = p_B_composite) +
    wrap_elements(full = p_C_composite) +
    plot_layout(design = design,
                widths  = c(a_w_mm, b_w_mm, c_w_mm),
                heights = c(header_h_mm, plasma_content_h_mm, legend_h_mm))
  
  # --- Header text placement (mm-from-top -> ggdraw y-fraction) ---
  y_frac <- function(mm_from_top) 1 - mm_from_top / comp_h_mm
  
  TAG_SZ <- round(scale_text(BASE_TAG, comp_w_mm))
  TTL_SZ <- TAG_SZ
  SUB_SZ <- round(TAG_SZ * 0.7)
  
  Y_TAG <- y_frac(3)
  Y_SUB <- y_frac(3 + SUB_SZ * 0.6)
  
  X_A <- 0.01
  X_B <- round(a_w_mm / comp_w_mm, 3) + 0.02
  X_C <- round((a_w_mm + b_w_mm) / comp_w_mm, 3) + 0.02
  X_TTL <- 0.03
  
  title_A <- "Differential Abundance \u2014 Volcano Plots"
  title_B <- "Heatmap \u2014 Muscle"
  title_C <- "Heatmap \u2014 Plasma"
  sub_BC  <- sprintf("%s and\n %s", unname(comp_labels[contrast_a]), unname(comp_labels[contrast_b]))
  
  fig_03_harmonized <- ggdraw(fig_raw) +
    draw_label("A", x = X_A, y = Y_TAG, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(title_A, x = X_A + X_TTL, y = Y_TAG, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label("B", x = X_B, y = Y_TAG, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(title_B, x = X_B + X_TTL, y = Y_TAG, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(sub_BC, x = X_B + X_TTL, y = Y_SUB, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40") +
    draw_label("C", x = X_C, y = Y_TAG, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(title_C, x = X_C + X_TTL, y = Y_TAG, size = TTL_SZ, fontface = "bold", hjust = 0, vjust = 1) +
    draw_label(sub_BC, x = X_C + X_TTL, y = Y_SUB, size = SUB_SZ, fontface = "bold.italic", hjust = 0, vjust = 1, colour = "grey40")
  
  # Only draw the legend if there's a real gap for it (see the
  # rows_plasma > rows_muscle fallback above).
  if (legend_h_mm > 0) {
    Y_LEG_BOTTOM <- y_frac(header_h_mm + plasma_content_h_mm + legend_h_mm)
    LEG_H_FRAC   <- legend_h_mm / comp_h_mm
    X_LEG        <- round((a_w_mm + b_w_mm) / comp_w_mm, 3)
    LEG_W_FRAC   <- c_w_mm / comp_w_mm
    
    fig_03_harmonized <- fig_03_harmonized +
      draw_plot(shared_legend_grob, x = X_LEG, y = Y_LEG_BOTTOM, width = LEG_W_FRAC, height = LEG_H_FRAC)
  }
  
  base_name <- sprintf("fig_03_composite_harmonized_%s_%s", contrast_a, contrast_b)
  save_scaled(fig_03_harmonized, fig_dir, base_name, width_mm = comp_w_mm, height_mm = comp_h_mm)
  
  list(composite = fig_03_harmonized, volcano = p_A, muscle_panel = p_B, plasma_panel = p_C,
       muscle_selected = muscle_data$selected, plasma_selected = plasma_data$selected)
}

fig_03_harmonized_result <- build_fig03_composite(df_vol)

message("Done. Harmonized (fig_05-style) composite written to: ", config$fig_dir)