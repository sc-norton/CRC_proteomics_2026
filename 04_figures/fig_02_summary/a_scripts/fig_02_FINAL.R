# Fig 02
#
# NOTE: GSEA data now comes from 02_gsea_pathway_consolidation.R's saved
# RDS rather than being loaded and recomputed from raw per-comparison RDS
# files here. This also changes what section 5's "Significant Pathways (n)"
# panel counts: it now excludes redundant pathways (see GSEA_COUNT_MODE
# below), so it's consistent with what the bubble plot / other GSEA figures
# actually show, rather than a raw per-database hit count that includes
# near-duplicate terms.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, patchwork, cowplot, scales, enrichplot, grid, ggnewscale,
               openxlsx)

# --- Configuration -----------------------------------------------------------
config <- list(
  limma_muscle = "03_analysis/c_data/muscle/limma",
  limma_plasma = "03_analysis/c_data/plasma/limma",
  gsea_data    = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  norm_muscle  = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS",
  norm_plasma  = "01_normalization/c_data/plasma/02_normalized_data_plasma.RDS",
  fig_dir      = "04_figures/fig_02_summary/b_figures",
  data_dir     = "04_figures/fig_02_summary/c_data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ----------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

tissue_colors <- c(Muscle = colors$pink_dark, Plasma = colors$aqua_dark)

# --- Shared constants --------------------------------------------------------
thresh      <- 0.05
adjp_thresh <- 0.05
pi_thresh   <- 0.05

# "representative": count only redundancy-filter survivors (excludes
#   near-duplicate pathways within a comparison x database, per the pooled
#   greedy filter in 02_gsea_pathway_consolidation.R). Matches what the
#   bubble plot and other GSEA figures show.
# "raw": count every significant pathway regardless of redundancy, i.e.
#   the original behavior of this panel. Kept as a one-line toggle in case
#   you want the two numbers side by side for a methods comparison.
GSEA_COUNT_MODE <- "representative"

BG_COLORS <- c(
  Muscle = scales::alpha(colors$pink,  0.2),
  Plasma = scales::alpha(colors$aqua, 0.2)
)

# Tissue label data for in-panel annotations (top-left corner of each facet)
tissue_label_df <- tibble(
  tissue = factor(c("Muscle", "Plasma"), levels = c("Muscle", "Plasma")),
  label  = c("Muscle", "Plasma")
)

gsea_group <- "SURV"  # matches the group column added by consolidation script

gsea_comps <- c("Cancer", "Residual", "Training")

comp_labels <- c(
  Cancer = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training     = "CRC+RET vs CRC",
  `Total Quantified` = "Total\nQuantified"  # new level added for panel A's reference bar (section 2b)
)

SET_LABELS <- c(
  Cancer = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training     = "CRC+RET vs CRC"
)

db_slots <- c(
  gsea_bp_clean = "GO:BP",
  gsea_cc_clean = "GO:CC",
  gsea_mf_clean = "GO:MF",
  gsea_bp_slim   = "GO Slim",
  gsea_kegg      = "KEGG",
  gsea_reactome  = "Reactome",
  gsea_hall      = "Hallmark"
)

db_colors <- c(
  "GO:BP"    = colors$red_light,
  "GO:CC"    = colors$yellow_light,
  "GO:MF"    = colors$blue_light,
  "GO Slim"  = colors$green_light,
  "KEGG"     = colors$purple_light,
  "Reactome" = colors$brown_light,
  "Hallmark" = colors$orange_light
)

DIR_COLORS <- c(Up = colors$red, Down = colors$blue, Mixed = colors$purple)

# --- 1: Load Data ------------------------------------------------------------
limma_results <- list(
  muscle = list(
    Cancer = read_excel(file.path(config$limma_muscle,
                                  "limma_results_all_muscle_cv.xlsx"), sheet = "Cancer"),
    Residual = read_excel(file.path(config$limma_muscle,
                                    "limma_results_all_muscle_cv.xlsx"), sheet = "Residual"),
    Training     = read_excel(file.path(config$limma_muscle,
                                        "limma_results_all_muscle_cv.xlsx"), sheet = "Training")
  ),
  plasma = list(
    Cancer = read_excel(file.path(config$limma_plasma,
                                  "limma_results_all_plasma.xlsx"), sheet = "Cancer"),
    Residual = read_excel(file.path(config$limma_plasma,
                                    "limma_results_all_plasma.xlsx"), sheet = "Residual"),
    Training     = read_excel(file.path(config$limma_plasma,
                                        "limma_results_all_plasma.xlsx"), sheet = "Training")
  )
)

# GSEA: read from the consolidated pathway data instead of raw per-
# comparison RDS files. gsea_supp_full carries every tested pathway across
# all db_slots databases (already tagged with tissue/comparison/database/
# significant); representative_lookup flags which (tissue, database, ID)
# combos survived the pooled redundancy filter.
if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run 02_gsea_pathway_consolidation.R first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

gsea_supp_full        <- gsea_pathway_data$gsea_supp_full        %>% filter(group == gsea_group)
representative_lookup <- gsea_pathway_data$representative_lookup %>% filter(group == gsea_group)

norm_m <- readRDS(file.path(config$norm_muscle))
norm_p <- readRDS(file.path(config$norm_plasma))

# --- 2: Data wrangling -------------------------------------------------------
make_da_long <- function(tissue) {
  bind_rows(lapply(names(limma_results[[tissue]]), function(comp) {
    df <- limma_results[[tissue]][[comp]]
    data.frame(
      comparison = comp,
      pval_n     = sum(df$P.Value   < thresh,    na.rm = TRUE),
      adjp_n     = sum(df$adj.P.Val < thresh,    na.rm = TRUE),
      pival_n    = sum(df$Pi.Val    < pi_thresh, na.rm = TRUE)
    )
  })) |>
    pivot_longer(cols = c(pval_n, adjp_n, pival_n),
                 names_to = "method", values_to = "n") |>
    mutate(
      method     = factor(method,
                          levels = c("pval_n", "adjp_n", "pival_n"),
                          labels = c("p-value", "FDR", "\u03a0-value")),
      comparison = factor(comparison, levels = names(limma_results[[tissue]])),
      tissue     = str_to_title(tissue)
    )
}

da_long_combined <- bind_rows(
  make_da_long("muscle"),
  make_da_long("plasma")
) |>
  mutate(
    tissue     = factor(tissue, levels = c("Muscle", "Plasma")),
    fill_color = BG_COLORS[as.character(tissue)]
  )

# --- 2b: Total quantified proteins per tissue (reference bar) ---------------
# Universe size for each tissue -- row count is identical across all 3
# contrasts within a tissue (same limma model, only the contrast differs),
# so any one sheet's row count is that tissue's total quantified proteins.
# Using "Cancer" here since it's the one contrast name guaranteed present
# and consistently spelled across both tissues' limma_results lists
# (the old "Baseline_SURVvCTL" sheet name this used to reference doesn't
# exist anymore under the Cancer/Training/Residual contrast naming).
#
# Rendered as a SEPARATE geom_col() layer at its own x position ("Total
# Quantified", added as a new level to the comparison factor) rather than
# folded into the p-value/FDR/Pi-value dodge trio -- it's a different kind
# of quantity (universe size vs. a significant subset), and cramming it
# into that dodge scheme would put it in one of 3 slots designed for a
# different 3-way comparison rather than reading as a standalone reference.
#
# NOTE ON SCALE: total quantified will be in the thousands while DAP counts
# are typically tens-to-hundreds. On a shared linear y-axis (kept below,
# since that's the literal ask -- see it in context against the universe)
# the DAP bars will look small next to this one. If that's not what you
# want once you see it rendered, the one-line fix is adding
# trans = "sqrt" (or "log10") to the scale_y_continuous() call in section 3.
total_quantified <- c(
  Muscle = nrow(limma_results$muscle$Cancer),
  Plasma = nrow(limma_results$plasma$Cancer)
)

COMPARISON_LEVELS <- c("Total Quantified", levels(da_long_combined$comparison))

da_long_combined <- da_long_combined |>
  mutate(comparison = factor(as.character(comparison), levels = COMPARISON_LEVELS))

total_df <- tibble(
  comparison = factor("Total Quantified", levels = COMPARISON_LEVELS),
  n          = unname(total_quantified[c("Muscle", "Plasma")]),
  tissue     = factor(c("Muscle", "Plasma"), levels = c("Muscle", "Plasma"))
) |>
  mutate(fill_color = BG_COLORS[as.character(tissue)])

# --- 3: DA counts plot -------------------------------------------------------
n_methods   <- length(levels(da_long_combined$method))
dodge_width <- 0.7
bar_width   <- dodge_width / n_methods

text_data_da <- da_long_combined |>
  group_by(tissue, comparison) |>
  mutate(
    rank_in_full = as.integer(method),
    x_dodge      = as.numeric(comparison) +
      (rank_in_full - (n_methods + 1) / 2) * bar_width
  ) |>
  ungroup() |>
  group_by(tissue) |>
  mutate(offset = max(n) * 0.03) |>
  ungroup()

# Text label for the Total Quantified reference bar -- own offset (its own
# data's max n), independent of text_data_da's offset which is scaled to
# the much smaller DAP counts.
text_data_total <- total_df |>
  group_by(tissue) |>
  mutate(offset = n * 0.03) |>
  ungroup()

p_da_combined <- ggplot(da_long_combined, aes(x = comparison, y = n)) +
  geom_rect(
    data        = da_long_combined |> distinct(tissue, fill_color),
    aes(fill    = fill_color),
    xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
    inherit.aes = FALSE
  ) +
  scale_fill_identity() +
  new_scale_fill() +
  # Total Quantified reference bar -- literal (non-aes) fill so it isn't
  # pulled into the method legend/scale_fill_manual() below; drawn at
  # dodge_width so it visually matches the full cluster width of the 3
  # dodged bars at the other x positions rather than looking undersized.
  geom_col(
    data        = total_df,
    aes(x = comparison, y = n),
    fill        = "grey60", colour = "black", linewidth = 0.3,
    width       = dodge_width,
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = text_data_total,
    aes(x = comparison, y = n + offset, label = n),
    vjust       = 0, size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_col(aes(fill = method),
           position = position_dodge(width = dodge_width), width = 0.6) +
  geom_text(
    data        = filter(text_data_da, n > 0),
    aes(x = x_dodge, y = n + offset, label = n),
    vjust       = 0, size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  geom_text(
    data        = filter(text_data_da, n == 0),
    aes(x = x_dodge, y = offset, label = "0"),
    vjust       = 0, size = 3, fontface = "bold",
    inherit.aes = FALSE
  ) +
  # Tissue label in top-left of each panel
  geom_text(
    data        = tissue_label_df,
    aes(x = -Inf, y = Inf, label = label, colour = label),
    hjust       = -0.15, vjust = 1.4,
    size        = 4, fontface = "bold",
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = tissue_colors, guide = "none") +
  scale_fill_manual(values = c(colors$blue_light, colors$blue, colors$blue_dark)) +
  scale_x_discrete(labels = comp_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_wrap(~ tissue, ncol = 1) +
  labs(x = NULL, y = "DAPs (n)", title = "DAPs per Contrast", fill = NULL) +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(
    legend.position      = c(0.97, 0.97),
    legend.justification = c(1, 1),
    legend.background    = element_rect(fill = NA, linewidth = 0),
    strip.text           = element_blank(),
    strip.background     = element_blank(),
    panel.spacing        = unit(1, "lines")
  )

# --- 4: UpSet plot -----------------------------------------------------------

# Pulled to the top level (out of make_combined_upset()) so both the
# per-tissue plot builder and the Section 8 export can reuse it without
# duplicating the intersection logic.
get_upset_components <- function(tissue, limma_results, pi_thresh, SET_LABELS) {
  results   <- limma_results[[tissue]]
  CONTRASTS <- names(results)
  
  sig_sets <- list()
  dir_map  <- list()
  
  for (ctr in CONTRASTS) {
    df      <- results[[ctr]]
    is_sig  <- !is.na(df$Pi.Val) & df$Pi.Val < pi_thresh
    sig_sets[[ctr]] <- df$uniprot[is_sig]
    dir_map[[ctr]]  <- setNames(
      ifelse(df$logFC[is_sig] > 0, "Up", "Down"),
      df$uniprot[is_sig]
    )
  }
  
  comb_grid <- expand.grid(
    Cancer = c(TRUE, FALSE),
    Residual = c(TRUE, FALSE),
    Training     = c(TRUE, FALSE)
  ) |>
    filter(rowSums(across(everything())) > 0)
  
  intersection_data <- bind_rows(lapply(seq_len(nrow(comb_grid)), function(i) {
    in_sets  <- CONTRASTS[which(as.logical(comb_grid[i, ]))]
    out_sets <- CONTRASTS[which(!as.logical(comb_grid[i, ]))]
    members  <- Reduce(intersect, sig_sets[in_sets])
    if (length(out_sets) > 0)
      members <- setdiff(members, unique(na.omit(unlist(sig_sets[out_sets]))))
    if (length(members) == 0) return(NULL)
    dirs <- sapply(members, function(prot) {
      d <- sapply(in_sets, \(ctr) dir_map[[ctr]][prot])
      d <- d[!is.na(d)]
      if (all(d == "Up")) "Up" else if (all(d == "Down")) "Down" else "Mixed"
    })
    tibble(
      in_sets   = list(in_sets), degree = length(in_sets),
      up        = sum(dirs == "Up"),
      down      = sum(dirs == "Down"),
      mixed     = sum(dirs == "Mixed"),
      total     = length(members),
      # protein-level membership + direction for this intersection, kept
      # as list-columns so identity-level exports (Section 8) can unnest
      # them without recomputing the intersection logic
      members    = list(unname(members)),
      directions = list(unname(dirs))
    )
  })) |>
    filter(total > 0) |>
    arrange(desc(total)) |>
    mutate(x = row_number())
  
  # long-format protein-level table: one row per (intersection, protein),
  # with its direction and which contrasts define that intersection
  members_long <- intersection_data |>
    dplyr::select(x, in_sets, members, directions) |>
    mutate(uniprot = members, direction = directions) |>
    dplyr::select(-members, -directions) |>
    unnest(c(uniprot, direction))
  
  bar_long <- intersection_data |>
    mutate(mixed_only = up == 0 & down == 0 & mixed > 0) |>
    dplyr::select(x, up, down, mixed_only, total) |>
    pivot_longer(c(up, down), names_to = "direction", values_to = "count") |>
    mutate(
      direction = factor(str_to_title(direction), levels = c("Down", "Up")),
      count     = ifelse(mixed_only, 0L, count)
    )
  
  mixed_bar <- intersection_data |>
    filter(up == 0 & down == 0 & mixed > 0) |>
    dplyr::select(x, total)
  
  dot_df <- intersection_data |>
    dplyr::select(x, in_sets) |>
    mutate(set = map(in_sets, \(s) tibble(set = CONTRASTS, active = CONTRASTS %in% s))) |>
    unnest(set) |>
    mutate(
      set_label = factor(SET_LABELS[set], levels = rev(unname(SET_LABELS))),
      ynum      = as.numeric(set_label)
    )
  
  seg_df <- dot_df |>
    filter(active) |>
    group_by(x) |>
    filter(n() > 1) |>
    summarise(ymin = min(ynum), ymax = max(ynum), .groups = "drop")
  
  list(intersection_data = intersection_data, bar_long = bar_long,
       mixed_bar = mixed_bar, dot_df = dot_df, seg_df = seg_df,
       members_long = members_long, CONTRASTS = CONTRASTS)
}

make_upset_plot <- function(tissue, title = NULL) {
  components        <- get_upset_components(tissue, limma_results, pi_thresh, SET_LABELS)
  intersection_data <- components$intersection_data
  bar_long          <- components$bar_long
  mixed_bar         <- components$mixed_bar
  dot_df            <- components$dot_df
  seg_df            <- components$seg_df
  CONTRASTS         <- components$CONTRASTS
  
  y_max           <- ceiling(max(c(bar_long$count, mixed_bar$total), na.rm = TRUE) * 1.15 / 10) * 10
  max_label_chars <- max(nchar(unname(SET_LABELS)))
  left_pad        <- max_label_chars * 5.5
  n_dirs          <- 2
  bar_width       <- 0.7 / n_dirs
  
  text_data_bars <- bar_long |>
    filter(count > 0) |>
    group_by(x) |>
    mutate(
      rank_in_full = as.integer(direction),
      x_dodge      = x + (rank_in_full - (n_dirs + 1) / 2) * bar_width,
      offset       = y_max * 0.03
    ) |>
    ungroup()
  
  p_bars <- ggplot(bar_long, aes(x = x, y = count, fill = direction)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6, linewidth = 0.3) +
    geom_col(
      data        = mixed_bar,
      aes(x = x, y = total, fill = "Mixed"),
      width       = 0.6,
      inherit.aes = FALSE
    ) +
    geom_text(
      data     = mixed_bar,
      aes(x = x, y = total + y_max * 0.03, label = total),
      vjust    = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = text_data_bars,
      aes(x = x_dodge, y = count + offset, label = count),
      vjust = 0, size = 3, inherit.aes = FALSE
    ) +
    scale_fill_manual(values = DIR_COLORS) +
    scale_x_continuous(breaks = intersection_data$x,
                       expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0)),
                       breaks = scales::breaks_pretty(n = 4)) +
    labs(y = "Intersection size", x = NULL, fill = NULL, title = title) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.position      = c(0.97, 0.97),
      legend.justification = c(1, 1),
      plot.margin          = margin(4, 4, 0, left_pad)
    )
  
  n_sets <- length(CONTRASTS)
  
  p_dots <- ggplot() +
    annotate("rect",
             xmin = 0.5, xmax = max(intersection_data$x) + 0.5,
             ymin = seq_len(n_sets) - 0.45, ymax = seq_len(n_sets) + 0.45,
             fill = rep(c("grey96", "grey90"), length.out = n_sets)) +
    {if (nrow(seg_df) > 0)
      geom_segment(data = seg_df,
                   aes(x = x, xend = x, y = ymin, yend = ymax),
                   linewidth = 0.7, color = "grey25")} +
    geom_point(data = filter(dot_df, !active),
               aes(x = x, y = ynum), color = "grey82", size = 2.5) +
    geom_point(data = filter(dot_df, active),
               aes(x = x, y = ynum), color = "grey15", size = 2.5) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = seq_len(n_sets),
                       labels = levels(dot_df$set_label),
                       expand = expansion(add = 0.4)) +
    labs(x = NULL, y = NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x  = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
      axis.text.y  = element_text(size = 9, face = "bold"),
      plot.margin  = margin(0, 4, 4, 4)
    )
  
  (p_bars / p_dots) + plot_layout(heights = c(0.75, 0.25))
}

p_upset_m <- make_upset_plot("muscle", "Contrast Overlap (UpSet) \u2014 Muscle")
p_upset_p <- make_upset_plot("plasma",  "Contrast Overlap (UpSet) \u2014 Plasma")

# --- Combined UpSet ----------------------------------------------------------
make_combined_upset <- function(limma_results, pi_thresh, SET_LABELS, DIR_COLORS) {
  
  m      <- get_upset_components("muscle", limma_results, pi_thresh, SET_LABELS)
  p      <- get_upset_components("plasma", limma_results, pi_thresh, SET_LABELS)
  offset <- max(m$intersection_data$x)
  
  shift_x    <- function(df) df |> mutate(x = x + offset)
  shift_x_mb <- function(df) if (nrow(df) > 0) df |> mutate(x = x + offset) else df
  
  p$intersection_data <- shift_x(p$intersection_data)
  p$bar_long          <- shift_x(p$bar_long)
  p$dot_df            <- shift_x(p$dot_df)
  p$seg_df            <- shift_x(p$seg_df)
  p$mixed_bar         <- shift_x_mb(p$mixed_bar)
  
  bar_all <- bind_rows(
    m$bar_long |> mutate(source = "Muscle"),
    p$bar_long |> mutate(source = "Plasma")
  )
  mixed_bar_all <- bind_rows(
    m$mixed_bar |> mutate(source = "Muscle"),
    p$mixed_bar |> mutate(source = "Plasma")
  )
  dot_all <- bind_rows(m$dot_df, p$dot_df)
  seg_all <- bind_rows(m$seg_df, p$seg_df)
  
  y_max <- ceiling(
    max(c(bar_all$count, mixed_bar_all$total), na.rm = TRUE) * 1.15 / 10
  ) * 10
  
  bg_df <- tibble(
    xmin  = c(0.5,          offset + 0.5),
    xmax  = c(offset + 0.5, max(bar_all$x) + 0.5),
    color = unname(BG_COLORS[c("Muscle", "Plasma")])
  )
  
  n_dirs    <- 2
  bar_width <- 0.7 / n_dirs
  
  text_df <- bar_all |>
    filter(count > 0) |>
    group_by(x) |>
    mutate(
      rank  = as.integer(direction),
      x_pos = x + (rank - (n_dirs + 1) / 2) * bar_width,
      y_pos = count + y_max * 0.02
    ) |>
    ungroup()
  
  p_bars <- ggplot(bar_all, aes(x = x, y = count, fill = direction)) +
    geom_rect(
      data        = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf),
      fill        = bg_df$color,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = offset + 0.5, linewidth = 0.6, color = "grey40") +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    {if (nrow(mixed_bar_all) > 0) list(
      geom_col(
        data        = mixed_bar_all,
        aes(x = x, y = total, fill = "Mixed"),
        width       = 0.6,
        inherit.aes = FALSE
      ),
      geom_text(
        data        = mixed_bar_all,
        aes(x = x, y = total + y_max * 0.03, label = total),
        vjust = 0, size = 3, fontface = "bold",
        inherit.aes = FALSE
      )
    )} +
    geom_text(
      data        = text_df,
      aes(x = x_pos, y = y_pos, label = count),
      size        = 3, vjust = 0, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = DIR_COLORS) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0))) +
    labs(y = "Intersection Size", x = NULL, fill = NULL,
         title = "Contrast Overlaps (UpSet)") +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x          = element_blank(),
      axis.ticks.x         = element_blank(),
      panel.grid.major.x   = element_blank(),
      legend.position      = c(0.90, 0.98),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = scales::alpha("white", 0.0),
                                          color = NA)
    )
  
  n_sets <- length(m$CONTRASTS)
  
  p_dots <- ggplot() +
    geom_rect(
      data        = bg_df,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill        = bg_df$color,
      inherit.aes = FALSE
    ) +
    geom_vline(xintercept = offset + 0.5, linewidth = 0.6, color = "grey40") +
    geom_segment(data = seg_all,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.7, color = "grey25") +
    geom_point(data = filter(dot_all, !active),
               aes(x = x, y = ynum), color = "grey82", size = 2.5) +
    geom_point(data = filter(dot_all, active),
               aes(x = x, y = ynum), color = "grey15", size = 2.5) +
    annotate("text",
             x = (offset + 0.5) / 2, y = 0.1, label = "Muscle",
             fontface = "bold", size = 3.5,
             color = tissue_colors["Muscle"], hjust = 0.5, vjust = 1) +
    annotate("text",
             x = (offset + 0.5 + max(bar_all$x) + 0.5) / 2,
             y = 0.1, label = "Plasma",
             fontface = "bold", size = 3.5,
             color = tissue_colors["Plasma"], hjust = 0.5, vjust = 1) +
    scale_x_continuous(expand = expansion(add = 0.5)) +
    scale_y_continuous(breaks = seq_len(n_sets),
                       labels = levels(dot_all$set_label),
                       expand = expansion(add = 0.4)) +
    labs(x = NULL, y = NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x  = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks   = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(0, 4, 14, 4)
    )
  
  (p_bars / p_dots) + plot_layout(heights = c(0.75, 0.25))
}

p_upset_c <- make_combined_upset(
  limma_results = limma_results,
  pi_thresh     = pi_thresh,
  SET_LABELS    = SET_LABELS,
  DIR_COLORS    = DIR_COLORS
)

# --- 5: GSEA counts plot -----------------------------------------------------
# Now sourced from gsea_supp_full (consolidated) instead of raw gsea_muscle/
# gsea_plasma RDS lists. Counting logic is controlled by GSEA_COUNT_MODE:
#   "representative" -- a pathway counts only if it's both significant in
#     THIS comparison AND flagged representative for its (tissue, database)
#     by the pooled redundancy filter (representative_lookup). The
#     redundancy flag itself is pooled across all three gsea_comps (see
#     02_gsea_pathway_consolidation.R), but whether it counts toward THIS
#     bar is still gated on significance in this specific comparison -- so
#     a pathway significant in two comparisons still contributes to both
#     bars, it just won't be double-counted against a near-duplicate
#     pathway within the same bar.
#   "raw" -- original behavior: every significant pathway counts,
#     redundant or not.
#
# Split into a data-prep half (build_gsea_counts) and a plot-building half
# (make_gsea_counts_combined), matching this project's separation
# convention -- also lets Section 8 reuse the same counts table for export.
build_gsea_counts <- function() {
  base <- gsea_supp_full %>%
    filter(comparison %in% gsea_comps, significant)
  
  if (GSEA_COUNT_MODE == "representative") {
    base <- base %>%
      semi_join(representative_lookup, by = c("tissue", "database", "ID"))
  } else if (GSEA_COUNT_MODE != "raw") {
    stop("GSEA_COUNT_MODE must be 'representative' or 'raw'")
  }
  
  base %>%
    mutate(tissue = as.character(tissue), database = as.character(database)) %>%
    group_by(tissue, comparison, database) %>%
    summarise(n = n(), .groups = "drop") %>%
    # ensure every tissue x comparison x database combo appears (0 counts
    # included), same as the old per-slot loop always emitting a row
    tidyr::complete(
      tissue     = c("Muscle", "Plasma"),
      comparison = gsea_comps,
      database   = unname(db_slots),
      fill       = list(n = 0L)
    ) |>
    mutate(
      comparison = factor(comparison, levels = gsea_comps),
      database   = factor(database,   levels = unname(db_slots)),
      tissue     = factor(tissue,     levels = c("Muscle", "Plasma")),
      fill_color = BG_COLORS[as.character(tissue)]
    )
}

make_gsea_counts_combined <- function(counts) {
  n_databases <- length(levels(counts$database))
  bar_width   <- 0.7 / n_databases
  
  text_data_gsea <- counts |>
    group_by(tissue, comparison) |>
    mutate(
      rank_in_full = as.integer(database),
      x_dodge      = as.numeric(comparison) +
        (rank_in_full - (n_databases + 1) / 2) * bar_width
    ) |>
    ungroup() |>
    group_by(tissue) |>
    mutate(offset = max(n) * 0.03) |>
    ungroup()
  
  ggplot(counts, aes(x = comparison, y = n)) +
    geom_rect(
      data        = counts |> distinct(tissue, fill_color),
      aes(fill    = fill_color),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_col(aes(fill = database),
             position  = position_dodge(width = 0.7), width = 0.6,
             color     = "black", linewidth = 0.3) +
    geom_text(
      data        = filter(text_data_gsea, n > 0),
      aes(x = x_dodge, y = n + offset, label = n),
      vjust       = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    geom_text(
      data        = filter(text_data_gsea, n == 0),
      aes(x = x_dodge, y = offset, label = "0"),
      vjust       = 0, size = 3, fontface = "bold",
      inherit.aes = FALSE
    ) +
    # Tissue label in top-left of each panel
    geom_text(
      data        = tissue_label_df,
      aes(x = -Inf, y = Inf, label = label, colour = label),
      hjust       = -0.15, vjust = 1.4,
      size        = 4, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = tissue_colors, guide = "none") +
    scale_fill_manual(values = db_colors) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    facet_wrap(~ tissue, ncol = 1) +
    labs(x = NULL, y = "Significant Pathways (n)",
         title = "Pathway Enrichment", fill = NULL,
         subtitle = if (GSEA_COUNT_MODE == "representative")
           "Redundant pathways (Jaccard/overlap-coefficient filtered) excluded"
         else NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      panel.grid.major.x   = element_blank(),
      panel.grid.minor     = element_blank(),
      legend.key.size      = unit(0.4, "cm"),
      legend.text          = element_text(size = 7),
      legend.position      = c(0.97, 0.46),
      legend.justification = c(1, 1),
      legend.background    = element_rect(fill = NA, linewidth = 0.3),
      strip.text           = element_blank(),
      strip.background     = element_blank(),
      panel.spacing        = unit(1, "lines"),
      plot.subtitle         = element_text(size = 8, face = "plain",
                                           hjust = 0, color = colors$gray_dark)
    )
}

gsea_counts     <- build_gsea_counts()
p_gsea_combined <- make_gsea_counts_combined(gsea_counts)

# --- 6: CV Plot --------------------------------------------------------------

build_cv_data <- function(data, annot, meta) {
  
  group_map <- bind_rows(
    meta %>% filter(cancer == "CTL")         %>% mutate(group = "CTL"),
    meta %>% filter(cancer_time == "CRC_T1") %>% mutate(group = "CRC"),
    meta %>% filter(cancer_time == "CRC_T2") %>% mutate(group = "CRC+RET"),
    meta %>% filter(supp_time  == "CRE_T1")  %>% mutate(group = "CRE_PRE"),
    meta %>% filter(supp_time  == "CRE_T2")  %>% mutate(group = "CRE_POST"),
    meta %>% filter(supp_time  == "PLA_T1")  %>% mutate(group = "PLA_PRE"),
    meta %>% filter(supp_time  == "PLA_T2")  %>% mutate(group = "PLA_POST")
  ) %>%
    dplyr::select(sample_id, group) %>%
    filter(!is.na(group))
  
  rn   <- rownames(data)
  data <- as.data.frame(base::apply(data, 2, as.numeric))
  rownames(data) <- rn
  
  data %>%
    rownames_to_column("uniprot_id") %>%
    left_join(annot %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
    filter(!is.na(gene)) %>%
    dplyr::select(gene, any_of(group_map$sample_id)) %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
    left_join(group_map, by = "sample_id", relationship = "many-to-many") %>%
    filter(!is.na(group), !is.na(intensity)) %>%
    group_by(gene, group) %>%
    summarise(
      mean_int = mean(intensity, na.rm = TRUE),
      sd_int   = sd(intensity,   na.rm = TRUE),
      n        = sum(!is.na(intensity)),
      cv       = (sd_int / mean_int) * 100,
      .groups  = "drop"
    ) %>%
    filter(mean_int > 0, n >= 3)
}

CRCvCTL_GROUPS <- c("CTL", "CRC", "CRC+RET")

GROUP_COLOURS <- c(
  CTL      = colors$green,
  CRC  = colors$blue,
  `CRC+RET` = colors$red
)

cv_df_m <- build_cv_data(
  data  = norm_m$data,
  annot = norm_m$annotation,
  meta  = norm_m$metadata
) %>% filter(group %in% CRCvCTL_GROUPS)

cv_df_p <- build_cv_data(
  data  = norm_p$data,
  annot = norm_p$annotation,
  meta  = norm_p$metadata
) %>% filter(group %in% CRCvCTL_GROUPS)

# --- 6b: Bootstrap median CV with 95% CI ------------------------------------
bootstrap_median_cv <- function(cv_df, groups, n_boot = 2000, seed = 42) {
  set.seed(seed)
  
  bind_rows(lapply(groups, function(grp) {
    vals <- cv_df$cv[cv_df$group == grp]
    vals <- vals[!is.na(vals)]
    
    boot_medians <- replicate(n_boot, median(sample(vals, replace = TRUE)))
    
    tibble(
      group  = grp,
      median = median(vals),
      ci_lo  = quantile(boot_medians, 0.025),
      ci_hi  = quantile(boot_medians, 0.975),
      n      = length(vals)
    )
  }))
}

boot_m <- bootstrap_median_cv(cv_df_m, CRCvCTL_GROUPS)
boot_p <- bootstrap_median_cv(cv_df_p, CRCvCTL_GROUPS)

cat("\nMuscle — bootstrap median CV (95% CI):\n")
print(boot_m)
cat("\nPlasma — bootstrap median CV (95% CI):\n")
print(boot_p)

# --- 6c: Brown-Forsythe global test (for methods reporting) ------------------
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
library(car)

bf_m <- car::leveneTest(
  cv ~ factor(group, levels = CRCvCTL_GROUPS),
  data   = cv_df_m |> filter(group %in% CRCvCTL_GROUPS),
  center = median
)

bf_p <- car::leveneTest(
  cv ~ factor(group, levels = CRCvCTL_GROUPS),
  data   = cv_df_p |> filter(group %in% CRCvCTL_GROUPS),
  center = median
)

cat("\nMuscle — Brown-Forsythe test (equal CV spread across groups):\n")
print(bf_m)
cat("\nPlasma — Brown-Forsythe test (equal CV spread across groups):\n")
print(bf_p)

for (tissue in c("Muscle", "Plasma")) {
  bf    <- if (tissue == "Muscle") bf_m else bf_p
  fstat <- round(bf$`F value`[1], 2)
  pval  <- bf$`Pr(>F)`[1]
  psym  <- ifelse(pval < 0.001, "p < 0.001",
                  ifelse(pval < 0.01,  sprintf("p = %.3f", pval),
                         sprintf("p = %.3f", pval)))
  cat(sprintf("%s Brown-Forsythe: F = %s, %s\n", tissue, fstat, psym))
}

# --- 6d: Violin fill palette — tissue-specific shades per group --------------
make_shades <- function(base_color, n = 3, lightest = 0.3, darkest = 0.85) {
  colorRampPalette(c("white", base_color))(100)[
    round(seq(lightest, darkest, length.out = n) * 100)
  ]
}

muscle_shades <- make_shades(colors$pink, n = 3)
plasma_shades <- make_shades(colors$aqua, n = 3)

# violin_fills <- c(
#   Muscle_CTL      = muscle_shades[1],
#   Muscle_CRC  = muscle_shades[2],
#   `Muscle_CRC+RET` = muscle_shades[3],
#   Plasma_CTL      = plasma_shades[1],
#   Plasma_CRC  = plasma_shades[2],
#   `Plasma_CRC+RET` = plasma_shades[3]
# )

violin_fills <- c(
  Muscle_CTL      = colors$green_light,
  Muscle_CRC  = colors$blue_light,
  `Muscle_CRC+RET` = colors$red_light,
  Plasma_CTL      = colors$green_light,
  Plasma_CRC  = colors$blue_light,
  `Plasma_CRC+RET` = colors$red_light
)

# --- 6e: CV violin plot ------------------------------------------------------
make_cv_violin_combined <- function(cv_df_m, cv_df_p, boot_m, boot_p) {
  
  cv_combined <- bind_rows(
    cv_df_m |> mutate(tissue = "Muscle"),
    cv_df_p |> mutate(tissue = "Plasma")
  ) |>
    mutate(
      tissue       = factor(tissue, levels = c("Muscle", "Plasma")),
      fill_color   = BG_COLORS[as.character(tissue)],
      group        = factor(group, levels = CRCvCTL_GROUPS),
      group_tissue = paste0(tissue, "_", group)
    )
  
  boot_combined <- bind_rows(
    boot_m |> mutate(tissue = "Muscle"),
    boot_p |> mutate(tissue = "Plasma")
  ) |>
    mutate(
      tissue = factor(tissue, levels = c("Muscle", "Plasma")),
      group  = factor(group,  levels = CRCvCTL_GROUPS)
    )
  
  # label_y positioned at a uniform height above each tissue panel's max
  # CV value (y_max + range * 0.08) -- same approach as make_cv_violin()
  # in the normalization reports scripts, rather than hovering each label
  # over its own violin's 95th percentile.
  violin_range <- cv_combined |>
    filter(group %in% CRCvCTL_GROUPS) |>
    group_by(tissue) |>
    summarise(y_max   = max(cv, na.rm = TRUE),
              y_range = diff(range(cv, na.rm = TRUE)),
              .groups = "drop")
  
  label_df <- boot_combined |>
    left_join(violin_range, by = "tissue") |>
    mutate(
      label_y   = y_max + y_range * 0.08,
      label_txt = sprintf("%.1f [%.1f\u2013%.1f]", median, ci_lo, ci_hi)
    )
  
  ggplot(cv_combined, aes(x = group, y = cv)) +
    geom_rect(
      data        = cv_combined |> distinct(tissue, fill_color),
      aes(fill    = fill_color),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    new_scale_fill() +
    geom_violin(aes(fill = group_tissue), alpha = 0.75,
                linewidth = 0.4, trim = FALSE) +
    geom_boxplot(aes(fill = group_tissue), width = 0.12, alpha = 0.95,
                 linewidth = 0.4, outlier.shape = NA,
                 colour = colors$gray_dark) +
    geom_point(
      data        = boot_combined,
      aes(x = group, y = median),
      inherit.aes = FALSE,
      size        = 2.5,
      shape       = 18,
      color       = colors$black
    ) +
    geom_errorbar(
      data        = boot_combined,
      aes(x = group, ymin = ci_lo, ymax = ci_hi),
      inherit.aes = FALSE,
      width       = 0.12,
      linewidth   = 0.8,
      color       = colors$black
    ) +
    geom_label(
      data        = label_df,
      aes(x = group, y = label_y, label = label_txt),
      inherit.aes = FALSE,
      size        = 2.8, lineheight = 0.9, fontface = "bold",
      color       = colors$gray_dark, fill = colors$white, label.size = 0.2
    ) +
    geom_text(
      data        = tissue_label_df,
      aes(x = -Inf, y = Inf, label = label, colour = label),
      hjust       = -0.15, vjust = 1.4,
      size        = 4, fontface = "bold",
      inherit.aes = FALSE
    ) +
    scale_colour_manual(values = tissue_colors, guide = "none") +
    scale_fill_manual(values = violin_fills, guide = "none") +
    scale_x_discrete(labels = c(CTL = "CTL", CRC = "CRC",
                                `CRC+RET` = "CRC+RET")) +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0, 0.15))) +
    facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
    labs(
      x       = NULL,
      y       = "CV (%)",
      title   = "Per-Group CV \u2014 Normalized",
      subtitle = "Median CV (95% bootstrap CI, n = 2,000 resamples)"
    ) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      axis.text.x      = element_text(face = "bold", size = 9),
      strip.text       = element_blank(),
      strip.background = element_blank(),
      panel.spacing    = unit(1, "lines"),
      plot.subtitle     = element_text(size = 8, face = "plain",
                                       hjust = 0, color = colors$gray_dark)
    )
}

p_cv_plot_combined <- make_cv_violin_combined(cv_df_m, cv_df_p, boot_m, boot_p)

# --- 7: Assemble & Save ------------------------------------------------------
f2 <- plot_grid(
  p_da_combined, p_gsea_combined, p_upset_c, p_cv_plot_combined,
  labels = "AUTO",
  nrow   = 2,
  ncol   = 2
)

ggsave(file.path(config$fig_dir, "fig_02.pdf"),
       f2, width = 16, height = 12, units = "in", device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02.png"),
       f2, width = 16, height = 12, units = "in", dpi = 300)

ggsave(file.path(config$fig_dir, "fig_02_A.pdf"),          p_da_combined,      width = 6,  height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_B.pdf"),          p_gsea_combined,    width = 6,  height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_C.pdf"),          p_upset_c,          width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_C_muscle.pdf"),   p_upset_m,          width = 8,  height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_C_plasma.pdf"),   p_upset_p,          width = 8,  height = 8, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_D.pdf"),          p_cv_plot_combined, width = 4,  height = 6, device = cairo_pdf)

ggsave(file.path(config$fig_dir, "fig_02_A.png"),          p_da_combined,      width = 6,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_B.png"),          p_gsea_combined,    width = 6,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_C.png"),          p_upset_c,          width = 10, height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_C_muscle.png"),   p_upset_m,          width = 8,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_C_plasma.png"),   p_upset_p,          width = 8,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_D.png"),          p_cv_plot_combined, width = 4,  height = 6,  dpi = 300)

# --- 8: Excel Data Export ----------------------------------------------------
# Companion workbook of the numbers behind each panel, one summary sheet
# per Panel x tissue (the plotted counts), plus an identity-level sheet
# for Panels A-C giving the actual protein/pathway names, direction, and
# supporting stats behind those counts. A README sheet describes all of it.

# --- Panel A: DA counts (per-method significant-protein counts) plus the
#     Total Quantified reference bar, by contrast
panel_a_export <- bind_rows(
  da_long_combined |>
    dplyr::select(tissue, comparison, method, n) |>
    mutate(method = as.character(method)),
  total_df |>
    dplyr::select(tissue, comparison, n) |>
    mutate(method = "Total Quantified")
) |>
  mutate(comparison = factor(as.character(comparison), levels = COMPARISON_LEVELS)) |>
  arrange(tissue, comparison, method)

panel_a_muscle <- panel_a_export |> filter(tissue == "Muscle") |> dplyr::select(-tissue)
panel_a_plasma <- panel_a_export |> filter(tissue == "Plasma") |> dplyr::select(-tissue)

# --- Panel A (identity-level): the significant proteins behind those counts.
#     A protein is included if it clears ANY of the three thresholds; the
#     sig_pvalue/sig_fdr/sig_pival flags show which one(s) it actually
#     cleared, and Pi.Val is this project's primary significance metric
#     (Xiao et al. 2014, p^|log2FC|) if you want to filter down further.
build_panel_a_proteins <- function(tissue) {
  bind_rows(lapply(names(limma_results[[tissue]]), function(comp) {
    df <- limma_results[[tissue]][[comp]]
    df |>
      mutate(
        comparison = comp,
        direction  = ifelse(logFC > 0, "Up", "Down"),
        sig_pvalue = P.Value   < thresh,
        sig_fdr    = adj.P.Val < adjp_thresh,
        sig_pival  = Pi.Val    < pi_thresh
      ) |>
      filter(sig_pvalue | sig_fdr | sig_pival) |>
      dplyr::select(any_of(c("uniprot", "gene")), comparison, direction,
                    logFC, P.Value, adj.P.Val, Pi.Val,
                    sig_pvalue, sig_fdr, sig_pival)
  })) |>
    mutate(comparison = factor(comparison, levels = names(limma_results[[tissue]]))) |>
    arrange(comparison, Pi.Val)
}

panel_a_muscle_sig <- build_panel_a_proteins("muscle")
panel_a_plasma_sig <- build_panel_a_proteins("plasma")

# --- Panel B: GSEA pathway counts per comparison x database
panel_b_muscle <- gsea_counts |>
  filter(tissue == "Muscle") |>
  dplyr::select(comparison, database, n) |>
  arrange(comparison, database)

panel_b_plasma <- gsea_counts |>
  filter(tissue == "Plasma") |>
  dplyr::select(comparison, database, n) |>
  arrange(comparison, database)

# --- Panel B (identity-level): the significant pathways behind those
#     counts, respecting GSEA_COUNT_MODE (redundant near-duplicates
#     excluded when "representative", same as the bars).
#
#     NOTE: column names below (ID/Description/NES/pvalue/p.adjust) are
#     assumed from a typical clusterProfiler/fgsea-style GSEA output --
#     VERIFY against str(gsea_supp_full) and adjust the any_of() list if
#     your actual column names differ (e.g. "pathway" vs "Description",
#     "padj" vs "p.adjust"). any_of() means this won't error on a
#     mismatch, it will just silently omit that column, so it's worth
#     checking the output has what you expect.
build_panel_b_pathways <- function(tissue_name) {
  base <- gsea_supp_full |>
    filter(tissue == tissue_name, comparison %in% gsea_comps, significant)
  
  if (GSEA_COUNT_MODE == "representative") {
    base <- base |>
      semi_join(representative_lookup, by = c("tissue", "database", "ID"))
  }
  
  if ("NES" %in% names(base)) {
    base <- base |> mutate(direction = ifelse(NES > 0, "Up", "Down"))
  } else {
    base <- base |> mutate(direction = NA_character_)
  }
  
  base |>
    mutate(comparison = factor(comparison, levels = gsea_comps)) |>
    dplyr::select(comparison, database,
                  any_of(c("ID", "Description", "pathway", "NES",
                           "pvalue", "p.adjust", "padj")),
                  direction) |>
    arrange(comparison, database)
}

panel_b_muscle_sig <- build_panel_b_pathways("Muscle")
panel_b_plasma_sig <- build_panel_b_pathways("Plasma")

# --- Panel C: UpSet intersection sizes (protein-level contrast overlaps)
upset_m_data <- get_upset_components("muscle", limma_results, pi_thresh, SET_LABELS)
upset_p_data <- get_upset_components("plasma", limma_results, pi_thresh, SET_LABELS)

format_intersection <- function(upset_components) {
  upset_components$intersection_data |>
    mutate(contrasts = map_chr(in_sets, ~ paste(SET_LABELS[.x], collapse = " & "))) |>
    dplyr::select(contrasts, degree, up, down, mixed, total) |>
    arrange(desc(total))
}

panel_c_muscle <- format_intersection(upset_m_data)
panel_c_plasma <- format_intersection(upset_p_data)

# --- Panel C (identity-level): the actual proteins in each contrast
#     overlap, with their direction (Up/Down/Mixed across the contrasts
#     that define that overlap). Gene symbol joined in from the same
#     normalized-data annotation table used for the CV panel; if a
#     uniprot ID has no gene mapping, `gene` comes back NA rather than
#     dropping the row.
format_intersection_proteins <- function(upset_components, annot, SET_LABELS) {
  upset_components$members_long |>
    mutate(contrasts = map_chr(in_sets, ~ paste(SET_LABELS[.x], collapse = " & "))) |>
    left_join(annot |> dplyr::select(uniprot_id, gene),
              by = c("uniprot" = "uniprot_id")) |>
    dplyr::select(contrasts, uniprot, gene, direction) |>
    arrange(contrasts, direction, gene)
}

panel_c_muscle_proteins <- format_intersection_proteins(upset_m_data, norm_m$annotation, SET_LABELS)
panel_c_plasma_proteins <- format_intersection_proteins(upset_p_data, norm_p$annotation, SET_LABELS)

# --- Panel D: per-gene CV values, per tissue; bootstrap median (95% CI)
#     summary goes in a companion sheet since it's a different grain
#     (per-group, not per-gene) but underlies the same panel
panel_d_muscle      <- cv_df_m |> dplyr::select(gene, group, mean_int, sd_int, n, cv)
panel_d_plasma      <- cv_df_p |> dplyr::select(gene, group, mean_int, sd_int, n, cv)
panel_d_muscle_boot <- boot_m
panel_d_plasma_boot <- boot_p

# --- README sheet
readme_df <- tibble(
  Sheet = c(
    "Panel_A_Muscle", "Panel_A_Plasma",
    "Panel_A_Muscle_Significant", "Panel_A_Plasma_Significant",
    "Panel_B_Muscle", "Panel_B_Plasma",
    "Panel_B_Muscle_Significant", "Panel_B_Plasma_Significant",
    "Panel_C_Muscle", "Panel_C_Plasma",
    "Panel_C_Muscle_Proteins", "Panel_C_Plasma_Proteins",
    "Panel_D_Muscle", "Panel_D_Plasma",
    "Panel_D_Muscle_Bootstrap", "Panel_D_Plasma_Bootstrap"
  ),
  Description = c(
    "Panel A (DAPs per Contrast), Muscle: significant-protein counts by method (p-value/FDR/Pi-value) per contrast, plus the Total Quantified reference row.",
    "Panel A (DAPs per Contrast), Plasma: same as Panel_A_Muscle.",
    "Panel A identity-level, Muscle: the individual proteins behind the counts -- uniprot/gene, contrast, direction (Up/Down), logFC, and the raw values + pass/fail flags for all three significance methods. Included if significant by ANY method.",
    "Panel A identity-level, Plasma: same as Panel_A_Muscle_Significant.",
    sprintf("Panel B (Pathway Enrichment), Muscle: significant-pathway counts per contrast x database. GSEA_COUNT_MODE = \"%s\" (%s).",
            GSEA_COUNT_MODE,
            if (GSEA_COUNT_MODE == "representative")
              "redundant near-duplicate pathways excluded, per the pooled greedy filter"
            else "every significant pathway counted, redundant or not"),
    "Panel B (Pathway Enrichment), Plasma: same as Panel_B_Muscle.",
    "Panel B identity-level, Muscle: the individual significant pathways behind the counts -- ID/name, database, contrast, direction (NES sign), and available stats. Same GSEA_COUNT_MODE filtering as the bars. Column names are assumed clusterProfiler/fgsea-style -- verify against gsea_supp_full.",
    "Panel B identity-level, Plasma: same as Panel_B_Muscle_Significant.",
    "Panel C (Contrast Overlap UpSet), Muscle: protein-level intersection sizes across the three contrasts (Cancer/Training/Residual), by contrast combination, split Up/Down/Mixed.",
    "Panel C (Contrast Overlap UpSet), Plasma: same as Panel_C_Muscle.",
    "Panel C identity-level, Muscle: the individual proteins in each contrast overlap -- uniprot, gene, and direction (Up/Down/Mixed across the contrasts defining that overlap).",
    "Panel C identity-level, Plasma: same as Panel_C_Muscle_Proteins.",
    "Panel D (Per-Group CV), Muscle: per-gene coefficient of variation (%) within each group (CTL/CRC/CRC+RET), from normalized data.",
    "Panel D (Per-Group CV), Plasma: same as Panel_D_Muscle.",
    "Panel D bootstrap summary, Muscle: median CV and 95% bootstrap CI per group (n = 2,000 resamples) underlying the diamond/error-bar annotation.",
    "Panel D bootstrap summary, Plasma: same as Panel_D_Muscle_Bootstrap."
  )
)

# --- Assemble workbook
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "README")
openxlsx::writeData(wb, "README", readme_df)
openxlsx::setColWidths(wb, "README", cols = 1:2, widths = c(26, 95))

sheet_data <- list(
  Panel_A_Muscle              = panel_a_muscle,
  Panel_A_Plasma              = panel_a_plasma,
  Panel_A_Muscle_Significant  = panel_a_muscle_sig,
  Panel_A_Plasma_Significant  = panel_a_plasma_sig,
  Panel_B_Muscle              = panel_b_muscle,
  Panel_B_Plasma              = panel_b_plasma,
  Panel_B_Muscle_Significant  = panel_b_muscle_sig,
  Panel_B_Plasma_Significant  = panel_b_plasma_sig,
  Panel_C_Muscle              = panel_c_muscle,
  Panel_C_Plasma              = panel_c_plasma,
  Panel_C_Muscle_Proteins     = panel_c_muscle_proteins,
  Panel_C_Plasma_Proteins     = panel_c_plasma_proteins,
  Panel_D_Muscle              = panel_d_muscle,
  Panel_D_Plasma              = panel_d_plasma,
  Panel_D_Muscle_Bootstrap    = panel_d_muscle_boot,
  Panel_D_Plasma_Bootstrap    = panel_d_plasma_boot
)

for (nm in names(sheet_data)) {
  openxlsx::addWorksheet(wb, nm)
  openxlsx::writeData(wb, nm, sheet_data[[nm]])
}

openxlsx::saveWorkbook(
  wb,
  file.path(config$data_dir, "fig_02_data.xlsx"),
  overwrite = TRUE
)