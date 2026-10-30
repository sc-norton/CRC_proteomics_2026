# Fig 02 - Supplement (CRE vs PLA) contrasts
#
# Mirrors fig_02.R structure exactly, but built around the 5 supplement-arm
# contrasts (Baseline_supp, Post_training_supp, Training_CRE, Training_PLA,
# Interaction_supp / gsea_group "CRE") instead of the 3 cancer-status
# contrasts (Cancer, Training, Residual / gsea_group "SURV").
#
# ONE SUBSTANTIVE DEPARTURE FROM A PURE 1:1 MIRROR: Section 6's CV panel
# now groups samples by CRE_PRE/CRE_POST/PLA_PRE/PLA_POST instead of
# CTL/CRC/CRC+RET, since that's the grouping the supp contrasts actually
# test (confirmed explicitly, not assumed). Everything else is a direct
# swap of the SURV-family contrast set for the CRE-family one.
#
# NOTE: GSEA data comes from 02_gsea_pathway_consolidation.R's saved RDS
# (same as fig_02.R). Section 5's "Significant Pathways (n)" panel excludes
# redundant pathways per GSEA_COUNT_MODE, consistent with other GSEA
# figures.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, patchwork, cowplot, scales, enrichplot, grid, ggnewscale)

# --- Configuration -----------------------------------------------------------
config <- list(
  limma_muscle = "03_analysis/c_data/muscle/limma",
  limma_plasma = "03_analysis/c_data/plasma/limma",
  gsea_data    = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data_consolidated.rds",
  norm_muscle  = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS",
  norm_plasma  = "01_normalization/c_data/plasma/02_normalized_data_plasma.RDS",
  fig_dir      = "04_figures/fig_02_supp",
  data_dir     = "04_figures/fig_02_supp/data"
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

gsea_group <- "CRE"  # matches the group column added by consolidation script

gsea_comps <- c("Baseline_supp", "Post_training_supp", "Training_CRE",
                "Training_PLA", "Interaction_supp")

# Display labels reused verbatim from 02_gsea_pathway_consolidation.R's
# comp_labels for the CRE-family comparisons, so wording stays consistent
# across the consolidation script, this figure, and the bubble/other GSEA
# figures.
comp_labels <- c(
  Baseline_supp       = "Baseline CRE vs PLA",
  Post_training_supp  = "Training CRE vs PLA",
  Training_CRE        = "POST vs PRE (CRE)",
  Training_PLA         = "POST vs PRE (PLA)",
  Interaction_supp     = "Supplement x Time interaction",
  `Total Quantified`   = "Total\nQuantified"  # new level added for panel A's reference bar (section 2b)
)

SET_LABELS <- c(
  Baseline_supp       = "Baseline CRE vs PLA",
  Post_training_supp  = "Training CRE vs PLA",
  Training_CRE        = "POST vs PRE (CRE)",
  Training_PLA         = "POST vs PRE (PLA)",
  Interaction_supp     = "Supplement x Time interaction"
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
    Baseline_supp = read_excel(file.path(config$limma_muscle,
                                         "limma_results_all_muscle_cv.xlsx"), sheet = "Baseline_supp"),
    Post_training_supp = read_excel(file.path(config$limma_muscle,
                                              "limma_results_all_muscle_cv.xlsx"), sheet = "Post_training_supp"),
    Training_CRE = read_excel(file.path(config$limma_muscle,
                                        "limma_results_all_muscle_cv.xlsx"), sheet = "Training_CRE"),
    Training_PLA = read_excel(file.path(config$limma_muscle,
                                        "limma_results_all_muscle_cv.xlsx"), sheet = "Training_PLA"),
    Interaction_supp = read_excel(file.path(config$limma_muscle,
                                            "limma_results_all_muscle_cv.xlsx"), sheet = "Interaction_supp")
  ),
  plasma = list(
    Baseline_supp = read_excel(file.path(config$limma_plasma,
                                         "limma_results_all_plasma.xlsx"), sheet = "Baseline_supp"),
    Post_training_supp = read_excel(file.path(config$limma_plasma,
                                              "limma_results_all_plasma.xlsx"), sheet = "Post_training_supp"),
    Training_CRE = read_excel(file.path(config$limma_plasma,
                                        "limma_results_all_plasma.xlsx"), sheet = "Training_CRE"),
    Training_PLA = read_excel(file.path(config$limma_plasma,
                                        "limma_results_all_plasma.xlsx"), sheet = "Training_PLA"),
    Interaction_supp = read_excel(file.path(config$limma_plasma,
                                            "limma_results_all_plasma.xlsx"), sheet = "Interaction_supp")
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
# Universe size for each tissue -- row count is identical across all 5
# supp-family contrasts within a tissue (same limma model, only the
# contrast differs), so any one sheet's row count is that tissue's total
# quantified proteins. Using "Baseline_supp" since it's guaranteed present
# and consistently spelled across both tissues' limma_results lists.
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
  Muscle = nrow(limma_results$muscle$Baseline_supp),
  Plasma = nrow(limma_results$plasma$Baseline_supp)
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
    panel.spacing        = unit(1, "lines"),
    axis.text.x          = element_text(size = 7)
  )

# --- 4: UpSet plot -----------------------------------------------------------
make_upset_plot <- function(tissue, title = NULL) {
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
  
  # 5 contrasts -> 2^5 - 1 = 31 possible non-empty combinations (vs 7 for
  # the 3-contrast SURV-family figure). Combinations with zero members are
  # dropped below, so the actually-displayed column count depends on your
  # data, but expect a visibly busier UpSet than fig_02.R's.
  comb_grid <- expand.grid(
    Baseline_supp       = c(TRUE, FALSE),
    Post_training_supp  = c(TRUE, FALSE),
    Training_CRE        = c(TRUE, FALSE),
    Training_PLA         = c(TRUE, FALSE),
    Interaction_supp     = c(TRUE, FALSE)
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
      in_sets = list(in_sets), degree = length(in_sets),
      up      = sum(dirs == "Up"),
      down    = sum(dirs == "Down"),
      mixed   = sum(dirs == "Mixed"),
      total   = length(members)
    )
  })) |>
    filter(total > 0) |>
    arrange(desc(total)) |>
    mutate(x = row_number())
  
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
  
  get_upset_components <- function(tissue) {
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
      Baseline_supp       = c(TRUE, FALSE),
      Post_training_supp  = c(TRUE, FALSE),
      Training_CRE        = c(TRUE, FALSE),
      Training_PLA         = c(TRUE, FALSE),
      Interaction_supp     = c(TRUE, FALSE)
    ) |> filter(rowSums(across(everything())) > 0)
    
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
        in_sets = list(in_sets),
        up      = sum(dirs == "Up"),
        down    = sum(dirs == "Down"),
        mixed   = sum(dirs == "Mixed"),
        total   = length(members)
      )
    })) |>
      filter(total > 0) |>
      arrange(desc(total)) |>
      mutate(x = row_number())
    
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
         CONTRASTS = CONTRASTS)
  }
  
  m      <- get_upset_components("muscle")
  p      <- get_upset_components("plasma")
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
# Sourced from gsea_supp_full (consolidated), filtered to group == "CRE"
# in Section 1. Counting logic is controlled by GSEA_COUNT_MODE:
#   "representative" -- a pathway counts only if it's both significant in
#     THIS comparison AND flagged representative for its (tissue, database)
#     by the pooled redundancy filter (representative_lookup). The
#     redundancy flag itself is pooled across all five gsea_comps (see
#     02_gsea_pathway_consolidation.R), but whether it counts toward THIS
#     bar is still gated on significance in this specific comparison -- so
#     a pathway significant in two comparisons still contributes to both
#     bars, it just won't be double-counted against a near-duplicate
#     pathway within the same bar.
#   "raw" -- original behavior: every significant pathway counts,
#     redundant or not.
make_gsea_counts_combined <- function() {
  base <- gsea_supp_full %>%
    filter(comparison %in% gsea_comps, significant)
  
  if (GSEA_COUNT_MODE == "representative") {
    base <- base %>%
      semi_join(representative_lookup, by = c("tissue", "database", "ID"))
  } else if (GSEA_COUNT_MODE != "raw") {
    stop("GSEA_COUNT_MODE must be 'representative' or 'raw'")
  }
  
  counts <- base %>%
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
                                           hjust = 0, color = colors$gray_dark),
      axis.text.x           = element_text(size = 7)
    )
}

p_gsea_combined <- make_gsea_counts_combined()

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

# Supp-arm groups (CRE vs PLA, pre/post) -- this is what the 5 supp
# contrasts actually test, unlike the CTL/CRC/CRC+RET grouping fig_02.R
# uses for the SURV-family contrasts.
CRE_PLA_GROUPS <- c("CRE_PRE", "CRE_POST", "PLA_PRE", "PLA_POST")

cv_df_m <- build_cv_data(
  data  = norm_m$data,
  annot = norm_m$annotation,
  meta  = norm_m$metadata
) %>% filter(group %in% CRE_PLA_GROUPS)

cv_df_p <- build_cv_data(
  data  = norm_p$data,
  annot = norm_p$annotation,
  meta  = norm_p$metadata
) %>% filter(group %in% CRE_PLA_GROUPS)

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

boot_m <- bootstrap_median_cv(cv_df_m, CRE_PLA_GROUPS)
boot_p <- bootstrap_median_cv(cv_df_p, CRE_PLA_GROUPS)

cat("\nMuscle — bootstrap median CV (95% CI):\n")
print(boot_m)
cat("\nPlasma — bootstrap median CV (95% CI):\n")
print(boot_p)

# --- 6c: Brown-Forsythe global test (for methods reporting) ------------------
if (!requireNamespace("car", quietly = TRUE)) install.packages("car")
library(car)

bf_m <- car::leveneTest(
  cv ~ factor(group, levels = CRE_PLA_GROUPS),
  data   = cv_df_m |> filter(group %in% CRE_PLA_GROUPS),
  center = median
)

bf_p <- car::leveneTest(
  cv ~ factor(group, levels = CRE_PLA_GROUPS),
  data   = cv_df_p |> filter(group %in% CRE_PLA_GROUPS),
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

# --- 6d: Violin fill palette — arm-by-timepoint shading -----------------
# Hue encodes timepoint (purple = PRE, orange = POST), light/dark shade
# encodes arm (light = CRE, dark = PLA) -- reusing the same convention
# already established in the normalization reports scripts' GROUP_COLOURS,
# so this figure's palette reads consistently with earlier QC figures
# rather than introducing a new scheme just for this panel. Same fill is
# used for both tissues (tissue is distinguished by the BG_COLORS panel
# background instead), matching fig_02.R's approach.
violin_fills <- c(
  Muscle_CRE_PRE   = colors$purple_light,
  Muscle_CRE_POST  = colors$orange_light,
  Muscle_PLA_PRE   = colors$purple_dark,
  Muscle_PLA_POST  = colors$orange_dark,
  Plasma_CRE_PRE   = colors$purple_light,
  Plasma_CRE_POST  = colors$orange_light,
  Plasma_PLA_PRE   = colors$purple_dark,
  Plasma_PLA_POST  = colors$orange_dark
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
      group        = factor(group, levels = CRE_PLA_GROUPS),
      group_tissue = paste0(tissue, "_", group)
    )
  
  boot_combined <- bind_rows(
    boot_m |> mutate(tissue = "Muscle"),
    boot_p |> mutate(tissue = "Plasma")
  ) |>
    mutate(
      tissue = factor(tissue, levels = c("Muscle", "Plasma")),
      group  = factor(group,  levels = CRE_PLA_GROUPS)
    )
  
  # label_y positioned at a uniform height above each tissue panel's max
  # CV value (y_max + range * 0.08) -- matches make_cv_violin() in the
  # normalization reports scripts.
  violin_range <- cv_combined |>
    filter(group %in% CRE_PLA_GROUPS) |>
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
    scale_x_discrete(labels = c(CRE_PRE = "CRE\nPRE", CRE_POST = "CRE\nPOST",
                                PLA_PRE = "PLA\nPRE", PLA_POST = "PLA\nPOST")) +
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
      axis.text.x      = element_text(face = "bold", size = 8),
      strip.text       = element_blank(),
      strip.background = element_blank(),
      panel.spacing    = unit(1, "lines"),
      plot.subtitle     = element_text(size = 8, face = "plain",
                                       hjust = 0, color = colors$gray_dark)
    )
}

p_cv_plot_combined <- make_cv_violin_combined(cv_df_m, cv_df_p, boot_m, boot_p)

# --- 7: Assemble & Save ------------------------------------------------------
# NOTE ON SIZING: widened relative to fig_02.R throughout -- panel A/B now
# have 6 x-axis categories instead of 4, the UpSet has 5 sets (and longer
# labels) instead of 3, and the CV panel has 4 groups instead of 3. These
# are starting-point estimates; re-check once rendered and adjust further
# if anything looks cramped.
f2_supp <- plot_grid(
  p_da_combined, p_gsea_combined, p_upset_c, p_cv_plot_combined,
  labels = "AUTO",
  nrow   = 2,
  ncol   = 2
)

ggsave(file.path(config$fig_dir, "fig_02_supp.pdf"),
       f2_supp, width = 18, height = 13, units = "in", device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp.png"),
       f2_supp, width = 18, height = 13, units = "in", dpi = 300)

ggsave(file.path(config$fig_dir, "fig_02_supp_A.pdf"),          p_da_combined,      width = 8,  height = 8,  device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp_B.pdf"),          p_gsea_combined,    width = 8,  height = 8,  device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp_C.pdf"),          p_upset_c,          width = 14, height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp_C_muscle.pdf"),   p_upset_m,          width = 10, height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp_C_plasma.pdf"),   p_upset_p,          width = 10, height = 10, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_02_supp_D.pdf"),          p_cv_plot_combined, width = 5,  height = 6,  device = cairo_pdf)

ggsave(file.path(config$fig_dir, "fig_02_supp_A.png"),          p_da_combined,      width = 8,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_supp_B.png"),          p_gsea_combined,    width = 8,  height = 8,  dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_supp_C.png"),          p_upset_c,          width = 14, height = 10, dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_supp_C_muscle.png"),   p_upset_m,          width = 10, height = 10, dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_supp_C_plasma.png"),   p_upset_p,          width = 10, height = 10, dpi = 300)
ggsave(file.path(config$fig_dir, "fig_02_supp_D.png"),          p_cv_plot_combined, width = 5,  height = 6,  dpi = 300)