# =============================================================================
# fig_05_panel_A_protein.R
# Panel A: Protein reversal scatter + REVERSAL/EXACERBATION ORA flanking bars.
# Depends on fig_05_style.R being sourced first.
# =============================================================================

pacman::p_load(dplyr, tidyr, tibble, ggplot2, ggrepel, patchwork)

# Multiplier applied on top of BASE_PATHWAY, just for the ORA flanking-bar
# labels in Panel A -- bump this if you want them bigger/smaller without
# touching BASE_PATHWAY globally.
PATHWAY_LABEL_SCALE <- 1.45

# Fraction of headroom added to the scatter's axis limits beyond the raw
# data-derived max(|logFC|). Without this, coord_cartesian's xlim/ylim sit
# EXACTLY at the most extreme point, so any protein at that max value (or
# its repelled label) lands flush against the panel edge and gets clipped.
# This only pads the coord_cartesian() view window -- `lim` itself is left
# untouched, since axis tick positions, the axis-title annotations, and the
# quadrant labels below are all computed relative to the raw `lim`. Raise
# this (e.g. 1.15) if points are still getting cut off; lower it (e.g. 1.03)
# if there's now too much empty space around the edges.
AXIS_PAD <- 1.08

ora_label_map <- c(
  "L-amino acid catabolic process"                                  = "L-AA Catabolism",
  "branched-chain amino acid catabolic process"                     = "BCAA Catabolism",
  "quinone metabolic process"                                       = "Quinone Metabolism",
  "pyruvate family amino acid catabolic process"                    = "Pyruvate AA Catabolism",
  "spermatid differentiation"                                       = "Spermatid Diff.",
  "calcium ion transmembrane transport"                             = "Calcium Transport",
  "actin cytoskeleton organization"                                 = "Actin Cytoskeleton Org.",
  "cell adhesion"                                                   = "Cell Adhesion",
  "vesicle-mediated transport"                                      = "Vesicle Transport",
  "cell motility"                                                   = "Cell Motility",
  "cell-substrate adhesion"                                         = "Cell-Substrate\nAdhesion",
  "regulation of actin cytoskeleton organization"                   = "Regulation of\nActin Cytoskeleton",
  "response to toxic substance"                                     = "Toxic Substance Response",
  "glutathione metabolic process"                                   = "Glutathione Metabolism",
  "Ras protein signal transduction"                                 = "Ras Signaling",
  "aldehyde metabolic process"                                      = "Aldehyde Metabolism",
  "regulation of proteolysis involved in protein catabolic process" = "Regulation of\nProteolysis",
  "proteolysis"                                                     = "Proteolysis",
  "proton transmembrane transport"                                  = "Proton Transport",
  "energy derivation by oxidation of organic compounds"             = "Oxidative Energy Deriv.",
  "mitochondrial electron transport, NADH to ubiquinone"            = "Mitochondrial\nElectron Trans.",
  "nucleoside phosphate biosynthetic process"                       = "Nucleoside-P\nBiosynthesis",
  "translation"                                                     = "Translation",
  "NADH dehydrogenase complex assembly"                             = "NADH\nDehydrogenase Assembly"
)

apply_labels <- function(df, map) {
  df |> mutate(label = if_else(Description %in% names(map), map[Description], Description))
}

# Kept independent of the pooled-Jaccard filter used for GSEA pathways --
# REVERSAL/EXACERBATION ORA is protein-level ORA, not one of the pooled
# GSEA contrasts.
jaccard_filter_ora <- function(res, threshold = 0.5,
                               padj_col = "p.adjust", id_col = "ID", gene_col = "geneID") {
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  df <- as.data.frame(res)
  df <- df[grepl("^GO:", df[[id_col]]), ]
  df <- df[order(df[[padj_col]]), ]
  sets <- strsplit(df[[gene_col]], "/")
  n    <- nrow(df); keep <- rep(TRUE, n)
  for (i in seq_len(n - 1)) {
    if (!keep[i]) next
    for (j in seq(i + 1, n)) {
      if (!keep[j]) next
      if (jaccard(sets[[i]], sets[[j]]) >= threshold) keep[j] <- FALSE
    }
  }
  df[keep, ]
}

make_ora_right <- function(df, fill_low, fill_high, panel_width_mm, show_x_title = FALSE,
                           outside_labels = character(0)) {
  df <- df |>
    mutate(neg_log10_p = -log10(p.adjust),
           y_order = reorder(label, neg_log10_p),
           sig_label = case_when(p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**",
                                 p.adjust < 0.05 ~ "*", TRUE ~ "ns"),
           is_outside = label %in% outside_labels)
  
  df_in  <- filter(df, !is_outside)
  df_out <- filter(df, is_outside)
  
  ggplot(df, aes(x = neg_log10_p, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    # in-bar labels: white text, anchored near the base of the bar
    geom_text(data = df_in, aes(label = label), x = 0.05, hjust = 0,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "white", fontface = "bold") +
    # flagged labels: dark text, anchored just past the bar's tip
    geom_text(data = df_out, aes(x = neg_log10_p, label = label), hjust = -0.05,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "grey15", fontface = "bold") +
    geom_text(aes(x = neg_log10_p, label = sig_label), hjust = -0.2, vjust = 0.75,
              size = scale_text(BASE_COUNT, panel_width_mm), colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_low, high = fill_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, if (nrow(df_out) > 0) 0.32 else 0.18)), position = "bottom") +
    scale_y_discrete(position = "left") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL) +
    fig_theme_scaled(panel_width_mm) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), legend.position = "none",
      plot.margin = margin(2, 10, 2, 0)
    )
}

make_ora_left <- function(df, fill_low, fill_high, panel_width_mm, show_x_title = FALSE,
                          outside_labels = character(0)) {
  df <- df |>
    mutate(neg_log10_p = -log10(p.adjust), neg_log10_p_neg = -neg_log10_p,
           y_order = reorder(label, neg_log10_p),
           sig_label = case_when(p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**",
                                 p.adjust < 0.05 ~ "*", TRUE ~ "ns"),
           is_outside = label %in% outside_labels)
  
  df_in  <- filter(df, !is_outside)
  df_out <- filter(df, is_outside)
  
  ggplot(df, aes(x = neg_log10_p_neg, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    # in-bar labels: white text, anchored near the base of the bar
    geom_text(data = df_in, aes(x = neg_log10_p_neg * 0.97, label = label), hjust = 0,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "white", fontface = "bold") +
    # flagged labels: dark text, anchored just past the bar's tip (further left)
    geom_text(data = df_out, aes(x = neg_log10_p_neg, label = label), hjust = 1,
              size = scale_text(BASE_PATHWAY * PATHWAY_LABEL_SCALE, panel_width_mm),
              colour = "grey15", fontface = "bold") +
    geom_text(aes(x = neg_log10_p_neg, label = sig_label), hjust = 1.4, vjust = 0.75,
              size = scale_text(BASE_COUNT, panel_width_mm), colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_low, high = fill_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(if (nrow(df_out) > 0) 0.32 else 0.18, 0)),
                       labels = function(x) abs(x), position = "bottom") +
    scale_y_discrete(position = "right") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL) +
    fig_theme_scaled(panel_width_mm) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      panel.grid = element_blank(), legend.position = "none",
      plot.margin = margin(2, 0, 2, 10)
    )
}

#' Build Panel A: protein reversal scatter + ORA flanking bars
#'
#' @param limma_baseline Raw Baseline_SURVvCTL sheet (unselected columns OK)
#' @param limma_training_surv Raw Training_SURV sheet
#' @param gsea_muscle Named list from load_gsea_rds(), must contain $REVERSAL / $EXACERBATION
#' @param panel_width_mm Total rendered width of Panel A in the final composite (mm)
#' @param pal Named color vector for the 4 reversal categories (default pal_shared)
build_panel_A <- function(limma_baseline, limma_training_surv, gsea_muscle,
                          panel_width_mm = PANEL_MD, pal = pal_shared,
                          ora_outside = list()) {
  
  df_PRE <- limma_baseline |>
    dplyr::select(uniprot, gene, logFC_CRC = logFC, Pi.Val_CRC = Pi.Val, adj.P.Val_CRC = adj.P.Val)
  df_RET <- limma_training_surv |>
    dplyr::select(uniprot, gene, logFC_RET = logFC, Pi.Val_RET = Pi.Val, adj.P.Val_RET = adj.P.Val)
  
  df_rev <- inner_join(df_PRE, df_RET, by = c("uniprot", "gene")) |>
    mutate(sig = case_when(
      Pi.Val_CRC < 0.05 & Pi.Val_RET < 0.05 ~ "Both",
      Pi.Val_CRC < 0.05                     ~ "CRC v CTL only",
      Pi.Val_RET < 0.05                     ~ "CRC+RET v CRC only",
      TRUE                                  ~ "Neither"
    ))
  
  lim <- max(abs(c(df_rev$logFC_CRC, df_rev$logFC_RET)), na.rm = TRUE)
  r_spear <- suppressWarnings(cor.test(df_rev$logFC_CRC, df_rev$logFC_RET, method = "spearman", exact = FALSE))
  
  protein_label_df <- df_rev |>
    filter(sig != "Neither", abs(logFC_CRC) > 1 | abs(logFC_RET) > 1) |>
    mutate(
      Pi.Val_RET = pmax(Pi.Val_RET, 1e-300), Pi.Val_CRC = pmax(Pi.Val_CRC, 1e-300),
      direction  = if_else(abs(logFC_CRC) >= abs(logFC_RET),
                           if_else(logFC_CRC > 0, "Up", "Down"),
                           if_else(logFC_RET > 0, "Up", "Down")),
      score = pmax(-log10(Pi.Val_RET), -log10(Pi.Val_CRC))
    ) |>
    group_by(sig, direction) |> slice_max(score, n = 5, with_ties = FALSE) |> ungroup()
  
  quad_counts <- df_rev |>
    mutate(q = case_when(
      logFC_CRC < 0 & logFC_RET > 0 ~ "TL", logFC_CRC > 0 & logFC_RET > 0 ~ "TR",
      logFC_CRC > 0 & logFC_RET < 0 ~ "BR", logFC_CRC < 0 & logFC_RET < 0 ~ "BL"
    )) |>
    filter(!is.na(q)) |> group_by(q) |>
    summarise(n_total = n(), n_sig = sum(sig != "Neither", na.rm = TRUE), .groups = "drop")
  
  quad_ann <- build_quad_label_df(
    counts = quad_counts, quad_ids = c("TL", "TR", "BR", "BL"),
    titles = c("Reversed (Up RET)", "Exacerbated (Up)",
               "Reversed (Down RET)", "Exacerbated (Down)"),
    label_fn = function(title, row) sprintf("%s\n%d / %d", title, row$n_sig, row$n_total),
    blue_quads = c("TL", "BR"), red_quads = c("TR", "BL")
  )
  
  th <- fig_theme_scaled(panel_width_mm)
  
  p_scatter <- ggplot(df_rev, aes(x = logFC_CRC, y = logFC_RET, color = sig, fill = sig)) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = colors$blue, alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = colors$blue, alpha = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
    geom_point(data = filter(df_rev, sig == "Neither"), size = 1.2, alpha = 0.5, shape = 16) +
    geom_point(data = filter(df_rev, sig != "Neither"), size = 2.2, alpha = 0.85, shape = 16) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    geom_quadrant_label(quad_ann, size = scale_text(BASE_QUADRANT, panel_width_mm)) +
    geom_label_repel(
      data = protein_label_df, aes(label = gene, fill = sig),
      color = colors$white, size = scale_text(BASE_GENE, panel_width_mm),
      max.overlaps = 20, box.padding = 0.3, label.padding = unit(0.15, "lines"),
      label.size = NA, segment.color = "grey60", segment.size = 0.3, show.legend = FALSE,
      xlim = c(-lim * 0.85, lim * 0.85), ylim = c(-lim * 0.85, lim * 0.85)
    ) +
    annotate("text", x = -lim * 0.80, y = 0.01, label = expression(log[2]*FC~"(CRC vs CTL)"),
             hjust = 0.5, vjust = -0.5, size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey20") +
    annotate("text", x = -0.01, y = -lim * 0.80, label = expression(log[2]*FC~"(CRC+RET vs CRC)"),
             hjust = 0.5, vjust = 1.3, angle = 90, size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey20") +
    coord_cartesian(xlim = c(-lim * AXIS_PAD, lim * AXIS_PAD),
                    ylim = c(-lim * AXIS_PAD, lim * AXIS_PAD), clip = "off") +
    th +
    theme(legend.position = "bottom", axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank())
  
  axis_breaks <- pretty(c(-lim, lim))
  axis_breaks <- axis_breaks[axis_breaks != 0 & abs(axis_breaks) < lim]
  p_scatter <- p_scatter +
    annotate("segment", x = axis_breaks, xend = axis_breaks, y = -lim * 0.015, yend = lim * 0.015,
             color = "grey30", linewidth = 0.4) +
    annotate("text", x = axis_breaks, y = -lim * 0.03, label = axis_breaks,
             vjust = 1, size = scale_text(BASE_STAT * 0.85, panel_width_mm), color = "grey20") +
    annotate("segment", y = axis_breaks, yend = axis_breaks, x = -lim * 0.015, xend = lim * 0.015,
             color = "grey30", linewidth = 0.4) +
    annotate("text", y = axis_breaks, x = -lim * 0.03, label = axis_breaks,
             hjust = 1, size = scale_text(BASE_STAT * 0.85, panel_width_mm), color = "grey20")
  
  # --- ORA flanking bars -------------------------------------------------
  ora_rev_up   <- jaccard_filter_ora(gsea_muscle$REVERSAL$ora_bp_up@result)   |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_rev_down <- jaccard_filter_ora(gsea_muscle$REVERSAL$ora_bp_down@result) |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_ex_up    <- jaccard_filter_ora(gsea_muscle$EXACERBATION$ora_bp_up@result)   |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  ora_ex_down  <- jaccard_filter_ora(gsea_muscle$EXACERBATION$ora_bp_down@result) |> arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
  
  n_enrich <- sum(gsea_muscle$REVERSAL$ora_bp_up@result$p.adjust < 0.05, na.rm = TRUE) +
    sum(gsea_muscle$REVERSAL$ora_bp_down@result$p.adjust < 0.05, na.rm = TRUE) +
    sum(gsea_muscle$EXACERBATION$ora_bp_up@result$p.adjust < 0.05, na.rm = TRUE) +
    sum(gsea_muscle$EXACERBATION$ora_bp_down@result$p.adjust < 0.05, na.rm = TRUE)
  
  # ORA flank panels sit at ~1.1/4.2 = ~0.26x the total composite's rendered
  # width (see plot_layout(widths) below: 1.1 : 2 : 1.1), so their own text
  # scale target is narrower than the scatter's.
  flank_width_mm <- panel_width_mm * (1.1 / 4.2)
  
  get_outside <- function(nm) if (is.null(ora_outside[[nm]])) character(0) else ora_outside[[nm]]
  
  left_col  <- make_ora_left(ora_rev_up,  colors$blue_dark, colors$blue_dark, flank_width_mm, FALSE,
                             outside_labels = get_outside("rev_up")) /
    make_ora_left(ora_ex_down, colors$red_dark,  colors$red_dark,  flank_width_mm, TRUE,
                  outside_labels = get_outside("ex_down"))
  right_col <- make_ora_right(ora_ex_up,    colors$red_dark,  colors$red_dark,  flank_width_mm, FALSE,
                              outside_labels = get_outside("ex_up")) /
    make_ora_right(ora_rev_down, colors$blue_dark, colors$blue_dark, flank_width_mm, TRUE,
                   outside_labels = get_outside("rev_down"))
  
  fig_compound <- (left_col | p_scatter | right_col) +
    plot_layout(widths = c(1.1, 2, 1.1)) &
    theme(plot.margin = margin(2, 0, 2, 0))
  
  ttl <- "Protein Reversal"
  sub <- sprintf("N = %d | %d sig (\u03a0) | %d ORA terms enriched (FDR) | \u03c1 = %.2f",
                 nrow(df_rev), sum(df_rev$sig != "Neither"), n_enrich, r_spear$estimate)
  
  list(
    plot           = fig_compound,                       # for standalone save (own legend)
    plot_composite = strip_for_composite(fig_compound),   # no legend/title, for the 380mm grid
    title = ttl, subtitle = sub,
    df_rev = df_rev
  )
}