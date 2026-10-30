# --- 0: Palette for comp1_grouped ------------------------
# 16-level qualitative palette for subcellular localization groups.
# "Other"/"Unknown" pinned to neutral grays so they recede visually.
colors$comp1_grouped <- c(
  "Cytoplasm/Cytosol"          = "#4E79A7",
  "Cytoskeleton/Contractile"   = "#F28E2B",
  "Membrane/Cell Surface"      = "#E15759",
  "Mitochondria"                = "#59A14F",
  "Protein Complex"             = "#B07AA1",
  "Extracellular/Secreted"      = "#EDC948",
  "Nucleus/Chromatin"            = "#76B7B2",
  "ER/Golgi/Vesicle"             = "#FF9DA7",
  "Cell Junction/Adhesion"       = "#9C755F",
  "Cilia/Centrosome/Neuron"      = "#BAB0AC",
  "Granule/Lysosome/Vesicle"     = "#D37295",
  "Blood Microparticle"          = "#8CD17D",
  "Sarcoplasmic/Myofibrillar"    = "#B6992D",
  "Acrosome"                     = "#499894",
  "Other"                        = "grey65",
  "Unknown"                      = "grey85"
)

# --- 3: Volcano plot function (comp1_grouped color) ---------------
# --- 3: Volcano plot function (comp1_grouped color, sig only) ---------------
make_volcano_plot <- function(df_vol, selected_comparisons,
                              n_labels = 5, pval_thresh = 0.05, lfc_thresh = 0.6) {
  df <- df_vol %>%
    filter(comparison %in% selected_comparisons) %>%
    mutate(
      comparison    = factor(comparison, levels = selected_comparisons),
      comp1_grouped = factor(comp1_grouped, levels = names(colors$comp1_grouped))
    )
  
  df_ns  <- df %>% filter(DE.pi == "NDE")
  df_sig <- df %>% filter(DE.pi != "NDE")
  
  df_label <- df_sig %>%
    mutate(nudge_x = ifelse(DE.pi == "UP", 2.5, -2.5), nudge_y = 0.5) %>%
    group_by(comparison, tissue, DE.pi) %>%
    slice_max(abs(logFC), n = n_labels) %>%
    ungroup()
  
  # Count UP and DOWN per facet -- still direction-based, independent of point color
  df_counts <- df_sig %>%
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
  
  ggplot(df, aes(logFC, -log10(P.Value))) +
    geom_rect(
      data        = tissue_bg,
      aes(fill    = fill_color),
      xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
      alpha       = 0.15,
      inherit.aes = FALSE
    ) +
    scale_fill_identity() +
    geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = colors$gray) +
    # Non-significant points: flat neutral background layer, drawn first
    geom_point(data = df_ns, color = "grey80", size = 1, alpha = 0.5) +
    # Significant points: colored by subcellular localization, drawn on top
    geom_point(data = df_sig, aes(color = comp1_grouped), size = 1.2, alpha = 0.85) +
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
    scale_color_manual(values   = colors$comp1_grouped,
                       na.value = "grey85", drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2)) +
    coord_cartesian(xlim = c(-4, 4)) +
    facet_grid2(tissue ~ comparison, strip = strip,
                labeller = labeller(comparison = comp_labels)) +
    labs(x = bquote(bold("log2(Fold Change)")), y = bquote(bold("-log10(p-value)")),
         color = "Subcellular\nLocalization\n(significant only)") +
    theme_cr(base_size = 10)
}

volcano_all <- make_volcano_plot(df_vol, c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV"))
