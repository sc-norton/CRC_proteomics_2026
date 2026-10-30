# =============================================================================
# fig_05_panel_B_pathway.R
# Panel B: Pathway (GSEA NES) reversal scatter.
# Depends on fig_05_style.R being sourced first.
# =============================================================================

pacman::p_load(dplyr, ggplot2, ggrepel, patchwork)

PANEL_B_TEXT_SCALE <- 1.3

pathway_label_map <- c(
  "cell adhesion"                                                                                                        = "Cell Adhesion",
  "cell motility"                                                                                                        = "Cell Motility",
  "cytoskeleton organization"                                                                                            = "Cytoskeleton Org.",
  "wound healing"                                                                                                        = "Wound Healing",
  "cell differentiation"                                                                                                 = "Cell Differentiation",
  "mitochondrial gene expression"                                                                                        = "Mito. Gene Expr.",
  "amino acid metabolic process"                                                                                         = "AA Metabolism",
  "tRNA metabolic process"                                                                                               = "tRNA Metabolism",
  "vesicle-mediated transport"                                                                                           = "Vesicle Transport",
  "mitochondrion organization"                                                                                           = "Mito. Organization",
  "extracellular matrix organization"                                                                                    = "ECM Organization",
  "generation of precursor metabolites and energy"                                                                       = "Precursor Metab. & Energy",
  "cell junction organization"                                                                                           = "Cell Junction Org.",
  "immune system process"                                                                                                = "Immune System",
  "transmembrane transport"                                                                                              = "Transmembrane Transport",
  "cilium organization"                                                                                                  = "Cilium Org.",
  "inflammatory response"                                                                                                = "Inflammation",
  "regulation of DNA-templated transcription"                                                                            = "Transcription Reg.",
  "detoxification"                                                                                                       = "Detoxification",
  "cytoplasmic translation"                                                                                              = "Cytoplasmic Translation",
  "DNA-templated transcription"                                                                                          = "Transcription",
  "reproductive process"                                                                                                 = "Reproductive Process",
  "sulfur compound metabolic process"                                                                                    = "Sulfur Metabolism",
  "modified amino acid metabolic process"                                                                                = "Modified AA Metabolism",
  "Genes encoding proteins involved in oxidative phosphorylation."                                                       = "OXPHOS",
  "Genes involved in metabolism of heme (a cofactor consisting of iron and porphyrin) and erythroblast differentiation." = "Heme Metabolism",
  "Genes defining epithelial-mesenchymal transition, as in wound healing, fibrosis and metastasis."                      = "EMT",
  "Genes encoding components of the complement system, which is part of the innate immune system."                       = "Complement",
  "Genes encoding components of apical junction complex."                                                                = "Apical Junction",
  "Genes encoding components of blood coagulation system; also up-regulated in platelets."                               = "Coagulation",
  "Genes encoding proteins involved in processing of drugs and other xenobiotics."                                       = "Xenobiotic Metab.",
  "Genes defining late response to estrogen."                                                                            = "Estrogen Response",
  "Genes up-regulated by STAT5 in response to IL2 stimulation."                                                          = "IL2-STAT5 Signaling",
  "Genes mediating programmed cell death (apoptosis) by activation of caspases."                                         = "Apoptosis",
  "Genes up-regulated by reactive oxigen species (ROS)."                                                                 = "ROS Pathway",
  "Genes down-regulated in response to ultraviolet (UV) radiation."                                                      = "UV Response",
  "Genes up-regulated in response to IFNG [GeneID=3458]."                                                                = "IFN\u03b3 Response",
  "Genes up-regulated in response to low oxygen levels (hypoxia)."                                                       = "Hypoxia",
  "Genes important for mitotic spindle assembly."                                                                        = "Mitotic Spindle"
)

#' Build Panel B: pathway (NES) reversal scatter
#'
#' @param gsea_pathway_data Consolidated GSEA RDS ($gsea_supp_full, $representative_lookup)
#' @param panel_width_mm Rendered width of Panel B in the final composite (mm)
#' @param pal Named color vector for the 4 reversal categories (default pal_shared)
build_panel_B <- function(gsea_pathway_data, panel_width_mm = PANEL_MD, pal = pal_shared) {
  
  reversal_group <- "SURV"
  reversal_dbs   <- c("GO Slim", "Hallmark")
  
  nes_pre <- gsea_pathway_data$gsea_supp_full %>%
    filter(group == reversal_group, tissue == "Muscle",
           database %in% reversal_dbs, comparison == "Cancer") %>%
    transmute(ID, Description, database = as.character(database), setSize,
              NES_PRE = NES, padj_PRE = p.adjust)
  
  nes_ret <- gsea_pathway_data$gsea_supp_full %>%
    filter(group == reversal_group, tissue == "Muscle",
           database %in% reversal_dbs, comparison == "Training") %>%
    transmute(ID, Description, database = as.character(database), NES_RET = NES, padj_RET = p.adjust)
  
  rep_lookup <- gsea_pathway_data$representative_lookup %>%
    filter(group == reversal_group, tissue == "Muscle", database %in% reversal_dbs) %>%
    mutate(database = as.character(database))
  
  nes_df <- full_join(nes_pre, nes_ret, by = c("ID", "Description", "database")) %>%
    filter(padj_PRE < 0.05 | padj_RET < 0.05) %>%
    semi_join(rep_lookup, by = c("database", "ID")) %>%
    mutate(
      sig_in = case_when(
        padj_PRE < 0.05 & padj_RET < 0.05 ~ "Both",
        padj_PRE < 0.05                    ~ "CRC v CTL only",
        padj_RET < 0.05                    ~ "CRC+RET v CRC only"
      ),
      go_set = database,
      label  = coalesce(pathway_label_map[Description], Description)
    )
  
  nes_rho <- suppressWarnings(cor.test(nes_df$NES_PRE, nes_df$NES_RET, method = "spearman", exact = FALSE))
  pw_rev_frac <- mean(sign(nes_df$NES_PRE) != sign(nes_df$NES_RET))
  
  nes_label_df <- nes_df |>
    mutate(
      quadrant = case_when(
        NES_PRE >= 0 & NES_RET >= 0 ~ "PP", NES_PRE < 0 & NES_RET < 0 ~ "NN",
        NES_PRE >= 0 & NES_RET < 0  ~ "PN", NES_PRE < 0 & NES_RET >= 0 ~ "NP"
      ),
      abs_nes_sum = abs(NES_PRE) + abs(NES_RET)
    ) |>
    group_by(quadrant) |> slice_max(abs_nes_sum, n = 10, with_ties = FALSE) |> ungroup()
  
  nes_lim <- max(abs(c(nes_df$NES_PRE, nes_df$NES_RET)), na.rm = TRUE)
  
  nes_quad_counts <- nes_df |>
    mutate(q = case_when(
      NES_PRE < 0 & NES_RET >= 0 ~ "TL", NES_PRE >= 0 & NES_RET >= 0 ~ "TR",
      NES_PRE >= 0 & NES_RET < 0 ~ "BR", NES_PRE < 0 & NES_RET < 0  ~ "BL"
    )) |> group_by(q) |>
    summarise(n_total = n(), n_sig = sum(sig_in != "Neither", na.rm = TRUE), .groups = "drop")
  
  nes_quad_labels <- build_quad_label_df(
    counts = nes_quad_counts, quad_ids = c("TL", "TR", "BR", "BL"),
    titles = c("Reversed (Up RET)", "Exacerbated (Up)", "Reversed (Down RET)", "Exacerbated (Down)"),
    label_fn = function(title, row) sprintf("%s\nn = %d (%s%%)", title, row$n_total,
                                            round(row$n_total / sum(nes_quad_counts$n_total) * 100, 1)),
    blue_quads = c("TL", "BR"), red_quads = c("TR", "BL")
  )
  
  th <- fig_theme_scaled(panel_width_mm)
  
  p <- ggplot(nes_df, aes(x = NES_PRE, y = NES_RET, color = sig_in, fill = sig_in)) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0, fill = colors$red, alpha = 0.2) +
    annotate("rect", xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0, fill = colors$blue, alpha = 0.2) +
    annotate("rect", xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf, fill = colors$blue, alpha = 0.2) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
    geom_point(data = filter(nes_df, sig_in == "Neither"), aes(shape = go_set, size = setSize), alpha = 0.8) +
    geom_point(data = filter(nes_df, sig_in != "Neither"), aes(shape = go_set, size = setSize), alpha = 0.85) +
    scale_color_manual(values = pal, name = NULL) +
    scale_fill_manual(values = pal, guide = "none") +
    scale_shape_manual(values = c("GO Slim" = 16, "Hallmark" = 17), name = "Database") +
    scale_size_continuous(name = "Set size", range = c(1, 6)) +
    geom_quadrant_label(nes_quad_labels, size = scale_text(BASE_QUADRANT * PANEL_B_TEXT_SCALE, panel_width_mm)) +
    geom_label_repel(
      data = nes_label_df, aes(label = label, fill = sig_in), color = colors$white,
      size = scale_text(BASE_GENE * PANEL_B_TEXT_SCALE, panel_width_mm), max.overlaps = 20, box.padding = 0.3,
      label.padding = unit(0.15, "lines"), label.size = NA, segment.color = "grey60",
      segment.size = 0.3, show.legend = FALSE,
      xlim = c(-nes_lim * 0.85, nes_lim * 0.85), ylim = c(-nes_lim * 0.85, nes_lim * 0.85)
    ) +
    annotate("text", x = -Inf, y = 0, label = " NES~'(CRC vs CTL)'", hjust = 0, vjust = -0.3,
             size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey30", parse = TRUE) +
    annotate("text", x = 0, y = -Inf, label = " NES~'(CRC+RET vs CRC)'", hjust = 0, vjust = -0.3, angle = 90,
             size = scale_text(BASE_STAT, panel_width_mm), fontface = "bold", color = "grey30", parse = TRUE) +
    coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
    th +
    guides(color = "none",  # color legend is merged with Panel A at composite level
           shape = guide_legend(nrow = 1, order = 2),
           size  = guide_legend(nrow = 1, order = 3)) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(t = 0),
          axis.title = element_blank(),
          axis.text = element_blank())
  
  ttl <- "Pathway Reversal"
  sub <- sprintf("\u03c1 = %.2f (p = %.2f) | %.0f%% reversed",
                 nes_rho$estimate, nes_rho$p.value, pw_rev_frac * 100)
  
  list(
    plot           = p + labs(title = ttl, subtitle = sub) +
      guides(color = guide_legend(nrow = 1, title = NULL, order = 1)),  # standalone gets its own color legend back
    plot_composite = strip_for_composite(p),
    title = ttl, subtitle = sub,
    nes_df = nes_df
  )
}