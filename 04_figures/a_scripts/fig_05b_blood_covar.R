# =============================================================================
# Figure 05: Protein Reversal | Pathway Reversal | Fry | Heatmap-Sankey
# =============================================================================
# Lean, final-figure-only script. Everything not needed to build the 4 panels
# below (ORA bar plots, GSEA dot plots, GSEA overlay curves, sig-count
# diagnostics, shared-pathway table export, pathway core-enrichment scatter)
# has been moved to fig_05_extras.R -- see that script for anything not
# reproduced here.

# --- 0: Setup ----------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readxl, readr, dplyr, tidyr, tibble, stringr, purrr,
  ggplot2, ggrepel, patchwork, ggtext, cowplot, scales
)

# --- 0a: Directory Config ----------------------------------------------------
config <- list(
  gsea_muscle      = "03_analysis/b_reports/muscle/blood_covar/SURVvCTL/gsea",
  gsea_data        = "03_analysis/c_data/blood_covar/gsea_consolidated/gsea_pathway_data.rds",
  limma_muscle     = "03_analysis/b_reports/muscle/blood_covar/all_data/limma",
  fry_barcode_script = "04_figures/a_scripts/06b_fry_barcode_muscle_SURV.R",
  sankey_script      = "04_figures/a_scripts/heatmap_sankey_muscle_surv.R",
  fig_dir          = "04_figures/blood_covar/fig_05/new",
  data_dir         = "04_figures/blood_covar/data/fig_05/new"
)

dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- 0b: Color Palette + Theme + save_plot() ----------------------------------
source("04_figures/a_scripts/00_theme.R")

# device = cairo_pdf (not the base "pdf" device) -- required for Unicode
# glyphs (arrows in quadrant labels, Greek letters in pathway names like
# "IFN\u03b3 Response") to render instead of throwing mbcsToSbcs conversion
# warnings and silently dropping the character.
save_plot <- function(p, name, width = 7, height = 6) {
  ggsave(file.path(config$fig_dir, paste0(name, ".pdf")), p,
         device = cairo_pdf, width = width, height = height)
  ggsave(file.path(config$fig_dir, paste0(name, ".png")), p,
         device = "png", width = width, height = height, dpi = 300)
}

# Shared color scale for reversal plots (protein + pathway).
# Blue = concordant/reversed, Red = discordant/exacerbated, Gray = neither.
pal_shared <- c(
  "Both"     = "#4E8B6F",
  "PRE v CTL only" = "#2471A3",
  "POST v PRE only" = "#C0392B",
  "Neither"  = "#9E9E9E"
)
pal_protein <- pal_shared

# --- 0c: Shared quadrant-label styling ---------------------------------------
# Keeps the protein-reversal scatter and pathway-reversal scatter rendering
# identically when placed side by side as fig_05 panels A/B.
geom_quadrant_label <- function(df, size = 3) {
  geom_label(
    data          = df,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes   = FALSE,
    fill          = df$bg_color,
    color         = colors$white,
    alpha         = 0.7,
    size          = size,
    fontface      = "bold",
    lineheight    = 1.2,
    label.size    = NA,
    label.padding = unit(0.25, "lines"),
    show.legend   = FALSE
  )
}

build_quad_label_df <- function(counts, quad_ids, titles, label_fn,
                                blue_quads, red_quads) {
  tibble(
    q     = quad_ids,
    title = titles,
    x     = c(-Inf,  Inf,  Inf, -Inf),
    y     = c( Inf,  Inf, -Inf, -Inf),
    hjust = c(   0,    1,    1,    0),
    vjust = c(   1,    1,    0,    0)
  ) |>
    left_join(counts, by = "q") |>
    mutate(
      label    = label_fn(title, pick(everything())),
      bg_color = if_else(q %in% blue_quads, colors$blue_dark, colors$red_dark)
    )
}

# --- 0d: Comparison Labels ---------------------------------------------------
surv_labels <- c(
  "Baseline_SURVvCTL" = "SURV_T1 v CTL",
  "Training_SURVvCTL" = "SURV_T2 v CTL",
  "Training_SURV"     = "SURV_T2 v SURV_T1"
)

remap_label <- function(x) {
  x |>
    str_replace_all("SURV", "CRC") |>
    str_replace_all("_T1",  "_PRE") |>
    str_replace_all("_T2",  "_POST") |>
    str_replace_all("T1",   "PRE") |>
    str_replace_all("T2",   "POST")
}

surv_labels <- setNames(remap_label(surv_labels), names(surv_labels))

# =============================================================================
# 1: Load Data
# =============================================================================

if (!file.exists(config$gsea_data)) {
  stop("Consolidated GSEA data not found at ", config$gsea_data,
       " -- run 02_gsea_pathway_consolidation.R first.")
}
gsea_pathway_data <- readRDS(config$gsea_data)

limma_results <- list(
  Baseline_SURVvCTL = read_excel(
    file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 11),
  Training_SURV     = read_excel(
    file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 13)
)

load_gsea_rds <- function(dir) {
  files <- list.files(dir, pattern = "^gsea_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^gsea_results_|\\.rds$", "",
                list.files(dir, pattern = "^gsea_results_.*\\.rds$"))
  setNames(lapply(files, readRDS), nms)
}

# Needed for the REVERSAL/EXACERBATION ORA results flanking panel A --
# these are protein-level ORA, not one of the three contrasts
# 02_gsea_pathway_consolidation.R pools, so they're loaded raw here rather
# than via gsea_pathway_data.
gsea_muscle <- load_gsea_rds(config$gsea_muscle)

# =============================================================================
# 2: Pathway (GSEA / NES) Reversal Data  -->  Panel B
# =============================================================================
# Sourced from the consolidated GSEA pathway data, restricted to the pooled
# redundancy-filter survivors (representative_lookup) built in
# 02_gsea_pathway_consolidation.R. See fig_05_extras.R history for the two
# intentional behavior changes vs. the old ad hoc per-script filter
# (pooled-3-contrast significance, Hallmark now filtered too).
reversal_group <- "SURV"
reversal_dbs   <- c("GO Slim", "Hallmark")

nes_pre <- gsea_pathway_data$gsea_supp_full %>%
  filter(group == reversal_group, tissue == "Muscle",
         database %in% reversal_dbs, comparison == "Baseline_SURVvCTL") %>%
  transmute(ID, Description, database = as.character(database), setSize,
            NES_PRE = NES, padj_PRE = p.adjust)

nes_ret <- gsea_pathway_data$gsea_supp_full %>%
  filter(group == reversal_group, tissue == "Muscle",
         database %in% reversal_dbs, comparison == "Training_SURV") %>%
  transmute(ID, Description, database = as.character(database),
            NES_RET = NES, padj_RET = p.adjust)

rep_lookup_reversal <- gsea_pathway_data$representative_lookup %>%
  filter(group == reversal_group, tissue == "Muscle", database %in% reversal_dbs) %>%
  mutate(database = as.character(database))

nes_df <- full_join(nes_pre, nes_ret, by = c("ID", "Description", "database")) %>%
  filter(padj_PRE < 0.05 | padj_RET < 0.05) %>%
  semi_join(rep_lookup_reversal, by = c("database", "ID")) %>%
  mutate(
    trajectory = case_when(
      sign(NES_PRE) != sign(NES_RET)                                ~ "reversed",
      sign(NES_PRE) == sign(NES_RET) & abs(NES_RET) > abs(NES_PRE)  ~ "exacerbated",
      sign(NES_PRE) == sign(NES_RET) & abs(NES_RET) < abs(NES_PRE)  ~ "attenuated",
      TRUE                                                           ~ "unchanged"
    ),
    sig_in = case_when(
      padj_PRE < 0.05 & padj_RET < 0.05 ~ "Both",
      padj_PRE < 0.05                    ~ "PRE v CTL only",
      padj_RET < 0.05                    ~ "POST v PRE only"
    ),
    go_set = database
  )

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

nes_df <- nes_df |>
  mutate(label = coalesce(pathway_label_map[Description], Description))

# =============================================================================
# 3: Protein (Limma) Reversal Data  -->  Panel A
# =============================================================================

df_PRE <- limma_results$Baseline_SURVvCTL |>
  dplyr::select(uniprot_id, gene,
                logFC_CRC    = logFC,
                Pi.Val_CRC = Pi.Val,
                adj.P.Val_CRC  = adj.P.Val)

df_RET <- limma_results$Training_SURV |>
  dplyr::select(uniprot_id, gene,
                logFC_RET    = logFC,
                Pi.Val_RET = Pi.Val,
                adj.P.Val_RET  = adj.P.Val)

df_rev <- inner_join(df_PRE, df_RET, by = c("uniprot_id", "gene")) |>
  mutate(
    sig = case_when(
      Pi.Val_CRC < 0.05 & Pi.Val_RET < 0.05 ~ "Both",
      Pi.Val_CRC < 0.05                          ~ "PRE v CTL only",   
      Pi.Val_RET < 0.05                          ~ "POST v PRE only",
      TRUE                                           ~ "Neither"
    )
  )

lim <- max(abs(c(df_rev$logFC_CRC, df_rev$logFC_RET)), na.rm = TRUE)

protein_label_df <- df_rev |>
  filter(sig != "Neither", abs(logFC_CRC) > 1 | abs(logFC_RET) > 1) |>
  mutate(
    Pi.Val_RET = pmax(Pi.Val_RET, 1e-300),
    Pi.Val_CRC = pmax(Pi.Val_CRC, 1e-300),
    direction      = if_else(abs(logFC_CRC) >= abs(logFC_RET),
                             if_else(logFC_CRC > 0, "Up", "Down"),
                             if_else(logFC_RET > 0, "Up", "Down")),
    score = pmax(-log10(Pi.Val_RET), -log10(Pi.Val_CRC))
  ) |>
  group_by(sig, direction) |>
  slice_max(score, n = 5, with_ties = FALSE) |>
  ungroup()

# =============================================================================
# 4: Panel A -- Protein Reversal Scatter
# =============================================================================

quad_counts <- df_rev |>
  mutate(q = case_when(
    logFC_CRC < 0 & logFC_RET > 0 ~ "TL",
    logFC_CRC > 0 & logFC_RET > 0 ~ "TR",
    logFC_CRC > 0 & logFC_RET < 0 ~ "BR",
    logFC_CRC < 0 & logFC_RET < 0 ~ "BL"
  )) |>
  filter(!is.na(q)) |>
  group_by(q) |>
  summarise(n_total = n(), n_sig = sum(sig != "Neither", na.rm = TRUE), .groups = "drop")

quad_ann <- build_quad_label_df(
  counts     = quad_counts,
  quad_ids   = c("TL", "TR", "BR", "BL"),
  titles     = c("Reversed (\u2193CRC, \u2191RET)", "Exacerbated Up",
                 "Reversed (\u2191CRC, \u2193RET)", "Exacerbated Down"),
  label_fn   = function(title, row) sprintf("%s\n%d / %d", title, row$n_sig, row$n_total),
  blue_quads = c("TL", "BR"),
  red_quads  = c("TR", "BL")
)

p_scatter_core <- ggplot(df_rev, aes(x = logFC_CRC, y = logFC_RET,
                                     color = sig, fill = sig)) +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = 0,    ymax =  Inf, fill = colors$blue, alpha = 0.2) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin = 0,    ymax =  Inf, fill = colors$red,  alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax = 0,    ymin = -Inf, ymax =  0,   fill = colors$red,  alpha = 0.2) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin = -Inf, ymax =  0,   fill = colors$blue, alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.6, color = "grey30", linetype = "dashed") +
  geom_point(data = filter(df_rev, sig == "Neither"),
             size = 1.2, alpha = 0.5, shape = 16) +
  geom_point(data = filter(df_rev, sig != "Neither"),
             size = 2.2, alpha = 0.85, shape = 16) +
  scale_color_manual(values = pal_protein, name = NULL) +
  scale_fill_manual(values = pal_protein, guide = "none") +
  geom_quadrant_label(quad_ann, size = 2.8) +
  geom_label_repel(
    data          = protein_label_df,
    aes(label = gene, fill = sig),
    color         = colors$white,
    size          = 2.5,
    max.overlaps  = 20,
    box.padding   = 0.3,
    label.padding = unit(0.15, "lines"),
    label.size    = NA,
    segment.color = "grey60",
    segment.size  = 0.3,
    show.legend   = FALSE,
    xlim          = c(-lim * 0.85, lim * 0.85),
    ylim          = c(-lim * 0.85, lim * 0.85)
  ) +
  annotate("text",
           x = -lim * 0.80, y = 0.01,
           label = expression(log[2]*FC~"(PRE vs CTL)"),
           hjust = 0.5, vjust = -0.5,
           size = 3.2, fontface = "bold", color = "grey20") +
  annotate("text",
           x = -0.01, y = -lim * 0.80,
           label = expression(log[2]*FC~"(POST vs PRE)"),
           hjust = 0.5, vjust = 1.3,
           angle = 90,
           size = 3.2, fontface = "bold", color = "grey20") +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim), clip = "off") +
  theme_cr() +
  theme(
    legend.position = "bottom",
    axis.title      = element_blank(),
    axis.text       = element_blank(),
    axis.ticks      = element_blank(),
    axis.line       = element_blank()
  )

axis_breaks <- pretty(c(-lim, lim))
axis_breaks <- axis_breaks[axis_breaks != 0 & abs(axis_breaks) < lim]

p_scatter_core <- p_scatter_core +
  annotate("segment",
           x = axis_breaks, xend = axis_breaks,
           y = -lim * 0.015, yend = lim * 0.015,
           color = "grey30", linewidth = 0.4) +
  annotate("text",
           x = axis_breaks, y = -lim * 0.03,
           label = axis_breaks,
           vjust = 1, size = 2.4, color = "grey20") +
  annotate("segment",
           y = axis_breaks, yend = axis_breaks,
           x = -lim * 0.015, xend = lim * 0.015,
           color = "grey30", linewidth = 0.4) +
  annotate("text",
           y = axis_breaks, x = -lim * 0.03,
           label = axis_breaks,
           hjust = 1, size = 2.4, color = "grey20")

# =============================================================================
# 4b: Panel A -- ORA Flanking Bars for Protein Reversal Scatter
# =============================================================================

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
  "mitochondrial electron transport, NADH to ubiquinone"            = "Mitochondrial\nElectron Transport",
  "nucleoside phosphate biosynthetic process"                       = "Nucleoside-P\nBiosynthesis",
  "translation"                                                     = "Translation",
  "NADH dehydrogenase complex assembly"                             = "NADH\nDehydrogenase Assembly"
)

apply_labels <- function(df, map) {
  df |> mutate(label = if_else(Description %in% names(map), map[Description], Description))
}

# Kept separate from the pooled-Jaccard filter in 02_gsea_pathway_consolidation.R
# -- REVERSAL/EXACERBATION ORA results are protein-level ORA, not one of the
# three contrasts that script currently pools.
jaccard_filter_ora <- function(res,
                               threshold = 0.5,
                               padj_col  = "p.adjust",
                               id_col    = "ID",
                               gene_col  = "geneID") {
  jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
  df <- as.data.frame(res)
  df <- df[grepl("^GO:", df[[id_col]]), ]
  df <- df[order(df[[padj_col]]), ]
  sets <- strsplit(df[[gene_col]], "/")
  n    <- nrow(df)
  keep <- rep(TRUE, n)
  for (i in seq_len(n - 1)) {
    if (!keep[i]) next
    for (j in seq(i + 1, n)) {
      if (!keep[j]) next
      if (jaccard(sets[[i]], sets[[j]]) >= threshold) keep[j] <- FALSE
    }
  }
  df[keep, ]
}

ora_rev_up_top6 <- jaccard_filter_ora(gsea_muscle$REVERSAL$ora_bp_up@result) |>
  arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
ora_rev_down_top6 <- jaccard_filter_ora(gsea_muscle$REVERSAL$ora_bp_down@result) |>
  arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
ora_ex_up_top6 <- jaccard_filter_ora(gsea_muscle$EXACERBATION$ora_bp_up@result) |>
  arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)
ora_ex_down_top6 <- jaccard_filter_ora(gsea_muscle$EXACERBATION$ora_bp_down@result) |>
  arrange(p.adjust) |> slice_head(n = 6) |> apply_labels(ora_label_map)

make_ora_right <- function(df, fill_color_low, fill_color_high, title = "", show_x_title = FALSE) {
  df <- df |>
    mutate(
      neg_log10_p = -log10(p.adjust),
      y_order     = reorder(label, neg_log10_p),
      sig_label   = case_when(
        p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**", p.adjust < 0.05 ~ "*", TRUE ~ "ns"
      )
    )
  ggplot(df, aes(x = neg_log10_p, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    geom_text(aes(label = label), x = 0.05, hjust = 0, size = 3.5, colour = "white", fontface = "bold") +
    geom_text(aes(x = neg_log10_p, label = sig_label), hjust = -0.2, vjust = 0.75, size = 4, colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_color_low, high = fill_color_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18)), position = "bottom") +
    scale_y_discrete(position = "left") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL, title = NULL) +
    theme_cr() +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7), axis.ticks.x = element_line(),
      panel.grid = element_blank(), legend.position = "none",
      plot.title = element_text(size = 9, hjust = 0, color = fill_color_high, face = "bold"),
      plot.margin = margin(2, 10, 2, 0)
    )
}

make_ora_left <- function(df, fill_color_low, fill_color_high, title = "", show_x_title = FALSE) {
  df <- df |>
    mutate(
      neg_log10_p     = -log10(p.adjust),
      neg_log10_p_neg = -neg_log10_p,
      y_order         = reorder(label, neg_log10_p),
      sig_label       = case_when(
        p.adjust < 0.001 ~ "***", p.adjust < 0.01 ~ "**", p.adjust < 0.05 ~ "*", TRUE ~ "ns"
      )
    )
  ggplot(df, aes(x = neg_log10_p_neg, y = y_order)) +
    geom_col(aes(fill = neg_log10_p), width = 0.7) +
    geom_text(aes(x = neg_log10_p_neg * 0.97, label = label), hjust = 0, size = 3.5, colour = "white", fontface = "bold") +
    geom_text(aes(x = neg_log10_p_neg, label = sig_label), hjust = 1.4, vjust = 0.75, size = 4, colour = "black", fontface = "bold") +
    scale_fill_gradient(low = fill_color_low, high = fill_color_high, guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0.18, 0)), labels = function(x) abs(x), position = "bottom") +
    scale_y_discrete(position = "right") +
    coord_cartesian(clip = "off") +
    labs(x = if (show_x_title) expression(-log[10](p[adj])) else NULL, y = NULL, title = NULL) +
    theme_cr() +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 7), axis.ticks.x = element_line(),
      panel.grid = element_blank(), legend.position = "none",
      plot.title = element_text(size = 9, hjust = 1, color = fill_color_high, face = "bold"),
      plot.margin = margin(2, 0, 2, 10)
    )
}

left_col  <- make_ora_left(ora_rev_up_top6,   colors$blue_light, colors$blue_dark, "Reversed (\u2193CRC \u2191RET)", FALSE) /
  make_ora_left(ora_ex_down_top6,  colors$red_light,  colors$red_dark,  "Exacerbated Down", TRUE)
right_col <- make_ora_right(ora_ex_up_top6,   colors$red_light,  colors$red_dark,  "Exacerbated Up", FALSE) /
  make_ora_right(ora_rev_down_top6, colors$blue_light, colors$blue_dark, "Reversed (\u2191CRC \u2193RET)", TRUE)

# p_scatter_core here is the fully annotated version from Section 4 above
# (quadrant labels, gene labels, inner axis ticks) -- this is the center
# panel of the compound; nothing about it changes for the flanking bars.
fig_05_compound <- (left_col | p_scatter_core | right_col) +
  plot_layout(widths = c(1.1, 2, 1.1)) &
  theme(plot.margin = margin(2, 0, 2, 0))

fig_05_compound

# =============================================================================
# 5: Panel B -- Pathway (NES) Reversal Scatter
# =============================================================================

nes_rho <- cor.test(nes_df$NES_PRE, nes_df$NES_RET, method = "spearman", exact = FALSE)
rho_label <- sprintf("\u03c1 = %.3f (p = %.2f)", nes_rho$estimate, nes_rho$p.value)

nes_label_df <- nes_df |>
  mutate(
    quadrant    = case_when(
      NES_PRE >= 0 & NES_RET >= 0 ~ "PP",
      NES_PRE <  0 & NES_RET <  0 ~ "NN",
      NES_PRE >= 0 & NES_RET <  0 ~ "PN",
      NES_PRE <  0 & NES_RET >= 0 ~ "NP"
    ),
    abs_nes_sum = abs(NES_PRE) + abs(NES_RET)
  ) |>
  group_by(quadrant) |>
  slice_max(abs_nes_sum, n = 10, with_ties = FALSE) |>
  ungroup()

nes_lim <- max(abs(c(nes_df$NES_PRE, nes_df$NES_RET)), na.rm = TRUE)

nes_quad_counts <- nes_df |>
  mutate(q = case_when(
    NES_PRE <  0 & NES_RET >= 0 ~ "TL",
    NES_PRE >= 0 & NES_RET >= 0 ~ "TR",
    NES_PRE >= 0 & NES_RET <  0 ~ "BR",
    NES_PRE <  0 & NES_RET <  0 ~ "BL"
  )) |>
  group_by(q) |>
  summarise(n_total = n(), n_sig = sum(sig_in != "Neither", na.rm = TRUE), .groups = "drop")

nes_quad_labels <- build_quad_label_df(
  counts     = nes_quad_counts,
  quad_ids   = c("TL", "TR", "BR", "BL"),
  titles     = c("Reversed (Up RET)", "Exacerbated (Up)",
                 "Reversed (Down RET)", "Exacerbated (Down)"),
  label_fn   = function(title, row) sprintf("%s\nn = %d (%s%%)", title, row$n_total,
                                            round(row$n_total / sum(nes_quad_counts$n_total) * 100, 1)),
  blue_quads = c("TL", "BR"),
  red_quads  = c("TR", "BL")
)

nes1 <- ggplot(nes_df, aes(x = NES_PRE, y = NES_RET,
                           color = sig_in, fill = sig_in)) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin =    0, ymax =  Inf, fill = colors$red, alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin = -Inf, ymax =    0, fill = colors$red, alpha = 0.2) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin = -Inf, ymax =    0, fill = colors$blue,  alpha = 0.2) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin =    0, ymax =  Inf, fill = colors$blue,  alpha = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_point(data = filter(nes_df, sig_in == "Neither"),
             aes(shape = go_set, size = setSize), alpha = 0.8) +
  geom_point(data = filter(nes_df, sig_in != "Neither"),
             aes(shape = go_set, size = setSize), alpha = 0.85) +
  scale_color_manual(values = pal_shared, name = NULL) +
  scale_fill_manual(values  = pal_shared, guide = "none") +
  scale_shape_manual(values = c("GO Slim" = 16, "Hallmark" = 17), name = "Database") +
  scale_size_continuous(name = "Set size", range = c(1, 6)) +
  labs(x = NULL, y = NULL, title = "Pathway Reversal", subtitle = rho_label) +
  annotate("text", x = -Inf, y = 0, label = " NES~'(CRC_PRE vs CTL)'",
           hjust = 0, vjust = -0.3, size = 3, fontface = "bold", color = "grey30", parse = TRUE) +
  annotate("text", x = 0, y = -Inf, label = " NES~'(CRC_POST vs CRC_PRE)'",
           hjust = 0, vjust = -0.3, angle = 90,
           size = 3, fontface = "bold", color = "grey30", parse = TRUE) +
  coord_cartesian(xlim = c(-nes_lim, nes_lim), ylim = c(-nes_lim, nes_lim)) +
  theme_cr() +
  theme(
    plot.title    = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0, size = 9, color = colors$gray_dark)
  )

nes2 <- nes1 + geom_quadrant_label(nes_quad_labels)

nes3 <- nes2 +
  geom_label_repel(
    data          = nes_label_df,
    aes(label = label, fill = sig_in),
    color         = colors$white,
    size          = 2.5,
    max.overlaps  = 20,
    box.padding   = 0.3,
    label.padding = unit(0.15, "lines"),
    label.size    = NA,
    segment.color = "grey60",
    segment.size  = 0.3,
    show.legend   = FALSE,
    xlim          = c(-nes_lim * 0.85, nes_lim * 0.85),
    ylim          = c(-nes_lim * 0.85, nes_lim * 0.85)
  )

nes5 <- nes3 +
  guides(
    color = guide_legend(nrow = 1, title = NULL, order = 1),
    shape = guide_legend(nrow = 1, order = 2),
    size  = guide_legend(nrow = 1, order = 3)
  ) +
  theme(
    legend.position   = "bottom",
    legend.box        = "vertical",
    legend.margin     = margin(t = 0),
    legend.spacing.y  = unit(0.2, "cm"),
    legend.key.size   = unit(0.4, "cm"),
    legend.text       = element_text(size = 8)
  )

# =============================================================================
# 6: Panel C -- Fry Barcode Plot (sourced live from 06b_fry_barcode_muscle_SURV.R)
# =============================================================================
config$fry_barcode_script <- "04_figures/a_scripts/06b_fry_barcode_muscle_SURV.R"

if (!file.exists(config$fry_barcode_script)) {
  stop("Fry barcode script not found at ", config$fry_barcode_script)
}

fry_env <- new.env()
source(config$fry_barcode_script, local = fry_env)

if (!exists("fry_barcode_fig", envir = fry_env)) {
  stop("fry_barcode_fig not found after sourcing ", config$fry_barcode_script,
       " -- check the object name still matches.")
}
p_fry <- fry_env$fry_barcode_fig

# =============================================================================
# 7: Panel D -- Heatmap + Sankey (sourced live, PANEL_SOURCE_MODE)
# =============================================================================
config$sankey_script <- "04_figures/a_scripts/fig_heatmap_sankey_from_ora_2.R"

if (!file.exists(config$sankey_script)) {
  stop("Heatmap/sankey script not found at ", config$sankey_script)
}

sankey_env <- new.env()
sankey_env$PANEL_SOURCE_MODE <- TRUE
source(config$sankey_script, local = sankey_env)

if (!exists("fig_ribbon_panel", envir = sankey_env)) {
  stop("fig_ribbon_panel not found after sourcing ", config$sankey_script,
       " -- check the object name still matches.")
}
p_sankey <- sankey_env$fig_ribbon_panel

# =============================================================================
# 7b: Shared A/B Legend (built once, extracted, applied to both panels)
# =============================================================================
# fig_05_compound and nes5 both use the same 4-category color scale (pal_shared,
# now aliased by pal_protein too) but patchwork's guides="collect" can't merge
# them because fig_05_compound is wrap_elements()-wrapped (it's a nested
# patchwork itself, and wrapping renders it as one opaque block -- collect
# can't see inside). Building and placing the legend manually instead.

# No-legend copy of the compound for the composite (standalone version with
# its own legend is still saved separately below, unchanged).
fig_05_compound_nolegend <- fig_05_compound &
  theme(legend.position = "none")

# nes5 keeps its shape (database) + size (setSize) legend -- only the color
# scale is being merged with panel A, since shape/size don't apply there.
nes5_embed <- nes5 +
  guides(color = "none")

# Minimal dummy plot purely to extract a clean single-row color legend.
shared_legend_plot <- ggplot(nes_df, aes(NES_PRE, NES_RET, color = sig_in)) +
  geom_point() +
  scale_color_manual(values = pal_shared, name = NULL) +
  guides(color = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "bottom", legend.text = element_text(size = 10))

shared_legend_grob <- cowplot::get_legend(shared_legend_plot)

# =============================================================================
# 8: Compose Final 4-Panel Figure
# =============================================================================

# --- Target print size -- pick this FIRST, everything else derives from it.
# 180mm = standard double-column journal width. Adjust to your target journal.
COMP_W <- 180
COMP_H <- 170
PRINT_SCALE <- COMP_W / 90   # 90mm = the width each sub-panel script was
# effectively authored for as a quarter of a
# ~180mm 2x2 grid; adjust the reference if your
# sub-panels were sized differently
TAG_SZ <- round(9 * PRINT_SCALE * 0.85)

layout <- paste(
  "AAAAAAAABBBB",
  "AAAAAAAABBBB",
  "AAAAAAAABBBB",
  "############",
  "CCCCDDDDDDDD",
  "CCCCDDDDDDDD",
  sep = "\n"
)

fig_raw <- wrap_elements(full = fig_05_compound_nolegend) +
  wrap_elements(full = nes5_embed) +
  wrap_elements(full = p_fry) +
  wrap_elements(full = p_sankey) +
  plot_layout(design = layout,
              widths  = rep(1, 12),
              heights = c(rep(6, 3), 1, rep(6, 2)))

fig_05 <- ggdraw(fig_raw) +
  draw_plot(shared_legend_grob, x = 0.30, y = 0.475, width = 0.4, height = 0.03) +
  draw_label("A", x = 0.01, y = 0.99, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("B", x = 0.67, y = 0.99, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("C", x = 0.01, y = 0.44, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1) +
  draw_label("D", x = 0.34, y = 0.44, size = TAG_SZ, fontface = "bold", hjust = 0, vjust = 1)

fig_05

ggsave(file.path(config$fig_dir, "fig_05_composite.pdf"), fig_05,
       width = COMP_W, height = COMP_H, units = "mm", device = cairo_pdf)
ggsave(file.path(config$fig_dir, "fig_05_composite.png"), fig_05,
       width = COMP_W, height = COMP_H, units = "mm", dpi = 300)
