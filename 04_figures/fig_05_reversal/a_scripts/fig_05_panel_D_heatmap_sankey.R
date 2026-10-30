# =============================================================================
# GO ORA (BP + CC + MF) — Muscle SURV contrasts
# + Protein-to-Pathway Figure (two versions: individual flows + ribbon flows)
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readxl, dplyr, tidyr, tibble, stringr, purrr, scales,
  ggplot2, ggnewscale, cowplot,
  clusterProfiler, org.Hs.eg.db, AnnotationDbi
)

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

config <- list(
  report_dir = "03_analysis/b_reports/muscle",
  SURV_limma = "03_analysis/c_data/muscle/limma/limma_results_all_muscle_cv.xlsx",
  figure_dir = "04_figures/fig_05_reversal/c_data"
)

OUTPUT_DIR <- file.path(config$figure_dir)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

organism        <- org.Hs.eg.db
PI_THRESH       <- 0.05
ORA_PVAL_CUTOFF <- 0.05
ORA_QVAL_CUTOFF <- 0.20
ORA_MIN_GSET    <- 5
ORA_MAX_GSET    <- 500
SIMPLIFY_CUTOFF <- 0.7

CONTRASTS <- c("Cancer", "Residual", "Training")

CONTRAST_LABELS <- c(
  Cancer = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training = "CRC+RET vs CRC"
)

# Figure layout constants
X_SIG    <- -0.6
X_HEAT   <- 0:2
X_TRAJ   <- 2.6
X_FLOW_R <- 5.5
X_BAR    <- 5.5
BAR_W    <- 4.0

FIG_WIDTH  <- 12
FIG_HEIGHT <- 8

traj_pal <- c(
  Reversed     = "#2471A3",  # blue
  Attenuated   = "#4E8B6F",  # green
  Exacerbated  = "#C0392B",  # red
  Unclassified = "#9E9E9E"   # gray
)
traj_levels <- c("Reversed", "Attenuated", "Exacerbated", "Unclassified")

# Panel sizing -----------------------------------------------------------
# PANEL_WIDTH_MM (mm) -- set this in the parent script's env *before*
# source() to control embedded sizing; PANEL_SOURCE_MODE suppresses this
# script's own standalone ggsave/xlsx exports.
#
# PANEL_SCALE now uses the same sqrt-scaling law as every other panel in
# the composite (scale_text(1, PANEL_WIDTH_MM), from fig_05_style.R) rather
# than a linear PANEL_WIDTH/FIG_WIDTH ratio -- keeps this panel's text
# consistent with how Panels A/B/C scale when embedded at a non-standalone
# width. scale_text() must already be defined in this environment (the
# stitcher assigns it into sankey_env before sourcing this file); falls
# back to a plain linear ratio if sourced standalone without it.
PANEL_SOURCE_MODE <- exists("PANEL_SOURCE_MODE") && isTRUE(PANEL_SOURCE_MODE)
PANEL_WIDTH_MM <- if (exists("PANEL_WIDTH_MM")) PANEL_WIDTH_MM else 178
PANEL_SCALE    <- if (exists("scale_text")) scale_text(1, PANEL_WIDTH_MM) else (PANEL_WIDTH_MM / 25.4) / FIG_WIDTH
PANEL_WIDTH    <- PANEL_WIDTH_MM / 25.4            # inches, only used for the FIG_WIDTH ratio below
PANEL_HEIGHT   <- PANEL_WIDTH * (FIG_HEIGHT / FIG_WIDTH)

# Legend layout ---------------------------------------------------------------
# All three legends (Significant in, logFC, Trajectory) live OUTSIDE the
# plot's data space entirely (composed via cowplot::plot_grid AFTER the
# ggplot is built), stacked as a single column to the right of the plot.
# This replaces the earlier annotation_custom()-in-data-space approach,
# which broke down because the footer band (title_y..bar_axis_y) actually
# fell inside the heatmap/ribbon's real y-range (protein rows run all the
# way down to y_pos = 1), so the legends overlapped live data.
# Tune by eye once you see the composed panel -- LEGEND_COL_REL_HEIGHTS
# controls the vertical split among the three stacked legends (sig pattern
# tends to need the most room since it can wrap to several rows).
LEGEND_COL_REL_WIDTH   <- 0.16              # right-hand legend column, fraction of total width
LEGEND_COL_REL_HEIGHTS <- c(0.45, 0.30, 0.25)  # sig_pattern, logFC, trajectory -- must sum to 1

# =============================================================================
# 1. LOAD SIG PROTEINS
# =============================================================================

if (!exists("sig_proteins")) {
  
  CONTRAST_SHEETS <- list(
    Cancer = "Cancer",
    Residual = "Residual",
    Training     = "Training"
  )
  
  da_tbl <- imap(
    map(CONTRAST_SHEETS, ~ read_excel(config$SURV_limma, sheet = .x)),
    function(df, contrast) {
      df %>% as_tibble() %>%
        transmute(
          uniprot  = uniprot, gene = gene, ensembl = Ensembl,
          logFC = logFC, pi_val = Pi.Val,
          sig = Pi.Val < PI_THRESH, contrast = contrast
        )
    }
  ) %>% bind_rows()
  
  SIG_COLS <- paste0("sig_", CONTRASTS)
  
  sig_proteins <- da_tbl %>%
    dplyr::select(uniprot, gene, ensembl, contrast, logFC) %>%
    group_by(uniprot) %>%
    mutate(ensembl = dplyr::first(na.omit(ensembl)),
           gene    = dplyr::first(na.omit(gene))) %>%
    ungroup() %>%
    pivot_wider(names_from = contrast, values_from = logFC,
                id_cols = c(uniprot, gene, ensembl)) %>%
    left_join(
      da_tbl %>%
        dplyr::select(uniprot, contrast, sig) %>%
        pivot_wider(names_from = contrast, values_from = sig,
                    values_fill = list(sig = FALSE)) %>%
        rename_with(~ paste0("sig_", .), all_of(CONTRASTS)),
      by = "uniprot"
    ) %>%
    filter(if_any(all_of(SIG_COLS), ~ .x == TRUE))
  
  message(sprintf("Loaded %d significant proteins", nrow(sig_proteins)))
}

# =============================================================================
# 2. ENTREZ ID CONVERSION
# =============================================================================

sig_entrez <- bitr(unique(na.omit(sig_proteins$ensembl)),
                   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)

universe_entrez <- bitr(unique(na.omit(da_tbl$ensembl)),
                        fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)

# =============================================================================
# 3. RUN enrichGO FOR EACH ONTOLOGY + SIMPLIFY
# =============================================================================

run_go_ora <- function(ont) {
  message(sprintf("\nRunning GO ORA: %s", ont))
  res <- enrichGO(
    gene = sig_entrez$ENTREZID, universe = universe_entrez$ENTREZID,
    OrgDb = organism, ont = ont, pAdjustMethod = "BH",
    pvalueCutoff = ORA_PVAL_CUTOFF, qvalueCutoff = ORA_QVAL_CUTOFF,
    minGSSize = ORA_MIN_GSET, maxGSSize = ORA_MAX_GSET, readable = FALSE
  )
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message(sprintf("  No enriched terms for %s", ont)); return(NULL)
  }
  res_simple <- simplify(res, cutoff = SIMPLIFY_CUTOFF,
                         by = "p.adjust", select_fun = min)
  df <- as.data.frame(res_simple) %>% 
    mutate(ontology = ont, n_proteins = str_count(geneID, "/") + 1L)
  message(sprintf("  Enriched: %d -> after simplify: %d",
                  nrow(as.data.frame(res)), nrow(df)))
  saveRDS(res_simple, file.path(OUTPUT_DIR,
                                sprintf("ora_go%s_muscle_surv.rds", tolower(ont))))
  df
}

go_results      <- list(BP = run_go_ora("BP"), CC = run_go_ora("CC"), MF = run_go_ora("MF"))
ora_df_combined <- bind_rows(go_results) %>% arrange(ontology, p.adjust)

write.csv(ora_df_combined,
          file.path(OUTPUT_DIR, "ora_go_all_muscle_surv_full.csv"), row.names = FALSE)

# =============================================================================
# 4. PATHWAY LABEL MAP & GROUP ASSIGNMENTS
# =============================================================================

pathway_label_map <- c(
  # Cytoskeleton
  "GO:0097435" = "Supramolecular fiber organization",
  "GO:0030036" = "Actin cytoskeleton organization",
  "GO:0030029" = "Actin filament-based process",
  "GO:0045103" = "Intermediate filament-based process",
  "GO:0099512" = "Supramolecular fiber",
  "GO:0015629" = "Actin cytoskeleton",
  "GO:0003779" = "Actin binding",
  
  # Muscle Development
  "GO:0055001" = "Muscle cell development",
  "GO:0030239" = "Myofibril assembly",
  "GO:0061061" = "Muscle structure development",
  "GO:0030016" = "Myofibril",
  
  # Development & Morphogenesis
  "GO:0048513" = "Animal organ development",
  "GO:0009653" = "Anatomical structure morphogenesis",
  "GO:0010927" = "Component assembly in morphogenesis",
  "GO:0009888" = "Tissue development",
  "GO:0048589" = "Developmental growth",
  
  # ECM & Adhesion
  "GO:0007155" = "Cell adhesion",
  "GO:0007044" = "Cell-substrate junction assembly",
  "GO:0150115" = "Cell-substrate junction organization",
  "GO:0031012" = "Extracellular matrix",
  "GO:0070161" = "Anchoring junction",
  "GO:0005201" = "ECM structural constituent",
  "GO:0005539" = "Glycosaminoglycan binding",
  "GO:0048146" = "Fibroblast proliferation (positive reg.)",
  
  # Collagen Trimer
  "GO:0005581" = "Collagen trimer",
  
  # Integrin Signaling & Binding
  "GO:0005178" = "Integrin binding",
  
  # Cell Surface
  "GO:0005938" = "Cell cortex",
  "GO:0031143" = "Pseudopodium",
  "GO:0045177" = "Apical part of cell",
  "GO:0031252" = "Cell leading edge",
  
  # Cell Motility & Locomotion
  "GO:0040011" = "Locomotion",
  
  # Blood & Coagulation
  "GO:0031093" = "Platelet alpha granule lumen",
  
  # Neuronal
  "GO:0043025" = "Neuronal cell body",
  "GO:0036477" = "Somatodendritic compartment",
  "GO:0030424" = "Axon",
  
  # Golgi Lumen
  "GO:0005796" = "Golgi lumen",
  
  # Signaling
  "GO:0005509" = "Calcium ion binding",
  "GO:0061134" = "Peptidase regulator activity",
  "GO:0050839" = "Cell adhesion molecule binding",
  "GO:0140677" = "Molecular function activator activity",
  "GO:0044877" = "Protein complex binding",
  
  # System Process (generic parent term -- consider dropping)
  "GO:0003008" = "System process"
)

pathway_groups <- c(
  "Supramolecular fiber organization"   = "Cytoskeleton",
  "Actin cytoskeleton organization"     = "Cytoskeleton",
  "Actin filament-based process"        = "Cytoskeleton",
  "Intermediate filament-based process" = "Cytoskeleton",
  "Supramolecular fiber"                = "Cytoskeleton",
  "Actin cytoskeleton"                  = "Cytoskeleton",
  "Actin binding"                       = "Cytoskeleton",
  
  "Muscle cell development"             = "Muscle Development",
  "Myofibril assembly"                  = "Muscle Development",
  "Muscle structure development"        = "Muscle Development",
  "Myofibril"                           = "Muscle Development",
  
  "Animal organ development"            = "Development & Morphogenesis",
  "Anatomical structure morphogenesis"  = "Development & Morphogenesis",
  "Component assembly in morphogenesis" = "Development & Morphogenesis",
  "Tissue development"                  = "Development & Morphogenesis",
  "Developmental growth"                = "Development & Morphogenesis",
  
  "Cell adhesion"                       = "ECM & Adhesion",
  "Cell-substrate junction assembly"    = "ECM & Adhesion",
  "Cell-substrate junction organization" = "ECM & Adhesion",
  "Extracellular matrix"                = "ECM & Adhesion",
  "Anchoring junction"                  = "ECM & Adhesion",
  "ECM structural constituent"          = "ECM & Adhesion",
  "Glycosaminoglycan binding"           = "ECM & Adhesion",
  "Fibroblast proliferation (positive reg.)" = "ECM & Adhesion",
  
  "Collagen trimer"                     = "Collagen Trimer",
  
  "Integrin binding"                    = "Integrin Signaling & Binding",
  
  "Cell cortex"                         = "Cell Surface",
  "Pseudopodium"                        = "Cell Surface",
  "Apical part of cell"                 = "Cell Surface",
  "Cell leading edge"                   = "Cell Surface",
  
  "Locomotion"                          = "Cell Motility & Locomotion",
  
  "Platelet alpha granule lumen"        = "Blood & Coagulation",
  
  "Neuronal cell body"                  = "Neuronal",
  "Somatodendritic compartment"         = "Neuronal",
  "Axon"                                = "Neuronal",
  
  "Golgi lumen"                         = "Golgi Lumen",
  
  "Calcium ion binding"                 = "Signaling",
  "Peptidase regulator activity"        = "Signaling",
  "Cell adhesion molecule binding"      = "Signaling",
  "Molecular function activator activity" = "Signaling",
  "Protein complex binding"             = "Signaling",
  
  "System process"                      = "System Process"
)

group_levels <- c("Cytoskeleton", "Muscle Development", "Development & Morphogenesis",
                  "ECM & Adhesion", "Collagen Trimer", "Integrin Signaling & Binding",
                  "Cell Surface", "Cell Motility & Locomotion", "Blood & Coagulation",
                  "Neuronal", "Golgi Lumen", "Signaling", "System Process")

# Reused your existing house palette where the theme carries over (same
# hex per group name as the original script); repurposed the previously-
# unused Plasma Membrane/Transport/Wound Healing slots (aqua/green/red_light)
# for the three new groups so nothing is a brand-new color.
group_pal <- c(
  "Cytoskeleton"                 = "#8C5A3C",  # brown
  "Muscle Development"           = "#5DAFA3",  # aqua
  "Development & Morphogenesis"  = "#7B5EA7",  # purple
  "ECM & Adhesion"                = "#D4A017",  # yellow
  "Collagen Trimer"              = "#C7A189",  # brown_light
  "Integrin Signaling & Binding" = "#D98C9C",  # pink
  "Cell Surface"                 = "#2471A3",  # blue
  "Cell Motility & Locomotion"   = "#A9CCE3",  # blue_light
  "Blood & Coagulation"          = "#C0392B",  # red
  "Neuronal"                     = "#4E8B6F",  # green
  "Golgi Lumen"                  = "#9E9E9E",  # gray
  "Signaling"                    = "#C4693A",  # orange
  "System Process"               = "#E6A19A"   # red_light
)

# =============================================================================
# 5. APPLY LABEL MAP — collapse to pathway groups
# =============================================================================

ora_df_final <- ora_df_combined %>%
  dplyr::filter(ID %in% names(pathway_label_map)) %>%
  mutate(
    term_label    = pathway_label_map[ID],
    pathway_group = factor(pathway_groups[term_label], levels = group_levels)
  ) %>%
  dplyr::filter(!is.na(term_label)) %>%
  arrange(pathway_group)

ora_groups <- ora_df_final %>%
  group_by(pathway_group) %>%
  summarise(
    geneID   = paste(unique(unlist(strsplit(geneID, "/"))), collapse = "/"),
    p.adjust = min(p.adjust),
    .groups  = "drop"
  ) %>%
  mutate(term_label = as.character(pathway_group))

group_order <- group_levels[group_levels %in% ora_groups$pathway_group]
ora_groups  <- ora_groups %>%
  mutate(term_label = factor(term_label, levels = group_order)) %>%
  arrange(term_label)

entrez_to_ensembl <- sig_entrez %>%
  dplyr::rename(ensembl = ENSEMBL, entrez = ENTREZID)

saveRDS(entrez_to_ensembl, file.path(OUTPUT_DIR, "entrez_to_ensembl_muscle_surv.rds"))
saveRDS(ora_df_final,      file.path(OUTPUT_DIR, "ora_go_all_muscle_surv_final.rds"))
write.csv(ora_df_final,    file.path(OUTPUT_DIR, "ora_go_all_muscle_surv_final.csv"),
          row.names = FALSE)

# =============================================================================
# 6. TRAJECTORY & SIGNIFICANCE PATTERN
# =============================================================================

SIG_COLS <- paste0("sig_", CONTRASTS)

# Trajectory — only compute if not already present
if (!"trajectory" %in% names(sig_proteins)) {
  sig_proteins <- sig_proteins %>%
    mutate(
      sign_base  = sign(Cancer),
      sign_tsvsc = sign(Residual),
      sign_ts    = sign(Training),
      trajectory = case_when(
        sign_base == sign_tsvsc & sign_base == sign_ts                               ~ "Exacerbated",
        sign_base != sign_tsvsc & sign_base != sign_ts                               ~ "Reversed",
        sign_base == sign_tsvsc & abs(Residual) < abs(Cancer) &
          sign_base != sign_ts                                                       ~ "Attenuated",
        TRUE                                                                         ~ "Unclassified"
      ),
      trajectory = factor(trajectory, levels = traj_levels)
    )
}

# Significance pattern — always recompute using display labels
sig_proteins <- sig_proteins %>%
  rowwise() %>%
  mutate(sig_pattern = {
    hits <- unname(CONTRAST_LABELS[CONTRASTS[c_across(all_of(SIG_COLS))]])
    if (length(hits) == 0) "None" else paste(hits, collapse = "\n+ ")
  }) %>%
  ungroup()

patterns <- sort(unique(sig_proteins$sig_pattern))
pattern_pal_colors <- c(
  "#D4A017",  # yellow  — single contrast hits
  "#5DAFA3",  # aqua
  "#4E8B6F",  # green
  "#7B5EA7",  # purple
  "#2471A3",  # blue
  "#C4693A",  # orange
  "#C0392B"   # red     — most contrasts significant
)
pattern_pal <- setNames(
  pattern_pal_colors[seq_len(length(patterns))],
  patterns
)

# =============================================================================
# 7. ROW ORDER & Y POSITIONS
# =============================================================================

mat <- sig_proteins %>%
  column_to_rownames("uniprot") %>%
  dplyr::select(all_of(CONTRASTS)) %>%
  as.matrix()

mat_nona <- replace(mat, is.na(mat), 0)

ordered_uniprots <- unlist(lapply(traj_levels, function(traj) {
  idx <- which(sig_proteins$trajectory == traj)
  if (length(idx) < 2) return(sig_proteins$uniprot[idx])
  sub <- mat_nona[sig_proteins$uniprot[idx], , drop = FALSE]
  hc  <- hclust(dist(sub, method = "euclidean"), method = "ward.D2")
  sig_proteins$uniprot[idx][hc$order]
}))

n_prot <- length(ordered_uniprots)
sig_proteins <- sig_proteins %>%
  mutate(y_pos = setNames(seq(n_prot, 1L), ordered_uniprots)[uniprot])

lfc_max <- max(abs(mat), na.rm = TRUE)

# =============================================================================
# 8. PROTEIN -> PATHWAY GROUP EDGES
# =============================================================================

term_order <- levels(ora_groups$term_label)
n_terms    <- length(term_order)

sankey_edges <- ora_groups %>%
  dplyr::select(term_label, geneID) %>%
  mutate(entrez = strsplit(geneID, "/")) %>%
  unnest(entrez) %>%
  left_join(entrez_to_ensembl, by = "entrez") %>%
  dplyr::filter(!is.na(ensembl))

ensembl_meta <- sig_proteins %>%
  dplyr::select(uniprot, gene, ensembl, trajectory, y_pos) %>%
  dplyr::filter(!is.na(ensembl)) %>%
  dplyr::distinct()

edges <- sankey_edges %>%
  left_join(ensembl_meta, by = "ensembl") %>%
  dplyr::filter(!is.na(trajectory)) %>%
  dplyr::distinct(gene, y_pos, trajectory, term_label)

message(sprintf("Proteins mapping to >=1 curated group: %d / %d",
                n_distinct(edges$gene), n_prot))

# =============================================================================
# 9. PATHWAY Y POSITIONS
# =============================================================================

# =============================================================================
# 9. PATHWAY Y POSITIONS + BAR SIZING
# =============================================================================

# Tune these two by eye:
#   BAR_GAP_FRAC -- fraction of each slot's half-height left as whitespace
#                   between adjacent bars. Raise for more gap, lower to
#                   make bars thicker/tighter together.
#   MARGIN_FRAC   -- top/bottom margin (as a fraction of n_prot) before the
#                   first/last bar.
BAR_GAP_FRAC <- 0.35
MARGIN_FRAC  <- 0.05

y_top_margin <- n_prot * MARGIN_FRAC
y_top        <- n_prot - y_top_margin
y_bottom     <- y_top_margin

# Distance between adjacent bar centers, now purely a function of the
# available range and n_terms -- BAR_H is then sized as a fraction of half
# that spacing, which guarantees bars can never overlap regardless of how
# many groups you have or how few significant proteins survive.
slot_height <- (y_top - y_bottom) / (n_terms - 1)
BAR_H       <- (slot_height / 2) * (1 - BAR_GAP_FRAC)

pathway_y <- setNames(seq(y_top, y_bottom, length.out = n_terms), term_order)

# =============================================================================
# 10. BUILD SHARED PLOT DATA
# =============================================================================

# Heatmap tiles
heat_df <- sig_proteins %>%
  dplyr::select(uniprot, y_pos, all_of(CONTRASTS)) %>%
  pivot_longer(all_of(CONTRASTS), names_to = "contrast", values_to = "logFC") %>%
  mutate(x = match(contrast, CONTRASTS) - 1L)

# Significance bar
sig_bar_df <- sig_proteins %>%
  dplyr::select(y_pos, sig_pattern) %>%
  mutate(x = X_SIG)

# Stacked bar data
# path_total named at source to avoid downstream rename conflicts
pathway_protein_counts <- edges %>%
  dplyr::distinct(gene, term_label) %>%
  count(term_label, name = "path_total")

max_total <- max(pathway_protein_counts$path_total)

bar_counts <- edges %>%
  dplyr::distinct(gene, trajectory, term_label) %>%
  count(term_label, trajectory, name = "n") %>%
  mutate(trajectory = factor(trajectory, levels = traj_levels)) %>%
  left_join(pathway_protein_counts, by = "term_label") %>%
  left_join(
    tibble(term_label = names(pathway_y), y_path = pathway_y),
    by = "term_label"
  ) %>%
  group_by(term_label) %>%
  arrange(trajectory, .by_group = TRUE) %>%
  mutate(
    x_end   = X_BAR + (cumsum(n) / max_total) * BAR_W,
    x_start = X_BAR + ((cumsum(n) - n) / max_total) * BAR_W
  ) %>%
  ungroup() %>%
  group_by(term_label) %>%
  mutate(x_end_max = max(x_end)) %>%
  ungroup()

# --- Individual bezier flows ---
bezier_pts <- edges %>%
  mutate(
    y_path = pathway_y[as.character(term_label)],
    x0 = X_TRAJ + 0.05, x3 = X_FLOW_R,
    x1 = X_TRAJ + 0.05 + (X_FLOW_R - X_TRAJ - 0.05) * 0.4,
    x2 = X_TRAJ + 0.05 + (X_FLOW_R - X_TRAJ - 0.05) * 0.6
  ) %>%
  rowwise() %>%
  do({
    row <- .
    t   <- seq(0, 1, length.out = 20)
    bx  <- (1-t)^3*row$x0 + 3*(1-t)^2*t*row$x1 +
      3*(1-t)*t^2*row$x2 + t^3*row$x3
    by  <- (1-t)^3*row$y_pos + 3*(1-t)^2*t*row$y_pos +
      3*(1-t)*t^2*row$y_path + t^3*row$y_path
    tibble(bx = bx, by = by,
           trajectory = row$trajectory,
           edge_id    = paste(row$gene, row$term_label, sep = "___"))
  }) %>%
  ungroup()

# --- Aggregated ribbon flows ---
# Left edge: span matches actual trajectory bar height (1 unit per protein)
# Right edge: span matches BAR_H (pathway bar height)

# --- Aggregated ribbon flows ---
traj_spans <- sig_proteins %>%
  group_by(trajectory) %>%
  summarise(
    y_left_max = max(y_pos) + 0.5,
    y_left_min = min(y_pos) - 0.5,
    traj_total = n_distinct(uniprot),
    .groups = "drop"
  )

flow_agg <- edges %>%
  dplyr::distinct(gene, trajectory, term_label) %>%
  count(trajectory, term_label, name = "n_proteins") %>%
  mutate(
    y_path     = pathway_y[as.character(term_label)],
    trajectory = factor(trajectory, levels = traj_levels),
    term_label = factor(term_label, levels = term_order)
  ) %>%
  left_join(traj_spans, by = "trajectory") %>%
  left_join(pathway_protein_counts, by = "term_label") %>%
  group_by(trajectory) %>%
  arrange(term_label, .by_group = TRUE) %>%
  mutate(
    slice_total   = sum(n_proteins),          # total flow weight, not unique proteins
    cum_l         = cumsum(n_proteins),
    y_left_top    = y_left_max - (cum_l - n_proteins) / slice_total * (y_left_max - y_left_min),
    y_left_bottom = y_left_max - cum_l        / slice_total * (y_left_max - y_left_min)
  ) %>%
  ungroup() %>%
  group_by(term_label) %>%
  arrange(trajectory, .by_group = TRUE) %>%
  mutate(
    slice_total_r  = sum(n_proteins),
    cum_r          = cumsum(n_proteins),
    y_right_top    = y_path + BAR_H - (cum_r - n_proteins) / slice_total_r * BAR_H * 2,
    y_right_bottom = y_path + BAR_H - cum_r                / slice_total_r * BAR_H * 2
  ) %>%
  ungroup()

ribbon_pts <- flow_agg %>%
  rowwise() %>%
  do({
    row <- .
    t   <- seq(0, 1, length.out = 40)
    x0  <- X_TRAJ + 0.05; x3 <- X_FLOW_R
    x1  <- x0 + (x3 - x0) * 0.4
    x2  <- x0 + (x3 - x0) * 0.6
    bx  <- (1-t)^3*x0 + 3*(1-t)^2*t*x1 + 3*(1-t)*t^2*x2 + t^3*x3
    by_top <- (1-t)^3*row$y_left_top    + 3*(1-t)^2*t*row$y_left_top    +
      3*(1-t)*t^2*row$y_right_top    + t^3*row$y_right_top
    by_bot <- (1-t)^3*row$y_left_bottom + 3*(1-t)^2*t*row$y_left_bottom +
      3*(1-t)*t^2*row$y_right_bottom + t^3*row$y_right_bottom
    tibble(
      x          = c(bx, rev(bx)),
      y          = c(by_top, rev(by_bot)),
      trajectory = row$trajectory,
      ribbon_id  = paste(row$trajectory, row$term_label, sep = "___")
    )
  }) %>%
  ungroup()

# =============================================================================
# 11. SHARED PLOT BUILDER
# =============================================================================

bar_axis_y      <- y_bottom - BAR_H
tick_len        <- n_prot * 0.015
label_y         <- bar_axis_y - tick_len - n_prot * 0.015
title_y         <- bar_axis_y - tick_len - n_prot * 0.05
bar_axis_breaks <- c(0, 50, 100, 150, 200)
bar_axis_x      <- X_BAR + BAR_W * (bar_axis_breaks / max_total)
col_labels      <- c("CRC v CTL", "CRC+RET v CTL", "CRC+RET v CRC")

# Legend grobs ------------------------------------------------------------
# All three legends (logFC, Trajectory, Significant in) are extracted as
# standalone grobs here and composed OUTSIDE the main plot's data space via
# compose_panel_with_legends() (section 11b) -- see the note above
# SIG_LEGEND_REL_WIDTH/FOOTER_LEGEND_REL_HEIGHT for why.

# get_legend() moved/changed signature across cowplot versions tied to
# ggplot2 >=3.5's guide-extraction internals; this falls back to manual
# gtable extraction if cowplot::get_legend() errors, so a version mismatch
# fails loudly with a clear NULL/empty grob rather than a cryptic error.
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

build_logfc_legend_grob <- function(scale = 1) {
  p <- ggplot(heat_df, aes(x = x, y = y_pos, fill = logFC)) +
    geom_tile() +
    scale_fill_gradient2(
      low = "#2471A3", mid = "#FFFFFF", high = "#C0392B", midpoint = 0,
      limits = c(-lfc_max, lfc_max), name = "logFC",
      guide = guide_colorbar(barwidth = unit(0.22 * scale, "cm"),
                             barheight = unit(2.2 * scale, "cm"),
                             title.position = "top", direction = "vertical")
    ) +
    theme(
      legend.text       = element_text(size = 6.5 * scale),
      legend.title      = element_text(size = 7.5 * scale),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.margin     = margin(0, 0, 0, 0)
    )
  get_legend_safe(p)
}

build_trajectory_legend_grob <- function(scale = 1) {
  p <- ggplot(sig_proteins, aes(x = 1, y = 1, fill = trajectory)) +
    geom_tile() +
    scale_fill_manual(
      values = traj_pal, name = "Trajectory",
      guide = guide_legend(ncol = 1, direction = "vertical", title.position = "top")
    ) +
    theme(
      legend.text       = element_text(size = 6.5 * scale),
      legend.title      = element_text(size = 7.5 * scale),
      legend.key.size   = unit(0.3 * scale, "cm"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.margin     = margin(0, 0, 0, 0)
    )
  get_legend_safe(p)
}

build_sig_pattern_legend_grob <- function(scale = 1) {
  p <- ggplot(sig_bar_df, aes(x = x, y = y_pos, fill = sig_pattern)) +
    geom_tile() +
    scale_fill_manual(
      values = pattern_pal, name = "Significant in",
      guide = guide_legend(ncol = 1, title.position = "top")  # single column, stacked vertically for the right-hand strip
    ) +
    theme(
      legend.text       = element_text(size = 6.5 * scale),
      legend.title      = element_text(size = 7.5 * scale),
      legend.key.size   = unit(0.3 * scale, "cm"),
      legend.background = element_blank(),
      legend.key        = element_blank(),
      legend.margin     = margin(0, 0, 0, 0)
    )
  get_legend_safe(p)
}

build_base_plot <- function(scale = 1) {
  
  ggplot() +
    
    # A: significance pattern bar
    geom_tile(data = sig_bar_df,
              aes(x = x, y = y_pos, fill = sig_pattern),
              width = 0.1, height = 1.0, colour = NA) +
    scale_fill_manual(
      values = pattern_pal, name = "Significant in",
      guide = "none"   # drawn separately via build_sig_pattern_legend_grob() + plot_grid
    ) +
    
    # B: heatmap tiles
    ggnewscale::new_scale_fill() +
    geom_tile(data = heat_df,
              aes(x = x, y = y_pos, fill = logFC),
              width = 1.0, height = 1.0) +
    scale_fill_gradient2(
      low = "#2471A3", mid = "#FFFFFF", high = "#C0392B",
      midpoint = 0, limits = c(-lfc_max, lfc_max), name = "logFC",
      guide = "none"   # drawn separately via build_logfc_legend_grob() + plot_grid
    ) +
    
    # C: trajectory colour bar
    ggnewscale::new_scale_fill() +
    geom_tile(data = sig_proteins %>% mutate(x = X_TRAJ),
              aes(x = x, y = y_pos, fill = trajectory),
              width = 0.1, height = 1.0, colour = NA) +
    scale_fill_manual(
      values = traj_pal, name = "Trajectory",
      guide = "none"   # drawn separately via build_trajectory_legend_grob() + plot_grid
    ) +
    
    # E: stacked bars
    ggnewscale::new_scale_fill() +
    geom_rect(data = bar_counts,
              aes(xmin = x_start, xmax = x_end,
                  ymin = y_path - BAR_H, ymax = y_path + BAR_H,
                  fill = trajectory),
              colour = "white", linewidth = 0.2 * scale) +
    scale_fill_manual(values = traj_pal, guide = "none") +
    
    geom_text(data = bar_counts %>% dplyr::distinct(term_label, y_path, x_end_max),
              aes(x = x_end_max + 0.15, y = y_path, label = term_label),
              hjust = 0, size = 3.0 * scale, colour = "grey15", fontface = "bold") +
    
    annotate("segment", x = X_BAR, xend = max(bar_axis_x),
             y = bar_axis_y, yend = bar_axis_y,
             colour = "grey40", linewidth = 0.4 * scale) +
    annotate("segment", x = bar_axis_x, xend = bar_axis_x,
             y = bar_axis_y, yend = bar_axis_y - tick_len,
             colour = "grey40", linewidth = 0.4 * scale) +
    annotate("text", x = bar_axis_x, y = label_y,
             label = as.character(bar_axis_breaks),
             size = 2.8 * scale, colour = "grey30", vjust = 1) +
    annotate("text", x = X_BAR + BAR_W / 2, y = title_y,
             label = "Protein count", size = 3 * scale, colour = "grey30",
             fontface = "bold", vjust = 1) +
    
    # Column labels -- angled diagonally rather than stacked with \n, so
    # they read as "PRE v CTL" etc. on a single baseline.
    annotate("text", x = X_HEAT, y = 0.3, label = col_labels,
             size = 3 * scale, colour = "grey15", fontface = "bold",
             angle = 45, hjust = 1, vjust = 1) +
    
    scale_x_continuous(breaks = X_HEAT, labels = c("", "", ""),
                       expand = expansion(add = c(0.2, 0))) +
    scale_y_continuous(expand = c(0.01, 0)) +
    coord_cartesian(clip = "off") +
    # No title/subtitle here -- fig_05_stitch.R supplies the panel title and
    # subtitle for the composite, so setting them in this panel too would
    # just duplicate them (this was rendering twice: once from this labs()
    # call, once from the stitcher). x/y left NULL since axis text is
    # already blanked below.
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = max(10 * scale, 5)) +
    theme(
      axis.text.y      = element_blank(),
      axis.ticks.y     = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.text.x      = element_blank(),
      panel.grid       = element_blank(),
      legend.position   = "none",   # all three legends composed externally via compose_panel_with_legends()
      plot.title        = element_blank(),
      plot.subtitle     = element_blank(),
      # Bottom margin has to fit the diagonal column labels' full rotated
      # extent, not just a nominal gap -- at angle = 45 with hjust = 1,
      # vjust = 1 the anchor point (y = 0.3) sits at the label's top-right
      # corner, so the whole string hangs down-left from there. coord_
      # cartesian(clip = "off") above only stops THIS panel's own grob from
      # clipping it; it does nothing for the composite once fig_05_stitch.R
      # assembles panels together, since that's a separate crop step. If
      # labels are still cut off after this, the fix has to happen in the
      # stitcher too -- e.g. give this panel's row more height, or reduce
      # its rel_height slightly, in whatever plot_grid()/patchwork call
      # places Panel D into the 380mm composite.
      plot.margin      = margin(10 * scale, 40 * scale, 30 * scale, 10 * scale)
    )
}

# =============================================================================
# 11b. LEGEND COMPOSITION -- assembles the bare main_plot (no inline legends)
# with all three legends (Significant in, logFC, Trajectory) stacked as a
# single column to its right. All three legends live outside the plot's
# data coordinate space, so they can't overlap ribbons/heatmap tiles
# regardless of n_prot.
# =============================================================================

compose_panel_with_legends <- function(main_plot, scale = 1) {
  sig_grob   <- build_sig_pattern_legend_grob(scale)
  logfc_grob <- build_logfc_legend_grob(scale)
  traj_grob  <- build_trajectory_legend_grob(scale)
  
  legend_col <- cowplot::plot_grid(
    sig_grob, logfc_grob, traj_grob,
    ncol = 1, rel_heights = LEGEND_COL_REL_HEIGHTS
  )
  
  cowplot::plot_grid(
    main_plot, legend_col, ncol = 2,
    rel_widths = c(1 - LEGEND_COL_REL_WIDTH, LEGEND_COL_REL_WIDTH)
  )
}

# =============================================================================
# 12. FIGURE BUILDING
# =============================================================================

# Full-size (standalone) plots — always built, cheap
fig_individual <- compose_panel_with_legends(
  build_base_plot(scale = 1) +
    geom_path(data = bezier_pts,
              aes(x = bx, y = by, group = edge_id, colour = trajectory),
              alpha = 0.25, linewidth = 0.3) +
    scale_colour_manual(values = traj_pal, guide = "none"),
  scale = 1
)

fig_ribbon <- compose_panel_with_legends(
  build_base_plot(scale = 1) +
    ggnewscale::new_scale_fill() +
    geom_polygon(data = ribbon_pts,
                 aes(x = x, y = y, group = ribbon_id, fill = trajectory),
                 alpha = 0.4, colour = NA) +
    scale_fill_manual(values = traj_pal, guide = "none"),
  scale = 1
)

# Panel-scaled versions — what the parent figure script should pull in
fig_individual_panel <- compose_panel_with_legends(
  build_base_plot(scale = PANEL_SCALE) +
    geom_path(data = bezier_pts,
              aes(x = bx, y = by, group = edge_id, colour = trajectory),
              alpha = 0.25, linewidth = 0.3 * PANEL_SCALE) +
    scale_colour_manual(values = traj_pal, guide = "none"),
  scale = PANEL_SCALE
)

fig_ribbon_panel <- compose_panel_with_legends(
  build_base_plot(scale = PANEL_SCALE) +
    ggnewscale::new_scale_fill() +
    geom_polygon(data = ribbon_pts,
                 aes(x = x, y = y, group = ribbon_id, fill = trajectory),
                 alpha = 0.4, colour = NA) +
    scale_fill_manual(values = traj_pal, guide = "none"),
  scale = PANEL_SCALE
)

if (!PANEL_SOURCE_MODE) {
  ggsave(file.path(OUTPUT_DIR, "heatmap_sankey_individual_flows.pdf"),
         plot = fig_individual, width = FIG_WIDTH, height = FIG_HEIGHT, device = cairo_pdf)
  ggsave(file.path(OUTPUT_DIR, "heatmap_sankey_ribbon_flows.pdf"),
         plot = fig_ribbon, width = FIG_WIDTH, height = FIG_HEIGHT, device = cairo_pdf)
  message("Saved both figures to: ", OUTPUT_DIR)
}

# =============================================================================
# 15. BUILD PLOT DATA (protein-level + ORA-level export tables)
# Computed unconditionally -- fig_05_stitch_FINAL.R pulls sheet1/sheet2 out
# of sankey_env for its own combined workbook even in PANEL_SOURCE_MODE.
# Only the standalone ggsave()/write_xlsx() calls below are gated on
# PANEL_SOURCE_MODE.
# =============================================================================

{
  # --- Sheet 1: Protein-level data ---
  protein_pathway_membership <- edges %>%
    dplyr::distinct(gene, term_label) %>%
    group_by(gene) %>%
    summarise(
      pathway_groups   = paste(sort(unique(as.character(term_label))), collapse = "; "),
      n_pathway_groups = n_distinct(term_label),
      .groups = "drop"
    )
  
  # Individual GO term membership (from ora_df_final, before collapsing)
  protein_term_membership <- ora_df_final %>%
    dplyr::select(term_label, pathway_group, geneID) %>%
    mutate(entrez = strsplit(geneID, "/")) %>%
    unnest(entrez) %>%
    left_join(entrez_to_ensembl, by = "entrez") %>%
    dplyr::filter(!is.na(ensembl)) %>%
    left_join(sig_proteins %>% dplyr::select(gene, ensembl) %>% dplyr::distinct(),
              by = "ensembl") %>%
    dplyr::filter(!is.na(gene)) %>%
    group_by(gene) %>%
    summarise(
      go_terms       = paste(sort(unique(term_label)),          collapse = "; "),
      go_term_groups = paste(sort(unique(as.character(pathway_group))), collapse = "; "),
      n_go_terms     = n_distinct(term_label),
      .groups        = "drop"
    )
  
  sheet1 <- sig_proteins %>%
    dplyr::select(
      uniprot, gene, ensembl,
      logFC_PREvCTL  = Cancer,
      logFC_POSTvCTL = Residual,
      logFC_POSTvPRE = Training,
      trajectory,
      sig_PREvCTL   = sig_Cancer,
      sig_POSTvCTL   = sig_Residual,
      sig_POSTvPRE       = sig_Training,
      sig_pattern
    ) %>%
    left_join(protein_pathway_membership, by = "gene") %>%
    left_join(protein_term_membership,    by = "gene") %>%
    mutate(
      maps_to_pathway = !is.na(pathway_groups),
      trajectory      = as.character(trajectory)
    ) %>%
    arrange(desc(maps_to_pathway), trajectory, gene)
  
  # Build Entrez → gene symbol lookup
  entrez_to_gene <- sig_entrez %>%
    left_join(sig_proteins %>% dplyr::select(ensembl, gene) %>% dplyr::distinct(),
              by = c("ENSEMBL" = "ensembl")) %>%
    dplyr::filter(!is.na(gene)) %>%
    dplyr::select(entrez = ENTREZID, gene) %>%
    dplyr::distinct()
  
  entrez_to_symbol <- setNames(entrez_to_gene$gene, entrez_to_gene$entrez)
  
  # Convert geneID column in ora_df_final before building sheet2
  sheet2 <- ora_df_final %>%
    mutate(
      gene_symbols = sapply(strsplit(geneID, "/"), function(ids) {
        syms <- entrez_to_symbol[ids]
        syms <- syms[!is.na(syms)]
        paste(sort(syms), collapse = "; ")
      })
    ) %>%
    dplyr::select(
      go_id         = ID,
      go_term       = term_label,
      pathway_group,
      ontology,
      gene_ratio    = GeneRatio,
      bg_ratio      = BgRatio,
      pvalue,
      p_adjust      = p.adjust,
      q_value       = qvalue,
      gene_count    = Count,
      n_proteins,
      gene_symbols
    ) %>%
    arrange(pathway_group, p_adjust)
}

if (!PANEL_SOURCE_MODE) {
  pacman::p_load(writexl)

  writexl::write_xlsx(
    list(
      "Proteins"    = sheet1 %>% mutate(across(where(is.factor), as.character)),
      "ORA_Results" = sheet2 %>% mutate(across(where(is.factor), as.character))
    ),
    path = file.path(OUTPUT_DIR, "heatmap_sankey_data.xlsx")
  )

  message("Data saved to: ", file.path(OUTPUT_DIR, "heatmap_sankey_data.xlsx"))

}