# =============================================================================
# Heatmap-Sankey — GSEA edition (keyword-selected pathways, muscle + plasma)
#
# Protein -> pathway flows only, colored by PATHWAY (not trajectory).
# Trajectory classification/bar/legend and the aggregated-ribbon version are
# removed -- ribbon aggregation only made geometric sense with trajectory
# giving contiguous groups on the left edge; without it, individual
# protein->pathway lines are the natural representation.
#
# MODIFIED from the hand-picked-ID version: pathways are now selected by
# keyword match against Description_clean (three themes: mitochondrial/
# respiration/ATP, immune, structural/ECM), rather than a manually curated
# ID list. HAND_PICKED_PATHWAYS is left in place below for reference/manual
# override, but is no longer what drives the default run at the bottom.
#
# Row order falls back to a single hierarchical clustering across all
# significant proteins' logFC vectors (no trajectory pre-grouping).
#
# Reuses the plotting machinery structure from the ORA-based heatmap_sankey
# script (build_base_plot, legend composition), adapted throughout.
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
  SURV_limma_muscle = "03_analysis/c_data/muscle/limma_results_all_muscle.xlsx",
  SURV_limma_plasma = "03_analysis/c_data/plasma/limma_results_all_plasma.xlsx",
  gsea_data         = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data.rds",
  figure_dir        = "04_figures/fig_05/heatmap_sankey_gsea"
)

OUTPUT_DIR <- config$figure_dir
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

organism  <- org.Hs.eg.db
PI_THRESH <- 0.05

CONTRASTS <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")

CONTRAST_LABELS <- c(
  Baseline_SURVvCTL = "CRC vs CTL",
  Training_SURVvCTL = "CRC+RET vs CTL",
  Training_SURV     = "CRC+RET vs CRC"
)

# --- HAND-PICKED PATHWAYS (legacy / manual override) ------------------------
# Left in place for reference or one-off manual runs. Not used by the default
# keyword-driven run at the bottom of this script -- see PATHWAY_THEMES and
# select_pathways_by_keyword() below instead.
HAND_PICKED_PATHWAYS <- tibble::tribble(
  ~id,          ~label_override,
  "GO:0005743", NA_character_,
  "GO:0031966", NA_character_,
  "GO:0005740", NA_character_,
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION", NA_character_,
  "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT", NA_character_,
  "GO:0005759", NA_character_,
  "GO:0098803", NA_character_,
  "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION", NA_character_,
  "GO:0015078", NA_character_,
  "GO:1902600", NA_character_,
  "REACTOME_CITRIC_ACID_CYCLE_TCA_CYCLE", NA_character_,
  "GO:0015986", NA_character_,
  "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT", NA_character_,
  "GO:0042776", NA_character_,
  "GO:0006754", NA_character_,
  "GO:0140053", NA_character_,
  "GO:0098798", NA_character_,
  "GO:0045333", NA_character_,
  "GO:0006119", NA_character_,
  "GO:0007005", NA_character_,
  "hsa00190", NA_character_,
  "GO:0032543", NA_character_,
  "GO:0045271", NA_character_,
  "GO:0005741", NA_character_,
  "hsa00020", NA_character_,
  "REACTOME_COMPLEX_I_BIOGENESIS", NA_character_,
  "REACTOME_MITOCHONDRIAL_BIOGENESIS", NA_character_,
  "GO:0033108", NA_character_,
  "GO:0042773", NA_character_,
  "GO:0042775", NA_character_,
  "GO:0006099", NA_character_,
  "REACTOME_MITOCHONDRIAL_TRANSLATION", NA_character_
)

# --- KEYWORD-BASED PATHWAY THEMES --------------------------------------------
# Edit these lists to tune what gets pulled in per theme. Matching is
# case-insensitive substring search against Description_clean (see
# select_pathways_by_keyword() in section 1b) -- keep terms short/generic
# enough to catch variants ("mitochond" catches "mitochondrial",
# "mitochondrion", etc.) but specific enough not to bleed into other themes.
PATHWAY_THEMES <- list(
  mitochondrial = c(
    "mitochond", "respiration", "respiratory electron", "oxidative phosphorylation",
    "electron transport", "atp synthesis", "atp metabolic", "atp biosynthetic",
    "tca cycle", "citric acid cycle", "oxphos", "cellular respiration", "nadh"
  ),
  immune = c(
    "immune", "inflammat", "cytokine", "leukocyte", "lymphocyte", "neutrophil",
    "interferon", "interleukin", "complement activation", "antigen", "chemokine",
    "innate immune", "adaptive immune", "t cell", "b cell", "macrophage",
    "phagocyt", "degranulation"
  ),
  structural_ecm = c(
    "extracellular matrix", "collagen", "cytoskeleton", "cell adhesion",
    "basement membrane", "cell junction", "actin", "focal adhesion",
    "structural constituent", "matrix organization", "cell-matrix"
  )
)

# Figure layout constants
X_GENE   <- -0.7   # gene name labels -- leftmost, outside the significance bar
X_SIG    <- -0.6   # significance bar, between gene labels and the heatmap
X_HEAT   <- 0:2
X_FLOW_L <- 2.6   # flows now originate right at the heatmap's edge
X_FLOW_R <- 5.5
X_BAR    <- 5.5
BAR_W    <- 4.0
# BAR_H is no longer a fixed constant -- computed per tissue in
# run_heatmap_sankey() based on n_terms/n_prot so pathway bars don't
# overlap regardless of how many pathways are hand-picked.

FIG_WIDTH  <- 12
FIG_HEIGHT <- 8

PANEL_SOURCE_MODE <- exists("PANEL_SOURCE_MODE") && isTRUE(PANEL_SOURCE_MODE)
PANEL_WIDTH_MM <- if (exists("PANEL_WIDTH_MM")) PANEL_WIDTH_MM else 178
PANEL_SCALE    <- if (exists("scale_text")) scale_text(1, PANEL_WIDTH_MM) else (PANEL_WIDTH_MM / 25.4) / FIG_WIDTH
PANEL_WIDTH    <- PANEL_WIDTH_MM / 25.4
PANEL_HEIGHT   <- PANEL_WIDTH * (FIG_HEIGHT / FIG_WIDTH)

LEGEND_COL_REL_WIDTH   <- 0.16
LEGEND_COL_REL_HEIGHTS <- c(0.55, 0.45)  # sig_pattern, logFC

# =============================================================================
# 1. LOAD GSEA + LIMMA DATA
# =============================================================================

gsea_pathway_data <- readRDS(config$gsea_data)
gsea_supp_full    <- gsea_pathway_data$gsea_supp_full %>% filter(group == "SURV")

# =============================================================================
# 1a. BROWSE SIGNIFICANT PATHWAYS
# Still useful for spot-checking individual terms or eyeballing what's out
# there before trusting the keyword match below.
# =============================================================================

pathway_browser <- gsea_supp_full %>%
  dplyr::filter(significant) %>%
  dplyr::select(tissue, comparison, database, ID, Description_clean, NES, p.adjust, setSize) %>%
  mutate(comparison = factor(comparison, levels = CONTRASTS)) %>%
  arrange(tissue, comparison, p.adjust)

write.csv(pathway_browser, file.path(OUTPUT_DIR, "significant_pathways_browser.csv"), row.names = FALSE)
message(sprintf("Wrote %d significant pathway rows to: %s",
                nrow(pathway_browser), file.path(OUTPUT_DIR, "significant_pathways_browser.csv")))

if (interactive()) View(pathway_browser)

#' Quick filtered browser -- narrow by tissue/comparison/database, or search
#' Description_clean for a keyword, then read off the ID column for
#' HAND_PICKED_PATHWAYS. Examples:
#'   browse_pathways(tissue_filter = "Muscle", comparison_filter = "Training_SURV")
#'   browse_pathways(search = "collagen")
#'   browse_pathways(database_filter = "KEGG", n = 50)
browse_pathways <- function(tissue_filter = NULL, comparison_filter = NULL,
                            database_filter = NULL, search = NULL, n = 30) {
  df <- pathway_browser
  if (!is.null(tissue_filter))     df <- df %>% dplyr::filter(tissue == tissue_filter)
  if (!is.null(comparison_filter)) df <- df %>% dplyr::filter(comparison == comparison_filter)
  if (!is.null(database_filter))   df <- df %>% dplyr::filter(database == database_filter)
  if (!is.null(search))            df <- df %>% dplyr::filter(str_detect(tolower(Description_clean), tolower(search)))
  df %>% arrange(p.adjust) %>% dplyr::slice_head(n = n) %>% print(n = n)
}

# =============================================================================
# 1b. KEYWORD-BASED PATHWAY SELECTION
#
# Pulls every DISTINCT significant pathway (pooled across tissue x contrast,
# since a term only needs to be significant SOMEWHERE to be worth including
# in the figure -- build_pathway_gene_pool() re-checks significance per
# tissue later anyway) whose Description_clean matches any keyword in a
# theme's list. Returns a tibble in the same {id, label_override} shape as
# HAND_PICKED_PATHWAYS, so it's a drop-in replacement for pathway_spec.
# =============================================================================

select_pathways_by_keyword <- function(keywords, exclude = NULL, tissues = NULL) {
  df <- gsea_supp_full %>% dplyr::filter(significant)
  if (!is.null(tissues)) df <- df %>% dplyr::filter(tissue %in% tissues)
  
  pattern <- paste(tolower(keywords), collapse = "|")
  df <- df %>% dplyr::filter(str_detect(tolower(Description_clean), pattern))
  
  if (!is.null(exclude)) {
    excl_pattern <- paste(tolower(exclude), collapse = "|")
    df <- df %>% dplyr::filter(!str_detect(tolower(Description_clean), excl_pattern))
  }
  
  df %>%
    dplyr::distinct(ID, Description_clean) %>%
    dplyr::arrange(Description_clean) %>%
    transmute(id = ID, label_override = NA_character_)
}

pathway_specs <- list(
  mitochondrial  = select_pathways_by_keyword(PATHWAY_THEMES$mitochondrial),
  immune         = select_pathways_by_keyword(PATHWAY_THEMES$immune),
  structural_ecm = select_pathways_by_keyword(PATHWAY_THEMES$structural_ecm)
)

# Preview counts + write each theme's matched pathway list to CSV before
# spending time on the full figure build -- same "check before you commit"
# pattern as the pathway_browser export above. Skim these; if a theme pulled
# in something off-topic, add it to that theme's `exclude` (rebuild the
# corresponding pathway_specs[[...]] entry with the exclude= argument) or
# tighten/loosen the keyword list and re-run this section.
iwalk(pathway_specs, function(spec, theme_name) {
  message(sprintf("Theme '%s': %d pathways matched", theme_name, nrow(spec)))
  preview <- gsea_supp_full %>%
    dplyr::filter(ID %in% spec$id) %>%
    dplyr::distinct(ID, Description_clean, database) %>%
    dplyr::arrange(Description_clean)
  write.csv(preview, file.path(OUTPUT_DIR, sprintf("theme_%s_matched_pathways.csv", theme_name)),
            row.names = FALSE)
})

CONTRAST_SHEETS <- list(
  Baseline_SURVvCTL = 11,
  Training_SURVvCTL = 12,
  Training_SURV     = 13
)

#' Build the full "significant in any of the 3 SURV contrasts" protein
#' table for one tissue. Handles two workbook schemas: muscle carries an
#' `Ensembl` column directly; plasma doesn't (biomaRt was removed from that
#' pipeline), so Ensembl IDs are derived via org.Hs.eg.db from `uniprot_id`,
#' which both workbooks do have.
load_sig_proteins <- function(limma_path) {
  sheets <- map(CONTRAST_SHEETS, ~ read_excel(limma_path, sheet = .x))
  
  has_ensembl <- "Ensembl" %in% names(sheets[[1]])
  
  if (!has_ensembl) {
    all_uniprot_ids <- unique(sheets[[1]]$uniprot_id)
    ens_map <- AnnotationDbi::select(
      organism, keys = all_uniprot_ids, keytype = "UNIPROT", columns = "ENSEMBL"
    ) %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::distinct(UNIPROT, .keep_all = TRUE) %>%   # collapse 1-to-many UniProt->Ensembl
      dplyr::rename(uniprot_id = UNIPROT, Ensembl = ENSEMBL)
    
    message(sprintf(
      "No Ensembl column in workbook -- derived via org.Hs.eg.db for %d / %d proteins",
      nrow(ens_map), length(all_uniprot_ids)
    ))
    
    sheets <- map(sheets, ~ left_join(.x, ens_map, by = "uniprot_id"))
  }
  
  da_tbl <- imap(sheets, function(df, contrast) {
    df %>% as_tibble() %>%
      transmute(
        uniprot = uniprot, gene = gene, ensembl = Ensembl,
        logFC = logFC, pi_val = Pi.Val,
        sig = Pi.Val < PI_THRESH, contrast = contrast
      )
  }) %>% bind_rows()
  
  SIG_COLS <- paste0("sig_", CONTRASTS)
  
  sp <- da_tbl %>%
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
    filter(if_any(all_of(SIG_COLS), ~ .x == TRUE)) %>%
    rowwise() %>%
    mutate(sig_pattern = {
      hits <- unname(CONTRAST_LABELS[CONTRASTS[c_across(all_of(SIG_COLS))]])
      if (length(hits) == 0) "None" else paste(hits, collapse = "\n+ ")
    }) %>%
    ungroup()
  
  message(sprintf("Loaded %d significant proteins from %s", nrow(sp), basename(limma_path)))
  sp
}

sig_proteins_by_tissue <- list(
  Muscle = load_sig_proteins(config$SURV_limma_muscle),
  Plasma = load_sig_proteins(config$SURV_limma_plasma)
)

patterns <- sort(unique(unlist(lapply(sig_proteins_by_tissue, function(x) x$sig_pattern))))
pattern_pal_colors <- c(
  "#D4A017", "#5DAFA3", "#4E8B6F", "#7B5EA7", "#2471A3", "#C4693A", "#C0392B"
)
pattern_pal <- setNames(pattern_pal_colors[seq_len(length(patterns))], patterns)

# =============================================================================
# 2. HARMONIZE PATHWAY MEMBERSHIP TO ENSEMBL
# =============================================================================
# core_enrichment for GO/Reactome/Hallmark is already ENSG-format; KEGG
# uses Entrez per your existing ID-harmonization convention and needs bitr().

harmonize_gene_ids <- function(id_vec, database) {
  if (database == "KEGG") {
    map_df <- bitr(id_vec, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = organism)
    unique(na.omit(map_df$ENSEMBL))
  } else {
    unique(id_vec)
  }
}

#' Pool core_enrichment ENSEMBL IDs for one pathway across every SURV
#' comparison where it was significant, for one tissue.
build_pathway_gene_pool <- function(pathway_id, tissue_name) {
  rows <- gsea_supp_full %>%
    filter(ID == pathway_id, tissue == tissue_name, significant)
  
  if (nrow(rows) == 0) {
    message(sprintf("[%s] %s: not significant in any SURV contrast -- skipping.",
                    tissue_name, pathway_id))
    return(NULL)
  }
  
  db          <- as.character(rows$database[1])
  raw_ids     <- unlist(strsplit(rows$core_enrichment, "/"))
  ensembl_ids <- harmonize_gene_ids(raw_ids, db)
  
  list(
    id               = pathway_id,
    label            = rows$Description_clean[1],
    database         = db,
    ensembl          = ensembl_ids,
    comparisons_sig  = unique(rows$comparison),
    nes_by_contrast  = setNames(rows$NES, rows$comparison),
    padj_by_contrast = setNames(rows$p.adjust, rows$comparison)
  )
}

# =============================================================================
# 3. BUILD PER-TISSUE SANKEY INPUTS (edges, pathway groups)
# =============================================================================

build_sankey_edges <- function(tissue_name, sig_proteins, pathway_spec = HAND_PICKED_PATHWAYS) {
  
  pools <- pmap(pathway_spec, function(id, label_override) {
    build_pathway_gene_pool(id, tissue_name)
  })
  pools <- compact(pools)
  
  if (length(pools) == 0) {
    warning(sprintf("[%s] None of the selected pathways were significant here.", tissue_name))
    return(NULL)
  }
  
  term_order <- map_chr(seq_along(pools), function(i) {
    override <- pathway_spec$label_override[pathway_spec$id == pools[[i]]$id]
    if (length(override) > 0 && !is.na(override)) override else pools[[i]]$label
  })
  
  # Disambiguate duplicate labels -- two different pathway IDs (often
  # across databases, e.g. a GO term and a Hallmark set) can clean down to
  # the same display string. factor() requires unique levels, and without
  # this, two distinct gene sets would silently merge into one bar/color.
  dup <- duplicated(term_order) | duplicated(term_order, fromLast = TRUE)
  if (any(dup)) {
    dup_ids <- map_chr(pools, "id")[dup]
    message("Duplicate pathway labels detected, disambiguating with ID suffix: ",
            paste(unique(term_order[dup]), collapse = ", "))
    term_order[dup] <- paste0(term_order[dup], " (", dup_ids, ")")
  }
  
  ora_groups <- tibble(
    term_label = factor(term_order, levels = term_order),
    geneID     = map_chr(pools, ~ paste(.x$ensembl, collapse = "/"))
  )
  
  edges <- ora_groups %>%
    dplyr::select(term_label, geneID) %>%
    mutate(ensembl = strsplit(geneID, "/")) %>%
    unnest(ensembl) %>%
    inner_join(
      sig_proteins %>% dplyr::select(gene, ensembl),
      by = "ensembl"
    ) %>%
    dplyr::distinct(gene, term_label)
  
  message(sprintf("[%s] %d / %d significant proteins map to >=1 selected pathway",
                  tissue_name, n_distinct(edges$gene), nrow(sig_proteins)))
  
  list(ora_groups = ora_groups, edges = edges, term_order = term_order,
       n_terms = length(term_order), pools = pools)
}

# =============================================================================
# 4-15. PER-TISSUE BUILD + PLOT
# =============================================================================

# theme_name is now threaded through purely for output pathing/filenames --
# it doesn't change any plotting logic, just keeps the three themes' figures
# and data exports from overwriting each other on disk.
run_heatmap_sankey <- function(tissue_name, pathway_spec = HAND_PICKED_PATHWAYS, theme_name = "custom") {
  
  sig_proteins <- sig_proteins_by_tissue[[tissue_name]]
  sankey_in    <- build_sankey_edges(tissue_name, sig_proteins, pathway_spec)
  if (is.null(sankey_in)) return(invisible(NULL))
  
  edges      <- sankey_in$edges
  term_order <- sankey_in$term_order
  n_terms    <- sankey_in$n_terms
  
  # Restrict to proteins that are BOTH significant AND map to >=1 selected
  # pathway -- edges was already built by inner-joining sig_proteins against
  # pathway membership, so edges$gene is exactly that set. Everything below
  # (row order, heatmap, sig bar) is now scoped to this subset only, rather
  # than showing every significant protein regardless of pathway membership.
  mapped_genes <- unique(edges$gene)
  sig_proteins <- sig_proteins %>% dplyr::filter(gene %in% mapped_genes)
  message(sprintf("[%s/%s] Figure scoped to %d proteins mapping to the selected pathways",
                  theme_name, tissue_name, nrow(sig_proteins)))
  
  # Pathway color palette -- generated fresh each run since the pathway list
  # can change; swap in a curated palette here if you want fixed colors.
  pathway_pal <- setNames(scales::hue_pal()(n_terms), term_order)
  
  tissue_dir <- file.path(OUTPUT_DIR, theme_name, tolower(tissue_name))
  dir.create(tissue_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Row order: single hierarchical clustering across ALL significant
  # proteins' logFC vectors (no trajectory pre-grouping) -----------------
  mat <- sig_proteins %>%
    column_to_rownames("uniprot") %>%
    dplyr::select(all_of(CONTRASTS)) %>%
    as.matrix()
  mat_nona <- replace(mat, is.na(mat), 0)
  
  hc <- hclust(dist(mat_nona, method = "euclidean"), method = "ward.D2")
  ordered_uniprots <- sig_proteins$uniprot[hc$order]
  
  n_prot <- length(ordered_uniprots)
  sig_proteins <- sig_proteins %>%
    mutate(y_pos = setNames(seq(n_prot, 1L), ordered_uniprots)[uniprot])
  
  edges <- edges %>%
    left_join(sig_proteins %>% dplyr::select(uniprot, gene, y_pos),
              by = "gene", relationship = "many-to-many") %>%
    dplyr::filter(!is.na(y_pos)) %>%
    dplyr::distinct(gene, y_pos, term_label)
  
  lfc_max <- max(abs(mat), na.rm = TRUE)
  
  # --- Pathway y positions ---------------------------------------------
  # Bar height/padding scale with n_terms and n_prot rather than the fixed
  # BAR_H constant -- a static height was fine for a small curated list
  # (the ORA-groups version topped out around a dozen redundancy-filtered
  # groups) but visibly overlaps once you have dozens of keyword-matched
  # pathways, especially ones with heavy pairwise gene overlap (e.g. several
  # mitochondrial/OXPHOS terms sharing most of the same core genes).
  PATHWAY_BAR_GAP_FRAC <- 0.12   # gap between adjacent bars, as a fraction of slot height
  top_bottom_pad <- max(n_prot * 0.03, 3)
  y_top    <- n_prot - top_bottom_pad
  y_bottom <- top_bottom_pad
  slot_height <- if (n_terms > 1) (y_top - y_bottom) / (n_terms - 1) else (y_top - y_bottom)
  BAR_H <- (slot_height * (1 - PATHWAY_BAR_GAP_FRAC)) / 2
  pathway_y <- if (n_terms > 1) {
    setNames(seq(y_top, y_bottom, length.out = n_terms), term_order)
  } else {
    setNames((y_top + y_bottom) / 2, term_order)
  }
  
  # --- Shared plot data ------------------------------------------------------
  heat_df <- sig_proteins %>%
    dplyr::select(uniprot, y_pos, all_of(CONTRASTS)) %>%
    pivot_longer(all_of(CONTRASTS), names_to = "contrast", values_to = "logFC") %>%
    mutate(x = match(contrast, CONTRASTS) - 1L)
  
  sig_bar_df <- sig_proteins %>%
    dplyr::select(y_pos, sig_pattern) %>%
    mutate(x = X_SIG)
  
  # One bar per pathway, single-colored by pathway (no trajectory split)
  bar_counts <- edges %>%
    dplyr::distinct(gene, term_label) %>%
    count(term_label, name = "n") %>%
    left_join(tibble(term_label = names(pathway_y), y_path = pathway_y), by = "term_label") %>%
    mutate(
      x_start = X_BAR,
      x_end   = X_BAR + (n / max(n)) * BAR_W
    )
  
  # --- Individual bezier flows, colored by PATHWAY -------------------------
  bezier_pts <- edges %>%
    mutate(
      y_path = pathway_y[as.character(term_label)],
      x0 = X_FLOW_L + 0.05, x3 = X_FLOW_R,
      x1 = X_FLOW_L + 0.05 + (X_FLOW_R - X_FLOW_L - 0.05) * 0.4,
      x2 = X_FLOW_L + 0.05 + (X_FLOW_R - X_FLOW_L - 0.05) * 0.6
    ) %>%
    rowwise() %>%
    do({
      row <- .
      t   <- seq(0, 1, length.out = 20)
      bx  <- (1-t)^3*row$x0 + 3*(1-t)^2*t*row$x1 + 3*(1-t)*t^2*row$x2 + t^3*row$x3
      by  <- (1-t)^3*row$y_pos + 3*(1-t)^2*t*row$y_pos + 3*(1-t)*t^2*row$y_path + t^3*row$y_path
      tibble(bx = bx, by = by, term_label = row$term_label,
             edge_id = paste(row$gene, row$term_label, sep = "___"))
    }) %>%
    ungroup() %>%
    mutate(term_label = factor(term_label, levels = term_order))
  
  # --- Plot builder ------------------------------------------------------
  bar_axis_y      <- y_bottom - BAR_H
  tick_len        <- n_prot * 0.015
  label_y         <- bar_axis_y - tick_len - n_prot * 0.015
  title_y         <- bar_axis_y - tick_len - n_prot * 0.05
  bar_axis_breaks <- pretty(c(0, max(bar_counts$n)))
  bar_axis_x      <- X_BAR + BAR_W * (bar_axis_breaks / max(bar_counts$n))
  col_labels      <- c("CRC v CTL", "CRC+RET v CTL", "CRC+RET v CRC")
  
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
      theme(legend.text = element_text(size = 6.5 * scale),
            legend.title = element_text(size = 7.5 * scale),
            legend.background = element_blank(), legend.key = element_blank(),
            legend.margin = margin(0, 0, 0, 0))
    get_legend_safe(p)
  }
  
  build_sig_pattern_legend_grob <- function(scale = 1) {
    p <- ggplot(sig_bar_df, aes(x = x, y = y_pos, fill = sig_pattern)) +
      geom_tile() +
      scale_fill_manual(values = pattern_pal, name = "Significant in",
                        guide = guide_legend(ncol = 1, title.position = "top")) +
      theme(legend.text = element_text(size = 6.5 * scale),
            legend.title = element_text(size = 7.5 * scale),
            legend.key.size = unit(0.3 * scale, "cm"),
            legend.background = element_blank(), legend.key = element_blank(),
            legend.margin = margin(0, 0, 0, 0))
    get_legend_safe(p)
  }
  
  build_base_plot <- function(scale = 1) {
    ggplot() +
      geom_tile(data = sig_bar_df, aes(x = x, y = y_pos, fill = sig_pattern),
                width = 0.1, height = 1.0, colour = NA) +
      scale_fill_manual(values = pattern_pal, name = "Significant in", guide = "none") +
      
      ggnewscale::new_scale_fill() +
      geom_tile(data = heat_df, aes(x = x, y = y_pos, fill = logFC), width = 1.0, height = 1.0) +
      scale_fill_gradient2(low = "#2471A3", mid = "#FFFFFF", high = "#C0392B",
                           midpoint = 0, limits = c(-lfc_max, lfc_max), name = "logFC", guide = "none") +
      
      geom_text(data = sig_proteins, aes(x = X_GENE, y = y_pos, label = gene),
                hjust = 1, size = 2.6 * scale, fontface = "italic", colour = "grey15") +
      
      ggnewscale::new_scale_fill() +
      geom_rect(data = bar_counts,
                aes(xmin = x_start, xmax = x_end, ymin = y_path - BAR_H, ymax = y_path + BAR_H, fill = term_label),
                colour = "white", linewidth = 0.2 * scale) +
      scale_fill_manual(values = pathway_pal, guide = "none") +
      
      geom_text(data = bar_counts, aes(x = x_end + 0.15, y = y_path, label = term_label),
                hjust = 0, size = 3.0 * scale, colour = "grey15", fontface = "bold") +
      
      annotate("segment", x = X_BAR, xend = max(bar_axis_x), y = bar_axis_y, yend = bar_axis_y,
               colour = "grey40", linewidth = 0.4 * scale) +
      annotate("segment", x = bar_axis_x, xend = bar_axis_x, y = bar_axis_y, yend = bar_axis_y - tick_len,
               colour = "grey40", linewidth = 0.4 * scale) +
      annotate("text", x = bar_axis_x, y = label_y, label = as.character(bar_axis_breaks),
               size = 2.8 * scale, colour = "grey30", vjust = 1) +
      annotate("text", x = X_BAR + BAR_W / 2, y = title_y, label = "Protein count",
               size = 3 * scale, colour = "grey30", fontface = "bold", vjust = 1) +
      
      annotate("text", x = X_HEAT, y = 0.3, label = col_labels,
               size = 3 * scale, colour = "grey15", fontface = "bold",
               angle = 45, hjust = 1, vjust = 1) +
      
      scale_x_continuous(breaks = X_HEAT, labels = c("", "", ""), expand = expansion(add = c(0.2, 0))) +
      scale_y_continuous(expand = c(0.01, 0)) +
      coord_cartesian(clip = "off") +
      labs(x = NULL, y = NULL) +
      theme_minimal(base_size = max(10 * scale, 5)) +
      theme(
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        panel.grid = element_blank(), legend.position = "none",
        plot.title = element_blank(), plot.subtitle = element_blank(),
        plot.margin = margin(10 * scale, 40 * scale, 30 * scale, 30 * scale)
      )
  }
  
  compose_panel_with_legends <- function(main_plot, scale = 1) {
    sig_grob   <- build_sig_pattern_legend_grob(scale)
    logfc_grob <- build_logfc_legend_grob(scale)
    legend_col <- cowplot::plot_grid(sig_grob, logfc_grob,
                                     ncol = 1, rel_heights = LEGEND_COL_REL_HEIGHTS)
    cowplot::plot_grid(main_plot, legend_col, ncol = 2,
                       rel_widths = c(1 - LEGEND_COL_REL_WIDTH, LEGEND_COL_REL_WIDTH))
  }
  
  # --- Figure building -------------------------------------------------------
  fig_sankey <- compose_panel_with_legends(
    build_base_plot(scale = 1) +
      geom_path(data = bezier_pts, aes(x = bx, y = by, group = edge_id, colour = term_label),
                alpha = 0.3, linewidth = 0.3) +
      scale_colour_manual(values = pathway_pal, guide = "none"),
    scale = 1
  )
  
  if (!PANEL_SOURCE_MODE) {
    ggsave(file.path(tissue_dir, sprintf("heatmap_sankey_gsea_%s_%s.pdf", theme_name, tolower(tissue_name))),
           plot = fig_sankey, width = FIG_WIDTH, height = FIG_HEIGHT, device = cairo_pdf)
    message(sprintf("[%s/%s] Saved figure to: %s", theme_name, tissue_name, tissue_dir))
  }
  
  # --- Data export -------------------------------------------------------
  if (!PANEL_SOURCE_MODE) {
    pacman::p_load(writexl)
    
    protein_pathway_membership <- edges %>%
      dplyr::distinct(gene, term_label) %>%
      group_by(gene) %>%
      summarise(
        pathways   = paste(sort(unique(as.character(term_label))), collapse = "; "),
        n_pathways = n_distinct(term_label),
        .groups = "drop"
      )
    
    sheet1 <- sig_proteins %>%
      dplyr::select(
        uniprot, gene, ensembl,
        logFC_PREvCTL  = Baseline_SURVvCTL,
        logFC_POSTvCTL = Training_SURVvCTL,
        logFC_POSTvPRE = Training_SURV,
        sig_PREvCTL  = sig_Baseline_SURVvCTL,
        sig_POSTvCTL = sig_Training_SURVvCTL,
        sig_POSTvPRE = sig_Training_SURV,
        sig_pattern
      ) %>%
      left_join(protein_pathway_membership, by = "gene") %>%
      mutate(maps_to_pathway = !is.na(pathways)) %>%
      arrange(desc(maps_to_pathway), gene)
    
    sheet2 <- bind_rows(lapply(sankey_in$pools, function(p) {
      tibble(
        pathway_id              = p$id,
        pathway_label           = p$label,
        database                = p$database,
        comparisons_significant = paste(unname(CONTRAST_LABELS[p$comparisons_sig]), collapse = "; "),
        n_ensembl_pooled        = length(p$ensembl),
        n_proteins_also_sig     = n_distinct(edges$gene[edges$term_label == p$label])
      )
    }))
    
    writexl::write_xlsx(
      list(
        "Proteins" = sheet1 %>% mutate(across(where(is.factor), as.character)),
        "Pathways" = sheet2
      ),
      path = file.path(tissue_dir, sprintf("heatmap_sankey_gsea_%s_%s_data.xlsx", theme_name, tolower(tissue_name)))
    )
    message(sprintf("[%s/%s] Data saved to: %s", theme_name, tissue_name,
                    file.path(tissue_dir, sprintf("heatmap_sankey_gsea_%s_%s_data.xlsx", theme_name, tolower(tissue_name)))))
  }
  
  invisible(list(fig_sankey = fig_sankey, sig_proteins = sig_proteins, edges = edges))
}

# =============================================================================
# 16. RUN FOR EACH THEME x TISSUE
# =============================================================================
# Produces 6 figures total: {mitochondrial, immune, structural_ecm} x
# {Muscle, Plasma}. Each theme's pathway_spec was built in section 1b from
# PATHWAY_THEMES -- edit the keyword lists there and re-run from section 1b
# down if a theme needs tuning.

sankey_results <- list()
for (theme_name in names(pathway_specs)) {
  spec <- pathway_specs[[theme_name]]
  sankey_results[[theme_name]] <- list(
    Muscle = run_heatmap_sankey("Muscle", pathway_spec = spec, theme_name = theme_name),
    Plasma = run_heatmap_sankey("Plasma", pathway_spec = spec, theme_name = theme_name)
  )
}
