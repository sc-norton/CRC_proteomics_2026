# =============================================================================
# Pathway-group bar chart -- mean logFC by contrast, with t-test significance
#
# Reuses the keyword-based pathway theme selection + protein/pathway
# harmonization machinery from the GSEA heatmap-sankey script
# (select_pathways_by_keyword(), build_pathway_gene_pool(),
# load_sig_proteins(), harmonize_gene_ids() -- copied near-verbatim below),
# but replaces the heatmap+sankey visualization entirely with a grouped bar
# chart: one bar pair (CRC vs CTL, CRC+RET vs CRC) per pathway theme
# (Metabolic / Immune / Structural), mean logFC across every significant
# protein mapping to that theme, error bars = SEM, and a one-sample t-test
# (H0: mean logFC = 0) per bar for significance asterisks.
#
# "Pathway group" here means a THEME (Metabolic/Immune/Structural), not an
# individual GO/KEGG/Hallmark term -- analogous to how the reference figure
# groups "Init factors" / "Elong factors" as x-axis categories, not
# individual genes or individual pathway IDs.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  readxl, dplyr, tidyr, tibble, stringr, purrr, scales,
  ggplot2, writexl,
  clusterProfiler, org.Hs.eg.db, AnnotationDbi
)

source("04_figures/a_scripts/00_theme.R")   # colors, theme_cr()

# =============================================================================
# 0. CONFIGURATION
# =============================================================================

config <- list(
  SURV_limma_muscle = "03_analysis/c_data/muscle/limma_results_all_muscle_cv.xlsx",
  SURV_limma_plasma = "03_analysis/c_data/plasma/limma_results_all_plasma.xlsx",
  gsea_data         = "03_analysis/c_data/gsea_consolidated/gsea_pathway_data.rds",
  fig_dir           = "04_figures/fig_pathway_bars",
  data_dir          = "04_figures/data/fig_pathway_bars"
)
dir.create(config$fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(config$data_dir, recursive = TRUE, showWarnings = FALSE)

organism  <- org.Hs.eg.db
PI_THRESH <- 0.05

# All 3 SURV contrasts are still needed to build load_sig_proteins()'s
# "significant in ANY contrast" protein pool (matches the original script's
# convention), but only these two go on the bar chart itself.
CONTRASTS      <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")
BAR_CONTRASTS  <- c("Baseline_SURVvCTL", "Training_SURV")   # CRC vs CTL, CRC+RET vs CRC

CONTRAST_LABELS <- c(
  Baseline_SURVvCTL = "CRC vs CTL",
  Training_SURVvCTL = "CRC+RET vs CTL",
  Training_SURV     = "CRC+RET vs CRC"
)

CONTRAST_SHEETS <- list(
  Baseline_SURVvCTL = 11,
  Training_SURVvCTL = 12,
  Training_SURV     = 13
)

# Same green/orange pair used for these two contrasts in the fig_03
# volcano+heatmap composite, so this figure reads consistently with it.
contrast_bar_colors <- c(
  Baseline_SURVvCTL = "#66C2A5",   # CRC vs CTL
  Training_SURV     = "#FC8D62"    # CRC+RET vs CRC
)

# --- KEYWORD-BASED PATHWAY THEMES --------------------------------------------
# Same keyword lists as the heatmap-sankey script's PATHWAY_THEMES (the
# first theme is renamed metabolic -> matches how you referred to it --
# same keywords, same underlying GO/KEGG/Hallmark/Reactome matches).
PATHWAY_THEMES <- list(
  metabolic = c(
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

THEME_LABELS <- c(
  metabolic      = "Metabolic",
  immune         = "Immune response",
  structural_ecm = "Extracellular matrix /\nstructural"
)

# =============================================================================
# 1. LOAD GSEA + LIMMA DATA
# =============================================================================

gsea_pathway_data <- readRDS(config$gsea_data)
gsea_supp_full    <- gsea_pathway_data$gsea_supp_full %>% filter(group == "SURV")

# =============================================================================
# 1a. KEYWORD-BASED PATHWAY SELECTION (verbatim from the heatmap-sankey script)
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

pathway_specs <- imap(PATHWAY_THEMES, function(keywords, theme_name) {
  spec <- select_pathways_by_keyword(keywords)
  message(sprintf("Theme '%s': %d pathways matched", theme_name, nrow(spec)))
  spec
})

# Preview export -- same "check before you commit" pattern as the sankey
# script, so you can skim what each theme actually pulled in.
iwalk(pathway_specs, function(spec, theme_name) {
  preview <- gsea_supp_full %>%
    dplyr::filter(ID %in% spec$id) %>%
    dplyr::distinct(ID, Description_clean, database) %>%
    dplyr::arrange(Description_clean)
  write.csv(preview, file.path(config$fig_dir, sprintf("theme_%s_matched_pathways.csv", theme_name)),
            row.names = FALSE)
})

# =============================================================================
# 1b. LOAD SIGNIFICANT PROTEIN TABLES (verbatim from the heatmap-sankey script)
# =============================================================================

load_sig_proteins <- function(limma_path) {
  sheets <- map(CONTRAST_SHEETS, ~ read_excel(limma_path, sheet = .x))
  
  has_ensembl <- "Ensembl" %in% names(sheets[[1]])
  
  if (!has_ensembl) {
    all_uniprot_ids <- unique(sheets[[1]]$uniprot_id)
    ens_map <- AnnotationDbi::select(
      organism, keys = all_uniprot_ids, keytype = "UNIPROT", columns = "ENSEMBL"
    ) %>%
      dplyr::filter(!is.na(ENSEMBL)) %>%
      dplyr::distinct(UNIPROT, .keep_all = TRUE) %>%
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

# =============================================================================
# 2. HARMONIZE PATHWAY MEMBERSHIP TO ENSEMBL (verbatim from the heatmap-sankey script)
# =============================================================================

harmonize_gene_ids <- function(id_vec, database) {
  if (database == "KEGG") {
    map_df <- bitr(id_vec, fromType = "ENTREZID", toType = "ENSEMBL", OrgDb = organism)
    unique(na.omit(map_df$ENSEMBL))
  } else {
    unique(id_vec)
  }
}

build_pathway_gene_pool <- function(pathway_id, tissue_name) {
  rows <- gsea_supp_full %>%
    filter(ID == pathway_id, tissue == tissue_name, significant)
  
  if (nrow(rows) == 0) return(NULL)
  
  db          <- as.character(rows$database[1])
  raw_ids     <- unlist(strsplit(rows$core_enrichment, "/"))
  ensembl_ids <- harmonize_gene_ids(raw_ids, db)
  
  list(id = pathway_id, label = rows$Description_clean[1], database = db, ensembl = ensembl_ids)
}

# =============================================================================
# 3. PROTEIN -> THEME MEMBERSHIP
# =============================================================================
# Union of Ensembl IDs across every one of a theme's keyword-matched
# pathways that's significant IN THIS TISSUE (not per-contrast -- a
# pathway only needs to be significant somewhere, matching
# build_pathway_gene_pool()'s existing per-tissue scoping).

build_theme_gene_set <- function(theme_name, tissue_name) {
  spec  <- pathway_specs[[theme_name]]
  pools <- map(spec$id, build_pathway_gene_pool, tissue_name = tissue_name)
  pools <- compact(pools)
  unique(unlist(map(pools, "ensembl")))
}

# =============================================================================
# 4. MEAN LOGFC + SEM + ONE-SAMPLE T-TEST PER THEME x CONTRAST
# =============================================================================
# For each theme, restrict to the tissue's significant-in-any-contrast
# proteins that ALSO map to >=1 of that theme's pathways. Bar height = mean
# logFC across those proteins for a given contrast (regardless of whether
# each individual protein was itself significant in THAT specific
# contrast -- same "group average" logic the reference figure uses).
# Significance = one-sample t-test against mu = 0.

get_stars <- function(p) {
  dplyr::case_when(
    is.na(p)     ~ "",
    p < 0.001    ~ "***",
    p < 0.01     ~ "**",
    p < 0.05     ~ "*",
    TRUE         ~ ""
  )
}

compute_theme_bar_stats <- function(tissue_name) {
  sig_proteins <- sig_proteins_by_tissue[[tissue_name]]
  
  map_dfr(names(PATHWAY_THEMES), function(theme_name) {
    ens_set <- build_theme_gene_set(theme_name, tissue_name)
    theme_proteins <- sig_proteins %>% dplyr::filter(ensembl %in% ens_set)
    
    if (nrow(theme_proteins) == 0) {
      warning(sprintf("[%s/%s] No proteins mapped to this theme -- skipping.", tissue_name, theme_name))
      return(NULL)
    }
    
    message(sprintf("[%s/%s] %d proteins mapped", tissue_name, theme_name, nrow(theme_proteins)))
    
    map_dfr(BAR_CONTRASTS, function(contrast) {
      vals <- theme_proteins[[contrast]]
      vals <- vals[!is.na(vals)]
      n <- length(vals)
      
      if (n < 2) {
        warning(sprintf("[%s/%s/%s] Only %d protein(s) with a logFC value -- t-test skipped.",
                        tissue_name, theme_name, contrast, n))
        return(tibble(theme = theme_name, contrast = contrast, n = n,
                      mean_logfc = if (n == 1) vals else NA_real_,
                      sem = NA_real_, p_value = NA_real_, stars = ""))
      }
      
      tt <- t.test(vals, mu = 0)
      tibble(
        theme = theme_name, contrast = contrast, n = n,
        mean_logfc = mean(vals), sem = sd(vals) / sqrt(n),
        p_value = tt$p.value, stars = get_stars(tt$p.value)
      )
    })
  }) %>%
    mutate(
      theme_label    = factor(unname(THEME_LABELS[theme]),
                              levels = unname(THEME_LABELS[names(PATHWAY_THEMES)])),
      contrast_label = factor(unname(CONTRAST_LABELS[contrast]),
                              levels = unname(CONTRAST_LABELS[BAR_CONTRASTS]))
    )
}

# =============================================================================
# 4b. PAIRED T-TEST -- CRC vs CTL logFC vs. CRC+RET vs CRC logFC, per theme
# =============================================================================
# Different question from compute_theme_bar_stats() above: not "does this
# contrast's mean logFC differ from 0", but "do the two contrasts differ
# from EACH OTHER" for the same set of proteins. Paired (not independent),
# since it's the same protein's logFC in both contrasts -- restricted to
# proteins with a non-NA value in BOTH contrasts so the pairing is exact.
compute_theme_comparison_stats <- function(tissue_name) {
  sig_proteins <- sig_proteins_by_tissue[[tissue_name]]
  contrast_a <- BAR_CONTRASTS[1]
  contrast_b <- BAR_CONTRASTS[2]
  
  map_dfr(names(PATHWAY_THEMES), function(theme_name) {
    ens_set <- build_theme_gene_set(theme_name, tissue_name)
    theme_proteins <- sig_proteins %>%
      dplyr::filter(ensembl %in% ens_set,
                    !is.na(.data[[contrast_a]]), !is.na(.data[[contrast_b]]))
    
    n <- nrow(theme_proteins)
    if (n < 2) {
      warning(sprintf("[%s/%s] Only %d paired protein(s) -- comparison t-test skipped.",
                      tissue_name, theme_name, n))
      return(tibble(theme = theme_name, n = n, mean_diff = NA_real_,
                    p_value = NA_real_, stars = ""))
    }
    
    vals_a <- theme_proteins[[contrast_a]]
    vals_b <- theme_proteins[[contrast_b]]
    tt <- t.test(vals_a, vals_b, paired = TRUE)
    
    tibble(theme = theme_name, n = n, mean_diff = mean(vals_a - vals_b),
           p_value = tt$p.value, stars = get_stars(tt$p.value))
  }) %>%
    mutate(theme_label = factor(unname(THEME_LABELS[theme]),
                                levels = unname(THEME_LABELS[names(PATHWAY_THEMES)])))
}

# =============================================================================
# 5. PLOT -- grouped bar chart, error bars = SEM, asterisks = t-test sig.,
# plus a bracket per theme for the paired contrast-vs-contrast comparison.
# =============================================================================

plot_theme_bars <- function(bar_stats, comparison_stats, tissue_name) {
  # Headroom above/below the tallest bar+error+asterisk, scaled to the
  # data range rather than a fixed constant.
  y_extent <- max(abs(bar_stats$mean_logfc) + bar_stats$sem, na.rm = TRUE)
  y_pad    <- y_extent * 0.35
  
  fill_values <- setNames(unname(contrast_bar_colors[BAR_CONTRASTS]),
                          unname(CONTRAST_LABELS[BAR_CONTRASTS]))
  
  # --- Bracket geometry, one per theme, spanning its two dodged bars ---
  # position_dodge(width = 0.8) with 2 groups centers bars at +-0.2 from
  # the theme's integer x position.
  # NOTE: bar_top below is a simplification (max(mean_logfc + sem) only,
  # not accounting for the case where a bar's SEM-adjusted low end is
  # further from 0 than its high end) -- fine for the small, mostly-
  # positive logFC values typical here; revisit if a theme's bars come
  # out strongly negative and the bracket ends up sitting oddly.
  tick_h <- y_pad * 0.12
  bracket_df <- comparison_stats %>%
    left_join(
      bar_stats %>% group_by(theme_label) %>%
        summarise(bar_top = max(mean_logfc + sem, na.rm = TRUE), .groups = "drop"),
      by = "theme_label"
    ) %>%
    mutate(
      x_left  = as.numeric(theme_label) - 0.2,
      x_right = as.numeric(theme_label) + 0.2,
      y       = bar_top + y_pad * 0.55,
      label   = dplyr::case_when(is.na(p_value) ~ "", p_value < 0.05 ~ stars, TRUE ~ "n.s.")
    )
  
  y_top <- max(y_extent + y_pad, max(bracket_df$y, na.rm = TRUE) + y_pad * 0.3)
  
  ggplot(bar_stats, aes(x = theme_label, y = mean_logfc, fill = contrast_label)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7,
             color = "black", linewidth = 0.3) +
    geom_errorbar(
      aes(ymin = mean_logfc - sem, ymax = mean_logfc + sem),
      position = position_dodge(width = 0.8), width = 0.2, linewidth = 0.4, color = "grey30"
    ) +
    geom_text(
      aes(y = mean_logfc + sign(mean_logfc) * (sem + y_pad * 0.15), label = stars),
      position = position_dodge(width = 0.8), size = 5, fontface = "bold"
    ) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
    # --- comparison bracket: horizontal bar + two downward end-ticks + label ---
    geom_segment(data = bracket_df, aes(x = x_left, xend = x_right, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket_df, aes(x = x_left, xend = x_left, y = y, yend = y - tick_h),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket_df, aes(x = x_right, xend = x_right, y = y, yend = y - tick_h),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = bracket_df, aes(x = (x_left + x_right) / 2, y = y + y_pad * 0.08, label = label),
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_fill_manual(values = fill_values, name = NULL) +
    coord_cartesian(ylim = c(-y_extent - y_pad, y_top)) +
    labs(x = NULL, y = "Log2 fold change",
         title = tissue_name,
         subtitle = "Mean logFC by pathway group (t-test vs. 0); bracket = paired t-test between contrasts") +
    theme_cr(base_size = 12) +
    theme(legend.position = "top",
          axis.text.x = element_text(face = "bold"))
}

# =============================================================================
# 5b. RAW PER-PROTEIN LOGFC (long format) -- for the boxplot+jitter version
# =============================================================================
# Same theme/protein membership as compute_theme_bar_stats(), but keeps
# every individual protein's logFC rather than collapsing to a mean --
# this is what the boxplot/jitter points are drawn from.
build_theme_protein_long <- function(tissue_name) {
  sig_proteins <- sig_proteins_by_tissue[[tissue_name]]
  
  map_dfr(names(PATHWAY_THEMES), function(theme_name) {
    ens_set <- build_theme_gene_set(theme_name, tissue_name)
    theme_proteins <- sig_proteins %>% dplyr::filter(ensembl %in% ens_set)
    
    if (nrow(theme_proteins) == 0) return(NULL)
    
    theme_proteins %>%
      dplyr::select(uniprot, gene, all_of(BAR_CONTRASTS)) %>%
      pivot_longer(all_of(BAR_CONTRASTS), names_to = "contrast", values_to = "logFC") %>%
      dplyr::filter(!is.na(logFC)) %>%
      mutate(theme = theme_name)
  }) %>%
    mutate(
      theme_label    = factor(unname(THEME_LABELS[theme]),
                              levels = unname(THEME_LABELS[names(PATHWAY_THEMES)])),
      contrast_label = factor(unname(CONTRAST_LABELS[contrast]),
                              levels = unname(CONTRAST_LABELS[BAR_CONTRASTS]))
    )
}

# =============================================================================
# 5c. PLOT -- boxplot + jitter version of the same figure
# =============================================================================
# Same t-test results as plot_theme_bars() (bar_stats, comparison_stats
# are passed in, not recomputed), but positions the significance
# annotations off the RAW data's max/min per group rather than mean+-SEM,
# since the boxplot+jitter spread can extend well past one SEM.
plot_theme_boxplot <- function(long_df, bar_stats, comparison_stats, tissue_name) {
  fill_values <- setNames(unname(contrast_bar_colors[BAR_CONTRASTS]),
                          unname(CONTRAST_LABELS[BAR_CONTRASTS]))
  
  data_extent <- max(abs(long_df$logFC), na.rm = TRUE)
  y_pad       <- data_extent * 0.35
  
  # Per theme x contrast top of the raw jitter spread, for star placement
  # -- dodge offset (-0.2 / +0.2) matches position_dodge(width = 0.8)
  # with 2 groups, same as the bar-chart version.
  star_pos <- long_df %>%
    group_by(theme_label, contrast_label) %>%
    summarise(top = max(logFC, na.rm = TRUE), .groups = "drop") %>%
    left_join(bar_stats %>% dplyr::select(theme_label, contrast_label, stars),
              by = c("theme_label", "contrast_label")) %>%
    mutate(x_dodge = as.numeric(theme_label) + if_else(as.numeric(contrast_label) == 1, -0.2, 0.2))
  
  # Comparison bracket -- same construction as plot_theme_bars(), but
  # bar_top now comes from the raw data's max per theme instead of
  # mean + sem.
  tick_h <- y_pad * 0.12
  bracket_df <- comparison_stats %>%
    left_join(
      long_df %>% group_by(theme_label) %>% summarise(bar_top = max(logFC, na.rm = TRUE), .groups = "drop"),
      by = "theme_label"
    ) %>%
    mutate(
      x_left  = as.numeric(theme_label) - 0.2,
      x_right = as.numeric(theme_label) + 0.2,
      y       = bar_top + y_pad * 0.55,
      label   = dplyr::case_when(is.na(p_value) ~ "", p_value < 0.05 ~ stars, TRUE ~ "n.s.")
    )
  
  y_top <- max(data_extent + y_pad, max(bracket_df$y, na.rm = TRUE) + y_pad * 0.3)
  
  ggplot(long_df, aes(x = theme_label, y = logFC, fill = contrast_label)) +
    geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
                 outlier.shape = NA, alpha = 0.7, linewidth = 0.4) +
    geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
               shape = 21, size = 1.4, alpha = 0.6, color = "black", stroke = 0.2) +
    geom_text(data = star_pos, aes(x = x_dodge, y = top + y_pad * 0.12, label = stars),
              inherit.aes = FALSE, size = 5, fontface = "bold") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black") +
    geom_segment(data = bracket_df, aes(x = x_left, xend = x_right, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket_df, aes(x = x_left, xend = x_left, y = y, yend = y - tick_h),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_segment(data = bracket_df, aes(x = x_right, xend = x_right, y = y, yend = y - tick_h),
                 inherit.aes = FALSE, linewidth = 0.4) +
    geom_text(data = bracket_df, aes(x = (x_left + x_right) / 2, y = y + y_pad * 0.08, label = label),
              inherit.aes = FALSE, size = 4, fontface = "bold") +
    scale_fill_manual(values = fill_values, name = NULL) +
    coord_cartesian(ylim = c(-data_extent - y_pad, y_top)) +
    labs(x = NULL, y = "Log2 fold change",
         title = tissue_name,
         subtitle = "Individual protein logFC by pathway group; asterisks = t-test vs. 0; bracket = paired t-test between contrasts") +
    theme_cr(base_size = 12) +
    theme(legend.position = "top",
          axis.text.x = element_text(face = "bold"))
}

# =============================================================================
# 6. RUN + SAVE
# =============================================================================

bar_results <- list()
box_results <- list()

for (tis in c("Muscle", "Plasma")) {
  message(sprintf("--- Building pathway-group bar chart: %s ---", tis))
  
  bar_stats        <- compute_theme_bar_stats(tis)
  comparison_stats <- compute_theme_comparison_stats(tis)
  p_bar <- plot_theme_bars(bar_stats, comparison_stats, tis)
  
  base_name <- sprintf("fig_pathway_bars_%s", tolower(tis))
  ggsave(file.path(config$fig_dir, paste0(base_name, ".pdf")), p_bar,
         width = 7, height = 6, device = cairo_pdf)
  ggsave(file.path(config$fig_dir, paste0(base_name, ".png")), p_bar,
         width = 7, height = 6, dpi = 300)
  write_xlsx(
    list(Group_vs_Zero = bar_stats, Contrast_Comparison = comparison_stats),
    file.path(config$data_dir, paste0(base_name, "_data.xlsx"))
  )
  
  bar_results[[tis]] <- list(stats = bar_stats, comparison_stats = comparison_stats, plot = p_bar)
  
  message(sprintf("--- Building pathway-group boxplot: %s ---", tis))
  
  long_df <- build_theme_protein_long(tis)
  p_box <- plot_theme_boxplot(long_df, bar_stats, comparison_stats, tis)
  
  box_base_name <- sprintf("fig_pathway_boxplot_%s", tolower(tis))
  ggsave(file.path(config$fig_dir, paste0(box_base_name, ".pdf")), p_box,
         width = 7, height = 6, device = cairo_pdf)
  ggsave(file.path(config$fig_dir, paste0(box_base_name, ".png")), p_box,
         width = 7, height = 6, dpi = 300)
  write_xlsx(long_df, file.path(config$data_dir, paste0(box_base_name, "_data.xlsx")))
  
  box_results[[tis]] <- list(long_df = long_df, plot = p_box)
}

message("Done. Pathway-group bar charts + boxplots written to: ", config$fig_dir)