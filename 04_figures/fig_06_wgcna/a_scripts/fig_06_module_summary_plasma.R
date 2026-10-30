# ============================================================
# fig_06 — data prep
# WGCNA co-expression modules: size, eigengene trajectory,
# top ORA pathways (multi-database)
#
# Reuses report_data.rds from 05_wgcna_muscle.R rather than
# recomputing the network / eigengenes / group stats.
# ============================================================

pacman::p_load(tidyverse, writexl, rprojroot, msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)

root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "04_figures", "a_scripts", "00_theme.R"))

config <- list(
  tissue          = "plasma",
  wgcna_bundle    = file.path(root, "05_wgcna", "report_data_plasma.rds"),
  n_top_pathways  = 5,
  fdr_method      = "BH",
  ora_pvalue_cutoff = 0.05,
  ora_qvalue_cutoff = 0.2,
  out_root        = file.path(root, "04_figures", "fig_06"),
  data_root       = file.path(root, "04_figures", "data", "fig_06")
)

build_output_path <- function(base, ...) {
  p <- file.path(base, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

rd <- readRDS(config$wgcna_bundle)
config$group_col <- rd$config$group_col  # e.g. "cancer_time" — carried forward for fig_06_composite.R

# ---- 1. Module sizes ----
module_sizes <- rd$module_membership %>%
  filter(module != "grey") %>%
  count(module, name = "n_proteins") %>%
  arrange(desc(n_proteins)) %>%
  mutate(module = factor(module, levels = module))  # size-descending order, reused everywhere below

module_order <- levels(module_sizes$module)

# ---- 2. Eigengene trajectory stats (mean +/- SE per module x group) ----
trajectory_stats <- rd$me_long %>%
  filter(module %in% module_order) %>%
  group_by(module, .data[[rd$config$group_col]]) %>%
  summarise(
    mean_eigengene = mean(eigengene, na.rm = TRUE),
    se_eigengene   = sd(eigengene, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(module = factor(module, levels = module_order))

# All 3 pairwise contrasts (CTL/CRC/CRC+RET), reused directly from
# rd$mc_long (05_wgcna_muscle.R section 6c) -- already sign-oriented to
# match the limma contrast convention, no string-matching against raw
# emmeans contrast labels needed (unlike the old single-contrast version).
#
# x1/x2 = positions on Panel B's x-axis (group_labels order: CTL=1,
# CRC=2, CRC+RET=3). tier = stacking level for bracket height in the
# composite script: the two adjacent short spans (CTL-CRC, CRC-CRC+RET)
# share tier 1; the wide CTL-to-CRC+RET span is tier 2, drawn higher so
# it clears the two shorter brackets below it.
contrast_bracket_spec <- tibble::tribble(
  ~named_contrast,      ~x1, ~x2, ~tier,
  "Baseline_SURVvCTL",   1,   2,   1,
  "Training_SURV",       2,   3,   1,
  "Training_SURVvCTL",   1,   3,   2
)

trajectory_brackets <- rd$mc_long %>%
  dplyr::inner_join(contrast_bracket_spec, by = "named_contrast") %>%
  dplyr::mutate(
    module = factor(module, levels = module_order),
    sig = p_adj < 0.05
  ) %>%
  dplyr::select(module, named_contrast, x1, x2, tier, estimate_oriented, p.value, p_adj, sig)

# ---- 3. Multi-database ORA per module ----
all_proteins <- rd$module_membership$protein
id_map <- bitr(all_proteins, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  distinct(UNIPROT, .keep_all = TRUE)
universe_entrez <- id_map$ENTREZID

hallmark_sets  <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene)
reactome_sets  <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene)

clean_term <- function(x, source) {
  x <- str_remove(x, "^HALLMARK_|^REACTOME_")
  x <- str_replace_all(x, "_", " ")
  str_to_title(x)
}

run_module_ora_multidb <- function(module_name) {
  mod_entrez <- id_map %>%
    filter(UNIPROT %in% (rd$module_membership %>% filter(module == module_name) %>% pull(protein))) %>%
    pull(ENTREZID)
  if (length(mod_entrez) < 5) return(NULL)
  
  go_bp_res <- tryCatch(
    enrichGO(mod_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
             pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
             qvalueCutoff = config$ora_qvalue_cutoff, readable = TRUE)@result %>%
      as_tibble() %>% mutate(source = "GO:BP"),
    error = function(e) NULL
  )
  go_cc_res <- tryCatch(
    enrichGO(mod_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "CC",
             pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
             qvalueCutoff = config$ora_qvalue_cutoff, readable = TRUE)@result %>%
      as_tibble() %>% mutate(source = "GO:CC"),
    error = function(e) NULL
  )
  go_mf_res <- tryCatch(
    enrichGO(mod_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db, ont = "MF",
             pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
             qvalueCutoff = config$ora_qvalue_cutoff, readable = TRUE)@result %>%
      as_tibble() %>% mutate(source = "GO:MF"),
    error = function(e) NULL
  )
  kegg_res <- tryCatch(
    enrichKEGG(mod_entrez, organism = "hsa", universe = universe_entrez,
               pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
               qvalueCutoff = config$ora_qvalue_cutoff)@result %>%
      as_tibble() %>% mutate(source = "KEGG"),
    error = function(e) NULL
  )
  hallmark_res <- tryCatch(
    enricher(mod_entrez, universe = universe_entrez, TERM2GENE = hallmark_sets,
             pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
             qvalueCutoff = config$ora_qvalue_cutoff)@result %>%
      as_tibble() %>% mutate(source = "Hallmark"),
    error = function(e) NULL
  )
  reactome_res <- tryCatch(
    enricher(mod_entrez, universe = universe_entrez, TERM2GENE = reactome_sets,
             pAdjustMethod = config$fdr_method, pvalueCutoff = config$ora_pvalue_cutoff,
             qvalueCutoff = config$ora_qvalue_cutoff)@result %>%
      as_tibble() %>% mutate(source = "Reactome"),
    error = function(e) NULL
  )
  
  bind_rows(go_bp_res, go_cc_res, go_mf_res, kegg_res, hallmark_res, reactome_res) %>%
    mutate(module = module_name,
           Description = if_else(source %in% c("Hallmark", "Reactome"),
                                 clean_term(Description, source), str_to_title(Description)))
}

ora_all <- map_dfr(module_order, run_module_ora_multidb)

top_ora <- ora_all %>%
  group_by(module) %>%
  slice_min(order_by = p.adjust, n = config$n_top_pathways, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    module = factor(module, levels = module_order),
    neg_log_fdr = -log10(p.adjust),
    bold_sig = p.adjust < 0.05
  ) %>%
  arrange(module, p.adjust) %>%
  group_by(module) %>%
  mutate(row_id = paste0(as.character(module), "__", row_number())) %>%
  ungroup() %>%
  mutate(row_id = factor(row_id, levels = rev(row_id)))
# row_id (not Description) is the y-aesthetic in the composite script — pathway
# names repeat across modules (e.g. "Muscle Contraction" in more than one),
# and a per-module factor() built on Description directly collides across
# groups, silently dropping rows whose label isn't in the winning group's
# level set. row_id is unique per row; Description is looked up as a label.

# ---- 3b. Hub proteins per module (Panel D) ----
# rd$hub_proteins (from 05_wgcna_muscle.R) already has the top N per module
# by kME_own_module -- this just adds a display gene symbol (falling back
# to the raw protein/UniProt ID if no symbol is found) and the same unique
# row_id pattern used for top_ora, since gene symbols can repeat across
# modules too (a protein isn't module-exclusive just because it's a hub in
# one; the same symbol could in principle appear as a hub elsewhere).
gene_symbol_map <- bitr(rd$hub_proteins$protein, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = org.Hs.eg.db) %>%
  distinct(UNIPROT, .keep_all = TRUE)

top_hubs <- rd$hub_proteins %>%
  left_join(gene_symbol_map, by = c("protein" = "UNIPROT")) %>%
  mutate(
    display_name = coalesce(SYMBOL, protein),
    module = factor(module, levels = module_order)
  ) %>%
  arrange(module, desc(kME_own_module)) %>%
  group_by(module) %>%
  mutate(row_id = paste0(as.character(module), "__", row_number())) %>%
  ungroup() %>%
  mutate(row_id = factor(row_id, levels = rev(row_id)))

n_no_symbol <- sum(is.na(top_hubs$SYMBOL))
if (n_no_symbol > 0) {
  message(sprintf("%d of %d hub proteins had no SYMBOL mapping -- falling back to raw protein ID for display.",
                  n_no_symbol, nrow(top_hubs)))
}

# ---- 3c. Module eigengene delta vs. phenotype delta (Panel E) ----
# rd$delta_trait_long is NULL if config$delta_pheno_cols wasn't filled in
# with real column names when 05_wgcna_muscle.R was run (see that script's
# config for the placeholder). delta_trait is left NULL here too in that
# case -- the composite script checks for this and skips Panel E rather
# than erroring, so the rest of the figure still builds.
if (is.null(rd$delta_trait_long)) {
  message("rd$delta_trait_long is NULL -- config$delta_pheno_cols wasn't set in 05_wgcna_muscle.R. ",
          "Panel E will be skipped until you fill that in and re-run from section 7b.")
  delta_trait <- NULL
} else {
  delta_trait <- rd$delta_trait_long %>%
    mutate(module = factor(module, levels = module_order)) %>%
    filter(!is.na(module))  # drop "grey" or anything outside module_order
}

# ---- 4. Save bundle + companion Excel ----
fig_06_data <- list(
  config              = config,
  module_sizes        = module_sizes,
  module_order        = module_order,
  trajectory_stats    = trajectory_stats,
  trajectory_brackets = trajectory_brackets,
  top_ora             = top_ora,
  ora_all             = ora_all,
  top_hubs            = top_hubs,
  delta_trait         = delta_trait
)
saveRDS(fig_06_data, build_output_path(config$data_root, "fig_06_data_plasma.rds"))

write_xlsx(
  c(list(module_sizes = module_sizes,
         trajectory_stats = trajectory_stats,
         trajectory_brackets = trajectory_brackets,
         top_ora = top_ora,
         top_hubs = top_hubs),
    if (!is.null(delta_trait)) list(delta_trait = delta_trait) else NULL),
  build_output_path(config$data_root, paste0("fig_06_data_", config$tissue, ".xlsx"))
)