# GSEA / ORA Script
# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, limma,
               writexl, clusterProfiler, org.Hs.eg.db, enrichplot, enrichR,
               msigdbr, ReactomePA, GSEABase, GO.db, AnnotationDbi)
# NOTE: biomaRt removed — GO Slim is now built entirely offline (see Section 0b).

# Configuration
config <- list(
  report_dir = "03_analysis/b_reports/muscle/blood_covar_v2/GSEA",
  data_dir   = "03_analysis/c_data/muscle/blood_covar_v2/GSEA",
  limma = "03_analysis/c_data/muscle/blood_covar_v2/limma/limma_results_all_muscle.xlsx",
  limma_cv = "03_analysis/c_data/muscle/blood_covar_v2/limma/limma_results_all_muscle_cv.xlsx",
  goslim_obo = "resources/goslim_generic.obo"  # see Section 0b for one-time download
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

organism <- org.Hs.eg.db

# --- 0b: Build GO Slim (BP) mapping — fully offline, no BioMart -------
# One-time step: download the generic GO slim OBO file if you don't have it yet.
# This is a static file from geneontology.org, not a live query service, so it
# won't hit the kind of outage BioMart/Ensembl does. Re-run only when you want
# to refresh to a newer GO release.
#
# dir.create("resources", showWarnings = FALSE)
# download.file(
#   "http://current.geneontology.org/ontology/subsets/goslim_generic.obo",
#   destfile = config$goslim_obo
# )

build_goslim_bp <- function(obo_path) {
  if (!file.exists(obo_path)) {
    stop(sprintf(
      "GO slim OBO file not found at '%s'. Download it once (see commented code above) before running GO Slim GSEA.",
      obo_path
    ))
  }
  
  cat(">> Building GO Slim (BP) term2gene map offline...\n")
  
  slim_collection <- GSEABase::getOBOCollection(obo_path)
  slim_ids <- GSEABase::ids(slim_collection)
  
  # Restrict to BP slim terms using GO.db's local ontology info
  slim_ont <- AnnotationDbi::select(GO.db, keys = slim_ids,
                                    columns = c("TERM", "ONTOLOGY"), keytype = "GOID")
  slim_bp_ids <- slim_ont %>% filter(ONTOLOGY == "BP") %>% pull(GOID)
  
  # Expand each slim term to itself + all descendant GO BP terms (local, via GO.db)
  offspring_list <- AnnotationDbi::mget(slim_bp_ids, GO.db::GOBPOFFSPRING, ifnotfound = NA)
  
  slim_expanded <- purrr::map2_dfr(names(offspring_list), offspring_list, function(term, offspring) {
    go_ids <- if (all(is.na(offspring))) term else c(term, offspring)
    data.frame(slim_term = term, GOALL = go_ids, stringsAsFactors = FALSE)
  })
  
  # Map all full GO BP terms to Ensembl genes (local, via org.Hs.eg.db)
  all_ensembl <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENSEMBL")
  full_bp_genes <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = all_ensembl,
    keytype = "ENSEMBL",
    columns = c("GOALL", "ONTOLOGYALL")
  ) %>%
    filter(ONTOLOGYALL == "BP") %>%
    dplyr::select(ENSEMBL, GOALL) %>%
    distinct()
  
  term2gene <- slim_expanded %>%
    inner_join(full_bp_genes, by = "GOALL") %>%
    distinct(slim_term, ENSEMBL) %>%
    dplyr::rename(term = slim_term, gene = ENSEMBL)
  
  term2name <- slim_ont %>%
    filter(GOID %in% slim_bp_ids) %>%
    dplyr::select(GOID, TERM) %>%
    dplyr::rename(term = GOID, name = TERM) %>%
    distinct()
  
  cat(sprintf("   %d slim BP terms, %d term-gene pairs\n", nrow(term2name), nrow(term2gene)))
  list(term2gene = term2gene, term2name = term2name)
}

goslim_available <- FALSE
goslim <- tryCatch(build_goslim_bp(config$goslim_obo), error = function(e) {
  cat("!! GO Slim setup failed — GO Slim GSEA will be skipped this run.\n")
  cat("   (", conditionMessage(e), ")\n")
  NULL
})

if (!is.null(goslim)) {
  slim_term2gene <- goslim$term2gene
  slim_term2name <- goslim$term2name
  goslim_available <- TRUE
}

# --- 1: Load limma Data for pathway analysis ----
sheet_names <- excel_sheets(config$limma_cv)
sheet_names <- sheet_names[-1]

limma_results <- sheet_names %>%
  set_names() %>%
  lapply(function(sh) read_excel(config$limma, sheet = sh))

# --- 2: Reversal Protein Data -------------------------
get_reversed_proteins <- function(df_cancer, df_training, lfc_thresh = 0, p_thresh = 0.05) {
  df_cancer   <- df_cancer   %>% dplyr::select(uniprot, gene, Ensembl, logFC, Pi.Val, adj.P.Val)
  df_training <- df_training %>% dplyr::select(uniprot, logFC, Pi.Val, adj.P.Val)
  
  colnames(df_cancer)[4:6]   <- c("logFC_cancer",   "Pi.Val_cancer",   "adj.P.Val_cancer")
  colnames(df_training)[2:4] <- c("logFC_training", "Pi.Val_training", "adj.P.Val_training")
  
  df <- inner_join(df_cancer, df_training, by = "uniprot")
  
  df %>%
    mutate(
      sig_cancer   = Pi.Val_cancer   < p_thresh,
      sig_training = Pi.Val_training < p_thresh,
      sig_cat = case_when(
        sig_cancer & sig_training   ~ "Both",
        sig_cancer & !sig_training  ~ "Cancer only",
        !sig_cancer & sig_training  ~ "Training only",
        TRUE                        ~ "Neither"
      ),
      direction_cancer   = sign(logFC_cancer),
      direction_training = sign(logFC_training),
      reversed    = direction_cancer != direction_training,
      exacerbated = direction_cancer == direction_training,
      concordance = case_when(
        reversed    & logFC_training > 0 ~ "Training_UP_Cancer_DOWN",
        reversed    & logFC_training < 0 ~ "Training_DOWN_Cancer_UP",
        exacerbated & logFC_training > 0 ~ "Training_UP_Cancer_UP",
        exacerbated & logFC_training < 0 ~ "Training_DOWN_Cancer_DOWN",
        TRUE ~ NA_character_
      )
    )
}

# Build concordance data frame
conc_df <- get_reversed_proteins(limma_results$Cancer, limma_results$Training)

# Build separate reversed and exacerbated data frames
rev_df  <- conc_df %>% filter(reversed)
ex_df   <- conc_df %>% filter(!reversed)

# --- 3: Identify DA Proteins (PRE/POST presence) ----
build_presence_sets <- function(df_pre, df_post) {
  df_pre$source  <- "PRE"
  df_post$source <- "POST"
  combined <- bind_rows(df_pre, df_post) %>%
    group_by(uniprot) %>%
    mutate(presence = case_when(
      "PRE" %in% source & "POST" %in% source ~ "Both",
      "PRE" %in% source                      ~ "PRE only",
      TRUE                                   ~ "POST only"
    )) %>%
    ungroup()
  list(
    pre_only  = combined %>% filter(presence == "PRE only"),
    post_only = combined %>% filter(presence == "POST only"),
    both      = combined %>% filter(presence == "Both", source == "PRE")
  )
}

up_sets   <- build_presence_sets(da_list$Baseline_SURVvCTL_UP,   da_list$Training_SURVvCTL_UP)
down_sets <- build_presence_sets(da_list$Baseline_SURVvCTL_DOWN, da_list$Training_SURVvCTL_DOWN)

da_only_list <- list(
  UP_PRE_only    = up_sets$pre_only,
  UP_POST_only   = up_sets$post_only,
  UP_both        = up_sets$both,
  DOWN_PRE_only  = down_sets$pre_only,
  DOWN_POST_only = down_sets$post_only,
  DOWN_both      = down_sets$both
)


# --- 4: GSEA Pipeline Function ---
run_gsea_pipeline <- function(df, label) {
  
  cat(sprintf("\n>> Running GSEA pipeline for: %s\n", label))
  set.seed(26)
  
  # --- Rank metric ---
  original_gene_list        <- df$t
  names(original_gene_list) <- df$Ensembl
  gene_list                 <- sort(na.omit(original_gene_list), decreasing = TRUE)
  
  # --- Hallmark GSEA ---
  hallmark <- msigdbr(collection = "H") %>% dplyr::select(gs_name, ensembl_gene)
  h2       <- msigdbr(collection = "H") %>% dplyr::select(gs_name, gs_description)
  gsea_hall <- GSEA(gene_list, TERM2GENE = hallmark, TERM2NAME = h2, eps = 0,
                    minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
                    pAdjustMethod = "BH", seed = 26, verbose = FALSE)
  
  # --- Reactome GSEA ---
  react  <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, ensembl_gene)
  react2 <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, gs_description)
  gsea_reactome <- GSEA(gene_list, TERM2GENE = react, TERM2NAME = react2, eps = 0,
                        minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
                        pAdjustMethod = "BH", seed = 26, verbose = FALSE)
  
  # --- GO GSEA (BP / CC / MF) ---
  run_go_gsea <- function(ont, pcut) {
    res <- gseGO(geneList = gene_list, ont = ont, keyType = "ENSEMBL", seed = 26,
                 eps = 0, minGSSize = 15, maxGSSize = 500, pvalueCutoff = pcut,
                 verbose = FALSE, OrgDb = organism, pAdjustMethod = "BH")
    clean        <- res
    clean@result <- res@result %>% filter(!is.na(NES), !is.na(pvalue))
    simple       <- simplify(clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
    list(clean = clean, simple = simple)
  }
  
  bp <- run_go_gsea("BP", 1)
  cc <- run_go_gsea("CC", 0.05)
  mf <- run_go_gsea("MF", 0.05)
  
  # --- GO Slim GSEA (BP) — offline mapping from Section 0b ---
  gsea_bp_slim <- if (goslim_available) {
    GSEA(gene_list, TERM2GENE = slim_term2gene, TERM2NAME = slim_term2name, eps = 0,
         minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
         pAdjustMethod = "BH", seed = 26, verbose = FALSE)
  } else NULL
  
  # --- KEGG GSEA ---
  ids       <- bitr(names(original_gene_list), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]
  df_kegg   <- df[df$Ensembl %in% dedup_ids$ENSEMBL, ]
  df_kegg$entrez <- dedup_ids$ENTREZID[match(df_kegg$Ensembl, dedup_ids$ENSEMBL)]
  
  entrez_gene_list        <- df_kegg$t
  names(entrez_gene_list) <- df_kegg$entrez
  entrez_gene_list        <- sort(na.omit(entrez_gene_list), decreasing = TRUE)
  
  kegg <- gseKEGG(geneList = entrez_gene_list, organism = "hsa", eps = 0,
                  minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, seed = 26,
                  pAdjustMethod = "BH", keyType = "ncbi-geneid")
  
  # --- GO ORA ---
  sig_genes      <- df %>% filter(Pi.Val < 0.05)
  sig_entrez_up   <- bitr(sig_genes %>% filter(logFC > 0) %>% pull(Ensembl), "ENSEMBL", "ENTREZID", organism)
  sig_entrez_down <- bitr(sig_genes %>% filter(logFC < 0) %>% pull(Ensembl), "ENSEMBL", "ENTREZID", organism)
  bg_entrez       <- bitr(df$Ensembl, "ENSEMBL", "ENTREZID", organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    if (nrow(entrez_ids) == 0) return(NULL)
    tryCatch(
      enrichGO(gene = entrez_ids$ENTREZID, universe = bg_entrez$ENTREZID, OrgDb = organism,
               ont = ont_type, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  ora <- list(
    bp_up   = run_go_ora(sig_entrez_up,   "BP", "UP"),
    cc_up   = run_go_ora(sig_entrez_up,   "CC", "UP"),
    mf_up   = run_go_ora(sig_entrez_up,   "MF", "UP"),
    bp_down = run_go_ora(sig_entrez_down, "BP", "DOWN"),
    cc_down = run_go_ora(sig_entrez_down, "CC", "DOWN"),
    mf_down = run_go_ora(sig_entrez_down, "MF", "DOWN")
  )
  
  results <- list(
    label = label,
    gsea_bp_clean = bp$clean, gsea_bp_simple = bp$simple, gsea_bp_slim = gsea_bp_slim,
    gsea_cc_clean = cc$clean, gsea_cc_simple = cc$simple,
    gsea_mf_clean = mf$clean, gsea_mf_simple = mf$simple,
    gsea_kegg = kegg, gsea_hall = gsea_hall, gsea_reactome = gsea_reactome,
    ora_bp_up = ora$bp_up, ora_cc_up = ora$cc_up, ora_mf_up = ora$mf_up,
    ora_bp_down = ora$bp_down, ora_cc_down = ora$cc_down, ora_mf_down = ora$mf_down
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("GSEA_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}


# --- 5: ORA-only Pipeline Function ---
run_ora_pipeline <- function(df, label, bg_df) {
  
  cat(sprintf("\n>> Running ORA pipeline for: %s\n", label))
  
  bg_entrez       <- bitr(bg_df$Ensembl, "ENSEMBL", "ENTREZID", organism)
  sig_entrez_up   <- bitr(df %>% filter(logFC > 0) %>% pull(Ensembl), "ENSEMBL", "ENTREZID", organism)
  sig_entrez_down <- bitr(df %>% filter(logFC < 0) %>% pull(Ensembl), "ENSEMBL", "ENTREZID", organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    if (nrow(entrez_ids) == 0) return(NULL)
    tryCatch(
      enrichGO(gene = entrez_ids$ENTREZID, universe = bg_entrez$ENTREZID, OrgDb = organism,
               ont = ont_type, pAdjustMethod = "BH", pvalueCutoff = 0.05,
               qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  results <- list(
    label = label,
    ora_bp_up   = run_go_ora(sig_entrez_up,   "BP", "UP"),
    ora_cc_up   = run_go_ora(sig_entrez_up,   "CC", "UP"),
    ora_mf_up   = run_go_ora(sig_entrez_up,   "MF", "UP"),
    ora_bp_down = run_go_ora(sig_entrez_down, "BP", "DOWN"),
    ora_cc_down = run_go_ora(sig_entrez_down, "CC", "DOWN"),
    ora_mf_down = run_go_ora(sig_entrez_down, "MF", "DOWN")
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("ORA_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}


# --- 6: Run GSEA for all limma contrasts ----
gsea_run_names <- c("Cancer", "Training", "Residual")

lapply(gsea_run_names, function(nm) {
  run_gsea_pipeline(df = limma_results[[nm]], label = nm)
})


# --- 7: Run ORA for PRE/POST-only presence sets + reversal/exacerbation ----
bg_full <- limma_results$Cancer

lapply(names(da_only_list), function(nm) {
  run_ora_pipeline(da_only_list[[nm]], nm, "cancer", bg_df = bg_full)
})

run_ora_pipeline(rev_df, "REVERSAL", bg_df = bg_full)
run_ora_pipeline(ex_df,  "EXACERBATION", bg_df = bg_full)


# --- 8: Load Results Helper Function ---
load_results <- function(label, folder) {
  path <- file.path(config$report_dir, folder, "gsea", sprintf("gsea_results_%s.rds", label))
  if (!file.exists(path)) {
    stop(sprintf("No saved results found for '%s'. Run run_gsea_pipeline() first.", label))
  }
  readRDS(path)
}