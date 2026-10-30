# GSEA / ORA Script - Plasma

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, limma,
               writexl, clusterProfiler, org.Hs.eg.db, enrichplot, enrichR,
               msigdbr, ReactomePA, GSEABase, GO.db, AnnotationDbi)
# NOTE: biomaRt removed from p_load — GO Slim mapping and UniProt->Ensembl
# mapping are now done locally via GSEABase/GO.db/org.Hs.eg.db, with no
# Ensembl web service dependency (see 2026-08 Ensembl mirror outage notes).

# Configuration
config <- list(
  report_dir  = "03_analysis/b_reports/plasma/GSEA",
  data_dir    = "03_analysis/c_data/plasma/GSEA",
  limma  = "03_analysis/c_data/plasma/limma/limma_results_all_plasma.xlsx",
  goslim_obo  = "00_input/resources/goslim_generic.obo"
  # one-time download, see setup note above:
  # download.file("https://current.geneontology.org/ontology/subsets/goslim_generic.obo",
  #               destfile = config$goslim_obo)
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# Set organism for GSEA
organism <- org.Hs.eg.db

# --- GO Slim setup (local OBO file, no biomaRt/Ensembl dependency) ----------
if (!file.exists(config$goslim_obo)) {
  stop(sprintf(
    "GO Slim OBO file not found at '%s'. Download it once with:\n  download.file(\"https://current.geneontology.org/ontology/subsets/goslim_generic.obo\", destfile = \"%s\")",
    config$goslim_obo, config$goslim_obo
  ))
}

slim_collection <- getOBOCollection(config$goslim_obo)
slim_all_ids    <- ids(slim_collection)

# Restrict to Biological Process slim terms (matches original goslim_goa BP scope)
slim_ontology <- suppressWarnings(AnnotationDbi::Ontology(GOTERM[slim_all_ids]))
slim_bp_ids   <- slim_all_ids[!is.na(slim_ontology) & slim_ontology == "BP"]

cat(sprintf(" GO Slim BP terms loaded: %d\n", length(slim_bp_ids)))

# GOALL keytype expands each slim (often ancestor-level) term to every gene
# annotated at or below it in the GO DAG — this reproduces what the old
# biomaRt goslim_goa mapping provided, purely from local annotation data.
slim_term2gene <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = slim_bp_ids,
  columns = "ENSEMBL",
  keytype = "GOALL"
) %>%
  filter(ONTOLOGYALL == "BP") %>%
  dplyr::select(GOALL, ENSEMBL) %>%
  dplyr::rename(term = GOALL, gene = ENSEMBL) %>%
  distinct()

slim_term2name <- tibble(
  term = slim_bp_ids,
  name = AnnotationDbi::Term(GOTERM[slim_bp_ids])
) %>%
  distinct()

cat(sprintf(" GO Slim term-to-gene mappings: %d rows across %d terms\n",
            nrow(slim_term2gene), n_distinct(slim_term2gene$term)))


# --- 1: Load limma Data ----
sheet_names <- excel_sheets(config$limma)
sheet_names <- sheet_names[-1]

limma_results <- sheet_names %>%
  set_names() %>%
  lapply(function(sh) read_excel(config$limma, sheet = sh))

# --- Fetch Ensembl IDs via UniProt mapping (local org.Hs.eg.db, no biomaRt) -
uniprot_ids <- unique(limma_results$Cancer$uniprot)

ensembl_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = uniprot_ids,
  keytype = "UNIPROT",
  columns = c("UNIPROT", "ENSEMBL")
) %>%
  filter(!is.na(UNIPROT), UNIPROT != "", !is.na(ENSEMBL), ENSEMBL != "") %>%
  distinct(UNIPROT, .keep_all = TRUE) %>%
  dplyr::rename(uniprot = UNIPROT, Ensembl = ENSEMBL)

limma_results <- lapply(limma_results, function(df) {
  left_join(df, ensembl_map, by = "uniprot")
})

# Check coverage — org.Hs.eg.db's UniProt keytype can map slightly fewer IDs
# than biomaRt's uniprotswissprot filter for some accessions; if coverage
# looks meaningfully worse than your prior biomaRt-based runs, flag it before
# proceeding rather than silently losing proteins from downstream GSEA/ORA.
n_mapped <- sum(!is.na(limma_results$Cancer$Ensembl))
n_total  <- nrow(limma_results$Cancer)
cat(sprintf(" UniProt -> Ensembl mapping coverage: %d / %d (%.1f%%)\n",
            n_mapped, n_total, 100 * n_mapped / n_total))


# --- 2: Reversal Protein Data ----
get_reversed_proteins <- function(df_cancer, df_training,
                                  lfc_thresh = 0,
                                  p_thresh   = 0.05) {
  df_cancer   <- df_cancer  %>% dplyr::select(uniprot, gene, Ensembl, logFC, Pi.Val, adj.P.Val)
  df_training <- df_training %>% dplyr::select(uniprot, logFC, Pi.Val, adj.P.Val)
  
  colnames(df_cancer)[4:6]   <- c("logFC_cancer",   "Pi.Val_cancer",   "adj.P.Val_cancer")
  colnames(df_training)[2:4] <- c("logFC_training",  "Pi.Val_training", "adj.P.Val_training")
  
  df <- inner_join(df_cancer, df_training, by = "uniprot")
  
  df %>%
    mutate(
      sig_cancer   = Pi.Val_cancer   < p_thresh,
      sig_training = Pi.Val_training < p_thresh,
      sig_cat = case_when(
        sig_cancer  & sig_training  ~ "Both",
        sig_cancer  & !sig_training ~ "Cancer only",
        !sig_cancer & sig_training  ~ "Training only",
        !sig_cancer & !sig_training ~ "Neither"
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
      ),
      logFC = logFC_training
    )
}

conc_df <- get_reversed_proteins(
  limma_results$Cancer,
  limma_results$Training
)

rev_df <- conc_df %>% filter(reversed == TRUE)
ex_df  <- conc_df %>% filter(reversed == FALSE)

# --- 4: GSEA Pipeline Function ----
run_gsea_pipeline <- function(df, label) {
  
  cat(sprintf("\n>> Running GSEA pipeline for: %s\n", label))
  set.seed(26)
  
  # --- Rank metric ---
  df$rank                   <- df$t
  original_gene_list        <- df$rank
  names(original_gene_list) <- df$Ensembl
  gene_list                 <- original_gene_list[!is.na(names(original_gene_list))]
  gene_list                 <- sort(gene_list, decreasing = TRUE)
  
  # --- Hallmark GSEA ---
  hallmark <- msigdbr(collection = "H") %>% dplyr::select(gs_name, ensembl_gene)
  h2       <- msigdbr(collection = "H") %>% dplyr::select(gs_name, gs_description)
  
  gsea_hall <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = hallmark,
    TERM2NAME     = h2,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    seed          = 26,
    verbose       = FALSE
  )
  
  # --- Reactome GSEA ---
  react  <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, ensembl_gene)
  react2 <- msigdbr(collection = "C2", subcollection = "CP:REACTOME") %>% dplyr::select(gs_name, gs_description)
  
  gsea_reactome <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = react,
    TERM2NAME     = react2,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    seed          = 26,
    verbose       = FALSE
  )
  
  # --- GO GSEA (BP) ---
  gsea_bp <- gseGO(
    geneList      = gene_list,
    ont           = "BP",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gsea_bp_clean        <- gsea_bp
  gsea_bp_clean@result <- gsea_bp@result %>% filter(!is.na(NES), !is.na(pvalue))
  gsea_bp_sig          <- gsea_bp_clean
  gsea_bp_sig@result   <- gsea_bp_clean@result %>% filter(p.adjust < 0.05)
  gsea_bp_simple       <- simplify(gsea_bp_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO GSEA (CC) ---
  gsea_cc <- gseGO(
    geneList      = gene_list,
    ont           = "CC",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gsea_cc_clean        <- gsea_cc
  gsea_cc_clean@result <- gsea_cc@result %>% filter(!is.na(NES), !is.na(pvalue))
  gsea_cc_sig          <- gsea_cc_clean
  gsea_cc_sig@result   <- gsea_cc_clean@result %>% filter(p.adjust < 0.05)
  gsea_cc_simple       <- simplify(gsea_cc_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO GSEA (MF) ---
  gsea_mf <- gseGO(
    geneList      = gene_list,
    ont           = "MF",
    keyType       = "ENSEMBL",
    seed          = 26,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 0.05,
    verbose       = FALSE,
    OrgDb         = organism,
    pAdjustMethod = "BH"
  )
  gsea_mf_clean        <- gsea_mf
  gsea_mf_clean@result <- gsea_mf@result %>% filter(!is.na(NES), !is.na(pvalue))
  gsea_mf_sig          <- gsea_mf_clean
  gsea_mf_sig@result   <- gsea_mf_clean@result %>% filter(p.adjust < 0.05)
  gsea_mf_simple       <- simplify(gsea_mf_clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # --- GO Slim GSEA (BP) ---
  gsea_bp_slim <- GSEA(
    geneList      = gene_list,
    TERM2GENE     = slim_term2gene,
    TERM2NAME     = slim_term2name,
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    seed          = 26,
    verbose       = FALSE
  )
  
  # --- KEGG GSEA ---
  ids       <- bitr(names(original_gene_list), fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = organism)
  dedup_ids <- ids[!duplicated(ids$ENSEMBL), ]
  df2       <- df[df$Ensembl %in% dedup_ids$ENSEMBL, ]
  df2$Y     <- dedup_ids$ENTREZID[match(df2$Ensembl, dedup_ids$ENSEMBL)]
  
  entrez_gene_list        <- df2$t
  names(entrez_gene_list) <- df2$Y
  entrez_gene_list        <- na.omit(entrez_gene_list)
  entrez_gene_list        <- sort(entrez_gene_list, decreasing = TRUE)
  
  kegg <- gseKEGG(
    geneList      = entrez_gene_list,
    organism      = "hsa",
    eps           = 0,
    minGSSize     = 15,
    maxGSSize     = 500,
    pvalueCutoff  = 1,
    seed          = 26,
    pAdjustMethod = "BH",
    keyType       = "ncbi-geneid"
  )
  
  # --- GO ORA ---
  sig_genes      <- df %>% filter(Pi.Val < 0.05)
  sig_genes_up   <- sig_genes %>% filter(logFC > 0)
  sig_genes_down <- sig_genes %>% filter(logFC < 0)
  
  sig_entrez_up   <- bitr(sig_genes_up$Ensembl,   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  sig_entrez_down <- bitr(sig_genes_down$Ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  bg_entrez       <- bitr(df$Ensembl,             fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    tryCatch(
      enrichGO(
        gene          = entrez_ids$ENTREZID,
        universe      = bg_entrez$ENTREZID,
        OrgDb         = organism,
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  ora_bp_up   <- run_go_ora(sig_entrez_up,   "BP", "UP")
  ora_cc_up   <- run_go_ora(sig_entrez_up,   "CC", "UP")
  ora_mf_up   <- run_go_ora(sig_entrez_up,   "MF", "UP")
  ora_bp_down <- run_go_ora(sig_entrez_down, "BP", "DOWN")
  ora_cc_down <- run_go_ora(sig_entrez_down, "CC", "DOWN")
  ora_mf_down <- run_go_ora(sig_entrez_down, "MF", "DOWN")
  
  # --- Save .rds ---
  results <- list(
    label          = label,
    gsea_bp_clean  = gsea_bp_clean,
    gsea_bp_simple = gsea_bp_simple,
    gsea_bp_slim   = gsea_bp_slim,
    gsea_cc_clean  = gsea_cc_clean,
    gsea_cc_simple = gsea_cc_simple,
    gsea_mf_clean  = gsea_mf_clean,
    gsea_mf_simple = gsea_mf_simple,
    gsea_kegg      = kegg,
    gsea_hall      = gsea_hall,
    gsea_reactome  = gsea_reactome,
    ora_bp_up      = ora_bp_up,
    ora_cc_up      = ora_cc_up,
    ora_mf_up      = ora_mf_up,
    ora_bp_down    = ora_bp_down,
    ora_cc_down    = ora_cc_down,
    ora_mf_down    = ora_mf_down
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}


# --- 5: ORA Only Pipeline Function ----
run_ora_pipeline <- function(df, label, bg_df) {
  
  cat(sprintf("\n>> Running ORA pipeline for: %s\n", label))
  
  bg_entrez <- bitr(bg_df$Ensembl, fromType = "ENSEMBL",
                    toType = "ENTREZID", OrgDb = organism)
  
  sig_entrez_up   <- bitr(df %>% filter(logFC > 0) %>% pull(Ensembl),
                          fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  sig_entrez_down <- bitr(df %>% filter(logFC < 0) %>% pull(Ensembl),
                          fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = organism)
  
  run_go_ora <- function(entrez_ids, ont_type, direction) {
    if (nrow(entrez_ids) == 0) return(NULL)
    tryCatch(
      enrichGO(
        gene          = entrez_ids$ENTREZID,
        universe      = bg_entrez$ENTREZID,
        OrgDb         = organism,
        ont           = ont_type,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      ),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n", ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  ora_bp_up   <- run_go_ora(sig_entrez_up,   "BP", "UP")
  ora_cc_up   <- run_go_ora(sig_entrez_up,   "CC", "UP")
  ora_mf_up   <- run_go_ora(sig_entrez_up,   "MF", "UP")
  ora_bp_down <- run_go_ora(sig_entrez_down, "BP", "DOWN")
  ora_cc_down <- run_go_ora(sig_entrez_down, "CC", "DOWN")
  ora_mf_down <- run_go_ora(sig_entrez_down, "MF", "DOWN")
  
  results <- list(
    label       = label,
    ora_bp_up   = ora_bp_up,
    ora_cc_up   = ora_cc_up,
    ora_mf_up   = ora_mf_up,
    ora_bp_down = ora_bp_down,
    ora_cc_down = ora_cc_down,
    ora_mf_down = ora_mf_down
  )
  
  gsea_dir <- file.path(config$data_dir)
  if (!dir.exists(gsea_dir)) dir.create(gsea_dir, recursive = TRUE)
  saveRDS(results, file.path(gsea_dir, sprintf("ORA_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  
  invisible(results)
}

# --- 6: Load Results Helper ----
load_results <- function(label, folder) {
  path <- file.path(config$report_dir, folder, "gsea",
                    sprintf("gsea_results_%s.rds", label))
  if (!file.exists(path)) {
    stop(sprintf("No saved results found for '%s'. Run run_gsea_pipeline() first.", label))
  }
  readRDS(path)
}


# --- 7: Run Pipelines ----
gsea_run_names <- c("Cancer", "Training", "Residual")
gsea_run_names <- c("Baseline_supp", "Post_training_supp", "Training_CRE",
                     "Training_PLA", "Interaction_supp")

lapply(gsea_run_names, function(nm) {
  run_gsea_pipeline(df = limma_results[[nm]], label = nm)
})


bg_full <- limma_results$Cancer

run_ora_pipeline(rev_df, "REVERSAL", bg_df = bg_full)
run_ora_pipeline(ex_df,  "EXACERBATION", bg_df = bg_full)


# --- 8: Inspect Results ----
res_PRE <- load_results("Baseline_SURVvCTL", "SURVvCTL")
res_POST <- load_results("Training_SURVvCTL", "SURVvCTL")
res_RET <- load_results("Training_SURV", "SURVvCTL")