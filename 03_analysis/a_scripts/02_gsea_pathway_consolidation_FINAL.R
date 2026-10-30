# GSEA Pathway Consolidation — one canonical pass over all tissues x databases
#
# PURPOSE
# Every downstream figure that touches GSEA results (bubble plot, heatmap-
# Sankey, ORA butterfly, reversal NES plots, etc.) needs the same three
# things: (1) the full tested-pathway table across every database, (2) the
# significance flag, and (3) which pathways survive redundancy filtering.
# Previously this lived inside fig_gsea_bubble.R (sections 1-3c) and was
# recomputed -- and could silently drift -- every time that script ran, and
# any other figure script that wanted redundancy-aware pathway sets had to
# either duplicate the logic or skip filtering entirely. This script runs
# that logic ONCE and writes a single RDS (+ companion xlsx) that every
# figure script loads. Nothing about the filtering math changes here --
# this is a straight extraction of fig_gsea_bubble.R sections 1-3c,
# generalized to loop over tissue x pooling-family instead of being
# hardcoded to Muscle/Plasma SURV.
#
# WHAT LIVES HERE vs WHAT STAYS IN FIGURE SCRIPTS
# - HERE: loading raw GSEA RDS, building the full tested-pathway table,
#   the significance flag, and the pooled greedy redundancy filter
#   (Jaccard/overlap coefficient on @geneSets). This is "what counts as a
#   real, non-redundant signal" -- an analysis decision, not a plotting one.
# - FIGURE SCRIPTS: TOP_N selection, per-contrast vs pooled ranking,
#   pathway ordering/clustering for display, and all ggplot calls. This is
#   "what subset do I show and how" -- a display decision that legitimately
#   differs per figure.
#
# REPORT WORKBOOKS (Section 5/6 below)
# In addition to the canonical RDS + the existing supplement xlsx, this
# script now also writes full human-readable GSEA report workbooks:
#   - one per tissue (b_reports/<tissue>/GSEA), with a README, a summary
#     sheet, an all-pathways sheet, and one sheet per database, and
#   - one combined workbook (c_data/gsea_consolidated) with a Muscle sheet
#     and a Plasma sheet of the same all-pathways data.
# All of these carry a "Proteins" column: the leading-edge genes
# (core_enrichment) that drove that pathway's enrichment score for that
# specific comparison. This is contrast-specific (not the full curated
# gene set), matching the core_enrichment field already extracted in
# Section 2 below.

# --- 0: Setup ------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, writexl)

# --- Configuration ---------------------------------------------------------
# gsea_dirs: tissue -> single directory of per-comparison GSEA_results_*.rds.
# Muscle and Plasma use the identical directory pattern -- there is no more
# SURV/CRE subfolder nesting; all comparisons (cancer-status family + supp
# family) live together in one directory per tissue.
config <- list(
  gsea_dirs = list(
    Muscle = "03_analysis/c_data/muscle/GSEA",
    Plasma = "03_analysis/c_data/plasma/GSEA"
  ),
  out_dir = "03_analysis/c_data/gsea_consolidated",
  # Per-tissue human-readable report workbooks (Section 5).
  report_dirs = list(
    Muscle = "03_analysis/b_reports/muscle/GSEA",
    Plasma = "03_analysis/b_reports/plasma/GSEA"
  )
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)
for (d in config$report_dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Comparisons pooled together per family -- same pooling-before-filtering
# rationale as before: redundancy is assessed on the UNION of significant
# pathways across all comps in the family, so a shared pathway is filtered
# consistently regardless of which comp(s) flagged it as significant.
#
# "group" below is now a pure pooling label (not a directory key) -- kept
# as SURV/CRE so it lines up with what downstream figure scripts already
# filter on. Both tissues carry all 8 comparisons, so the same comp_defs
# apply to Muscle and Plasma alike.
comp_defs <- list(
  SURV = c("Cancer", "Training", "Residual"),
  CRE  = c("Baseline_supp", "Post_training_supp", "Training_CRE",
           "Training_PLA", "Interaction_supp")
)

# NOTE: display labels for Cancer / Training / Residual are NOT confirmed --
# these don't have an obvious counterpart in the prior SURVvCTL-style naming
# (Baseline_SURVvCTL / Training_SURVvCTL / Training_SURV), and given how
# much the contrast definitions matter downstream (biopsy confound flags,
# causal exercise claims), guessing was worse than leaving them unmapped.
# comparison_label falls back to the raw name for anything not listed here
# (see Section 2) -- fill these in once confirmed.
#
# Baseline_supp/Post_training_supp ARE inferred as renames of the CRE-family
# comparisons formerly called Baseline_CREvPLA/Training_CREvPLA -- the other
# three CRE-family labels (Training_CRE, Training_PLA, Interaction_supp)
# matched exactly, which is what the inference rests on. Confirm before
# trusting these two in a figure/manuscript.
comp_labels <- c(
  Baseline_supp      = "Baseline CRE vs PLA",   # inferred rename of Baseline_CREvPLA -- confirm
  Post_training_supp = "Training CRE vs PLA",   # inferred rename of Training_CREvPLA -- confirm
  Training_CRE       = "POST vs PRE (CRE)",
  Training_PLA       = "POST vs PRE (PLA)",
  Interaction_supp   = "Supplement x Time interaction",
  Cancer = "CRC vs CTL",
  Training = "CRC+RET vs CRC",
  Residual = "CRC+RET vs CTL"
)

db_slots <- c(
  gsea_bp_clean  = "GO:BP",
  gsea_cc_clean  = "GO:CC",
  gsea_mf_clean  = "GO:MF",
  gsea_bp_slim   = "GO Slim",
  gsea_kegg      = "KEGG",
  gsea_reactome  = "Reactome",
  gsea_hall      = "Hallmark"
)
slot_for_database <- setNames(names(db_slots), unname(db_slots))

adjp_thresh <- 0.05

# --- Redundancy filter parameters -------------------------------------------
# Centralized here so every figure inherits the same definition of
# "redundant" instead of each figure script picking its own threshold.
JACCARD_THRESH <- 0.5
OVERLAP_THRESH <- 0.5
MAX_N_FOR_REDUNDANCY <- Inf  # cap per tissue/group/database if O(n^2) gets slow

# --- 1: Load ----------------------------------------------------------------
# One RDS per comparison, each holding a named list of gseaResult objects
# keyed by db_slots names, i.e. gsea_list[[comparison]][[slot]] -> gseaResult.
# Saved files are named "GSEA_results_<label>.rds" (capitalized) and share a
# directory with "ORA_results_*.rds" files, which this loader intentionally
# ignores (only the GSEA prefix is matched).
load_gsea_rds <- function(dir) {
  files <- list.files(dir, pattern = "^GSEA_results_.*\\.rds$",
                      full.names = TRUE, ignore.case = TRUE)
  nms   <- gsub("^GSEA_results_|\\.rds$", "",
                list.files(dir, pattern = "^GSEA_results_.*\\.rds$", ignore.case = TRUE),
                ignore.case = TRUE)
  setNames(lapply(files, readRDS), nms)
}

# tissue -> loaded gsea_list (flat -- no more group-level directory nesting).
gsea_lists <- list()
for (tissue_name in names(config$gsea_dirs)) {
  gsea_lists[[tissue_name]] <- load_gsea_rds(config$gsea_dirs[[tissue_name]])
}

# --- Diagnostic check --------------------------------------------------------
# Run before trusting anything downstream -- if names don't match
# comp_defs/db_slots, fix the mismatch here, not by patching symptoms below.
all_expected_comps <- unique(unlist(comp_defs, use.names = FALSE))

for (tissue_name in names(gsea_lists)) {
  gl <- gsea_lists[[tissue_name]]
  cat("=== ", tissue_name, " ===\n", sep = "")
  cat("Comparisons found:", paste(names(gl), collapse = ", "), "\n")
  cat("Slots in first comparison:", paste(names(gl[[1]]), collapse = ", "), "\n")
  
  if (!all(all_expected_comps %in% names(gl))) {
    stop(tissue_name,
         " is missing comparisons: ",
         paste(setdiff(all_expected_comps, names(gl)), collapse = ", "))
  }
  if (!all(names(db_slots) %in% names(gl[[1]]))) {
    stop(tissue_name,
         " is missing db_slots entries: ",
         paste(setdiff(names(db_slots), names(gl[[1]])), collapse = ", "))
  }
}

# --- 2: Data wrangling -------------------------------------------------------
# Full tested-pathway table across every database in db_slots (not just
# whatever a given figure plots), with core_enrichment retained for the
# redundancy-filter fallback below.
make_gsea_full_data <- function(gsea_list, tissue_name, group_name, comps) {
  bind_rows(lapply(comps, function(comp) {
    bind_rows(lapply(names(db_slots), function(slot) {
      res <- gsea_list[[comp]][[slot]]
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      if (!"core_enrichment" %in% colnames(df)) df$core_enrichment <- NA_character_
      df %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize, core_enrichment) %>%
        dplyr::mutate(
          comparison = comp,
          database   = db_slots[[slot]],
          tissue     = tissue_name,
          group      = group_name
        )
    }))
  }))
}

# tissue x pooling-family grid -- both families apply to both tissues now,
# since there's no directory-level restriction on which comps a tissue has.
tissue_group_grid <- do.call(rbind, lapply(names(gsea_lists), function(t) {
  data.frame(tissue = t, group = names(comp_defs), stringsAsFactors = FALSE)
}))

gsea_supp_full <- purrr::pmap_dfr(
  tissue_group_grid,
  function(tissue, group) {
    make_gsea_full_data(gsea_lists[[tissue]], tissue, group, comp_defs[[group]])
  }
) %>%
  mutate(
    comparison_label  = dplyr::coalesce(comp_labels[as.character(comparison)], as.character(comparison)),
    database          = factor(database, levels = unname(db_slots)),
    tissue            = factor(tissue,   levels = names(config$gsea_dirs)),
    Description_clean = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    significant = p.adjust < adjp_thresh
  )

cat("\n=== Combined data check ===\n")
cat("Rows:", nrow(gsea_supp_full), "\n")
print(table(gsea_supp_full$database, gsea_supp_full$tissue, gsea_supp_full$group))

# format_proteins(): core_enrichment stores leading-edge genes "/"-delimited
# (clusterProfiler convention) -- reformat to a comma-separated list for
# human-readable xlsx output. Vectorized so it can be used inside mutate()
# directly instead of purrr::map_chr + a per-row is.na() check.
format_proteins <- function(x) {
  ifelse(is.na(x) | x == "", NA_character_, str_replace_all(x, "/", ", "))
}

gsea_supp_raw <- gsea_supp_full %>%
  mutate(Proteins = format_proteins(core_enrichment)) %>%
  dplyr::select(tissue, group, database,
                comparison = comparison_label,
                ID, Description = Description_clean,
                NES, pvalue, p.adjust, setSize, significant, Proteins) %>%
  arrange(tissue, group, database, comparison, p.adjust)

gsea_supp_summary <- gsea_supp_raw %>%
  group_by(tissue, group, database) %>%
  summarise(
    n_tested      = n(),
    n_significant = sum(significant, na.rm = TRUE),
    .groups = "drop"
  )

cat("\n=== Summary sheet check (pre-redundancy) ===\n")
print(gsea_supp_summary)

# --- 3: Pooled redundancy filtering (Jaccard + overlap coefficient) ---------
# Greedy significance-ordered FILTER, not clustering: pathways ranked by
# significance, then a pathway is dropped if it's redundant (Jaccard OR
# overlap coefficient >= threshold) with an already-kept, more-significant
# pathway. No transitive grouping -- A and C are never both dropped/merged
# just because each is separately similar to some B.
#
# Filtering is done ONCE on the pooled (union) set of significant pathways
# across all comps in a pooling family, per tissue x group x database --
# NOT independently per contrast (see rationale in header comment above).
#
# Similarity is computed on FULL pathway gene-set membership (gseaResult's
# @geneSets slot), not core_enrichment/leading-edge genes, since leading
# edge is contrast-specific. core_enrichment is used only as a fallback if
# @geneSets isn't available, and flagged as such in the output.
parse_gene_set <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  str_split(x, "/")[[1]]
}

jaccard_index <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / length(union(a, b))
}

overlap_coefficient <- function(a, b) {
  if (length(a) == 0 || length(b) == 0) return(0)
  length(intersect(a, b)) / min(length(a), length(b))
}

get_reference_gene_sets <- function(gsea_list, comp_priority, slot) {
  for (comp in comp_priority) {
    res <- gsea_list[[comp]][[slot]]
    if (is.null(res)) next
    gs <- tryCatch(res@geneSets, error = function(e) NULL)
    if (!is.null(gs) && length(gs) > 0) return(gs)
  }
  NULL
}

empty_redundancy_result <- tibble(
  ID = character(), Description_clean = character(), min_padj = double(),
  tissue = character(), group = character(), database = character(),
  is_representative = logical(), redundant_with = character(), gene_source = character()
)

build_redundancy_mapping <- function(tissue_name, group_name, database_name) {
  slot <- slot_for_database[[database_name]]
  gene_sets_ref <- get_reference_gene_sets(
    gsea_lists[[tissue_name]], comp_defs[[group_name]], slot
  )
  used_fallback <- is.null(gene_sets_ref)
  
  sig_df <- gsea_supp_full %>%
    filter(tissue == tissue_name, group == group_name,
           database == database_name, significant) %>%
    group_by(ID, Description_clean) %>%
    summarise(
      min_padj            = min(p.adjust, na.rm = TRUE),
      core_enrichment_any = dplyr::first(na.omit(core_enrichment)),
      .groups = "drop"
    ) %>%
    arrange(min_padj)
  
  n <- nrow(sig_df)
  if (n == 0) return(empty_redundancy_result)
  
  if (is.finite(MAX_N_FOR_REDUNDANCY) && n > MAX_N_FOR_REDUNDANCY) {
    cat("Note:", tissue_name, "/", group_name, "/", database_name, "--", n,
        "significant pathways exceeds MAX_N_FOR_REDUNDANCY; filtering only the top",
        MAX_N_FOR_REDUNDANCY, "by p.adjust.\n")
    sig_df <- sig_df %>% slice_head(n = MAX_N_FOR_REDUNDANCY)
    n <- nrow(sig_df)
  }
  
  gene_source <- if (used_fallback) "core_enrichment" else "full_gene_set"
  if (used_fallback) {
    cat("Note:", tissue_name, "/", group_name, "/", database_name,
        "-- no @geneSets found on the stored gseaResult; falling back to",
        "core_enrichment (leading edge), which is contrast-specific and a",
        "less exact measure of pathway overlap.\n")
    gene_lists <- lapply(sig_df$core_enrichment_any, parse_gene_set)
  } else {
    gene_lists <- gene_sets_ref[sig_df$ID]
  }
  
  if (n == 1) {
    return(sig_df %>%
             dplyr::select(-core_enrichment_any) %>%
             mutate(tissue = tissue_name, group = group_name, database = database_name,
                    is_representative = TRUE, redundant_with = NA_character_,
                    gene_source = gene_source))
  }
  
  keep           <- rep(TRUE, n)
  redundant_with <- rep(NA_character_, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (!keep[j]) next
      a <- gene_lists[[i]]; b <- gene_lists[[j]]
      if (is.null(a) || is.null(b) || length(a) == 0 || length(b) == 0) next
      if (jaccard_index(a, b) >= JACCARD_THRESH || overlap_coefficient(a, b) >= OVERLAP_THRESH) {
        keep[i] <- FALSE
        redundant_with[i] <- sig_df$ID[j]
        break
      }
    }
  }
  
  sig_df %>%
    dplyr::select(-core_enrichment_any) %>%
    mutate(tissue = tissue_name, group = group_name, database = database_name,
           is_representative = keep, redundant_with = redundant_with,
           gene_source = gene_source)
}

redundancy_grid <- gsea_supp_full %>%
  distinct(tissue, group, database) %>%
  mutate(tissue = as.character(tissue), database = as.character(database)) %>%
  arrange(tissue, group, database)

pathway_redundancy_mapping <- purrr::pmap_dfr(
  redundancy_grid,
  ~ build_redundancy_mapping(tissue_name = ..1, group_name = ..2, database_name = ..3)
) %>%
  relocate(tissue, group, database, .before = ID) %>%
  arrange(tissue, group, database, min_padj)

pathway_redundancy_summary <- redundancy_grid %>%
  left_join(
    pathway_redundancy_mapping %>%
      group_by(tissue, group, database) %>%
      summarise(n_sig_unique = n(), n_representative = sum(is_representative), .groups = "drop"),
    by = c("tissue", "group", "database")
  ) %>%
  mutate(
    n_sig_unique     = coalesce(n_sig_unique, 0L),
    n_representative = coalesce(n_representative, 0L)
  )

cat("\n=== Pooled redundancy filtering (Jaccard >=", JACCARD_THRESH,
    "or Overlap Coefficient >=", OVERLAP_THRESH, ") ===\n")
print(pathway_redundancy_summary)

# representative_lookup: THE lookup every figure script should use to
# restrict candidate pathways to redundancy-filter survivors before doing
# any figure-specific TOP_N/ranking.
representative_lookup <- pathway_redundancy_mapping %>%
  filter(is_representative) %>%
  distinct(tissue, group, database, ID)

# Fold unique/representative counts into the same summary as n_tested/n_significant
gsea_supp_summary <- gsea_supp_summary %>%
  mutate(tissue = as.character(tissue), database = as.character(database)) %>%
  left_join(pathway_redundancy_summary, by = c("tissue", "group", "database"))

cat("\n=== Summary sheet check (with redundancy columns) ===\n")
print(gsea_supp_summary)

# --- 4: Save -----------------------------------------------------------------
# ONE canonical object. Every figure script loads this instead of
# recomputing sections 1-3 above.
gsea_pathway_data <- list(
  gsea_supp_full             = gsea_supp_full,
  gsea_supp_raw              = gsea_supp_raw,
  gsea_supp_summary          = gsea_supp_summary,
  pathway_redundancy_mapping = pathway_redundancy_mapping,
  representative_lookup      = representative_lookup,
  pathway_redundancy_summary = pathway_redundancy_summary,
  params = list(
    adjp_thresh          = adjp_thresh,
    jaccard_thresh       = JACCARD_THRESH,
    overlap_thresh       = OVERLAP_THRESH,
    max_n_for_redundancy = MAX_N_FOR_REDUNDANCY,
    comp_defs            = comp_defs,
    comp_labels           = comp_labels,
    db_slots              = db_slots,
    generated_at          = Sys.time()
  )
)

saveRDS(gsea_pathway_data, file.path(config$out_dir, "gsea_pathway_data_consolidated.rds"))

# gsea_supp_raw now carries the Proteins column (added in Section 2), so the
# Muscle/Plasma sheets below get leading-edge protein mappings for free.
# Redundancy_Mapping is left as-is: it's a filter decision table (kept/
# dropped/redundant_with), not a per-pathway results table, and it already
# records which gene source (full_gene_set vs core_enrichment) drove each
# filtering decision via gene_source.
writexl::write_xlsx(
  list(
    Summary            = gsea_supp_summary,
    Redundancy_Mapping = pathway_redundancy_mapping,
    Muscle             = filter(gsea_supp_raw, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma             = filter(gsea_supp_raw, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$out_dir, "gsea_pathway_data_supplement.xlsx")
)

cat("\nSaved consolidated GSEA pathway data to:\n  ",
    file.path(config$out_dir, "gsea_pathway_data_consolidated.rds"), "\n  ",
    file.path(config$out_dir, "gsea_pathway_data_supplement.xlsx"), "\n", sep = "")

# --- 5: Per-tissue GSEA report workbooks -------------------------------------
# Full, human-readable report per tissue -- every tested pathway, across
# every database, with proteins mapped in. Structure:
#   README        -- outline of the sheets below
#   Summary       -- tested/significant pathway counts per group x database
#   All_Pathways  -- every tested pathway across every database, one table
#   <database>    -- one sheet per database (GO:BP, GO:CC, ... Hallmark),
#                     same columns as All_Pathways, for quick lookup
#
# sanitize_sheet_name(): Excel sheet names can't contain : \ / ? * [ ] and
# are capped at 31 chars. Only the GO:BP/CC/MF labels trip this (the colon).
sanitize_sheet_name <- function(x) {
  x <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", x)
  substr(x, 1, 31)
}

databases_present <- unname(db_slots)  # preserves db_slots' declared order

build_tissue_report <- function(tissue_name) {
  full_data <- gsea_supp_full %>%
    filter(tissue == tissue_name) %>%
    mutate(Proteins = format_proteins(core_enrichment),
           database = as.character(database)) %>%
    dplyr::select(group, comparison = comparison_label, database,
                  ID, Description = Description_clean,
                  NES, pvalue, p.adjust, setSize, significant, Proteins) %>%
    arrange(group, database, comparison, p.adjust)
  
  summary_tbl <- gsea_supp_summary %>%
    filter(tissue == tissue_name) %>%
    dplyr::select(-tissue)
  
  readme_tbl <- tibble::tibble(
    Sheet = c("Summary", "All_Pathways", sanitize_sheet_name(databases_present)),
    Description = c(
      paste0("Tested/significant pathway counts (p.adjust < ", adjp_thresh,
             ") per comparison-family (group) and database, plus redundancy-",
             "filter counts (n_sig_unique / n_representative)."),
      "Every tested pathway across every database and comparison, with NES, p-values, and mapped proteins.",
      paste0("Every tested pathway from the ", databases_present, " database only (same columns as All_Pathways).")
    )
  )
  readme_tbl <- dplyr::bind_rows(
    tibble::tibble(Sheet = "README", Description = "This sheet -- what each sheet in this workbook contains."),
    readme_tbl
  )
  readme_tbl$Note <- c(
    paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), " from the consolidated GSEA pathway data."),
    "Group SURV = cancer-status family (Cancer/Training/Residual); group CRE = supplement family.",
    "Proteins = leading-edge genes (core_enrichment) driving that pathway's enrichment score for that specific comparison -- not the full curated gene set.",
    rep(NA_character_, length(databases_present))
  )
  
  db_sheets <- setNames(
    lapply(databases_present, function(db) {
      full_data %>% filter(database == db) %>% dplyr::select(-database)
    }),
    sanitize_sheet_name(databases_present)
  )
  
  sheets <- c(
    list(README = readme_tbl, Summary = summary_tbl, All_Pathways = full_data),
    db_sheets
  )
  
  out_path <- file.path(config$report_dirs[[tissue_name]],
                        sprintf("GSEA_full_results_%s.xlsx", tolower(tissue_name)))
  writexl::write_xlsx(sheets, path = out_path)
  out_path
}

tissue_report_paths <- sapply(names(config$report_dirs), build_tissue_report)

cat("\nSaved per-tissue GSEA report workbooks to:\n  ",
    paste(tissue_report_paths, collapse = "\n  "), "\n", sep = "")

# --- 6: Combined muscle + plasma workbook ------------------------------------
# Same all-pathways data as Section 5's All_Pathways sheets, but both
# tissues in one workbook (one sheet each) for side-by-side reference.
# Lives alongside the RDS/supplement in gsea_consolidated as a distinct
# file -- NOT a replacement for gsea_pathway_data_supplement.xlsx.
build_combined_tissue_sheet <- function(tissue_name) {
  gsea_supp_full %>%
    filter(tissue == tissue_name) %>%
    mutate(Proteins = format_proteins(core_enrichment),
           database = as.character(database)) %>%
    dplyr::select(group, comparison = comparison_label, database,
                  ID, Description = Description_clean,
                  NES, pvalue, p.adjust, setSize, significant, Proteins) %>%
    arrange(group, database, comparison, p.adjust)
}

combined_path <- file.path(config$out_dir, "gsea_full_data_by_tissue.xlsx")
writexl::write_xlsx(
  list(
    Muscle = build_combined_tissue_sheet("Muscle"),
    Plasma = build_combined_tissue_sheet("Plasma")
  ),
  path = combined_path
)

cat("\nSaved combined muscle + plasma GSEA workbook to:\n  ", combined_path, "\n", sep = "")