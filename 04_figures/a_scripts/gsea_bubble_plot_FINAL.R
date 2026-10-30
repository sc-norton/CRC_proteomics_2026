# Fig GSEA Bubble — Hallmark + GO:BP, pathway x contrast, faceted by database

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               patchwork, cowplot, scales, enrichplot, writexl)

# --- Configuration -----------------------------------------------------------
# Same fields as fig_02.R's config -- only fig_dir/data_dir differ.
config <- list(
  gsea_muscle = "03_analysis/c_data/muscle/GSEA",
  gsea_plasma = "03_analysis/c_data/plasma/GSEA",
  fig_dir     = "04_figures/fig_gsea_bubble",
  data_dir    = "04_figures/fig_gsea_bubble/data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ------------------------------------------------
source("04_figures/a_scripts/00_theme.R")

tissue_colors <- c(Muscle = colors$pink_dark, Plasma = colors$aqua_dark)

# --- Shared constants --------------------------------------------------------
adjp_thresh <- 0.05

BG_COLORS <- c(
  Muscle = scales::alpha(colors$pink,  0.2),
  Plasma = scales::alpha(colors$aqua, 0.2)
)

gsea_comps <- c("Cancer", "Residual", "Training")

comp_labels <- c(
  Cancer = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training     = "CRC+RET vs CRC"
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

# Panels for THIS figure -- Hallmark + GO:BP only. Add slots here later
# (e.g. gsea_reactome) to extend to a supplementary all-database version.
# Order here sets facet row order (top to bottom).
bubble_slots <- c(gsea_hall = "Hallmark", gsea_bp_clean = "GO:BP")

# NES fill uses the same Up/Down red-blue convention as DIR_COLORS in fig_02.R
NES_LOW  <- colors$blue   # negative NES
NES_HIGH <- colors$red    # positive NES

# --- 1: Load Data ------------------------------------------------------------
# Reuses fig_02.R's loader: reads one RDS per comparison from a directory,
# each RDS holding a named list of gseaResult objects keyed by db_slots names,
# i.e. gsea_muscle[[comparison]][[slot]] -> gseaResult
load_gsea_rds_p <- function(dir) {
  files <- list.files(dir, pattern = "^gsea_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^gsea_results_|\\.rds$", "",
                list.files(dir, pattern = "^gsea_results_.*\\.rds$"))
  setNames(lapply(files, readRDS), nms)
}

load_gsea_rds_m <- function(dir) {
  files <- list.files(dir, pattern = "^GSEA_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^GSEA_results_|\\.rds$", "",
                list.files(dir, pattern = "^GSEA_results_.*\\.rds$"))
  setNames(lapply(files, readRDS), nms)
}

gsea_muscle <- load_gsea_rds_m(config$gsea_muscle)
gsea_plasma <- load_gsea_rds_p(config$gsea_plasma)

# --- Diagnostic check --------------------------------------------------------
# Run this before trusting anything downstream. If comparison names or slot
# names below don't match gsea_comps / bubble_slots, fix the mismatch here
# (not by patching symptoms in the wrangling code further down).
cat("=== gsea_muscle structure ===\n")
cat("Comparisons found:", paste(names(gsea_muscle), collapse = ", "), "\n")
cat("Slots in first comparison:",
    paste(names(gsea_muscle[[1]]), collapse = ", "), "\n")
cat("Class of gsea_hall result:", class(gsea_muscle[[1]][["gsea_hall"]]), "\n")
cat("Columns in gsea_hall @result:\n")
print(colnames(as.data.frame(gsea_muscle[[1]][["gsea_hall"]])))

stopifnot(
  "gsea_muscle comparison names must match gsea_comps" =
    all(gsea_comps %in% names(gsea_muscle)),
  "gsea_muscle must contain all bubble_slots" =
    all(names(bubble_slots) %in% names(gsea_muscle[[1]]))
)

# --- 2: Data wrangling -------------------------------------------------------
# Mirrors make_gsea_counts_combined()'s traversal of gsea_muscle/gsea_plasma,
# but pulls full result rows instead of just significant counts.
make_gsea_bubble_data <- function(gsea_list, tissue_name, slots) {
  bind_rows(lapply(gsea_comps, function(comp) {
    bind_rows(lapply(names(slots), function(slot) {
      res <- gsea_list[[comp]][[slot]]
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      df %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize) %>%
        dplyr::mutate(
          comparison = comp,
          database   = slots[[slot]],
          tissue     = tissue_name
        )
    }))
  }))
}

# Same traversal as make_gsea_bubble_data, but also retains core_enrichment
# (the leading-edge gene string) -- needed for the redundancy-filtering
# fallback in section 3c. Kept as a separate function rather than adding
# core_enrichment to make_gsea_bubble_data's output so the main bubble-plot
# data stays exactly as before.
make_gsea_bubble_data_full <- function(gsea_list, tissue_name, slots) {
  bind_rows(lapply(gsea_comps, function(comp) {
    bind_rows(lapply(names(slots), function(slot) {
      res <- gsea_list[[comp]][[slot]]
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      if (!"core_enrichment" %in% colnames(df)) df$core_enrichment <- NA_character_
      df %>%
        dplyr::select(ID, Description, NES, pvalue, p.adjust, setSize, core_enrichment) %>%
        dplyr::mutate(
          comparison = comp,
          database   = slots[[slot]],
          tissue     = tissue_name
        )
    }))
  }))
}

gsea_bubble_raw <- bind_rows(
  make_gsea_bubble_data(gsea_muscle, "Muscle", bubble_slots),
  make_gsea_bubble_data(gsea_plasma, "Plasma", bubble_slots)
) %>%
  mutate(
    comparison = factor(comparison, levels = gsea_comps),
    database   = factor(database,   levels = unname(bubble_slots)),
    tissue     = factor(tissue,     levels = c("Muscle", "Plasma"))
  )

cat("\n=== Combined bubble data check ===\n")
cat("Rows:", nrow(gsea_bubble_raw), "\n")
print(table(gsea_bubble_raw$database, gsea_bubble_raw$comparison, gsea_bubble_raw$tissue))

# --- 3: Clean pathway labels --------------------------------------------------
gsea_bubble_raw <- gsea_bubble_raw %>%
  mutate(
    Description_clean = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence()
  )

# --- 3b: Full-database export (supplement) -----------------------------------
# Same extraction as the plot data, but across every database in db_slots
# (not just bubble_slots' Hallmark + GO:BP), and with no TOP_N / significance
# filtering -- every tested pathway from every db_slots entry is retained.
# gsea_supp_full additionally carries core_enrichment (needed by section 3c);
# gsea_supp_raw is the same slimmed-down sheet as before for the Muscle/Plasma
# supplement tabs.
gsea_supp_full <- bind_rows(
  make_gsea_bubble_data_full(gsea_muscle, "Muscle", db_slots),
  make_gsea_bubble_data_full(gsea_plasma, "Plasma", db_slots)
) %>%
  mutate(
    comparison         = factor(comparison, levels = gsea_comps),
    comparison_label   = comp_labels[as.character(comparison)],
    database           = factor(database, levels = unname(db_slots)),
    tissue             = factor(tissue,   levels = c("Muscle", "Plasma")),
    Description_clean  = Description %>%
      str_remove("^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_sentence(),
    significant         = p.adjust < adjp_thresh
  )

gsea_supp_raw <- gsea_supp_full %>%
  dplyr::select(tissue, database,
                comparison = comparison_label,
                ID, Description = Description_clean,
                NES, pvalue, p.adjust, setSize, significant) %>%
  arrange(tissue, database, comparison, p.adjust)

cat("\n=== Full supplement data check ===\n")
cat("Rows:", nrow(gsea_supp_raw), "\n")
print(table(gsea_supp_raw$database, gsea_supp_raw$tissue))

# --- Summary sheet: pathway counts per database x tissue ---------------------
# Total tested pathways (mirrors the table() check above) plus how many
# clear the significance threshold, so reviewers can see both the universe
# size and the "hit rate" per database/tissue without opening the full sheet.
# n_significant here is total (contrast, pathway) hits -- i.e. a pathway
# significant in all 3 contrasts counts 3 times. Section 3c below adds the
# unique-pathway and redundancy-collapsed counts alongside this.
gsea_supp_summary <- gsea_supp_raw %>%
  group_by(tissue, database) %>%
  summarise(
    n_tested      = n(),
    n_significant = sum(significant, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from  = tissue,
    values_from = c(n_tested, n_significant),
    names_glue  = "{tissue}_{.value}"
  ) %>%
  dplyr::select(database,
                Muscle_n_tested, Muscle_n_significant,
                Plasma_n_tested, Plasma_n_significant) %>%
  arrange(match(database, unname(db_slots)))

cat("\n=== Summary sheet check (pre-redundancy) ===\n")
print(gsea_supp_summary)

# --- 3c: Pooled redundancy filtering (Jaccard + overlap coefficient) ----------
# Answers "how many pathways should I report" alongside n_significant above.
# This is a greedy significance-ordered FILTER, not clustering: pathways are
# ranked by significance, then a pathway is dropped if it's redundant
# (Jaccard OR overlap coefficient >= threshold) with an already-kept,
# more-significant pathway. No transitive grouping -- A and C are never
# both dropped/merged just because each is separately similar to some B.
#
# Filtering is done ONCE on the pooled (union) set of significant pathways
# across all three gsea_comps per database x tissue -- NOT independently
# per contrast. Per-contrast filtering would let a pathway shared between
# two contrasts' significant lists get filtered against DIFFERENT
# more-significant neighbors in each (since what's "already kept" at that
# point depends on what else is in that contrast's list), which silently
# breaks any later cross-contrast overlap comparison. Pooling first
# guarantees a shared term is filtered consistently regardless of which
# contrast(s) flagged it as significant.
#
# Similarity is computed on FULL pathway gene-set membership (gseaResult's
# @geneSets slot), not core_enrichment/leading-edge genes -- leading edge is
# contrast-specific (it depends on that contrast's ranked gene list), so
# using it here would make "how similar are pathway A and B" itself
# contrast-dependent, defeating the point of pooling. core_enrichment is
# used only as a fallback if @geneSets isn't available on the stored object,
# and is flagged as such in the output.
JACCARD_THRESH <- 0.5
OVERLAP_THRESH <- 0.5

# Optional safety valve: if a database/tissue has a very large number of
# significant pathways, the O(n^2) pairwise comparison below can get slow.
# Set to a finite number (e.g. 500) to cap filtering to the top-N most
# significant pathways (by min p.adjust across contrasts) if needed; Inf
# disables the cap.
MAX_N_FOR_REDUNDANCY <- Inf

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

# Pulls the full pathway gene-set list off whichever contrast's gseaResult
# has it first (should be identical across contrasts since all three were
# run against the same reference database -- this just guards against one
# contrast's object being malformed/incomplete).
get_reference_gene_sets <- function(gsea_list, comp_priority, slot) {
  for (comp in comp_priority) {
    res <- gsea_list[[comp]][[slot]]
    if (is.null(res)) next
    gs <- tryCatch(res@geneSets, error = function(e) NULL)
    if (!is.null(gs) && length(gs) > 0) return(gs)
  }
  NULL
}

tissue_gsea_lists <- list(Muscle = gsea_muscle, Plasma = gsea_plasma)
slot_for_database  <- setNames(names(db_slots), unname(db_slots))

empty_redundancy_result <- tibble(
  ID = character(), Description_clean = character(), min_padj = double(),
  tissue = character(), database = character(),
  is_representative = logical(), redundant_with = character(), gene_source = character()
)

# Greedy significance-ordered filter (same principle as
# clusterProfiler::simplify(), generalized to work by gene-set overlap
# rather than GO semantic similarity): rank pathways by significance, then
# walk down the list and drop a pathway if it's redundant (Jaccard OR
# overlap coefficient >= threshold) with an ALREADY-KEPT, more significant
# pathway. Each dropped pathway's redundant_with points to the specific
# more-significant pathway that subsumed it.
#
# Deliberately NOT connected-components clustering: clustering would let A
# and C get grouped together (and only one survive) purely because both
# happen to be similar to some intermediate B, even if A and C aren't
# themselves similar. The greedy filter only ever drops a pathway against
# something concrete and still present in the list, so what "redundant
# with" means stays interpretable and doesn't chain arbitrarily.
build_redundancy_mapping <- function(tissue_name, database_name) {
  slot <- slot_for_database[[database_name]]
  gene_sets_ref <- get_reference_gene_sets(tissue_gsea_lists[[tissue_name]], gsea_comps, slot)
  used_fallback <- is.null(gene_sets_ref)
  
  sig_df <- gsea_supp_full %>%
    filter(tissue == tissue_name, database == database_name, significant) %>%
    group_by(ID, Description_clean) %>%
    summarise(
      min_padj            = min(p.adjust, na.rm = TRUE),
      core_enrichment_any = dplyr::first(na.omit(core_enrichment)),
      .groups = "drop"
    ) %>%
    arrange(min_padj)  # most significant first -- these are preferentially kept
  
  n <- nrow(sig_df)
  if (n == 0) return(empty_redundancy_result)
  
  if (is.finite(MAX_N_FOR_REDUNDANCY) && n > MAX_N_FOR_REDUNDANCY) {
    cat("Note:", tissue_name, "/", database_name, "--", n,
        "significant pathways exceeds MAX_N_FOR_REDUNDANCY; filtering only the top",
        MAX_N_FOR_REDUNDANCY, "by p.adjust.\n")
    sig_df <- sig_df %>% slice_head(n = MAX_N_FOR_REDUNDANCY)
    n <- nrow(sig_df)
  }
  
  gene_source <- if (used_fallback) "core_enrichment" else "full_gene_set"
  if (used_fallback) {
    cat("Note:", tissue_name, "/", database_name,
        "-- no @geneSets found on the stored gseaResult; falling back to",
        "core_enrichment (leading edge) for the redundancy filter, which is",
        "contrast-specific and a less exact measure of pathway overlap.\n")
    gene_lists <- lapply(sig_df$core_enrichment_any, parse_gene_set)
  } else {
    gene_lists <- gene_sets_ref[sig_df$ID]
  }
  
  if (n == 1) {
    return(sig_df %>%
             dplyr::select(-core_enrichment_any) %>%
             mutate(tissue = tissue_name, database = database_name,
                    is_representative = TRUE, redundant_with = NA_character_,
                    gene_source = gene_source))
  }
  
  keep           <- rep(TRUE, n)
  redundant_with <- rep(NA_character_, n)
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      if (!keep[j]) next  # only compare against pathways still kept (more significant)
      a <- gene_lists[[i]]; b <- gene_lists[[j]]
      if (is.null(a) || is.null(b) || length(a) == 0 || length(b) == 0) next
      if (jaccard_index(a, b) >= JACCARD_THRESH || overlap_coefficient(a, b) >= OVERLAP_THRESH) {
        keep[i] <- FALSE
        redundant_with[i] <- sig_df$ID[j]
        break  # redundant with at least one kept pathway -- drop and move on
      }
    }
  }
  
  sig_df %>%
    dplyr::select(-core_enrichment_any) %>%
    mutate(tissue = tissue_name, database = database_name,
           is_representative = keep, redundant_with = redundant_with,
           gene_source = gene_source)
}

redundancy_grid <- gsea_supp_full %>% distinct(tissue, database) %>% arrange(tissue, database)

pathway_redundancy_mapping <- purrr::pmap_dfr(
  redundancy_grid,
  ~ build_redundancy_mapping(tissue_name = ..1, database_name = ..2)
) %>%
  relocate(tissue, database, .before = ID) %>%
  arrange(tissue, database, min_padj)

pathway_redundancy_summary <- pathway_redundancy_mapping %>%
  group_by(tissue, database) %>%
  summarise(n_sig_unique = n(), n_representative = sum(is_representative), .groups = "drop")

# Lookup used by section 4 to actually restrict which pathways CAN be
# plotted to survivors of the greedy redundancy filter -- this is what
# removes near-duplicate pathways (e.g. "cell projection organization" vs
# "plasma membrane bounded cell projection organization") from the figure
# itself, as opposed to pathway_redundancy_summary above which only reports
# counts. tissue/database forced to character since pmap_dfr can leave
# factor-typed grouping columns behind depending on redundancy_grid's types.
representative_lookup <- pathway_redundancy_mapping %>%
  filter(is_representative) %>%
  mutate(tissue = as.character(tissue), database = as.character(database)) %>%
  distinct(tissue, database, ID)

# Ensure every tissue x database combo appears even if it had 0 significant
# pathways (left_join + coalesce rather than relying on group_by to surface
# empty groups, which it won't).
pathway_redundancy_summary <- redundancy_grid %>%
  left_join(pathway_redundancy_summary, by = c("tissue", "database")) %>%
  mutate(
    n_sig_unique     = coalesce(n_sig_unique, 0L),
    n_representative = coalesce(n_representative, 0L)
  )

cat("\n=== Pooled redundancy filtering (Jaccard >=", JACCARD_THRESH,
    "or Overlap Coefficient >=", OVERLAP_THRESH, ") ===\n")
print(pathway_redundancy_summary)

# Fold the unique/representative counts into the same summary table as
# n_tested/n_significant, so the full progression --
#   tested -> significant (raw hits) -> significant (unique pathways) ->
#   representative (redundancy-collapsed)
# -- is visible in one place per database x tissue.
gsea_supp_summary <- gsea_supp_summary %>%
  left_join(
    pathway_redundancy_summary %>%
      pivot_wider(names_from = tissue, values_from = c(n_sig_unique, n_representative),
                  names_glue = "{tissue}_{.value}"),
    by = "database"
  ) %>%
  dplyr::select(database,
                Muscle_n_tested, Muscle_n_significant, Muscle_n_sig_unique, Muscle_n_representative,
                Plasma_n_tested, Plasma_n_significant, Plasma_n_sig_unique, Plasma_n_representative)

cat("\n=== Summary sheet check (with redundancy columns) ===\n")
print(gsea_supp_summary)

# One sheet per tissue, all databases/contrasts included, sorted by adj p
# within each, plus the redundancy mapping for transparency (which
# pathways were filtered out and which more-significant pathway each was
# filtered against). Written to data_dir alongside the plot's underlying data.
writexl::write_xlsx(
  list(
    Summary             = gsea_supp_summary,
    Redundancy_Mapping  = pathway_redundancy_mapping,
    Muscle              = filter(gsea_supp_raw, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma              = filter(gsea_supp_raw, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$data_dir, "gsea_results_all_databases_supplement.xlsx")
)

cat("\nSupplement data saved to:",
    file.path(config$data_dir, "gsea_results_all_databases_supplement.xlsx"), "\n")

# --- 4: Select pathways to plot ----------------------------------------------
# gsea_hall / gsea_bp_clean were run with pvalueCutoff = 1 (per pipeline
# notes) so every tested pathway is present here. Cap what's plotted:
#   adjp_thresh (shared constant, 0.05)  : "significant in >=1 contrast" gate
#   representative_lookup (section 3c)  : drops non-representative pathways
#                                          from redundant clusters BEFORE ranking
#   TOP_N_HALLMARK / TOP_N_GOBP          : cap per database, ranked by best p.adjust
# The redundancy filter and TOP_N cap now work together: candidates are
# reduced by the greedy significance-ordered Jaccard/overlap-coefficient
# filter first (see select_pathways below), then the top N of what's left
# (by p.adjust) are plotted. This is what actually removes near-duplicate
# pairs like "cell projection organization" / "plasma membrane bounded cell
# projection organization" from the figure, not just from a reported count.
TOP_N_HALLMARK <- 30
TOP_N_GOBP     <- 30
TOP_N <- c(Hallmark = TOP_N_HALLMARK, `GO:BP` = TOP_N_GOBP)

select_pathways <- function(df, top_n, tissue_name, database_name) {
  ranked <- df %>%
    group_by(ID, Description_clean) %>%
    summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
    filter(min_padj < adjp_thresh) %>%
    arrange(min_padj)
  
  if (nrow(ranked) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh,
            " -- falling back to top ", top_n, " by p.adjust regardless of significance.")
    ranked <- df %>%
      group_by(ID, Description_clean) %>%
      summarise(min_padj = min(p.adjust, na.rm = TRUE), .groups = "drop") %>%
      arrange(min_padj)
  } else {
    # Restrict candidates to survivors of the pooled greedy redundancy
    # filter (section 3c) BEFORE capping to top_n -- this is the actual
    # fix: it removes near-duplicate/parent-child pathways from plot
    # candidacy rather than just reporting their count separately. Only
    # applied when representative info exists for this tissue/database (it
    # always will here since `ranked` being non-empty means there WERE
    # significant pathways for section 3c to filter).
    rep_ids <- representative_lookup %>%
      filter(tissue == tissue_name, database == database_name) %>%
      pull(ID)
    if (length(rep_ids) > 0) {
      n_before <- nrow(ranked)
      ranked <- ranked %>% filter(ID %in% rep_ids)
      cat("  [", tissue_name, "/", database_name, "] redundancy filter: ",
          n_before, " -> ", nrow(ranked), " candidate pathways\n", sep = "")
    }
  }
  
  ranked %>% slice_head(n = top_n) %>% pull(Description_clean)
}

build_bubble_selection <- function(tissue_name) {
  df <- filter(gsea_bubble_raw, tissue == tissue_name)
  
  keep <- bind_rows(lapply(names(TOP_N), function(db) {
    tibble(database = db,
           Description_clean = select_pathways(filter(df, database == db), TOP_N[[db]],
                                               tissue_name = tissue_name, database_name = db))
  }))
  
  df %>% semi_join(keep, by = c("database", "Description_clean"))
}

# --- 5: Order pathways within each database (hierarchical clustering on NES) -
# facet_grid(scales = "free_y") lets each database panel keep its own
# category set, so a per-database cluster order concatenated across
# databases is sufficient -- no need to force a single global ordering.
order_pathways <- function(df, key_col = "Description_clean") {
  nes_mat <- df %>%
    dplyr::select(all_of(key_col), comparison, NES) %>%
    pivot_wider(names_from = comparison, values_from = NES) %>%
    column_to_rownames(key_col) %>%
    as.matrix()
  
  nes_mat[is.na(nes_mat)] <- 0
  
  if (nrow(nes_mat) > 2) {
    hc <- hclust(dist(nes_mat), method = "average")
    rownames(nes_mat)[hc$order]
  } else {
    rownames(nes_mat)
  }
}

# --- 6: Plot -------------------------------------------------------------
# Single ggplot, database as facet rows (facet_grid + space = "free_y" so
# panel height scales with pathway count per database), one shared legend.
make_gsea_bubble_plot <- function(tissue_name) {
  sel <- build_bubble_selection(tissue_name)
  
  # order within each database, then set as one factor spanning all panels
  ordered_levels <- unlist(lapply(names(TOP_N), function(db) {
    order_pathways(filter(sel, database == db))
  }))
  
  sel <- sel %>%
    mutate(
      Description_clean = factor(Description_clean, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  # Outline (color + stroke) flags significance; fill stays reserved for NES
  # and size for -log10(FDR), so the outline reads as a separate signal
  # rather than duplicating what fill/size already show.
  p <- ggplot(sel, aes(x = comparison, y = Description_clean)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)), 
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = paste0("GSEA \u2014 ", tissue_name), x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position   = "right",
      panel.grid.major  = element_line(color = colors$gray_light),
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y      = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

fig_muscle <- make_gsea_bubble_plot("Muscle")
fig_plasma <- make_gsea_bubble_plot("Plasma")  # likely too sparse to be
# useful in this format -- see note at bottom of script.

# --- 7: Assemble & Save -------------------------------------------------------
plot_height_muscle <- max(4, 0.22 * fig_muscle$n_pathways + 1.5)

ggsave(file.path(config$fig_dir, "gsea_bubble_muscle.pdf"),
       fig_muscle$plot, width = 12, height = plot_height_muscle,
       device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_bubble_muscle.png"),
       fig_muscle$plot, width = 12, height = plot_height_muscle,
       dpi = 300, bg = "white")

cat("\nSaved to:", config$fig_dir, "\n")

# --- 8: All-database top-N bubble plot -----------------------------------------
# Same visual format as sections 4-7, but pools ALL databases (db_slots, not
# just bubble_slots' Hallmark + GO:BP) and selects a single GLOBAL top N by
# p.adjust across every database combined -- database becomes purely a
# facet label rather than each database getting its own quota. A database
# with unusually strong/abundant signal (e.g. GO:BP, just by virtue of
# having many more granular terms than Hallmark) can end up supplying most
# or all of the 30 shown; that's the intended behavior here, not a bug --
# use the per-database TOP_N figure above instead if you want every
# database guaranteed some representation regardless of relative signal.
#
# Redundancy filtering (section 3c/representative_lookup) is still applied
# WITHIN each database before the pooled ranking -- comparing gene-set
# overlap between, say, a Hallmark term and a KEGG term isn't meaningful,
# so cross-database redundancy is out of scope here; only the final top-N
# selection is pooled across databases.
TOP_N_ALL_DB <- 30

# "padj"  -- rank significant survivors by p.adjust (most statistically
#            confident first). Safer default; doesn't privilege noisy
#            small gene sets that can inflate NES.
# "NES"   -- rank significant survivors by |NES| (strongest effect size
#            among pathways already cleared for significance). Reasonable
#            when the goal is "of what I'm confident is enriched, what
#            shifted the most" -- NOT reasonable if used to bypass the
#            significance gate entirely, since NES isn't adjusted for the
#            extra variance small gene sets get from the running-sum
#            statistic (a high |NES| on a small set can just be noise).
# The significance gate (min_padj < adjp_thresh) applies either way -- only
# the sort order among survivors changes.
RANK_METRIC <- "padj"

select_pathways_global <- function(tissue_name, top_n) {
  df   <- gsea_supp_full %>% filter(tissue == tissue_name)
  reps <- representative_lookup %>% filter(tissue == tissue_name)
  
  ranked <- df %>%
    group_by(database, ID, Description_clean) %>%
    summarise(
      min_padj    = min(p.adjust, na.rm = TRUE),
      max_abs_nes = max(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # Restrict to redundancy-filter survivors within each database (falls
    # through unfiltered for any database with 0 significant pathways --
    # same fallback logic as select_pathways() in section 4).
    dplyr::group_by(database) %>%
    dplyr::group_modify(~ {
      db_name <- as.character(.y$database)
      rep_ids <- reps %>% filter(database == db_name) %>% pull(ID)
      if (length(rep_ids) == 0) return(.x)
      dplyr::filter(.x, ID %in% rep_ids)
    }) %>%
    dplyr::ungroup() %>%
    filter(min_padj < adjp_thresh)
  
  ranked <- if (RANK_METRIC == "NES") {
    ranked %>% arrange(desc(max_abs_nes))
  } else {
    ranked %>% arrange(min_padj)
  }
  
  if (nrow(ranked) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh,
            " in any database for ", tissue_name)
    return(ranked %>% dplyr::select(database, Description_clean))
  }
  
  ranked %>% slice_head(n = top_n) %>% dplyr::select(database, Description_clean)
}

build_bubble_selection_global <- function(tissue_name) {
  keep <- select_pathways_global(tissue_name, TOP_N_ALL_DB)
  gsea_supp_full %>%
    filter(tissue == tissue_name) %>%
    semi_join(keep, by = c("database", "Description_clean"))
}

make_gsea_bubble_plot_alldb <- function(tissue_name) {
  sel <- build_bubble_selection_global(tissue_name)
  
  # Facet order follows db_slots' canonical order, restricted to whichever
  # databases actually placed a pathway in the top N -- empty databases
  # don't get a blank panel.
  present_dbs <- unname(db_slots)[unname(db_slots) %in% unique(as.character(sel$database))]
  
  # Pooling across databases means two DIFFERENT pathways (from two
  # different databases) can share the same cleaned display text purely by
  # coincidence (e.g. some GO:BP term and some Reactome/KEGG term both
  # cleaning to an identical string) -- harmless as data, but factor()
  # requires unique level strings, so building the y-axis directly off
  # Description_clean fails as soon as that happens (even though the two
  # occurrences land in different facet panels). pathway_key (database::ID)
  # is unique by construction; Description_clean is kept only as the
  # printed label via scale_y_discrete(labels = ...) below.
  sel <- sel %>%
    mutate(pathway_key = paste(database, ID, sep = "::"))
  
  ordered_levels <- unlist(lapply(present_dbs, function(db) {
    order_pathways(filter(sel, database == db), key_col = "pathway_key")
  }))
  
  label_lookup <- setNames(sel$Description_clean, sel$pathway_key)
  
  sel <- sel %>%
    mutate(
      database    = factor(as.character(database), levels = present_dbs),
      pathway_key = factor(pathway_key, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  p <- ggplot(sel, aes(x = comparison, y = pathway_key)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)),
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_discrete(labels = label_lookup) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = paste0("GSEA \u2014 ", tissue_name,
                        " (Top ", TOP_N_ALL_DB, ", pooled across all databases)"),
         x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position  = "right",
      panel.grid.major = element_line(color = colors$gray_light),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y     = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

fig_muscle_alldb <- make_gsea_bubble_plot_alldb("Muscle")
fig_plasma_alldb <- make_gsea_bubble_plot_alldb("Plasma")  # pooling across
# all 7 databases may rescue plasma from the sparsity noted below, or may
# not -- worth checking n_pathways before assuming it's useful.

plot_height_muscle_alldb <- max(4, 0.22 * fig_muscle_alldb$n_pathways + 1.5)
plot_height_plasma_alldb <- max(4, 0.22 * fig_plasma_alldb$n_pathways + 1.5)

ggsave(file.path(config$fig_dir, "gsea_bubble_muscle_top30_alldb.pdf"),
       fig_muscle_alldb$plot, width = 12, height = plot_height_muscle_alldb,
       device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_bubble_muscle_top30_alldb.png"),
       fig_muscle_alldb$plot, width = 12, height = plot_height_muscle_alldb,
       dpi = 300, bg = "white")

if (fig_plasma_alldb$n_pathways > 0) {
  ggsave(file.path(config$fig_dir, "gsea_bubble_plasma_top30_alldb.pdf"),
         fig_plasma_alldb$plot, width = 12, height = plot_height_plasma_alldb,
         device = cairo_pdf)
  ggsave(file.path(config$fig_dir, "gsea_bubble_plasma_top30_alldb.png"),
         fig_plasma_alldb$plot, width = 12, height = plot_height_plasma_alldb,
         dpi = 300, bg = "white")
} else {
  cat("\nPlasma all-database top-30 plot skipped -- 0 pathways passed",
      "adjp_thresh in any database.\n")
}

cat("\nSaved all-database top-", TOP_N_ALL_DB, " figure(s) to: ", config$fig_dir, "\n", sep = "")

# --- 9: Per-Contrast Top-N Bubble Plots (pooled across all databases) --------
# Same redundancy-filtered candidate pool as section 8 (representative_lookup,
# built from ALL db_slots databases pooled per section 3c), but selection is
# now driven by EACH CONTRAST'S OWN p.adjust rather than pooling significance
# across gsea_comps -- i.e. "top 30 pathways diagnostic for POST vs PRE
# specifically" rather than "top 30 pathways significant anywhere."
#
# The resulting plot still shows all three comparison columns (not just the
# selecting contrast) so the selected pathways' behavior across the other two
# contrasts remains visible -- this is deliberate: given the linear
# dependency between the three contrasts, seeing a pathway's NES pattern
# across all three columns is often the point (e.g. confirming a
# POST-vs-PRE-driven pathway isn't just baseline noise carried through).
TOP_N_PER_CONTRAST <- 30

# Restricted database set for the GO:BP/CC/MF-only variant added in section
# 9c below. Kept as a named constant (not inlined) so it's a single edit
# point if you want a different subset later.
GO_DBS <- c("GO:BP", "GO:CC", "GO:MF")

# `databases = NULL` (default) reproduces the exact original behavior --
# pooled across every db_slots database, unchanged from before this
# parameter was added. Pass a character vector (e.g. GO_DBS) to restrict
# candidacy to just those databases; this filter is applied BEFORE the
# per-database redundancy restriction and the top_n cap, so "top N" means
# top N within the restricted set, not top N overall then filtered down.
select_pathways_per_contrast <- function(tissue_name, contrast_name, top_n, databases = NULL) {
  reps <- representative_lookup %>% filter(tissue == tissue_name)
  
  df <- gsea_supp_full %>%
    filter(tissue == tissue_name, comparison == contrast_name)
  
  if (!is.null(databases)) df <- df %>% filter(as.character(database) %in% databases)
  
  df <- df %>%
    # Restrict to redundancy-filter survivors within each database, same
    # per-database restriction logic as select_pathways_global() in section 8
    # (redundancy status itself is unchanged -- still pooled across contrasts
    # per section 3c's rationale; only the RANKING here is contrast-specific).
    dplyr::group_by(database) %>%
    dplyr::group_modify(~ {
      db_name <- as.character(.y$database)
      rep_ids <- reps %>% filter(database == db_name) %>% pull(ID)
      if (length(rep_ids) == 0) return(.x)
      dplyr::filter(.x, ID %in% rep_ids)
    }) %>%
    dplyr::ungroup() %>%
    filter(p.adjust < adjp_thresh) %>%
    arrange(p.adjust)
  
  if (nrow(df) == 0) {
    warning("No pathways pass adjp_thresh = ", adjp_thresh,
            " in any database for ", tissue_name, " / ", contrast_name)
    return(df %>% dplyr::select(database, Description_clean))
  }
  
  df %>% slice_head(n = top_n) %>% dplyr::select(database, Description_clean)
}

# Pulls ALL THREE contrasts' rows for the selected pathway set, not just the
# selecting contrast's row -- this is what keeps the cross-contrast columns
# populated in the resulting plot. `databases` is threaded through to the
# selection step AND to the base gsea_supp_full filter below, so the plot
# never shows a non-GO database row even if a GO-selected pathway happened
# to share a cleaned label with something outside the restricted set.
build_bubble_selection_per_contrast <- function(tissue_name, contrast_name, databases = NULL) {
  keep <- select_pathways_per_contrast(tissue_name, contrast_name, TOP_N_PER_CONTRAST, databases = databases)
  df <- gsea_supp_full %>% filter(tissue == tissue_name)
  if (!is.null(databases)) df <- df %>% filter(as.character(database) %in% databases)
  df %>% semi_join(keep, by = c("database", "Description_clean"))
}

# Visually identical to make_gsea_bubble_plot_alldb() (same pathway_key
# scheme to avoid duplicate-label factor collisions across databases, same
# aesthetics/legend), differing only in title and which selection function
# feeds it. `databases`/`title_suffix` let section 9c reuse this unchanged
# for the GO-only variant rather than duplicating the whole plot function.
make_gsea_bubble_plot_percontrast <- function(tissue_name, contrast_name, databases = NULL, title_suffix = "") {
  sel <- build_bubble_selection_per_contrast(tissue_name, contrast_name, databases = databases)
  
  present_dbs <- unname(db_slots)[unname(db_slots) %in% unique(as.character(sel$database))]
  
  sel <- sel %>%
    mutate(pathway_key = paste(database, ID, sep = "::"))
  
  ordered_levels <- unlist(lapply(present_dbs, function(db) {
    order_pathways(filter(sel, database == db), key_col = "pathway_key")
  }))
  
  label_lookup <- setNames(sel$Description_clean, sel$pathway_key)
  
  sel <- sel %>%
    mutate(
      database    = factor(as.character(database), levels = present_dbs),
      pathway_key = factor(pathway_key, levels = ordered_levels),
      significant = p.adjust < adjp_thresh
    )
  
  n_pathways <- length(ordered_levels)
  
  p <- ggplot(sel, aes(x = comparison, y = pathway_key)) +
    geom_point(aes(size = -log10(p.adjust), fill = NES,
                   color = significant, stroke = significant),
               shape = 21) +
    scale_fill_gradient2(
      low = NES_LOW, mid = colors$white, high = NES_HIGH, midpoint = 0,
      name = "NES"
    ) +
    scale_color_manual(
      values = c(`TRUE` = colors$black, `FALSE` = colors$gray_light),
      labels = c(`TRUE` = paste0("FDR < ", adjp_thresh),
                 `FALSE` = paste0("FDR \u2265 ", adjp_thresh)),
      name = "Significance"
    ) +
    scale_discrete_manual(
      aesthetics = "stroke",
      values = c(`TRUE` = 1, `FALSE` = 0.3),
      guide = "none"
    ) +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1, 7)) +
    scale_x_discrete(labels = comp_labels) +
    scale_y_discrete(labels = label_lookup) +
    facet_grid(rows = vars(database), scales = "free_y", space = "free_y") +
    labs(title = paste0("GSEA \u2014 ", tissue_name, " \u2014 Top ", TOP_N_PER_CONTRAST,
                        " for ", comp_labels[[contrast_name]],
                        " (deduplicated, pooled across databases)", title_suffix),
         x = NULL, y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(
      axis.text.y      = element_text(size = 7, face = "plain"),
      axis.text.x      = element_text(size = 8, face = "bold"),
      legend.position  = "right",
      panel.grid.major = element_line(color = colors$gray_light),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = colors$gray_light, color = NA),
      strip.text.y     = element_text(angle = 0)
    )
  
  list(plot = p, n_pathways = n_pathways)
}

# --- 9b: Build, save, and export companion data -------------------------------
percontrast_selection_data <- purrr::map_dfr(c("Muscle", "Plasma"), function(tissue_name) {
  purrr::map_dfr(gsea_comps, function(contrast_name) {
    df <- build_bubble_selection_per_contrast(tissue_name, contrast_name)
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(selected_for = comp_labels[[contrast_name]]) %>%
      dplyr::select(selected_for, tissue, database, ID,
                    Description = Description_clean,
                    comparison = comparison_label,
                    NES, pvalue, p.adjust, setSize, significant)
  })
})

writexl::write_xlsx(
  list(
    Muscle = filter(percontrast_selection_data, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma = filter(percontrast_selection_data, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection.xlsx",
                                            TOP_N_PER_CONTRAST))
)
cat("\nPer-contrast selection data saved to:",
    file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection.xlsx",
                                       TOP_N_PER_CONTRAST)), "\n")

walk(c("Muscle", "Plasma"), function(tissue_name) {
  walk(gsea_comps, function(contrast_name) {
    res <- make_gsea_bubble_plot_percontrast(tissue_name, contrast_name)
    
    if (res$n_pathways == 0) {
      cat("\n", tissue_name, "/", contrast_name,
          " per-contrast top-", TOP_N_PER_CONTRAST,
          " plot skipped -- 0 pathways passed adjp_thresh.\n", sep = "")
      return(invisible(NULL))
    }
    
    plot_height <- max(4, 0.22 * res$n_pathways + 1.5)
    file_stub <- sprintf("gsea_bubble_%s_%s_top%d_percontrast",
                         tolower(tissue_name), contrast_name, TOP_N_PER_CONTRAST)
    
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
           res$plot, width = 12, height = plot_height, device = cairo_pdf)
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
           res$plot, width = 12, height = plot_height, dpi = 300, bg = "white")
    
    cat(" Saved: ", file_stub, ".pdf/.png (", res$n_pathways, " pathways)\n", sep = "")
  })
})

# --- 9c: Per-Contrast Top-N Bubble Plots, GO:BP/CC/MF only -------------------
# Same logic as section 9, restricted to GO_DBS (GO:BP, GO:CC, GO:MF) via
# the `databases` argument added to the functions above. A separate file
# stub/companion workbook so these don't overwrite the all-database
# outputs from 9b -- both variants are kept side by side.
percontrast_selection_data_go <- purrr::map_dfr(c("Muscle", "Plasma"), function(tissue_name) {
  purrr::map_dfr(gsea_comps, function(contrast_name) {
    df <- build_bubble_selection_per_contrast(tissue_name, contrast_name, databases = GO_DBS)
    if (nrow(df) == 0) return(NULL)
    df %>%
      mutate(selected_for = comp_labels[[contrast_name]]) %>%
      dplyr::select(selected_for, tissue, database, ID,
                    Description = Description_clean,
                    comparison = comparison_label,
                    NES, pvalue, p.adjust, setSize, significant)
  })
})

writexl::write_xlsx(
  list(
    Muscle = filter(percontrast_selection_data_go, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma = filter(percontrast_selection_data_go, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection_GOonly.xlsx",
                                            TOP_N_PER_CONTRAST))
)
cat("\nPer-contrast GO:BP/CC/MF-only selection data saved to:",
    file.path(config$data_dir, sprintf("gsea_bubble_percontrast_top%d_selection_GOonly.xlsx",
                                       TOP_N_PER_CONTRAST)), "\n")

walk(c("Muscle", "Plasma"), function(tissue_name) {
  walk(gsea_comps, function(contrast_name) {
    res <- make_gsea_bubble_plot_percontrast(tissue_name, contrast_name,
                                             databases = GO_DBS,
                                             title_suffix = " [GO:BP/CC/MF only]")
    
    if (res$n_pathways == 0) {
      cat("\n", tissue_name, "/", contrast_name,
          " GO-only per-contrast top-", TOP_N_PER_CONTRAST,
          " plot skipped -- 0 pathways passed adjp_thresh in GO:BP/CC/MF.\n", sep = "")
      return(invisible(NULL))
    }
    
    plot_height <- max(4, 0.22 * res$n_pathways + 1.5)
    file_stub <- sprintf("gsea_bubble_%s_%s_top%d_percontrast_GOonly",
                         tolower(tissue_name), contrast_name, TOP_N_PER_CONTRAST)
    
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".pdf")),
           res$plot, width = 12, height = plot_height, device = cairo_pdf)
    ggsave(file.path(config$fig_dir, paste0(file_stub, ".png")),
           res$plot, width = 12, height = plot_height, dpi = 300, bg = "white")
    
    cat(" Saved: ", file_stub, ".pdf/.png (", res$n_pathways, " pathways)\n", sep = "")
  })
})

# --- Note on plasma ---------------------------------------------------------
# Given the sparse plasma GSEA results discussed earlier, this bubble format
# (which selects on adjp_thresh significance) may return very few or zero
# pathways for plasma using the same thresholds. Two options if/when you
# extend this to plasma:
#   1. Relax TOP_N to show top-ranked pathways by NES regardless of
#      significance, clearly labeled as such (matches earlier framing of
#      sparse plasma results as a finding rather than noise).
#   2. Keep plasma out of the bubble plot and report it as a compact
#      supplementary table instead (all databases, all contrasts).
#   3. Or, add tissue as a second facet dimension
#      (facet_grid(rows = vars(database), cols = vars(tissue))) once plasma
#      has enough pathways to be worth showing side-by-side with muscle.

# --- 10: 3x3 Composite — GSEA lollipop + ORA up/down bars, per contrast -----
# One 3-column x 3-row figure per tissue:
#   column = contrast (gsea_comps order: Cancer, Residual,
#            Training)
#   row 1  = GSEA lollipop, top 8 pathways pooled across ALL db_slots
#            databases, ranked by p.adjust regardless of NES sign (up and
#            down compete for the same 8 slots, not 8+8)
#   row 2  = ORA bar chart, top 8 pooled across GO:BP/CC/MF, UP-regulated
#            protein set
#   row 3  = same, DOWN-regulated protein set
# Per your call: raw top 8 by p.adjust, no Jaccard/overlap redundancy
# filtering (section 3c) applied here — that filter stays specific to the
# bubble plots in sections 4-9.

# --- 10a: Load & pool ORA data ------------------------------------------------
ora_up_slots   <- c(ora_bp_up = "GO:BP", ora_cc_up = "GO:CC", ora_mf_up = "GO:MF")
ora_down_slots <- c(ora_bp_down = "GO:BP", ora_cc_down = "GO:CC", ora_mf_down = "GO:MF")

# Diagnostic — run once before trusting the extraction below.
cat("=== ORA diagnostic ===\n")
cat("Class of ora_bp_up:", class(gsea_muscle[[1]][["ora_bp_up"]]), "\n")
cat("Columns in ora_bp_up @result:\n")
print(colnames(as.data.frame(gsea_muscle[[1]][["ora_bp_up"]])))

make_ora_data <- function(gsea_list, tissue_name, slots) {
  bind_rows(lapply(gsea_comps, function(comp) {
    bind_rows(lapply(names(slots), function(slot) {
      res <- gsea_list[[comp]][[slot]]
      if (is.null(res)) return(NULL)
      df <- as.data.frame(res)
      if (nrow(df) == 0) return(NULL)
      df %>%
        dplyr::select(ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, Count) %>%
        dplyr::mutate(comparison = comp, database = slots[[slot]], tissue = tissue_name)
    }))
  }))
}

ora_supp_full <- bind_rows(
  make_ora_data(gsea_muscle, "Muscle", ora_up_slots)   %>% mutate(direction = "Up"),
  make_ora_data(gsea_muscle, "Muscle", ora_down_slots) %>% mutate(direction = "Down"),
  make_ora_data(gsea_plasma, "Plasma", ora_up_slots)   %>% mutate(direction = "Up"),
  make_ora_data(gsea_plasma, "Plasma", ora_down_slots) %>% mutate(direction = "Down")
) %>%
  mutate(
    comparison         = factor(comparison, levels = gsea_comps),
    database           = factor(database, levels = c("GO:BP", "GO:CC", "GO:MF")),
    tissue              = factor(tissue, levels = c("Muscle", "Plasma")),
    Description_clean  = str_to_sentence(Description),
    significant         = p.adjust < adjp_thresh
  )

cat("\n=== ORA pooled data check ===\n")
cat("Rows:", nrow(ora_supp_full), "\n")
print(table(ora_supp_full$database, ora_supp_full$direction, ora_supp_full$tissue))

# --- 10b: Top-8 selection (raw p.adjust, no redundancy filter) ---------------
# Same significance-gate-with-fallback pattern as select_pathways() /
# select_pathways_global() elsewhere in this script, for consistency: if
# nothing clears adjp_thresh, fall back to top 8 by p.adjust regardless
# (with a warning) rather than silently returning an empty panel.
select_gsea_top8 <- function(tissue_name, contrast_name) {
  df <- gsea_supp_full %>% filter(tissue == tissue_name, comparison == contrast_name)
  ranked <- df %>% filter(p.adjust < adjp_thresh) %>% arrange(p.adjust)
  if (nrow(ranked) == 0) {
    warning("No GSEA pathways pass adjp_thresh for ", tissue_name, "/", contrast_name,
            " -- falling back to top 8 by p.adjust regardless of significance.")
    ranked <- df %>% arrange(p.adjust)
  }
  ranked %>% slice_head(n = 8)
}

select_ora_top8 <- function(tissue_name, contrast_name, direction_name) {
  df <- ora_supp_full %>%
    filter(tissue == tissue_name, comparison == contrast_name, direction == direction_name)
  ranked <- df %>% filter(p.adjust < adjp_thresh) %>% arrange(p.adjust)
  if (nrow(ranked) == 0) {
    warning("No ORA (", direction_name, ") pathways pass adjp_thresh for ",
            tissue_name, "/", contrast_name,
            " -- falling back to top 8 by p.adjust regardless of significance.")
    ranked <- df %>% arrange(p.adjust)
  }
  ranked %>% slice_head(n = 8)
}

# Pads a selection out to exactly n rows with blank/NA filler so every panel
# in the 3x3 grid reserves the same vertical space -- otherwise a contrast
# with only 3-4 significant ORA-down hits produces a visibly shorter panel
# than its neighbors and the grid stops reading as a clean 3x3. Filler rows
# have NA numeric fields (geoms drop them via na.rm = TRUE below) and a
# blank display label.
pad_to_n <- function(df, label_col, n = 8) {
  cur <- nrow(df)
  if (cur >= n) return(df)
  filler <- df[rep(1, n - cur), ]
  filler[[label_col]] <- paste0("__blank_", seq_len(n - cur))
  num_cols <- names(filler)[vapply(filler, is.numeric, logical(1))]
  filler[num_cols] <- NA_real_
  bind_rows(df, filler)
}

blank_labels <- function(x) ifelse(grepl("^__blank_", x), "", x)

# --- 10c: Panel builders ------------------------------------------------------
make_gsea_lollipop <- function(tissue_name, contrast_name) {
  sel <- select_gsea_top8(tissue_name, contrast_name) %>%
    pad_to_n("Description_clean", n = 8) %>%
    arrange(desc(p.adjust)) %>%   # least-significant first -> most significant plots on top
    mutate(
      Description_clean = factor(Description_clean, levels = unique(Description_clean)),
      direction          = dplyr::case_when(is.na(NES) ~ NA_character_, NES >= 0 ~ "Up", TRUE ~ "Down")
    )
  
  label_lookup <- setNames(blank_labels(as.character(sel$Description_clean)),
                           as.character(sel$Description_clean))
  
  ggplot(sel, aes(x = NES, y = Description_clean)) +
    geom_vline(xintercept = 0, color = colors$gray_light, linewidth = 0.4) +
    geom_segment(aes(x = 0, xend = NES, yend = Description_clean, color = direction),
                 linewidth = 0.7, na.rm = TRUE) +
    geom_point(aes(fill = direction, size = -log10(p.adjust)), shape = 21,
               color = colors$black, na.rm = TRUE) +
    scale_fill_manual(values = c(Up = NES_HIGH, Down = NES_LOW), na.translate = FALSE, name = NULL) +
    scale_color_manual(values = c(Up = NES_HIGH, Down = NES_LOW), na.translate = FALSE, guide = "none") +
    scale_size_continuous(name = expression(-log[10]~"(FDR)"), range = c(1.5, 5)) +
    scale_y_discrete(labels = label_lookup) +
    labs(title = paste0("GSEA (all DBs) \u2014 ", comp_labels[[contrast_name]]), x = "NES", y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(axis.text.y = element_text(size = 7), legend.position = "right")
}

make_ora_bar <- function(tissue_name, contrast_name, direction_name) {
  bar_color <- if (direction_name == "Up") NES_HIGH else NES_LOW
  
  sel <- select_ora_top8(tissue_name, contrast_name, direction_name) %>%
    pad_to_n("Description_clean", n = 8) %>%
    arrange(desc(p.adjust)) %>%
    mutate(Description_clean = factor(Description_clean, levels = unique(Description_clean)))
  
  label_lookup <- setNames(blank_labels(as.character(sel$Description_clean)),
                           as.character(sel$Description_clean))
  
  ggplot(sel, aes(x = -log10(p.adjust), y = Description_clean)) +
    geom_col(fill = bar_color, na.rm = TRUE) +
    scale_y_discrete(labels = label_lookup) +
    labs(title = paste0("ORA ", direction_name, " (GO:BP/CC/MF) \u2014 ", comp_labels[[contrast_name]]),
         x = expression(-log[10]~"(FDR)"), y = NULL) +
    theme_cr(base_size = 9, title_hjust = 0) +
    theme(axis.text.y = element_text(size = 7))
}

# --- 10d: Assemble 3x3 composite, per tissue ----------------------------------
build_composite_figure <- function(tissue_name) {
  columns <- lapply(gsea_comps, function(contrast_name) {
    make_gsea_lollipop(tissue_name, contrast_name) /
      make_ora_bar(tissue_name, contrast_name, "Up") /
      make_ora_bar(tissue_name, contrast_name, "Down")
  })
  
  wrap_plots(columns, ncol = 3) +
    plot_annotation(title = paste0("GSEA + ORA summary \u2014 ", tissue_name)) &
    theme(plot.title = element_text(face = "bold"))
}

fig_muscle_composite <- build_composite_figure("Muscle")
fig_plasma_composite <- build_composite_figure("Plasma")

ggsave(file.path(config$fig_dir, "gsea_ora_composite_muscle_3x3.pdf"),
       fig_muscle_composite, width = 16, height = 12, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_ora_composite_muscle_3x3.png"),
       fig_muscle_composite, width = 16, height = 12, dpi = 300, bg = "white")

ggsave(file.path(config$fig_dir, "gsea_ora_composite_plasma_3x3.pdf"),
       fig_plasma_composite, width = 16, height = 12, device = cairo_pdf)
ggsave(file.path(config$fig_dir, "gsea_ora_composite_plasma_3x3.png"),
       fig_plasma_composite, width = 16, height = 12, dpi = 300, bg = "white")

cat("\nSaved 3x3 composite figures to:", config$fig_dir, "\n")

# --- 10e: Companion data export -----------------------------------------------
composite_selection_data <- purrr::map_dfr(c("Muscle", "Plasma"), function(tissue_name) {
  purrr::map_dfr(gsea_comps, function(contrast_name) {
    bind_rows(
      select_gsea_top8(tissue_name, contrast_name) %>%
        transmute(tissue = tissue_name, comparison = comp_labels[[contrast_name]],
                  panel = "GSEA (all DBs)", database, ID, Description = Description_clean,
                  NES, pvalue, p.adjust, setSize),
      select_ora_top8(tissue_name, contrast_name, "Up") %>%
        transmute(tissue = tissue_name, comparison = comp_labels[[contrast_name]],
                  panel = "ORA Up", database, ID, Description = Description_clean,
                  GeneRatio, pvalue, p.adjust, Count),
      select_ora_top8(tissue_name, contrast_name, "Down") %>%
        transmute(tissue = tissue_name, comparison = comp_labels[[contrast_name]],
                  panel = "ORA Down", database, ID, Description = Description_clean,
                  GeneRatio, pvalue, p.adjust, Count)
    )
  })
})

writexl::write_xlsx(
  list(
    Muscle = filter(composite_selection_data, tissue == "Muscle") %>% dplyr::select(-tissue),
    Plasma = filter(composite_selection_data, tissue == "Plasma") %>% dplyr::select(-tissue)
  ),
  path = file.path(config$data_dir, "gsea_ora_composite_3x3_selection.xlsx")
)
cat("Composite selection data saved to:",
    file.path(config$data_dir, "gsea_ora_composite_3x3_selection.xlsx"), "\n")
