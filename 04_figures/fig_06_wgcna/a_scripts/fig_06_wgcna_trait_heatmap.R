# ============================================================
# fig_wgcna_trait_heatmap.R
#
# Standalone module eigengene - trait correlation heatmap, separate
# from the main fig_06 composite. Four trait families in one figure,
# grouped as facet columns:
#
#   1. Baseline (PRE) ME vs baseline (pre_) traits -- CRC
#      -- Uses each subject's CRC_T1 (PRE-timepoint) module eigengene
#      value specifically -- not a delta, not the full mixed sample set
#      -- correlated against the pre_<name> metadata columns (subject-
#      level constants, confirmed identical across a subject's two rows).
#      CRC-only because the pre_/post_ column structure only exists for
#      the two-visit CRC arm.
#
#   1b. Delta ME vs delta (post - pre) traits -- CRC
#      -- Computed directly from the SAME resolved trait set as the
#      Baseline (CRC) panel above (not reused from
#      rd$delta_trait_long/05_wgcna_muscle.R section 7b, which was built
#      against a separate, smaller delta_pheno_cols list defined in that
#      other script -- that mismatch is why this panel used to only ever
#      show 4 variables regardless of config$baseline_pheno_cols here).
#      Computing it locally, from the same resolved column list, means
#      the Delta and Baseline (CRC) panels always cover the exact same
#      trait set by construction. A trait needs BOTH a pre_<name> and a
#      post_<name> column with non-NA data to appear here -- if
#      post_<name> is missing for a variable that did have pre_<name>
#      (so it appears in Baseline but not here), that's flagged with a
#      warning rather than silently dropped.
#
#   2. Baseline ME vs baseline traits -- CTL
#      -- CTL subjects have a single visit, so there's no pre_/post_
#      pair to speak of, but they DO have their own baseline eigengene +
#      trait values, which are worth correlating on their own: this tells
#      you whether an ME-trait association is specific to the cancer/
#      training context or shows up in healthy controls too.
#      ASSUMPTION (verify against your metadata): CTL traits are assumed
#      to live under the unprefixed column name (e.g. "fat_kg"); the
#      resolver below falls back to "pre_<name>" if that's actually
#      where a CTL subject's value is, and only uses a candidate column
#      if it has non-NA data for CTL subjects.
#
#   3. Cancer-specific outcomes (stage, yrs_since_dx, yrs_since_chemo)
#      vs BASELINE (CRC_T1) ME only -- these are subject-level constants
#      duplicated across a subject's CRC_T1/CRC_T2 rows, so including
#      CRC_T2 would pseudoreplicate each subject and mix a fixed
#      baseline characteristic with an already-training-affected POST
#      eigengene. CTL rows are structurally excluded here (these
#      variables are undefined, not just missing, for someone never
#      diagnosed) -- not a choice, a fact about what the variables mean.
#      Uses SPEARMAN correlation (stage is ordinal; yrs_since_dx/chemo
#      are typically right-skewed) rather than Pearson, unlike panels
#      1/1b/2. The Student's-t p-value from corPvalueStudent() is an
#      approximation when applied to a rank correlation -- treat panel 3
#      p-values as more approximate than the others, especially at low n.
#
# All correlations use WGCNA::cor(use="pairwise.complete.obs") + PER-
# COLUMN nSamples for corPvalueStudent (not a single matrix-wide
# nSamples), so a trait with more missingness than another doesn't
# borrow the full matrix's sample size for its p-value.
#
# FDR is computed ONCE, after all four panels are combined, with the
# scope controlled by config$fdr_scope:
#   "per_category" -- one BH correction per panel (closest to the
#     original behavior of this script, now made explicit/configurable
#     instead of being buried three different ways across panels).
#   "global" -- one BH correction across every test in the whole
#     heatmap, so a cell's significance doesn't depend on which panel it
#     happened to land in. Recommended now that there are four panels of
#     different sizes.
#
# Cells built on fewer than config$min_n pairwise-complete observations
# get a "\u2020" marker appended to their label (in addition to the FDR
# outline) -- a correlation from n=4 shouldn't read identically to one
# from n=40 just because both happen to be significant.
# ============================================================

pacman::p_load(tidyverse, WGCNA, writexl, rprojroot, ggnewscale, ggh4x)

root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "04_figures", "a_scripts", "00_theme.R"))

# Display names for the four facet panels -- a NAMED vector so every
# block-building call below references CATEGORY_LEVELS[["delta"]] etc.
# instead of retyping the literal string. Rename a panel here and it
# propagates everywhere (block category tags, the facet strips, and the
# x-axis trait_key construction) automatically -- previously these
# strings were duplicated in five separate places and had to be kept in
# sync by hand, which is exactly the kind of thing that silently breaks
# (a mismatched string here makes bind_rows -> factor(category, levels=)
# turn that whole panel's rows into NA, and ggplot drops NA facets
# without erroring).
CATEGORY_LEVELS <- c(
  module_swatch = "",   # blank facet strip -- see module_color_df below
  delta         = "\u0394 Trait (POST - PRE)",
  baseline_crc  = "Baseline Traits (CRC)",
  baseline_ctl  = "Baseline Traits (CRC & CTL)",
  cancer        = "Cancer Outcomes"
)

config <- list(
  tissue       = "muscle",
  wgcna_bundle = file.path(root, "05_wgcna", "report_data_muscle.rds"),
  out_dir      = file.path(root, "04_figures", "fig_06_wgcna"),
  
  # Baseline (PRE) trait base names -- also drives the Delta panel now
  # (see section 1b), since Delta reuses this exact resolved list rather
  # than an independently-defined delta_pheno_cols. Also used for the
  # CTL baseline panel, since CTL should be compared on the same trait
  # set as CRC's baseline panel to be a meaningful contrast.
  baseline_pheno_cols = c("fat_kg", "fat_pct", "ALM_kg", "ALM_index", "BMI", 
                          "leg_ext_lbs", "chest_press_lbs", "grip_lbs"),
  
  # Cancer-specific outcome variables, raw values.
  cancer_outcome_cols  = c("stage", "yrs_since_dx", "yrs_since_chemo"),
  
  # Label used for the CTL group in group_col. Verify this matches your
  # actual metadata coding -- if rd$me_long[[group_col]] doesn't contain
  # this value, the CTL baseline panel will warn and come back empty
  # rather than silently mis-grouping subjects.
  ctl_group_label = "CTL_T1",
  
  # "per_category" (one BH correction per panel, closest to this
  # script's original per-block behavior) or "global" (one BH correction
  # across the whole combined heatmap -- recommended).
  fdr_scope  = "per_category",
  fdr_method = "BH",
  
  # Cells with fewer than this many pairwise-complete observations get a
  # "\u2020" marker on their label, regardless of significance.
  min_n = 5,
  
  # Display labels for the x-axis, keyed on the BASE trait name (i.e.
  # without the "\u0394 " prefix used in the delta panel -- that panel
  # adds "\u0394 " in front of whatever this map returns, so one entry
  # covers a trait in every panel it appears in). Any trait not listed
  # here falls back to its raw name with underscores turned to spaces
  # and title-cased (see pretty_trait_label() below) -- so adding a new
  # baseline_pheno_col or cancer_outcome_col without updating this map
  # won't break anything, it'll just look less polished until you add
  # an entry for it.
  trait_label_map = c(
    fat_kg          = "Fat Mass (kg)",
    fat_pct         = "Body Fat %",
    ALM_kg          = "ALM (kg)",
    ALM_index       = "ALM Index",
    BMI             = "BMI",
    leg_ext_lbs     = "Leg Extension (lbs)",
    chest_press_lbs = "Chest Press (lbs)",
    grip_lbs        = "Grip Strength (lbs)",
    stage           = "Cancer Stage",
    yrs_since_dx    = "Years Since Diagnosis",
    yrs_since_chemo = "Years Since Chemo"
  )
)
dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)
# NOTE: added -- pdf_path/ggsave/write_xlsx below write into
# file.path(config$out_dir, "b_figures"/"c_data", ...), but only
# out_dir itself was being created. Without these, cairo_pdf() will
# error on a fresh run since it doesn't create missing directories.
dir.create(file.path(config$out_dir, "b_figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(config$out_dir, "c_data"),    recursive = TRUE, showWarnings = FALSE)

# Resolves a raw trait name (possibly "\u0394 "-prefixed) to a display
# label via config$trait_label_map, falling back to a title-cased,
# underscore-stripped version of the name if it isn't in the map.
pretty_trait_label <- function(trait, label_map = config$trait_label_map) {
  is_delta <- startsWith(trait, "\u0394 ")
  base     <- ifelse(is_delta, sub("^\u0394 ", "", trait), trait)
  fallback <- stringr::str_to_title(gsub("_", " ", base))
  pretty   <- dplyr::coalesce(unname(label_map[base]), fallback)
  ifelse(is_delta, paste0("\u0394 ", pretty), pretty)
}

# ---- 0. Load data ----
rd <- readRDS(config$wgcna_bundle)

subject_col <- rd$config$subject_col
group_col   <- rd$config$group_col

module_order <- rd$module_membership %>%
  dplyr::filter(module != "grey") %>%
  dplyr::count(module) %>%
  dplyr::arrange(dplyr::desc(n)) %>%
  dplyr::pull(module)

# me_long already carries every metadata column joined in (built in
# 05_wgcna_muscle.R as MEs -> left_join(meta) -> pivot_longer), so it
# doubles as the metadata source here -- no separate meta object needed.
meta_by_sample  <- rd$me_long %>% dplyr::distinct(sample_id, .keep_all = TRUE)
meta_by_subject <- rd$me_long %>% dplyr::distinct(.data[[subject_col]], .keep_all = TRUE)

# Shared helper: correlate an ME matrix against a trait matrix (same
# rownames on both), with PER-COLUMN nSamples for corPvalueStudent -- a
# trait column with more missingness than another must not borrow the
# full-matrix sample size for its p-value (see 05_wgcna_muscle.R sections
# 7/7b for the original bug this fixes).
#
# Returns r/p/n only -- NOT fdr/sig. FDR is computed once, after all
# panels are combined (see section 4), so the correction scope is a
# single deliberate choice (config$fdr_scope) instead of three-plus
# independent corrections quietly baked in at different stages.
correlate_me_traits <- function(me_mat, trait_mat, category_label, cor_method = "pearson") {
  cor_mat <- WGCNA::cor(me_mat, trait_mat, use = "pairwise.complete.obs", method = cor_method)
  n_per_trait <- sapply(trait_mat, function(x) sum(!is.na(x)))
  pval_mat <- cor_mat
  for (tr in colnames(cor_mat)) {
    pval_mat[, tr] <- WGCNA::corPvalueStudent(cor_mat[, tr], nSamples = n_per_trait[tr])
  }
  
  tibble::as_tibble(cor_mat, rownames = "module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "r") %>%
    dplyr::left_join(
      tibble::as_tibble(pval_mat, rownames = "module") %>%
        tidyr::pivot_longer(-module, names_to = "trait", values_to = "p"),
      by = c("module", "trait")
    ) %>%
    dplyr::mutate(n = n_per_trait[trait], category = category_label)
}

empty_block <- function() {
  tibble::tibble(module = character(0), trait = character(0), r = double(0),
                 p = double(0), n = integer(0), category = character(0))
}

# ============================================================
# 1. Baseline (PRE) ME vs baseline (pre_) traits -- CRC
# ============================================================
me_baseline_wide <- rd$me_long %>%
  dplyr::filter(.data[[group_col]] == "CRC_T1", module %in% module_order) %>%
  dplyr::distinct(.data[[subject_col]], module, .keep_all = TRUE) %>%
  dplyr::select(subject = dplyr::all_of(subject_col), module, eigengene) %>%
  tidyr::pivot_wider(names_from = module, values_from = eigengene)

present_pre_cols <- paste0("pre_", config$baseline_pheno_cols)
missing_pre_cols <- setdiff(present_pre_cols, colnames(meta_by_subject))
if (length(missing_pre_cols) > 0) {
  warning("baseline_pheno_cols: couldn't find these columns in metadata, skipping: ",
          paste(missing_pre_cols, collapse = ", "))
}
present_pre_cols <- intersect(present_pre_cols, colnames(meta_by_subject))

if (length(present_pre_cols) == 0) {
  warning("No valid baseline_pheno_cols found -- CRC baseline block will be empty.")
  baseline_block       <- empty_block()
  baseline_trait_names  <- character(0)
  baseline_traits       <- tibble::tibble(subject = character(0))
} else {
  baseline_traits <- meta_by_subject %>%
    dplyr::select(subject = dplyr::all_of(subject_col), dplyr::all_of(present_pre_cols)) %>%
    dplyr::rename_with(~ gsub("^pre_", "", .x), .cols = dplyr::all_of(present_pre_cols))
  
  baseline_merged <- dplyr::inner_join(me_baseline_wide, baseline_traits, by = "subject")
  
  baseline_me_mat    <- baseline_merged %>% dplyr::select(dplyr::all_of(module_order)) %>% as.data.frame()
  baseline_trait_names <- setdiff(colnames(baseline_traits), "subject")
  baseline_trait_mat <- baseline_merged %>% dplyr::select(dplyr::all_of(baseline_trait_names)) %>% as.data.frame()
  rownames(baseline_me_mat)    <- baseline_merged$subject
  rownames(baseline_trait_mat) <- baseline_merged$subject
  
  baseline_block <- correlate_me_traits(baseline_me_mat, baseline_trait_mat, CATEGORY_LEVELS[["baseline_crc"]])
}

# ============================================================
# 1b. Delta ME vs delta (post - pre) traits -- CRC
# ============================================================
# Reuses baseline_trait_names/baseline_traits (the resolved pre_ trait
# set from section 1 above) rather than an independent variable list, so
# this panel always covers the exact same traits as Baseline (CRC) --
# see the module header for why that's the point of computing this
# locally instead of reusing rd$delta_trait_long.
if (length(baseline_trait_names) == 0) {
  warning("No baseline traits resolved -- delta block will be empty (it depends ",
          "on the same trait set as the Baseline (CRC) block above).")
  delta_block <- empty_block()
} else {
  present_post_cols <- paste0("post_", baseline_trait_names)
  missing_post_cols <- setdiff(present_post_cols, colnames(meta_by_subject))
  if (length(missing_post_cols) > 0) {
    warning("delta traits: couldn't find these post_ columns in metadata, so these ",
            "variables appear in the Baseline (CRC) panel but NOT the Delta panel: ",
            paste(gsub("^post_", "", missing_post_cols), collapse = ", "))
  }
  present_post_cols <- intersect(present_post_cols, colnames(meta_by_subject))
  
  if (length(present_post_cols) == 0) {
    warning("No usable post_<name> trait columns found -- delta block will be empty.")
    delta_block <- empty_block()
  } else {
    delta_trait_names <- gsub("^post_", "", present_post_cols)
    
    # post_<name> (renamed to <name>_post) joined against the matching
    # pre_<name> values already resolved in baseline_traits (renamed
    # <name>_pre via the join suffix), then subtracted per trait.
    delta_traits <- meta_by_subject %>%
      dplyr::select(subject = dplyr::all_of(subject_col), dplyr::all_of(present_post_cols)) %>%
      dplyr::rename_with(~ gsub("^post_", "", .x), .cols = dplyr::all_of(present_post_cols)) %>%
      dplyr::inner_join(
        baseline_traits %>% dplyr::select(subject, dplyr::all_of(delta_trait_names)),
        by = "subject", suffix = c("_post", "_pre")
      )
    
    for (nm in delta_trait_names) {
      delta_traits[[nm]] <- delta_traits[[paste0(nm, "_post")]] - delta_traits[[paste0(nm, "_pre")]]
    }
    delta_trait_mat <- delta_traits %>% dplyr::select(dplyr::all_of(delta_trait_names)) %>% as.data.frame()
    rownames(delta_trait_mat) <- delta_traits$subject
    
    me_post_wide <- rd$me_long %>%
      dplyr::filter(.data[[group_col]] == "CRC_T2", module %in% module_order) %>%
      dplyr::distinct(.data[[subject_col]], module, .keep_all = TRUE) %>%
      dplyr::select(subject = dplyr::all_of(subject_col), module, eigengene) %>%
      tidyr::pivot_wider(names_from = module, values_from = eigengene)
    
    delta_me_merged <- dplyr::inner_join(
      me_baseline_wide %>% dplyr::select(subject, dplyr::all_of(module_order)),
      me_post_wide      %>% dplyr::select(subject, dplyr::all_of(module_order)),
      by = "subject", suffix = c("_pre", "_post")
    ) %>%
      dplyr::filter(subject %in% rownames(delta_trait_mat))
    
    # Selected by exact expected name (module + suffix) rather than
    # ends_with(), so the two matrices' columns are guaranteed to line
    # up module-for-module regardless of how dplyr happened to order
    # the joined columns.
    post_mat <- as.matrix(delta_me_merged[paste0(module_order, "_post")])
    pre_mat  <- as.matrix(delta_me_merged[paste0(module_order, "_pre")])
    colnames(post_mat) <- module_order
    colnames(pre_mat)  <- module_order
    
    delta_me_mat <- as.data.frame(post_mat - pre_mat)
    rownames(delta_me_mat) <- delta_me_merged$subject
    
    delta_trait_mat <- delta_trait_mat[rownames(delta_me_mat), , drop = FALSE]
    
    delta_block <- correlate_me_traits(delta_me_mat, delta_trait_mat, CATEGORY_LEVELS[["delta"]])
  }
}

# ============================================================
# 2. Baseline ME vs baseline traits -- CTL
# ============================================================
if (!(config$ctl_group_label %in% unique(rd$me_long[[group_col]]))) {
  warning("No \"", config$ctl_group_label, "\" level found in rd$me_long[[group_col]] -- ",
          "CTL baseline block will be empty. Check config$ctl_group_label against ",
          "the actual values: ", paste(unique(rd$me_long[[group_col]]), collapse = ", "))
  ctl_baseline_block <- empty_block()
} else {
  ctl_subject_ids <- rd$me_long %>%
    dplyr::filter(.data[[group_col]] == config$ctl_group_label) %>%
    dplyr::pull(.data[[subject_col]]) %>%
    unique()
  
  # Tries the unprefixed base name first (assumed CTL convention), falls
  # back to pre_<name> -- only returns a candidate if it actually has
  # non-NA values among CTL subjects. Returns NA if neither works, and
  # that trait is silently dropped from the CTL panel (with a warning
  # below if that empties the whole panel).
  resolve_ctl_trait_col <- function(base_name) {
    for (cand in c(base_name, paste0("pre_", base_name))) {
      if (cand %in% colnames(meta_by_subject)) {
        vals <- meta_by_subject %>%
          dplyr::filter(.data[[subject_col]] %in% ctl_subject_ids) %>%
          dplyr::pull(.data[[cand]])
        if (any(!is.na(vals))) return(cand)
      }
    }
    NA_character_
  }
  
  ctl_col_map <- purrr::set_names(
    purrr::map_chr(config$baseline_pheno_cols, resolve_ctl_trait_col),
    config$baseline_pheno_cols
  )
  ctl_col_map <- ctl_col_map[!is.na(ctl_col_map)]
  
  if (length(ctl_col_map) == 0) {
    warning("No usable CTL baseline trait columns resolved (tried unprefixed and ",
            "pre_-prefixed names) -- CTL baseline block will be empty. Verify how ",
            "CTL subjects' baseline traits are actually named in your metadata.")
    ctl_baseline_block <- empty_block()
  } else {
    me_ctl_wide <- rd$me_long %>%
      dplyr::filter(.data[[group_col]] == config$ctl_group_label, module %in% module_order) %>%
      dplyr::distinct(.data[[subject_col]], module, .keep_all = TRUE) %>%
      dplyr::select(subject = dplyr::all_of(subject_col), module, eigengene) %>%
      tidyr::pivot_wider(names_from = module, values_from = eigengene)
    
    ctl_traits <- meta_by_subject %>%
      dplyr::filter(.data[[subject_col]] %in% ctl_subject_ids) %>%
      dplyr::select(subject = dplyr::all_of(subject_col), dplyr::all_of(unname(ctl_col_map))) %>%
      dplyr::rename_with(~ names(ctl_col_map)[match(.x, ctl_col_map)],
                         .cols = dplyr::all_of(unname(ctl_col_map)))
    
    ctl_merged <- dplyr::inner_join(me_ctl_wide, ctl_traits, by = "subject")
    
    ctl_me_mat      <- ctl_merged %>% dplyr::select(dplyr::all_of(module_order)) %>% as.data.frame()
    ctl_trait_names <- names(ctl_col_map)
    ctl_trait_mat   <- ctl_merged %>% dplyr::select(dplyr::all_of(ctl_trait_names)) %>% as.data.frame()
    rownames(ctl_me_mat)    <- ctl_merged$subject
    rownames(ctl_trait_mat) <- ctl_merged$subject
    
    ctl_baseline_block <- correlate_me_traits(ctl_me_mat, ctl_trait_mat, CATEGORY_LEVELS[["baseline_ctl"]])
  }
}

# ============================================================
# 3. Cancer-specific outcomes vs baseline (CRC_T1) ME
# ============================================================
# stage/yrs_since_dx/yrs_since_chemo are subject-level constants,
# duplicated across a subject's CRC_T1 and CRC_T2 rows -- same pattern
# already confirmed for the pre_/post_ columns. Restricted to CRC_T1
# only, same as the Baseline block above: including CRC_T2 would
# pseudoreplicate each subject (same trait value paired against two
# different ME values) and confound a fixed baseline characteristic with
# a POST-training eigengene that's already been shaped by the
# intervention. CTL is structurally excluded: these variables are
# undefined (not just missing) for a subject who was never diagnosed.
#
# Uses SPEARMAN (cor_method below), not Pearson: stage is ordinal and
# yrs_since_dx/yrs_since_chemo are typically right-skewed. See the
# module-level header note on why panel 3's p-values should be read as
# more approximate than the other panels'.
present_cancer_cols <- intersect(config$cancer_outcome_cols, colnames(meta_by_sample))
missing_cancer_cols <- setdiff(config$cancer_outcome_cols, colnames(meta_by_sample))
if (length(missing_cancer_cols) > 0) {
  warning("cancer_outcome_cols: couldn't find these columns in metadata, skipping: ",
          paste(missing_cancer_cols, collapse = ", "))
}

if (length(present_cancer_cols) == 0) {
  warning("No valid cancer_outcome_cols found -- cancer outcomes block will be empty.")
  cancer_block <- empty_block()
} else {
  baseline_sample_ids <- meta_by_sample %>%
    dplyr::filter(.data[[group_col]] == "CRC_T1") %>%
    dplyr::pull(sample_id)
  
  me_wide_baseline <- rd$MEs[rownames(rd$MEs) %in% baseline_sample_ids, module_order, drop = FALSE]
  
  cancer_trait_mat <- meta_by_sample %>%
    dplyr::filter(sample_id %in% rownames(me_wide_baseline)) %>%
    dplyr::arrange(match(sample_id, rownames(me_wide_baseline))) %>%
    dplyr::select(dplyr::all_of(present_cancer_cols)) %>%
    as.data.frame()
  rownames(cancer_trait_mat) <- rownames(me_wide_baseline)
  
  cancer_block <- correlate_me_traits(as.data.frame(me_wide_baseline), cancer_trait_mat,
                                      CATEGORY_LEVELS[["cancer"]], cor_method = "spearman")
}

# ============================================================
# 4. Combine + FDR + plot
# ============================================================

# Module color swatch: a one-column facet panel placed FIRST (leftmost)
# via CATEGORY_LEVELS[["module_swatch"]] = "", so it renders immediately
# to the right of the shared y-axis module labels -- facet_grid only
# draws axis text on the leftmost panel here since y isn't a free scale.
# ASSUMPTION: WGCNA module names ARE valid R color strings by convention
# (that's what labels2colors() assigns them for), so the module's own
# name doubles as its fill color with no separate color map needed. If
# any module names here were customized away from that convention,
# scale_fill_identity() below will error at render time naming the bad
# value, rather than silently drawing a wrong color.
module_color_df <- tibble::tibble(
  module       = factor(module_order, levels = rev(module_order)),
  module_color = module_order,
  category     = factor(CATEGORY_LEVELS[["module_swatch"]], levels = unname(CATEGORY_LEVELS)),
  trait_key    = paste(CATEGORY_LEVELS[["module_swatch"]], "swatch", sep = "::")
)

# trait_key (category::trait), not raw trait name, is the actual x
# factor -- the same variable name can legitimately appear in more than
# one category (baseline_pheno_cols is shared across the CRC and CTL
# baseline panels by design, so they're comparable), and facet_grid
# showing "fat_kg" as a valid tick in two different facet panels is
# fine, but a single global factor() on raw trait name isn't, since
# factor() requires unique level strings across the whole column even
# though only one category's facet ever displays a given tick at once.
trait_key_levels <- c(
  unique(module_color_df$trait_key),
  paste(CATEGORY_LEVELS[["delta"]],        unique(delta_block$trait),         sep = "::"),
  paste(CATEGORY_LEVELS[["baseline_crc"]], unique(baseline_block$trait),      sep = "::"),
  paste(CATEGORY_LEVELS[["baseline_ctl"]], unique(ctl_baseline_block$trait),  sep = "::"),
  paste(CATEGORY_LEVELS[["cancer"]],       unique(cancer_block$trait),        sep = "::")
)

all_blocks <- dplyr::bind_rows(delta_block, baseline_block, ctl_baseline_block, cancer_block) %>%
  dplyr::mutate(category = factor(category, levels = unname(CATEGORY_LEVELS)))

# FDR computed once, here, across whichever scope config$fdr_scope
# specifies -- see the module header for what each option means.
if (config$fdr_scope == "global") {
  all_blocks <- all_blocks %>%
    dplyr::mutate(fdr = p.adjust(p, method = config$fdr_method))
} else if (config$fdr_scope == "per_category") {
  all_blocks <- all_blocks %>%
    dplyr::group_by(category) %>%
    dplyr::mutate(fdr = p.adjust(p, method = config$fdr_method)) %>%
    dplyr::ungroup()
} else {
  stop("config$fdr_scope must be \"per_category\" or \"global\", got: ", config$fdr_scope)
}

all_blocks <- all_blocks %>%
  dplyr::mutate(
    sig       = p < 0.05,
    low_n     = n < config$min_n,
    module    = factor(module, levels = rev(module_order)),
    trait_key = paste(category, trait, sep = "::"),
    trait_key = factor(trait_key, levels = trait_key_levels),
    label     = paste0(sprintf("%.2f\n(%.3f)", r, p), ifelse(low_n, "\n\u2020", ""))
  )

# Delta-panel traits are stored under their BASE name (e.g. "fat_kg") so
# they share trait_label_map entries with the baseline panels, but what's
# actually correlated there is post - pre -- prefix "Δ " for that
# panel only, which pretty_trait_label() then carries through to the
# display label (e.g. "Δ Fat Mass (kg)").
trait_key_labels <- all_blocks %>%
  dplyr::distinct(trait_key, trait, category) %>%
  dplyr::mutate(
    label_input   = dplyr::if_else(category == CATEGORY_LEVELS[["delta"]],
                                   paste0("Δ ", trait), trait),
    display_label = pretty_trait_label(label_input)
  ) %>%
  dplyr::select(trait_key, display_label) %>%
  tibble::deframe()

# Blank tick label for the swatch column -- no trait name to show there.
trait_key_labels <- c(setNames("", unique(module_color_df$trait_key)), trait_key_labels)

message(sprintf(
  paste("Trait heatmap: %d delta rows, %d CRC-baseline rows, %d CTL-baseline rows,",
        "%d cancer-outcome rows (%d modules x %d trait slots total, fdr_scope = \"%s\")"),
  nrow(delta_block), nrow(baseline_block), nrow(ctl_baseline_block), nrow(cancer_block),
  length(module_order), length(trait_key_levels), config$fdr_scope
))

p_trait_heat <- ggplot(all_blocks, aes(x = trait_key, y = module, fill = r)) +
  geom_tile(color = "white") +
  geom_tile(data = dplyr::filter(all_blocks, sig), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = label, fontface = ifelse(sig, "bold", "plain")), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red,
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  # module color swatch -- separate fill scale (new_scale_fill()) so its
  # identity colors don't get folded into the r-value gradient above
  ggnewscale::new_scale_fill() +
  geom_tile(data = module_color_df, aes(x = trait_key, y = module, fill = module_color),
            color = "white", inherit.aes = FALSE) +
  scale_fill_identity(guide = "none") +
  scale_x_discrete(labels = trait_key_labels) +
  # scales/space = "free_x" alone can only make the swatch as narrow as
  # ONE data column -- to go thinner, panel widths are forced explicitly
  # below via ggh4x, so "space" here is just for independent per-panel
  # x-ranges/breaks, not sizing.
  facet_grid(. ~ category, scales = "free_x") +
  labs(x = NULL, y = NULL,
       title = paste("Module eigengene \u2013 trait correlations \u2014", config$tissue),
       subtitle = NULL) +
  theme_cr() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold")
  )

# Force explicit panel widths: the swatch gets SWATCH_WIDTH_FRACTION (a
# fraction of one normal data column's width -- lower it for an even
# thinner strip, raise it if the module colors feel cramped), and every
# other panel keeps the same relative width facet_grid's
# space = "free_x" was giving it before (proportional to its own number
# of trait columns).
#
# IMPORTANT: force_panelsizes() matches its `cols` vector to panels
# POSITIONALLY, so the vector must have exactly one entry per panel that
# actually gets drawn, in the same left-to-right order. Building that
# list from the block variables independently (assuming all five
# categories always render, in a fixed assumed order) is fragile: if any
# one block came back empty, facet_grid's default drop = TRUE means that
# panel doesn't exist at all, which shifts every width after it by one
# slot -- exactly the kind of thing that silently breaks a panel further
# down the line (e.g. "Cancer Outcomes" going missing) without erroring.
#
# So instead this counts trait_key columns directly from the data that's
# actually about to be plotted (all_blocks + module_color_df), grouped
# and ordered by the category factor's own level order -- which is the
# same order facet_grid uses for panels -- so the vector's length and
# order always match reality, however many categories actually ended up
# with data.
SWATCH_WIDTH_FRACTION <- 0.35

panel_trait_counts <- dplyr::bind_rows(
  all_blocks      %>% dplyr::distinct(category, trait_key),
  module_color_df %>% dplyr::distinct(category, trait_key)
) %>%
  dplyr::mutate(category = droplevels(category)) %>%
  dplyr::count(category, name = "n_cols") %>%
  dplyr::arrange(category)

message("Panels actually being drawn (left to right): ",
        paste(sprintf("%s (%d col%s)", panel_trait_counts$category, panel_trait_counts$n_cols,
                      ifelse(panel_trait_counts$n_cols == 1, "", "s")),
              collapse = ", "))

panel_widths <- ifelse(
  as.character(panel_trait_counts$category) == CATEGORY_LEVELS[["module_swatch"]],
  SWATCH_WIDTH_FRACTION,
  panel_trait_counts$n_cols
)

p_trait_heat <- p_trait_heat + ggh4x::force_panelsizes(cols = panel_widths)

pdf_path <- file.path(config$out_dir, "b_figures", paste0("wgcna_trait_heatmap_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = max(10, length(trait_key_levels) * 0.8), height = max(6, length(module_order) * 0.5))
print(p_trait_heat)
dev.off()

ggsave(file.path(config$out_dir, "b_figures", paste0("wgcna_trait_heatmap_", config$tissue, ".png")),
       p_trait_heat, width = max(10, length(trait_key_levels) * 0.8), height = max(6, length(module_order) * 0.5),
       dpi = 300, bg = "white")

write_xlsx(
  list(all_correlations = all_blocks %>% dplyr::select(-label, -trait_key)),
  file.path(config$out_dir, "c_data", paste0("wgcna_trait_heatmap_data_", config$tissue, ".xlsx"))
)

message("Done. Saved to: ", config$out_dir)