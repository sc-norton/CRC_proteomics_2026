# ============================================================
# fig_wgcna_trait_heatmap.R
#
# Standalone module eigengene - trait correlation heatmap, separate
# from the main fig_06 composite. Three trait families in one figure,
# grouped as facet columns:
#
#   1. Delta ME vs delta traits (post - pre change scores)
#      -- reuses rd$delta_trait_long from 05_wgcna_muscle.R section 7b.
#      Structurally CRC-only (only subjects with a SURV_T1+SURV_T2 pair
#      contribute).
#
#   2. Baseline (PRE) ME vs baseline (pre_) traits
#      -- NEW correlation computed here. Uses each subject's SURV_T1
#      (PRE-timepoint) module eigengene value specifically -- not a
#      delta, not the full mixed sample set -- correlated against the
#      pre_<name> metadata columns (which are subject-level constants,
#      confirmed identical across a subject's two rows in your actual
#      metadata). Structurally CRC-only: CTL subjects don't have this
#      pre_/post_ trial structure.
#
#   3. Cancer-specific outcomes (stage, yrs_since_dx, yrs_since_chemo)
#      vs BASELINE (SURV_T1) ME only -- these are subject-level constants
#      duplicated across a subject's SURV_T1/SURV_T2 rows, so including
#      SURV_T2 would pseudoreplicate each subject and mix a fixed
#      baseline characteristic with an already-training-affected POST
#      eigengene. Same pairwise.complete.obs + per-column nSamples
#      approach as 05_wgcna_muscle.R section 7. CTL rows drop out
#      structurally either way (NA for these variables).
#
# All three blocks use IDENTICAL statistical machinery
# (WGCNA::cor(use="pairwise.complete.obs") + per-column corPvalueStudent,
# not a single matrix-wide nSamples) so results are comparable across
# the three panels rather than three different methods bolted together.
# ============================================================

pacman::p_load(tidyverse, WGCNA, writexl, rprojroot)

root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "04_figures", "a_scripts", "00_theme.R"))

config <- list(
  tissue       = "muscle",
  wgcna_bundle = file.path(root, "05_wgcna", "report_data.rds"),
  out_dir      = file.path(root, "04_figures", "fig_wgcna_traits"),
  
  # Baseline (PRE) trait base names -- defaults to the same variables
  # used for delta_pheno_cols in 05_wgcna_muscle.R. Edit independently if
  # you want a different set of variables for baseline vs delta.
  baseline_pheno_cols = c("fat_kg", "fat_pct", "ALM_kg", "ALM_index", "BMI", 
                          "leg_ext_lbs", "chest_press_lbs", "sts_max_pwr",
                          "grip_lbs"),
  
  # Cancer-specific outcome variables, raw values.
  cancer_outcome_cols  = c("stage", "yrs_since_dx", "yrs_since_chemo"),
  
  fdr_method = "BH"
)
dir.create(config$out_dir, recursive = TRUE, showWarnings = FALSE)

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
correlate_me_traits <- function(me_mat, trait_mat, category_label) {
  cor_mat <- WGCNA::cor(me_mat, trait_mat, use = "pairwise.complete.obs")
  n_per_trait <- sapply(trait_mat, function(x) sum(!is.na(x)))
  pval_mat <- cor_mat
  for (tr in colnames(cor_mat)) {
    pval_mat[, tr] <- WGCNA::corPvalueStudent(cor_mat[, tr], nSamples = n_per_trait[tr])
  }
  fdr_mat <- matrix(
    p.adjust(pval_mat, method = config$fdr_method),
    nrow = nrow(pval_mat), dimnames = dimnames(pval_mat)
  )
  
  tibble::as_tibble(cor_mat, rownames = "module") %>%
    tidyr::pivot_longer(-module, names_to = "trait", values_to = "r") %>%
    dplyr::left_join(
      tibble::as_tibble(pval_mat, rownames = "module") %>%
        tidyr::pivot_longer(-module, names_to = "trait", values_to = "p"),
      by = c("module", "trait")
    ) %>%
    dplyr::left_join(
      tibble::as_tibble(fdr_mat, rownames = "module") %>%
        tidyr::pivot_longer(-module, names_to = "trait", values_to = "fdr"),
      by = c("module", "trait")
    ) %>%
    dplyr::mutate(sig = fdr < 0.05, n = n_per_trait[trait], category = category_label)
}

# ============================================================
# 1. Delta ME vs delta traits (reuse from 05_wgcna_muscle.R section 7b)
# ============================================================
if (is.null(rd$delta_trait_long)) {
  warning("rd$delta_trait_long is NULL -- delta block will be empty. ",
          "Make sure 05_wgcna_muscle.R section 7b ran with valid delta_pheno_cols.")
  delta_block <- tibble::tibble(module = character(0), trait = character(0), r = double(0),
                                p = double(0), fdr = double(0), sig = logical(0),
                                n = integer(0), category = character(0))
} else {
  delta_block <- rd$delta_trait_long %>%
    dplyr::rename(trait = variable) %>%
    dplyr::mutate(category = "\u0394 (post - pre)")
}

# ============================================================
# 2. Baseline (PRE) ME vs baseline (pre_) traits
# ============================================================
me_baseline_wide <- rd$me_long %>%
  dplyr::filter(.data[[group_col]] == "SURV_T1", module %in% module_order) %>%
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
  warning("No valid baseline_pheno_cols found -- baseline block will be empty.")
  baseline_block <- tibble::tibble(module = character(0), trait = character(0), r = double(0),
                                   p = double(0), fdr = double(0), sig = logical(0),
                                   n = integer(0), category = character(0))
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
  
  baseline_block <- correlate_me_traits(baseline_me_mat, baseline_trait_mat, "Baseline (PRE)")
}

# ============================================================
# 3. Cancer-specific outcomes vs baseline (SURV_T1) ME
# ============================================================
# stage/yrs_since_dx/yrs_since_chemo are subject-level constants,
# duplicated across a subject's SURV_T1 and SURV_T2 rows -- same pattern
# already confirmed for the pre_/post_ columns. Restricted to SURV_T1
# only, same as the Baseline block above: including SURV_T2 would
# pseudoreplicate each subject (same trait value paired against two
# different ME values) and confound a fixed baseline characteristic with
# a POST-training eigengene that's already been shaped by the
# intervention. CTL samples drop out structurally either way (these
# variables are NA for them), so restricting to SURV_T1 only changes
# which SURV rows are used, not whether CTL contributes.
present_cancer_cols <- intersect(config$cancer_outcome_cols, colnames(meta_by_sample))
missing_cancer_cols <- setdiff(config$cancer_outcome_cols, colnames(meta_by_sample))
if (length(missing_cancer_cols) > 0) {
  warning("cancer_outcome_cols: couldn't find these columns in metadata, skipping: ",
          paste(missing_cancer_cols, collapse = ", "))
}

if (length(present_cancer_cols) == 0) {
  warning("No valid cancer_outcome_cols found -- cancer outcomes block will be empty.")
  cancer_block <- tibble::tibble(module = character(0), trait = character(0), r = double(0),
                                 p = double(0), fdr = double(0), sig = logical(0),
                                 n = integer(0), category = character(0))
} else {
  baseline_sample_ids <- meta_by_sample %>%
    dplyr::filter(.data[[group_col]] == "SURV_T1") %>%
    dplyr::pull(sample_id)
  
  me_wide_baseline <- rd$MEs[rownames(rd$MEs) %in% baseline_sample_ids, module_order, drop = FALSE]
  
  cancer_trait_mat <- meta_by_sample %>%
    dplyr::filter(sample_id %in% rownames(me_wide_baseline)) %>%
    dplyr::arrange(match(sample_id, rownames(me_wide_baseline))) %>%
    dplyr::select(dplyr::all_of(present_cancer_cols)) %>%
    as.data.frame()
  rownames(cancer_trait_mat) <- rownames(me_wide_baseline)
  
  cancer_block <- correlate_me_traits(as.data.frame(me_wide_baseline), cancer_trait_mat, "Cancer outcomes")
}

# ============================================================
# 4. Combine + plot
# ============================================================
CATEGORY_LEVELS <- c("\u0394 (post - pre)", "Baseline (PRE)", "Cancer outcomes")

# trait_key (category::trait), not raw trait name, is the actual x
# factor -- the same variable name can legitimately appear in more than
# one category (baseline_pheno_cols defaults to the same names as
# delta_pheno_cols), and facet_grid showing "fat_kg" as a valid tick in
# two different facet panels is fine, but a single global factor() on
# raw trait name isn't, since factor() requires unique level strings
# across the whole column even though only one category's facet ever
# displays a given tick at once.
trait_key_levels <- c(
  paste("\u0394 (post - pre)", unique(delta_block$trait), sep = "::"),
  paste("Baseline (PRE)", unique(baseline_block$trait), sep = "::"),
  paste("Cancer outcomes", unique(cancer_block$trait), sep = "::")
)

all_blocks <- dplyr::bind_rows(
  delta_block    %>% dplyr::select(module, trait, r, p, fdr, sig, n, category),
  baseline_block %>% dplyr::select(module, trait, r, p, fdr, sig, n, category),
  cancer_block   %>% dplyr::select(module, trait, r, p, fdr, sig, n, category)
) %>%
  dplyr::mutate(
    module    = factor(module, levels = rev(module_order)),
    category  = factor(category, levels = CATEGORY_LEVELS),
    trait_key = paste(category, trait, sep = "::"),
    trait_key = factor(trait_key, levels = trait_key_levels),
    label     = sprintf("%.2f\n(%.3f)", r, p)
  )

trait_key_labels <- all_blocks %>%
  dplyr::distinct(trait_key, trait) %>%
  tibble::deframe()

message(sprintf(
  "Trait heatmap: %d delta rows, %d baseline rows, %d cancer-outcome rows (%d modules x %d trait slots total)",
  nrow(delta_block), nrow(baseline_block), nrow(cancer_block),
  length(module_order), length(trait_key_levels)
))

p_trait_heat <- ggplot(all_blocks, aes(x = trait_key, y = module, fill = r)) +
  geom_tile(color = "white") +
  geom_tile(data = dplyr::filter(all_blocks, sig), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red,
                       midpoint = 0, limits = c(-1, 1), name = "r") +
  scale_x_discrete(labels = trait_key_labels) +
  facet_grid(. ~ category, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = NULL,
       title = paste("Module eigengene \u2013 trait correlations \u2014", config$tissue),
       subtitle = paste(
         "\u0394 = post-pre change score (CRC only, subject-level);",
         "Baseline = PRE-timepoint ME (CRC only, subject-level);",
         "Cancer outcomes = baseline (SURV_T1) ME only.",
         "Black outline = FDR q < 0.05"
       )) +
  theme_cr() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold")
  )

pdf_path <- file.path(config$out_dir, paste0("wgcna_trait_heatmap_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = max(10, length(trait_key_levels) * 0.8), height = max(6, length(module_order) * 0.5))
print(p_trait_heat)
dev.off()

ggsave(file.path(config$out_dir, paste0("wgcna_trait_heatmap_", config$tissue, ".png")),
       p_trait_heat, width = max(10, length(trait_key_levels) * 0.8), height = max(6, length(module_order) * 0.5),
       dpi = 300, bg = "white")

write_xlsx(
  list(all_correlations = all_blocks %>% dplyr::select(-label, -trait_key)),
  file.path(config$out_dir, paste0("wgcna_trait_heatmap_data_", config$tissue, ".xlsx"))
)

message("Done. Saved to: ", config$out_dir)