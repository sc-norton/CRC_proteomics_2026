# =============================================================================
# fig_heatmap_prepost_zscore_muscle_plasma.R
# Heatmap: proteins significant in Training_SURV (POST vs PRE), CRC survivors
# only (sample_id starts with "CR"), z-scored per protein across samples.
# MUSCLE + PLASMA combined, faceted by tissue. Uses IMPUTED data for
# display/clustering; cells that were NA in the normalized (pre-imputation)
# data are bordered to flag them as imputed.
#
# Rows: significant proteins per tissue (gene symbol labels), clustered
#       WITHIN tissue via standard Euclidean distance -- muscle and plasma
#       do not share a protein list, so they cannot share a row order.
# Columns: CR samples for that tissue, ordered by timepoint then subject
#          (not clustered). Muscle and plasma samples are different
#          columns entirely, so each tissue panel gets its own x-axis via
#          free facet scales.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, stringr, tibble, purrr, ggplot2, ggnewscale, writexl)

source("04_figures/a_scripts/00_theme.R")

# --- 0: Config -------------------------------------------------------------
# One block per tissue -- everything tissue-specific lives here so sections
# 1-2 can be a single function (load_tissue) called once per tissue via map().
tissue_configs <- list(
  muscle = list(
    tissue_label    = "Muscle",
    dalist_norm_rds = "01_normalization/c_data/muscle/FINAL/SURV/02_DAList_normalized_SURV_muscle.RDS",  
    dalist_imp_rds  = "02_imputation/c_data/muscle_imputed_dalist.RDS",                                   
    limma_xlsx      = "03_analysis/c_data/muscle/reformat/limma_results_all_muscle.xlsx",                 
    limma_sheet     = "Training_SURV"                                                                     
  ),
  plasma = list(
    tissue_label    = "Plasma",
    dalist_norm_rds = "01_normalization/c_data/plasma/01_normalized_DAList_SURV.RDS", 
    dalist_imp_rds  = "02_imputation/c_data/plasma_imputed_dalist.RDS",                                   
    limma_xlsx      = "03_analysis/c_data/plasma/reformat/limma_results_all_plasma.xlsx",                 
    limma_sheet     = "Training_SURV"                                                                     
  )
)

config <- list(
  fig_dir    = "04_figures/fig_heatmap",
  data_dir   = "04_figures/data/fig_heatmap",
  sig_col    = "Pi.Val",
  sig_thresh = 0.05
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 10) {
  ggsave(file.path(config$fig_dir, paste0(name, ".pdf")), p, device = cairo_pdf, width = width, height = height)
  ggsave(file.path(config$fig_dir, paste0(name, ".png")), p, device = "png", width = width, height = height, dpi = 300)
}

# =============================================================================
# 1: Per-tissue loader -- identical logic to the original single-tissue
# script (load, subset, z-score, cluster, order), parameterized by tc.
# Returns the pieces section 2/3 need: long-format heat_df/ann_df, z_mat +
# na_mask (for export), sample_meta, and the clustered row order.
# =============================================================================

load_tissue <- function(tc) {
  
  dal_norm <- readRDS(tc$dalist_norm_rds)
  dal_imp  <- readRDS(tc$dalist_imp_rds)
  
  mat_norm <- dal_norm$data   # >>> VERIFY <<< proteins x samples, contains NA
  mat_imp  <- dal_imp$data    # >>> VERIFY <<< same shape, no NA -- confirm dimnames match mat_norm's
  
  limma_res <- read_excel(tc$limma_xlsx, sheet = tc$limma_sheet)
  
  sig_proteins <- limma_res %>%
    filter(.data[[config$sig_col]] < config$sig_thresh) %>%
    pull(uniprot_id)   # >>> VERIFY <<< column name
  
  cat(sprintf("[%s] Significant proteins (%s < %.2f): %d\n",
              tc$tissue_label, config$sig_col, config$sig_thresh, length(sig_proteins)))
  
  cr_samples <- colnames(mat_imp)[str_starts(colnames(mat_imp), "CR")]
  
  proteins_available <- intersect(sig_proteins, intersect(rownames(mat_norm), rownames(mat_imp)))
  missing_from_imp <- setdiff(sig_proteins, rownames(mat_imp))
  if (length(missing_from_imp) > 0) {
    warning(sprintf("[%s] %d significant protein(s) not found in the imputed matrix (dropped during imputation filtering?): %s",
                    tc$tissue_label, length(missing_from_imp), paste(head(missing_from_imp, 10), collapse = ", ")))
  }
  
  mat_sig_imp  <- mat_imp[proteins_available, cr_samples, drop = FALSE]
  mat_sig_norm <- mat_norm[proteins_available, cr_samples, drop = FALSE]
  
  na_mask <- is.na(mat_sig_norm)   # TRUE = was NA pre-imputation
  
  cat(sprintf("[%s] Matrix subset: %d proteins x %d CR samples | %d cells imputed (%.1f%%)\n",
              tc$tissue_label, nrow(mat_sig_imp), ncol(mat_sig_imp), sum(na_mask), 100 * mean(na_mask)))
  if (nrow(mat_sig_imp) == 0) stop(sprintf("[%s] No significant proteins matched -- check the uniprot_id/rowname join in Section 1.", tc$tissue_label))
  
  extract_timepoint <- function(ids) {
    case_when(
      str_detect(ids, regex("_?T1$|_?PRE$",  ignore_case = TRUE)) ~ "PRE",
      str_detect(ids, regex("_?T2$|_?POST$", ignore_case = TRUE)) ~ "POST",
      TRUE ~ NA_character_
    )
  }
  
  sample_meta <- tibble(sample_id = cr_samples) %>%
    mutate(
      timepoint  = extract_timepoint(sample_id),
      subject_id = str_remove(sample_id, "_?(T1|T2|PRE|POST)$")   # >>> VERIFY <<<
    )
  
  if (any(is.na(sample_meta$timepoint))) {
    warning(sprintf("[%s] %d sample(s) had no PRE/POST match:\n%s",
                    tc$tissue_label, sum(is.na(sample_meta$timepoint)),
                    paste(sample_meta$sample_id[is.na(sample_meta$timepoint)], collapse = ", ")))
  }
  
  # --- Row-wise z-score, on imputed values (no NA remaining) -----------------
  z_mat <- t(scale(t(mat_sig_imp)))
  
  # --- Row clustering, WITHIN this tissue -------------------------------------
  dist_rows <- dist(z_mat, method = "euclidean")
  hc_rows <- hclust(dist_rows, method = "ward.D2")
  row_order <- rownames(z_mat)[hc_rows$order]
  
  # --- Column ordering: timepoint, then subject, WITHIN this tissue ----------
  sample_meta <- sample_meta %>%
    mutate(timepoint = factor(timepoint, levels = c("PRE", "POST"))) %>%
    arrange(timepoint, subject_id)
  col_order <- sample_meta$sample_id
  
  # --- Gene-symbol labels ------------------------------------------------------
  gene_lookup <- limma_res %>%
    distinct(uniprot_id, gene) %>%   # >>> VERIFY <<<
    filter(uniprot_id %in% rownames(z_mat))
  
  dup_genes <- gene_lookup$gene[duplicated(gene_lookup$gene) | duplicated(gene_lookup$gene, fromLast = TRUE)]
  
  gene_lookup <- gene_lookup %>%
    mutate(gene_label = if_else(gene %in% dup_genes, paste0(gene, " (", uniprot_id, ")"), gene))
  
  row_order_genes <- gene_lookup$gene_label[match(row_order, gene_lookup$uniprot_id)]
  
  # --- Long format --------------------------------------------------------------
  heat_df <- as.data.frame(z_mat) %>%
    rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "sample_id", values_to = "zscore")
  
  mask_df <- as.data.frame(na_mask) %>%
    rownames_to_column("protein") %>%
    pivot_longer(-protein, names_to = "sample_id", values_to = "is_imputed")
  
  heat_df <- heat_df %>%
    left_join(mask_df, by = c("protein", "sample_id")) %>%
    left_join(gene_lookup %>% dplyr::select(uniprot_id, gene_label), by = c("protein" = "uniprot_id")) %>%
    left_join(sample_meta, by = "sample_id") %>%
    mutate(
      sample_id     = factor(sample_id, levels = col_order),
      imputed_label = factor(if_else(is_imputed, "Imputed", "Observed"), levels = c("Observed", "Imputed")),
      tissue        = tc$tissue_label
    )
  
  ann_df <- sample_meta %>%
    mutate(sample_id = factor(sample_id, levels = col_order), tissue = tc$tissue_label)
  
  list(
    heat_df         = heat_df,
    ann_df          = ann_df,
    z_mat           = z_mat,
    na_mask         = na_mask,
    sample_meta     = sample_meta,
    row_order_genes = row_order_genes
  )
}

tissue_data <- map(tissue_configs, load_tissue)

# =============================================================================
# 2: Combine across tissues.
# y-axis factor levels = muscle's clustered order followed by plasma's
# clustered order. With facet_grid(scales = "free"), each panel only draws
# the levels present in its own data and drops the rest, so this single
# combined factor preserves each tissue's within-cluster order independently
# -- it does NOT force muscle and plasma proteins onto a shared ordering.
# =============================================================================

combined_row_levels <- unlist(map(tissue_data, "row_order_genes"), use.names = FALSE)

heat_df_all <- map_dfr(tissue_data, "heat_df") %>%
  mutate(
    gene_label = factor(gene_label, levels = combined_row_levels),
    tissue     = factor(tissue, levels = c("Muscle", "Plasma"))
  )

ann_df_all <- map_dfr(tissue_data, "ann_df") %>%
  mutate(tissue = factor(tissue, levels = c("Muscle", "Plasma")))

# =============================================================================
# 3: Plot -- facet_grid(tissue ~ ., scales = "free", space = "free") gives
# each tissue its own protein rows AND its own sample columns (muscle/plasma
# are different subjects x different sample IDs, so they cannot share an
# x-axis). space = "free" also makes panel height proportional to protein
# count, so plasma won't get squashed if muscle has far more sig. proteins.
#
# NOTE ON THE "Timepoint" ANNOTATION STRIP: it's drawn as an extra row per
# tissue (same trick as the original script), placed inside each facet panel
# via ann_df_all's own `tissue` column. This is the part most worth checking
# first when you run it -- mixing a literal "Timepoint" y-value with the
# gene_label factor across two facets is a little fragile. If it doesn't sit
# at the top of each panel correctly, the more robust fix is to stop cramming
# it into the shared y-scale and instead build it as a separate strip via
# patchwork (annotation row + heatmap, per tissue, stacked with
# `patchwork::wrap_plots(ncol = 1)`), which sidesteps the shared-scale issue
# entirely and mirrors how your Figure 5 composite is already assembled.
# =============================================================================

fig_heatmap <- ggplot() +
  geom_tile(data = ann_df_all, aes(x = sample_id, y = "Timepoint", fill = timepoint)) +
  scale_fill_manual(values = c(PRE = colors$blue_dark, POST = colors$red_dark), name = "Timepoint") +
  ggnewscale::new_scale_fill() +
  geom_tile(data = heat_df_all,
            aes(x = sample_id, y = gene_label, fill = zscore, color = imputed_label),
            width = 0.92, height = 0.92, linewidth = 0.4) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red,
                       midpoint = 0, na.value = "grey85", name = "Z-score") +
  scale_color_manual(values = c(Observed = NA, Imputed = "black"), name = "Value",
                     guide = guide_legend(override.aes = list(fill = "grey90"))) +
  facet_grid(tissue ~ ., scales = "free", space = "free") +
  labs(x = NULL, y = NULL,
       title = "Significant Proteins (POST vs PRE), CRC Survivors",
       subtitle = sprintf("Muscle: %d proteins x %d samples | Plasma: %d proteins x %d samples | \u03A0 < %.2f | imputed cells bordered",
                          nrow(tissue_data$muscle$z_mat), ncol(tissue_data$muscle$z_mat),
                          nrow(tissue_data$plasma$z_mat), ncol(tissue_data$plasma$z_mat),
                          config$sig_thresh)) +
  theme_cr(base_size = 9) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    axis.text.y  = element_text(size = 5),
    axis.ticks.y = element_blank(),
    panel.grid   = element_blank(),
    plot.title    = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0, size = 8, color = "grey40"),
    strip.background = element_rect(fill = "grey90", color = NA),  # >>> VERIFY <<< swap fill per-tissue (e.g. colors$muscle/colors$plasma) if your other faceted figures color the strip itself rather than just label it
    strip.text        = element_text(face = "bold")
  )

fig_heatmap

total_rows <- sum(map_int(tissue_data, ~ nrow(.x$z_mat)))
save_plot(fig_heatmap, "fig_heatmap_prepost_zscore_muscle_plasma",
          width = 10, height = max(8, total_rows * 0.08))

# =============================================================================
# 4: Companion export -- z-scores/imputed-mask/metadata per tissue, kept as
# separate sheets (row sets and sample sets differ between tissues, so a
# single merged sheet would just be two block-diagonal chunks of NA).
# =============================================================================

write_xlsx(
  list(
    "Muscle Z-scores"  = as.data.frame(tissue_data$muscle$z_mat) %>% rownames_to_column("protein"),
    "Muscle Imputed"   = as.data.frame(tissue_data$muscle$na_mask) %>% rownames_to_column("protein"),
    "Muscle Metadata"  = tissue_data$muscle$sample_meta,
    "Plasma Z-scores"  = as.data.frame(tissue_data$plasma$z_mat) %>% rownames_to_column("protein"),
    "Plasma Imputed"   = as.data.frame(tissue_data$plasma$na_mask) %>% rownames_to_column("protein"),
    "Plasma Metadata"  = tissue_data$plasma$sample_meta
  ),
  path = file.path(config$data_dir, "fig_heatmap_prepost_zscore_muscle_plasma_data.xlsx")
)