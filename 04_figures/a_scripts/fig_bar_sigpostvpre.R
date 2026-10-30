# =============================================================================
# fig_bar_prepost_intensity_muscle.R
# Bar plot: normalized intensity (mean +/- SD) at PRE vs POST, per protein,
# for proteins significant in Training_SURV (POST vs PRE), CRC survivors
# only (sample_id starts with "CR"). One bar-pair per protein, faceted,
# paginated across multiple PDF pages since the significant-protein list is
# likely too large for one page of panels.
#
# Shares data prep with fig_heatmap_prepost_zscore_muscle.R (same DAList,
# same significance filter, same CR-sample subset, same gene-label join) --
# kept as a separate standalone script per project convention rather than
# sourcing the heatmap script, since the two produce genuinely different
# artifacts (raw mean+/-SD vs. row-wise z-score) and shouldn't share
# execution state.
# =============================================================================

setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, stringr, tibble, ggplot2, ggforce, writexl)

source("04_figures/a_scripts/00_theme.R")

# --- 0: Config -----------------------------------------------------------
config <- list(
  dalist_rds  = "01_normalization/c_data/muscle/FINAL/SURV/02_DAList_normalized_SURV_muscle.RDS",  
  limma_xlsx  = "03_analysis/c_data/muscle/reformat/limma_results_all_muscle.xlsx",
  limma_sheet = "Training_SURV",
  fig_dir     = "04_figures/fig_bar_intensity",
  data_dir    = "04_figures/data/fig_bar_intensity",
  sig_col     = "Pi.Val",
  sig_thresh  = 0.05,
  panels_per_page = 20,   # 4 rows x 5 cols -- adjust to taste
  ncol_page   = 5
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1: Load normalized DAList + limma results (same as heatmap script)
# =============================================================================

dal <- readRDS(config$dalist_rds)
mat_raw <- dal$data   # >>> VERIFY <<< proteins x samples

limma_res <- read_excel(config$limma_xlsx, sheet = config$limma_sheet)

sig_proteins <- limma_res %>%
  filter(.data[[config$sig_col]] < config$sig_thresh) %>%
  pull(uniprot_id)   # >>> VERIFY <<< column name

cat(sprintf("Significant proteins (%s < %.2f): %d\n", config$sig_col, config$sig_thresh, length(sig_proteins)))

# =============================================================================
# 2: Subset to significant proteins x CR samples
# =============================================================================

cr_samples <- colnames(mat_raw)[str_starts(colnames(mat_raw), "CR")]
mat_sig <- mat_raw[rownames(mat_raw) %in% sig_proteins, cr_samples, drop = FALSE]

cat(sprintf("Matrix subset: %d proteins x %d CR samples\n", nrow(mat_sig), ncol(mat_sig)))
if (nrow(mat_sig) == 0) stop("No significant proteins matched rownames(mat_raw) -- check the uniprot_id/rowname join in Section 1.")

# --- Timepoint extraction from sample_id -----------------------------------
# >>> KNOWN ISSUE, NOT YET FIXED <<< this is the same placeholder regex
# flagged as broken in fig_heatmap_prepost_zscore_muscle.R -- it was
# mis-splitting T1/T2 there. Apply the same fix here once real sample_id
# examples confirm the correct pattern; both scripts need the identical
# correction since they share this exact logic.
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
    subject_id = str_remove(sample_id, "_?(T1|T2|PRE|POST)$")   # >>> VERIFY <<< same known-fragile extraction
  )

if (any(is.na(sample_meta$timepoint))) {
  warning(sprintf(
    "%d sample(s) had no PRE/POST match:\n%s",
    sum(is.na(sample_meta$timepoint)),
    paste(sample_meta$sample_id[is.na(sample_meta$timepoint)], collapse = ", ")
  ))
}

# =============================================================================
# 3: Gene-symbol labels (same disambiguation logic as the heatmap script)
# =============================================================================

gene_lookup <- limma_res %>%
  distinct(uniprot_id, gene) %>%   # >>> VERIFY <<< column names
  filter(uniprot_id %in% rownames(mat_sig))

dup_genes <- gene_lookup$gene[duplicated(gene_lookup$gene) | duplicated(gene_lookup$gene, fromLast = TRUE)]

gene_lookup <- gene_lookup %>%
  mutate(gene_label = if_else(gene %in% dup_genes,
                              paste0(gene, " (", uniprot_id, ")"),
                              gene))

# =============================================================================
# 4: Long format -- raw normalized intensity per protein x sample, joined
# to timepoint + gene label. NOT z-scored -- this figure shows actual
# intensity scale, unlike the heatmap.
# =============================================================================

long_df <- as.data.frame(mat_sig) %>%
  rownames_to_column("uniprot_id") %>%
  pivot_longer(-uniprot_id, names_to = "sample_id", values_to = "intensity") %>%
  left_join(sample_meta, by = "sample_id") %>%
  left_join(gene_lookup %>% dplyr::select(uniprot_id, gene_label), by = "uniprot_id") %>%
  filter(!is.na(timepoint), !is.na(intensity)) %>%   # drop samples with no PRE/POST match + missing values
  mutate(timepoint = factor(timepoint, levels = c("PRE", "POST")))

# =============================================================================
# 5: Summarize -- mean +/- SD per gene x timepoint
# =============================================================================

summary_df <- long_df %>%
  group_by(gene_label, uniprot_id, timepoint) %>%
  summarise(mean_intensity = mean(intensity), sd_intensity = sd(intensity), n = n(), .groups = "drop")

# Order genes alphabetically for stable pagination -- swap to order by
# p-value/logFC if you'd rather the most significant proteins appear first.
gene_order <- sort(unique(summary_df$gene_label))
summary_df <- summary_df %>% mutate(gene_label = factor(gene_label, levels = gene_order))
long_df    <- long_df    %>% mutate(gene_label = factor(gene_label, levels = gene_order))

# =============================================================================
# 6: Single wide bar chart -- all significant proteins on x-axis, dodged by
# timepoint, mean +/- SD error bars.
# =============================================================================

n_proteins <- length(gene_order)

fig_bar <- ggplot(summary_df, aes(x = gene_label, y = mean_intensity, fill = timepoint)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, color = NA) +
  geom_errorbar(
    aes(ymin = mean_intensity - sd_intensity, ymax = mean_intensity + sd_intensity),
    position = position_dodge(width = 0.75), width = 0.25, linewidth = 0.4, color = "grey20"
  ) +
  scale_fill_manual(values = c(PRE = colors$blue, POST = colors$red), name = "Timepoint") +
  labs(x = NULL, y = "Normalized Intensity",
       title = "Significant Proteins (POST vs PRE), CRC Survivors",
       subtitle = sprintf("%d proteins (Π < %.2f)", n_proteins, config$sig_thresh)) +
  theme_cr(base_size = 10) +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
    legend.position = "top",
    plot.title    = element_text(hjust = 0),
    plot.subtitle = element_text(hjust = 0, size = 9, color = "grey40"),
    panel.grid.major.x = element_blank()
  )

fig_bar

# Width scales with protein count so labels don't collide -- ~0.18in per
# protein, floor of 12in for small significant-protein counts.
save_width <- max(12, n_proteins * 0.18)
ggsave(file.path(config$fig_dir, "fig_bar_prepost_intensity_muscle.pdf"), fig_bar,
       device = cairo_pdf, width = save_width, height = 7, limitsize = FALSE)
ggsave(file.path(config$fig_dir, "fig_bar_prepost_intensity_muscle.png"), fig_bar,
       device = "png", width = save_width, height = 7, dpi = 300, limitsize = FALSE)

# =============================================================================
# 7: Companion export
# =============================================================================

write_xlsx(
  list("Summary" = summary_df, "Per-sample" = long_df),
  path = file.path(config$data_dir, "fig_bar_prepost_intensity_muscle_data.xlsx")
)