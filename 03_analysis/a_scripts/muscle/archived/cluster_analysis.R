# =============================================================================
# Mfuzz Soft Clustering of Proteins Across 3 Contrasts
# Contrasts: PREvCTL, POSTvPRE, POSTvCTL
# Input: limma results RDS (named list of data frames with logFC column)
# =============================================================================

# --- 0: Setup --------------------------
# Set working directory
setwd(rprojroot::find_rstudio_root_file())
getwd()
# Load packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("Mfuzz")
pacman::p_load(Biobase, dplyr, ggplot2, tibble, patchwork, RColorBrewer)
library(Mfuzz)
library(readxl)

# Configuration
config <- list(
  # Directories
  report_dir = "03_analysis/b_reports/muscle",
  data_dir = "03_analysis/c_data/muscle",
  # Data files
  SURV_RDS = "01_normalization/c_data/muscle/SURV/stringent/BM15/02_DAList_normalized_SURV_stringent_BM15_muscle.RDS",
  CRE_RDS = "01_normalization/c_data/muscle/CRE/stringent/BM15/02_DAList_normalized_CRE_stringent_BM15_muscle.RDS"
)

dir.create(config$report_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

# =============================================================================
# 1. USER SETTINGS — adjust these to match your environment
# =============================================================================

rds_path      <- "03_analysis/c_data/muscle/limma_results_all_muscle.xlsx"  # RDS with named list of results
output_dir    <- "03_analysis/b_reports/muscle/SURVvCTL/cluster"                    # where to save plots/tables
contrast_names <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")  # must match names in RDS
logfc_col     <- "logFC"        # column name for log2 fold change in each result df
protein_col   <- "gene"           # set to column name if protein IDs aren't rownames; NULL if rownames
min_membership <- 0.5           # membership threshold for high-confidence assignments
k_range       <- seq(2, 12, 1)  # cluster numbers to test in Dmin
k_final       <- 2           # set integer after inspecting Dmin plot; NULL = auto (not recommended)
repeats_dmin  <- 3              # repeats for Dmin stability
standardize   <- TRUE           # TRUE = cluster by shape; FALSE = cluster by shape + magnitude

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 2. LOAD RESULTS AND BUILD INPUT MATRIX
# =============================================================================

results <- lapply(contrast_names, function(ct) {
  read_xlsx(rds_path, sheet = ct)  # assumes sheet names match contrast_names
})
names(results) <- contrast_names

# Extract logFC per contrast and join into wide matrix
fc_list <- lapply(contrast_names, function(ct) {
  df <- results[[ct]]
  if (!is.null(protein_col)) {
    ids <- df[[protein_col]]
  } else {
    ids <- rownames(df)
  }
  data.frame(protein = ids, logFC = df[[logfc_col]], stringsAsFactors = FALSE)
})
names(fc_list) <- contrast_names

fc_wide <- Reduce(function(a, b) full_join(a, b, by = "protein", suffix = c("", "_dup")), fc_list)
colnames(fc_wide) <- c("protein", contrast_names)

# Remove proteins with any NA across contrasts
n_before <- nrow(fc_wide)
fc_wide  <- fc_wide[complete.cases(fc_wide), ]
n_after  <- nrow(fc_wide)
cat(sprintf("Proteins retained after removing NAs: %d / %d\n", n_after, n_before))

rownames(fc_wide) <- fc_wide$protein
fc_mat <- as.matrix(fc_wide[, contrast_names])

cat("Input matrix dimensions:", nrow(fc_mat), "proteins x", ncol(fc_mat), "contrasts\n")

# =============================================================================
# 3. BUILD ExpressionSet AND STANDARDIZE
# =============================================================================

eset <- new("ExpressionSet", exprs = fc_mat)

if (standardize) {
  eset_use <- standardise(eset)
  cat("Standardization applied (clustering by trajectory shape).\n")
} else {
  eset_use <- eset
  cat("No standardization (clustering by shape + magnitude).\n")
}

# =============================================================================
# 4. ESTIMATE FUZZIFIER m
# =============================================================================

m_est <- mestimate(eset_use)
cat(sprintf("Estimated fuzzifier m: %.3f\n", m_est))

# =============================================================================
# 5. DETERMINE OPTIMAL CLUSTER NUMBER WITH Dmin
# =============================================================================

cat("Running Dmin to estimate optimal cluster number...\n")

pdf(file.path(output_dir, "dmin_plot.pdf"), width = 7, height = 5)
dmin_result <- Dmin(eset_use, m = m_est, crange = k_range, repeats = repeats_dmin, visu = TRUE)
dev.off()

cat("Dmin values:\n")
print(dmin_result)
cat("\nInspect dmin_plot.pdf and set k_final above, then re-run from section 6.\n")

# =============================================================================
# 6. RUN FUZZY CLUSTERING
# (Re-run from here after setting k_final)
# =============================================================================

if (is.null(k_final)) {
  stop("Set k_final based on the Dmin plot, then re-run from section 6.")
}

set.seed(26)
clust_result <- mfuzz(eset_use, c = k_final, m = m_est)

cat(sprintf("Clustering complete: %d clusters, m = %.3f\n", k_final, m_est))
cat("Cluster sizes (all proteins, hard-assigned to max membership):\n")
print(table(apply(clust_result$membership, 1, which.max)))

# =============================================================================
# 7. BUILT-IN Mfuzz TRAJECTORY PLOTS
# =============================================================================

n_cols <- min(3, k_final)
n_rows <- ceiling(k_final / n_cols)

pdf(file.path(output_dir, "mfuzz_trajectory_plots.pdf"),
    width = n_cols * 4, height = n_rows * 3.5)
mfuzz.plot2(
  eset_use,
  cl         = clust_result,
  mfrow      = c(n_rows, n_cols),
  time.labels = contrast_names,
  min.mem    = min_membership,
  x11        = FALSE
)
dev.off()

cat("Trajectory plots saved to mfuzz_trajectory_plots.pdf\n")

# Algebra check:
fc_wide %>%
  mutate(LFC_dif = round((Baseline_SURVvCTL + Training_SURV) - Training_SURVvCTL, 4))

# =============================================================================
# 8. EXTRACT MEMBERSHIP AND ASSIGNMENTS
# =============================================================================

memb <- as.data.frame(clust_result$membership)
colnames(memb) <- paste0("cluster_", seq_len(k_final), "_memb")
memb$protein   <- rownames(memb)

assignments <- memb %>%
  mutate(
    cluster  = apply(clust_result$membership, 1, which.max),
    max_memb = apply(clust_result$membership, 1, max)
  ) %>%
  select(protein, cluster, max_memb, everything()) %>%
  arrange(cluster, desc(max_memb))

high_conf <- assignments %>% filter(max_memb >= min_membership)

cat(sprintf("\nHigh-confidence assignments (membership >= %.2f): %d / %d proteins\n",
            min_membership, nrow(high_conf), nrow(assignments)))
cat("High-confidence cluster sizes:\n")
print(table(high_conf$cluster))

write.csv(assignments, file.path(output_dir, "mfuzz_all_assignments.csv"), row.names = FALSE)
write.csv(high_conf,   file.path(output_dir, "mfuzz_highconf_assignments.csv"), row.names = FALSE)

# =============================================================================
# 9. CLUSTER CENTROIDS — interpret trajectory patterns
# =============================================================================

centroids <- as.data.frame(clust_result$centers)
colnames(centroids) <- contrast_names
centroids$cluster   <- paste0("Cluster_", seq_len(k_final))
centroids <- centroids %>% select(cluster, everything())

cat("\nCluster centroids (standardised scale):\n")
print(round(centroids, 3))

# Annotate trajectory pattern from centroid signs
centroids$pattern <- apply(centroids[, contrast_names], 1, function(x) {
  signs <- ifelse(x > 0.3, "^", ifelse(x < -0.3, "v", "-"))
  paste(signs, collapse = " | ")
})

cat("\nCluster trajectory patterns:\n")
print(centroids[, c("cluster", "pattern", contrast_names)])

write.csv(centroids, file.path(output_dir, "mfuzz_centroids.csv"), row.names = FALSE)

# =============================================================================
# 10. GGPLOT2 TRAJECTORY PLOT (publication-quality)
# =============================================================================

# Prepare long-format data for high-confidence proteins
plot_data <- high_conf %>%
  select(protein, cluster, max_memb) %>%
  left_join(
    fc_wide,  # already has protein column, no rownames_to_column needed
    by = "protein"
  ) %>%
  tidyr::pivot_longer(
    cols = all_of(contrast_names),
    names_to  = "contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    contrast = factor(contrast, levels = contrast_names),
    cluster  = factor(paste0("Cluster ", cluster))
  )

# Centroid lines in long format
centroid_long <- centroids %>%
  tidyr::pivot_longer(
    cols = all_of(contrast_names),
    names_to  = "contrast",
    values_to = "logFC"
  ) %>%
  mutate(
    contrast = factor(contrast, levels = contrast_names),
    cluster  = factor(cluster)
  )

# Add cluster sizes to facet labels
cluster_sizes <- high_conf %>%
  count(cluster) %>%
  mutate(label = paste0("Cluster ", cluster, "\n(n = ", n, ")"))
size_lookup <- setNames(cluster_sizes$label, paste0("Cluster ", cluster_sizes$cluster))

p_traj <- ggplot(plot_data, aes(x = contrast, y = logFC, group = protein)) +
  geom_line(aes(alpha = max_memb), color = "steelblue3", linewidth = 0.3) +
  geom_line(
    data = centroid_long,
    aes(x = contrast, y = logFC, group = cluster),
    color = "firebrick", linewidth = 1.2, inherit.aes = FALSE
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.4) +
  scale_alpha_continuous(range = c(0.05, 0.5), guide = "none") +
  facet_wrap(~ cluster, labeller = labeller(cluster = size_lookup)) +
  labs(
    title    = "Mfuzz Soft Clustering — Protein Trajectories Across Contrasts",
    subtitle = paste0("k = ", k_final, ", m = ", round(m_est, 2),
                      ", membership threshold = ", min_membership),
    x        = "Contrast",
    y        = "log2 Fold Change",
    caption  = "Red line = cluster centroid; Blue lines = individual proteins (alpha = membership)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    axis.text.x      = element_text(angle = 30, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(output_dir, "mfuzz_trajectories_ggplot.pdf"),
  plot     = p_traj,
  width    = n_cols * 4,
  height   = n_rows * 3.5
)

cat("ggplot trajectory figure saved.\n")

# =============================================================================
# 11. FLAG REVERSAL CLUSTERS
# Flag clusters where PREvCTL and POSTvPRE have opposite signs
# (i.e., cancer dysregulation partially or fully reversed by exercise)
# =============================================================================

centroids <- centroids %>%
  mutate(
    reversal_flag = case_when(
      PREvCTL > 0.3  & POSTvPRE < -0.3 ~ "Reversal (UP in cancer → DN after exercise)",
      PREvCTL < -0.3 & POSTvPRE >  0.3 ~ "Reversal (DN in cancer → UP after exercise)",
      TRUE ~ "Not reversal"
    )
  )

cat("\nReversal cluster flags:\n")
print(centroids[, c("cluster", "pattern", "reversal_flag")])

# Save final annotated centroids
write.csv(centroids, file.path(output_dir, "mfuzz_centroids_annotated.csv"), row.names = FALSE)

# =============================================================================
# 12. SAVE CLUSTERING OBJECT
# =============================================================================

saveRDS(
  list(
    clust_result = clust_result,
    eset_use     = eset_use,
    assignments  = assignments,
    high_conf    = high_conf,
    centroids    = centroids,
    m_est        = m_est,
    k_final      = k_final
  ),
  file = file.path(output_dir, "mfuzz_results.rds")
)

cat("\nAll outputs saved to:", output_dir, "\n")
cat("Done.\n")