# Effect Size Distribution - SURV Contrasts
# Maps logFC distributions across the three SURV contrasts

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, ggplot2, tidyr, stringr, purrr)

# --- 1: Load Results -------------------
surv_rds <- readRDS("03_analysis/b_reports/muscle/SURVvCTL/limma/limma_out_SURV.RDS")
results  <- surv_rds$results_list

# --- 2: Combine into long format -------
df_long <- bind_rows(
  lapply(names(results), function(cname) {
    results[[cname]] |>
      select(uniprot, gene, logFC) |>
      mutate(contrast = cname)
  })
) |>
  mutate(contrast = factor(contrast, levels = c("Baseline_SURVvCTL",
                                                "Training_SURVvCTL",
                                                "Training_SURV")))

# --- 3: Compute summary stats ----------
summary_stats <- df_long |>
  group_by(contrast) |>
  summarise(
    n       = n(),
    mean    = mean(logFC, na.rm = TRUE),
    median  = median(logFC, na.rm = TRUE),
    sd      = sd(logFC, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)

# --- 4: Plot ---------------------------

contrast_labels <- c(
  "Baseline_SURVvCTL" = "CRC_PRE v CTL",
  "Training_SURVvCTL" = "CRC_POST vs CTL",
  "Training_SURV"     = "CRC_POST vs CRC_PRE"
)

contrast_colors <- c(
  "Baseline_SURVvCTL" = "#4C72B0",
  "Training_SURVvCTL" = "#DD8452",
  "Training_SURV"     = "#55A868"
)

# Density plot
p_density <- ggplot(df_long, aes(x = logFC, fill = contrast, color = contrast)) +
  geom_density(alpha = 0.3, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_vline(data = summary_stats,
             aes(xintercept = median, color = contrast),
             linetype = "dotted", linewidth = 0.8) +
  scale_fill_manual(values  = contrast_colors,
                    labels  = contrast_labels,
                    name    = NULL) +
  scale_color_manual(values = contrast_colors,
                     labels = contrast_labels,
                     name   = NULL) +
  labs(
    title    = "logFC Distribution Across SURV Contrasts",
    subtitle = "Dotted lines indicate median logFC per contrast",
    x        = expression(log[2]~"Fold Change"),
    y        = "Density"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "top",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold")
  )

# Histogram faceted
p_hist <- ggplot(df_long, aes(x = logFC, fill = contrast)) +
  geom_histogram(bins = 60, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
  geom_vline(data = summary_stats,
             aes(xintercept = median),
             linetype = "dotted", color = "black", linewidth = 0.8) +
  facet_wrap(~ contrast, ncol = 1,
             labeller = labeller(contrast = contrast_labels)) +
  scale_fill_manual(values = contrast_colors, guide = "none") +
  labs(
    title    = "logFC Distribution Across SURV Contrasts",
    subtitle = "Dotted lines indicate median logFC per contrast",
    x        = expression(log[2]~"Fold Change"),
    y        = "Count"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )

# --- 5: Save ---------------------------
output_dir <- "03_analysis/b_reports/muscle/SURVvCTL/limma"

ggsave(file.path(output_dir, "logFC_distribution_density.pdf"),
       p_density, width = 7, height = 5)

ggsave(file.path(output_dir, "logFC_distribution_histogram.pdf"),
       p_hist, width = 7, height = 8)

cat("Plots saved to:", output_dir, "\n")
