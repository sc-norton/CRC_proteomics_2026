# Contamination Assessment
setwd(rprojroot::find_rstudio_root_file())

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, org.Hs.eg.db, AnnotationDbi, tibble, ggpubr,
               Hmisc, ggnewscale)

# Configuration
config <- list(
  # Data file
  report_data = "01_normalization/c_data/muscle/00_report_data.rds",
  # Directories
  report_dir = "01_normalization/b_reports/muscle",
  data_dir = "01_normalization/c_data/muscle"
)

# Color pallet
colors <- list(
  # Standard color list
  red = "#C0392B", red_light  = "#E6A19A", red_dark   = "#7B241C",
  blue = "#2471A3", blue_light = "#A9CCE3", blue_dark  = "#154360",
  green = "#4E8B6F", green_light= "#A9C9B8", green_dark = "#2E5E4E",
  pink = "#D98C9C", pink_light = "#F2C6CF", pink_dark  = "#A65C6B",
  purple = "#7B5EA7", purple_light = "#C5B4DD", purple_dark  = "#4A3A73",
  yellow = "#D4A017", yellow_light = "#F1D77A", yellow_dark  = "#9A6F0A",
  orange = "#C4693A", orange_light = "#E5A07B", orange_dark  = "#8A4321",
  aqua = "#5DAFA3", aqua_light = "#A9D6CF", aqua_dark  = "#3A7F76",
  brown = "#8C5A3C", brown_light = "#C7A189", brown_dark  = "#5A3A28",
  red_pale = "#FADBD8", blue_pale  = "#D6EAF8",
  gray_light = "#F5F5F5", gray = "#9E9E9E", gray_dark = "#4A4A4A",
  black = "#000000", white = "#FFFFFF",
  auburn     = c(navy   = "#0D234A", auburn = "#E96A2C")
)

# Load data
report_data <- readRDS(file.path(config$report_data))
marker_pct_df <- report_data$marker_pct_df
dals <- report_data$dals

make_long_from_dal <- function(d, nm) {
  
  df <- bind_cols(
    d$annotation %>% dplyr::select(uniprot_id, gene, comp1_grouped),
    d$data
  )
  
  long <- df %>%
    pivot_longer(
      cols = -c(uniprot_id, gene, comp1_grouped),
      names_to = "sample_id",
      values_to = "intensity"
    ) %>%
    left_join(d$metadata, by = "sample_id")
  
  # Add parsed metadata from name
  long <- long %>%
    mutate(
      dataset = ifelse(grepl("^SURV", nm), "SURV", "CRE"),
      filter  = stringr::str_extract(nm, "none|moderate|stringent"),
      ratio   = stringr::str_extract(nm, "BM[^_]+") %>% str_remove("BM")
    )
  
  return(long)
}

all_long <- purrr::imap_dfr(dals, make_long_from_dal)

# Compartment summaries
comp_summary <- all_long %>%
  filter(!is.na(comp1_grouped),
         !comp1_grouped %in% c("Unknown", "Other")) %>%
  
  group_by(dataset, filter, ratio, sample_id, comp1_grouped) %>%
  summarise(
    comp_intensity = sum(intensity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  group_by(dataset, filter, ratio, sample_id) %>%
  mutate(
    total_intensity = sum(comp_intensity),
    comp_pct = comp_intensity / total_intensity
  ) %>%
  ungroup() %>%
  
  # bring in contamination metric
  left_join(
    all_long %>%
      distinct(sample_id, blood_pct),
    by = "sample_id"
  )

# Compute correlations
cor_results <- comp_summary %>%
  group_by(dataset, filter, ratio, comp1_grouped) %>%
  summarise(
    cor = cor(comp_pct, blood_pct,
              method = "spearman",
              use = "pairwise.complete.obs"),
    p = suppressWarnings(
      cor.test(comp_pct, blood_pct,
               method = "spearman")$p.value
    ),
    .groups = "drop"
  )

# Pipeline label for plotting
cor_results <- cor_results %>%
  mutate(pipeline = paste(dataset, filter, ratio, sep = "_"))

# Add p-value (BH_adj)
cor_results <- cor_results %>%
  group_by(dataset, filter, ratio) %>%  # adjust within each pipeline
  mutate(
    p_adj = p.adjust(p, method = "BH")
  ) %>%
  ungroup()

# P-value labels for heatmap
cor_results <- cor_results %>%
  mutate(
    sig_label = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Heatmap
bc_cor <- ggplot(cor_results,
       aes(x = comp1_grouped, y = pipeline, fill = cor)) +
  geom_tile() +
  geom_text(aes(label = sig_label), size = 3) +
  scale_fill_gradient2(
    low = colors$blue, mid = colors$white, high = colors$red, midpoint = 0
  ) +
  theme_minimal() +
  labs(
    x = "Cellular compartment",
    y = "Pipeline",
    fill = "Spearman r",
    title = "Blood contam. corr. with GO:CC"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(config$report_dir, "11_blood_contam_correlation.pdf"), bc_cor, device = "pdf")

comp_summary %>%
  filter(comp1_grouped == "Blood Microparticle",
        ratio == 15) %>%
  ggplot(aes(x = blood_pct, y = comp_pct, color = dataset)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_grid(filter ~ ratio) +
  theme_minimal() +
  labs(
    title = "Sarco/Myo vs contamination",
    x = "Contamination (%)",
    y = "Proportion"
  )

comp_summary %>%
  filter(dataset == "SURV", filter == "stringent") %>%
  ggplot(aes(x = blood_pct, y = comp_intensity)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    ggpubr::stat_cor(
      method = "spearman",
      aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),
      label.x.npc = "left",
      label.y.npc = "top"
      ) +
    facet_grid(comp1_grouped ~ dataset + filter + ratio) +
    theme_minimal()

# Rank Compartments
comp_rank <- cor_results %>%
  group_by(comp1_grouped) %>%
  summarise(
    mean_abs_cor = mean(abs(cor), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_abs_cor))

# Rank pipelines
pipeline_score <- cor_results %>%
  group_by(dataset, filter, ratio) %>%
  summarise(
    mean_abs_cor = mean(abs(cor), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_abs_cor)

pipeline_score <- cor_results %>%
  group_by(dataset, filter, ratio) %>%
  summarise(
    mean_abs_cor = mean(abs(cor), na.rm = TRUE),
    n_sig        = sum(p_adj < 0.05, na.rm = TRUE),
    frac_sig     = mean(p_adj < 0.05, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_abs_cor)

comp_sum_5 <- comp_summary %>%
  filter(ratio == 5, dataset =="SURV", comp1_grouped == "Mitochondria", filter == "stringent")
comp_sum_10 <- comp_summary %>%
  filter(ratio == 10, dataset =="SURV", comp1_grouped == "Mitochondria", filter == "stringent")
comp_sum_15 <- comp_summary %>%
  filter(ratio == 15, dataset =="SURV", comp1_grouped == "Mitochondria", filter == "stringent")
comp_sum_20 <- comp_summary %>%
  filter(ratio == 20, dataset =="SURV", comp1_grouped == "Mitochondria", filter == "stringent")

# --- PCA Correlations --------------

