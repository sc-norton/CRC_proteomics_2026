setwd(rprojroot::find_rstudio_root_file())
getwd()

# Load Packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(proteoDA, readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, org.Hs.eg.db, AnnotationDbi, tibble, ggpubr,
               Hmisc, ggnewscale)
# Configuration
config <- list(
  # Data
  report_data = "01_normalization/c_data/muscle/00_report_data.rds",
  # Directories
  report_dir = "01_normalization/b_reports/muscle",
  data_dir = "01_normalization/c_data/muscle"
)

# Color palette
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

# Load Data
report_data <- readRDS(file.path(config$report_data))
flagged_hpa <- report_data$flagged_hpa
flagged_blood_moderate <- report_data$flagged_blood_moderate
flagged_blood_stringent <- report_data$flagged_blood_stringent
marker_pct_df <- report_data$marker_pct_df
sample_summary <- report_data$sample_summary
outlier_diag <- report_data$outlier_diag
filter_logs <- report_data$filter_logs

build_output_path <- function(nm, dir) {
  
  dataset <- ifelse(grepl("^SURV", nm), "SURV", "CRE")
  fl      <- stringr::str_extract(nm, "none|moderate|stringent")
  bm      <- stringr::str_extract(nm, "BM[^_]+")
  
  out_dir <- file.path(dir, dataset, fl, bm)
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  return(out_dir)
}
# Flagged Contaminants
flagged_contaminants_all <- bind_rows(flagged_hpa, 
                                      flagged_blood_moderate,
                                      flagged_blood_stringent) %>%
  group_by(gene) %>%
  summarise(
    protein           = dplyr::first(protein),
    gene              = dplyr::first(gene),
    description       = dplyr::first(description),
    in_hpa_muscle     = dplyr::first(in_hpa_muscle),
    in_moderate = dplyr::first(in_moderate),
    in_stringent = dplyr::first(in_stringent),
    .groups = "drop"
  )

write.csv(flagged_hpa, 
          file.path(config$report_dir, "01_flagged_contaminants_hpa_muscle.csv"))
write.csv(flagged_blood_moderate,
          file.path(config$report_dir, "01_flagged_contaminants_blood_moderate.csv"))
write.csv(flagged_blood_stringent,
          file.path(config$report_dir, "01_flagged_contaminants_blood_stringent.csv"))
write.csv(flagged_contaminants_all, 
          file.path(config$report_dir, "01_flagged_contaminants_all_muscle.csv"))
cat("\n Saved flagged_contaminants.csv : \n", config$report_dir, "\n")

# Pre-compute label positions (top of each Blood stack)
blood_label_df <- marker_pct_df %>%
  dplyr::select(sample_id, ALB_pct, HBB_pct, HBA1_pct, B_M_ratio) %>%
  mutate(blood_total = ALB_pct + HBB_pct + HBA1_pct) %>%
  arrange(desc(blood_total)) %>%
  mutate(
    marker_type = "Blood",
    sample_id   = factor(sample_id, levels = unique(sample_id))
  )

blood_order <- levels(blood_label_df$sample_id)

p_bm <- marker_pct_df %>%
  dplyr::select(sample_id, ALB_pct, HBB_pct, HBA1_pct, MB_pct, CKM_pct) %>%
  pivot_longer(-sample_id, names_to = "marker", values_to = "pct") %>%
  mutate(
    marker      = gsub("_pct", "", marker),
    marker_type = case_when(
      marker %in% c("ALB", "HBB", "HBA1") ~ "Blood",
      marker %in% c("MB", "CKM")          ~ "Muscle"
    ),
    sample_id = factor(sample_id, levels = blood_order)
  ) %>%
  ggplot(aes(x = sample_id, y = pct)) +
  # --- Blood layer ---
  geom_col(
    data = ~ dplyr::filter(.x, marker_type == "Blood"),
    aes(fill = marker),
    width = 0.8
  ) +
  scale_fill_manual(
    name   = "Blood markers",
    values = c("ALB"  = colors$red_dark,
               "HBB"  = colors$red_light,
               "HBA1" = colors$red)
  ) +
  new_scale_fill() +
  # --- Muscle layer ---
  geom_col(
    data = ~ dplyr::filter(.x, marker_type == "Muscle"),
    aes(fill = marker),
    width = 0.8
  ) +
  scale_fill_manual(
    name   = "Muscle markers",
    values = c("MB"  = colors$aqua,
               "CKM" = colors$aqua_dark)
  ) +
  facet_wrap(~ marker_type, ncol = 1, scales = "free_y") +
  labs(
    x = "Sample",
    y = "% of total intensity",
    title = "Blood vs muscle marker % of total intensity per sample"
  ) +
  geom_text(
    data    = blood_label_df,
    aes(x = sample_id, y = blood_total, label = round(B_M_ratio, 2)),
    vjust   = -0.4,
    size    = 2.5,
    inherit.aes = FALSE
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    strip.text  = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(config$report_dir, "01_blood_v_muscle_intensity_plot_raw_muscle.pdf"),
       p_bm, width = 14, height = 9)

write.csv(marker_pct_df, file.path(config$data_dir, "01_blood_contamination.csv"))

# Save filter logs and filtered protein lists
walk(names(filter_logs), function(nm) {
  log <- filter_logs[[nm]]
  d   <- dals[[nm]]
  
  dir <- build_output_path(nm, config$data_dir)
  
  write_csv(log, file.path(report_dir, sprintf("02_filtering_effects_%s_muscle.csv", nm)))
  
  kept_uniprots    <- rownames(d$data)
  removed_uniprots <- setdiff(annotation$uniprot_id, kept_uniprots)
  
  filtered_proteins <- annotation %>%
    filter(uniprot_id %in% removed_uniprots) %>%
    dplyr::select(uniprot_id, gene, description, in_hpa_muscle, in_blood_none, in_blood_moderate, in_blood_stringent) %>%
    mutate(reason = if_else(!in_hpa_muscle | in_blood_moderate | in_blood_stringent, "Contaminant", "Missingness"))
  
  write_csv(filtered_proteins, file.path(report_dir, sprintf("03_filtered_proteins_%s_muscle.csv", nm)))
  cat(sprintf(" Saved: 03_filtering_effects_%s_muscle.csv\n", nm))
  cat(sprintf(" Saved: 03_filtered_proteins_%s_muscle.csv\n", nm))
})

# Build filter plot function
make_filter_plot <- function(log, label) {
  ggplot(mutate(log, step = factor(step, levels = step)), aes(x = step, y = n_after)) +
    geom_col(fill = colors$blue, width = 0.5) +
    geom_text(aes(label = n_after), vjust = -0.3, size = 3) +
    labs(x = NULL, y = "Number of Proteins",
         title = paste("Proteins Retained Through Filtering -", label)) +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

# Save filter plots
walk(names(filter_logs), function(nm) {
  p <- make_filter_plot(filter_logs[[nm]], nm)
  dir <- build_output_path(nm, config$report_dir)
  ggsave(file.path(dir, sprintf("04_filtering_plot_%s_muscle.pdf", nm)),
         p, width = 8, height = 5)
  cat(sprintf(" Saved: 04_filtering_plot_%s_muscle.pdf\n", nm))
})

# Contaminant plots
flagged_contaminants <- flagged_contaminants_all %>%
  mutate(
    contaminant_category = case_when(
      grepl("KERATIN|KRT",           description, ignore.case = TRUE) ~ "Keratins",
      grepl("SERUM|ALBUMIN",         description, ignore.case = TRUE) ~ "Serum Proteins",
      grepl("IMMUNOGLOBULIN|IGG",    description, ignore.case = TRUE) ~ "Immunoglobulins",
      grepl("TRYP",                  description, ignore.case = TRUE) ~ "Trypsins",
      grepl("HEMOGLOBIN",            description, ignore.case = TRUE) ~ "Hemoglobins",
      grepl("LIPOPROTEIN",           description, ignore.case = TRUE) ~ "Lipoproteins",
      TRUE                                                            ~ "Other"
    )
  )

contaminant_summary <- flagged_contaminants %>%
  count(contaminant_category) %>%
  mutate(percentage = round(n / sum(n) * 100, 1))

# Arrange to match plot
contaminant_summary <- contaminant_summary %>%
  arrange(contaminant_category)

p_contaminants <- ggplot(contaminant_summary,
                         aes(x = 2, y = n, fill = contaminant_category)) +
  geom_col(color = "white") +
  coord_polar(theta = "y") +
  xlim(c(0.5, 2.5)) +
  scale_fill_brewer(
    palette = "Set2",
    labels  = paste0(contaminant_summary$contaminant_category, " (", 
                     contaminant_summary$percentage, "%)")
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.text     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(title = "Contaminant Classification by Type",
       fill  = NULL)

ggsave(file.path(config$report_dir, "05_contaminant_classification_muscle.pdf"),
       p_contaminants, width = 8, height = 12)
ggsave(file.path(config$report_dir, "05_contaminant_classification_muscle.png"),
       p_contaminants, width = 8, height = 12)
cat(" Saved: 05_contaminant_classification_muscle.pdf\n")

# Save outlier data
walk(names(outlier_diag), function(nm) {
  dir <- build_output_path(nm, config$data_dir)
  write_csv(
    outlier_diag[[nm]]$outlier_diag,
    file.path(dir, sprintf("06_outlier_diagnostics_%s_muscle.csv", nm))
  )
})

# Make outlier plots
res <- outlier_diag[[nm]]
df  <- res$outlier_diag
thr <- res$thresholds
pca <- res$pca_res

# Diagnostic plots
p1 <- ggplot(outlier_diag %>% filter(!is.na(delta_missing)),
             aes(x = pct_missing, y = delta_missing)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  geom_hline(yintercept = c(delta_lower, delta_threshold),
             linetype = "dashed", color = "red", alpha = 0.5) +
  geom_vline(xintercept = miss_threshold,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "% Missing", y = "Delta Missing (T2 - T1)", title = "A: Paired Missingness") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")

var_explained <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)
pc_df <- as.data.frame(pca_res$x[, 1:2])
pc_df$sample_id <- rownames(pc_df)
pc_df <- left_join(pc_df, outlier_diag %>% dplyr::select(sample_id, consensus_outlier), by = "sample_id")

p2 <- ggplot(pc_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = consensus_outlier), size = 3) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = sprintf("PC1 (%.1f%%)", var_explained[1]),
       y = sprintf("PC2 (%.1f%%)", var_explained[2]),
       title = "B: PCA Outliers (Mahalanobis)") +
  theme_minimal(base_size = 11) + theme(legend.position = "none")

p3 <- ggplot(outlier_diag, aes(x = reorder(sample_id, sample_median), y = sample_median)) +
  geom_point(aes(color = consensus_outlier), size = 2.5) +
  geom_hline(yintercept = global_median, color = "black") +
  geom_hline(yintercept = global_median + c(-3, 3) * mad_val,
             linetype = "dashed", color = "red", alpha = 0.5) +
  scale_color_manual(values = c("FALSE" = "gray40", "TRUE" = "red")) +
  labs(x = "Sample", y = "Median log2 intensity", title = "C: MAD Intensity Outliers") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
        legend.position = "none")

ggsave(file.path(report_dir, sprintf("06_outlier_diagnostic_plots_%s_muscle.pdf", nm)),
       p1 + p2 + p3 + plot_layout(ncol = 1), width = 10, height = 12)
cat(sprintf("   Saved: 06_outlier_diagnostic_plots_%s.pdf\n", nm))

# --- 7: Normalized Data Assessment -------------------------------

# --- CV plots ---
# CV dataframe function
build_cv_data <- function(data, annot, meta) {
  
  group_map <- bind_rows(
    # CTL group
    meta %>%
      filter(cancer == "CTL") %>%
      mutate(group = "CTL"),
    # SURV by timepoint
    meta %>%
      filter(cancer_time == "SURV_T1") %>%
      mutate(group = "CRC_T1"),
    meta %>%
      filter(cancer_time == "SURV_T2") %>%
      mutate(group = "CRC_T2"),
    meta %>%
      filter(supp_time == "CRE_T1") %>%
      mutate(group = "CRE_T1"),
    meta %>%
      filter(supp_time == "CRE_T2") %>%
      mutate(group = "CRE_T2"),
    meta %>%
      filter(supp_time == "PLA_T1") %>%
      mutate(group = "PLA_T1"),
    meta %>%
      filter(supp_time == "PLA_T2") %>%
      mutate(group = "PLA_T2"),
  ) %>%
    dplyr::select(sample_id, group) %>%
    filter(!is.na(group))
  
  rn   <- rownames(data)
  data <- as.data.frame(base::apply(data, 2, as.numeric))
  rownames(data) <- rn
  
  data %>%
    rownames_to_column("uniprot_id") %>%
    left_join(annot %>% dplyr::select(uniprot_id, gene), by = "uniprot_id") %>%
    filter(!is.na(gene)) %>%
    dplyr::select(gene, any_of(group_map$sample_id)) %>%
    pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
    left_join(group_map, by = "sample_id", relationship = "many-to-many") %>%
    filter(!is.na(group), !is.na(intensity)) %>%
    group_by(gene, group) %>%
    summarise(
      mean_int = mean(intensity, na.rm = TRUE),
      sd_int   = sd(intensity,   na.rm = TRUE),
      n        = sum(!is.na(intensity)),
      cv       = (sd_int / mean_int) * 100,
      .groups  = "drop"
    ) %>%
    filter(mean_int > 0, n >= 3)
}

# Make CV Wide data function
make_cv_wide <- function(cv_df, groups) {
  cv_df %>%
    filter(group %in% groups) %>%
    dplyr::select(gene, group, cv) %>%
    pivot_wider(names_from = group, values_from = cv)
}

# Scatter plot function
make_cv_scatter <- function(cv_wide, x_group, y_group,
                            n_label   = 20,
                            delta_col = NULL,
                            label = NULL) {
  
  # Check groups exist
  missing <- setdiff(c(x_group, y_group), colnames(cv_wide))
  if (length(missing) > 0) stop("Groups not found in cv_wide: ", paste(missing, collapse = ", "))
  
  # Delta column — auto-derive if not supplied
  if (is.null(delta_col)) {
    delta_col <- paste0("cv_delta_", x_group, "v", y_group)
    cv_wide <- cv_wide %>%
      mutate(!!delta_col := .data[[y_group]] - .data[[x_group]])
  }
  
  df <- cv_wide %>%
    filter(!is.na(.data[[x_group]]), !is.na(.data[[y_group]]))
  
  r_val     <- cor(df[[x_group]], df[[y_group]], use = "pairwise.complete.obs")
  delta_lim <- max(abs(df[[delta_col]]), na.rm = TRUE)
  
  df_label <- df %>% slice_max(abs(.data[[delta_col]]), n = n_label)
  
  ggplot(df, aes(x = .data[[x_group]], y = .data[[y_group]],
                 colour = .data[[delta_col]])) +
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", colour = "grey40", linewidth = 0.6) +
    geom_point(size = 1.8, alpha = 0.8) +
    geom_smooth(method = "lm", colour = colors$black, fill = "grey80",
                alpha = 0.15, linewidth = 0.8) +
    geom_text_repel(
      data          = df_label,
      aes(label     = gene),
      colour        = "black",
      fontface      = "bold",
      size          = 2.8,
      max.overlaps  = Inf,
      box.padding   = 0.4,
      point.padding = 0.2,
      segment.color = "grey60",
      segment.size  = 0.3,
      show.legend   = FALSE,
      seed          = 42
    ) +
    scale_colour_gradientn(
      colours = c(colors$purple_dark, colors$purple, colors$orange,
                  colors$yellow,
                  colors$orange,    colors$purple, colors$purple_dark),
      values  = scales::rescale(c(-delta_lim, -delta_lim * 0.5, -delta_lim * 0.15,
                                  0,
                                  delta_lim * 0.15,  delta_lim * 0.5,  delta_lim)),
      limits  = c(-delta_lim, delta_lim),
      name    = sprintf("\u0394CV\n(%s \u2212 %s)", y_group, x_group)
    ) +
    annotate("text",
             x = Inf, y = -Inf, hjust = 1.1, vjust = -0.8,
             label    = sprintf("r = %.2f\n", r_val),
             size     = 4, fontface = "bold", colour = colors$gray_dark) +
    labs(
      x     = sprintf("CV \u2014 %s (%%)", x_group),
      y     = sprintf("CV \u2014 %s (%%)", y_group),
      title = sprintf("Per-protein CV: %s vs %s", x_group, y_group),
      subtitle = sprintf("Dataset: %s", label)
    ) +
    theme_minimal()
}

# ── 1. Build CV data ──────────────────────────────────────────────────────────

cv_df_raw <- build_cv_data(data  = intensity,
                           annot = annotation,
                           meta  = meta)

cv_df     <- build_cv_data(data  = dals$SURV_moderate$data,
                           annot = dals$SURV_moderate$annotation,
                           meta  = dals$SURV_moderate$metadata)

cv_df_str <- build_cv_data(data  = dals$SURV_stringent$data,
                           annot = dals$SURV_stringent$annotation,
                           meta  = dals$SURV_stringent$metadata)


# ── 2. Wide format + delta columns ───────────────────────────────────────────

ALL_GROUPS <- c("CTL", "CRC_T1", "CRC_T2",
                "CRE_T1", "CRE_T2", "PLA_T1", "PLA_T2")

build_cv_surv <- function(cv_df) {
  cv_df %>%
    filter(group %in% ALL_GROUPS) %>%
    dplyr::select(gene, group, cv) %>%
    pivot_wider(names_from = group, values_from = cv) %>%
    filter(!is.na(SURV_T1), !is.na(SURV_T2)) %>%
    mutate(
      cv_delta_t1vt2   = CRC_T2 - CRC_T1,
      cv_delta_t1vctl  = CRC_T1 - CTL,
      cv_delta_t2vctl  = CRC_T2 - CTL,
      cv_delta_ct1vct2 = CRE_T2  - CRE_T1,
      cv_delta_pt1vpt2 = PLA_T2  - PLA_T1,
      cv_delta_ct1vpt1 = CRE_T1  - PLA_T1,
      cv_delta_ct2vpt2 = CRE_T2  - PLA_T2
    )
}

cv_surv     <- build_cv_surv(cv_df)
cv_surv_str <- build_cv_surv(cv_df_str)
cv_surv_raw <- build_cv_surv(cv_df_raw)

# Wide format without delta columns (used by make_cv_scatter with auto-delta)
cv_wide_mod <- make_cv_wide(cv_df,     ALL_GROUPS)
cv_wide_str <- make_cv_wide(cv_df_str, ALL_GROUPS)
cv_wide_raw <- make_cv_wide(cv_df_raw, ALL_GROUPS)


# ── 3. CV violin plot ─────────────────────────────────────────────────────────

GROUP_COLOURS <- c(
  CTL     = colors$green,
  SURV_T1 = colors$purple,
  SURV_T2 = colors$orange,
  CRE_T1  = colors$purple_light,
  CRE_T2  = colors$orange_light,
  PLA_T1  = colors$purple_dark,
  PLA_T2  = colors$orange_dark
)

p_cv_violin <- ggplot(cv_df_raw, aes(x = group, y = cv, fill = group, color = group)) +
  geom_violin(alpha = 0.5, linewidth = 0.5, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.8, linewidth = 0.5,
               outlier.shape = NA, colour = colors$gray_dark) +
  scale_fill_manual(values   = GROUP_COLOURS, guide = "none") +
  scale_colour_manual(values = GROUP_COLOURS, guide = "none") +
  scale_y_continuous(limits = c(0, NA)) +
  labs(x = NULL, y = "CV (%)", title = "Per-protein CV by group") +
  theme_minimal()


# ── 4. CV scatter plots ───────────────────────────────────────────────────────

# Named list of comparisons: x_group, y_group, delta_col, cv_surv object, label
cv_scatter_specs <- list(
  list(data = cv_surv_str, x = "CRC_T1", y = "CTL",     delta = "cv_delta_t1vctl",  label = "str", n = 30),
  list(data = cv_surv_str, x = "CRC_T1", y = "CRC_T2", delta = "cv_delta_t1vt2",   label = "str", n = 30),
  list(data = cv_surv_raw, x = "CRC_T1", y = "CTL",     delta = "cv_delta_t1vctl",  label = "raw", n = 50),
  list(data = cv_surv_raw, x = "CRC_T1", y = "CRC_T2", delta = "cv_delta_t1vt2",   label = "raw", n = 50)
)

cv_scatter_plots <- lapply(cv_scatter_specs, function(s) {
  make_cv_scatter(s$data, s$x, s$y,
                  delta_col = s$delta,
                  label     = s$label,
                  n_label   = s$n)
})

cv_scatter_plots[[4]]

# Save all
lapply(cv_scatter_specs, function(s) {
  plot_name <- sprintf("cv_scatter_%s_v_%s_%s.pdf", s$x, s$y, s$label)
  ggsave(
    file.path(report_dir, plot_name),
    make_cv_scatter(s$data, s$x, s$y,
                    delta_col = s$delta,
                    label     = s$label,
                    n_label   = s$n),
    width = 7, height = 6, device = "pdf"
  )
})

# --- 8: Protein Assessments ------------

invisible(lapply(names(dals), function(nm) {
  d         <- dals[[nm]]
  miss_pct  <- rowMeans(is.na(d$data))
  
  report <- tibble(
    Total      = nrow(d$data),
    Complete   = sum(miss_pct == 0),
    Incomplete = sum(miss_pct > 0)
  )
  dir <- build_output_path(nm, config$data_dir)
  write.csv(report, file.path(dir, sprintf("11_protein_report_%s_muscle.csv", nm)))
  cat(sprintf("  [%s] Total: %d | Complete: %d | Incomplete: %d\n", 
              nm, report$Total, report$Complete, report$Incomplete))
}))

protein_summary <- map_dfr(names(dals), function(nm) {
  
  d <- dals[[nm]]
  miss_pct <- rowMeans(is.na(d$data))
  
  tibble(
    dataset = ifelse(str_detect(nm, "^SURV"), "SURV", "CRE"),
    filter  = str_extract(nm, "none|moderate|stringent"),
    ratio   = str_extract(nm, "BM[^_]+") |> str_remove("BM"),
    
    Total      = nrow(d$data),
    Complete   = sum(miss_pct == 0),
    Incomplete = sum(miss_pct > 0)
  )
})

# Optional: clean up ratio column
protein_summary <- protein_summary %>%
  mutate(ratio = ifelse(ratio == "none", "none", as.character(ratio)))

# Save single combined CSV
write.csv(protein_summary,
          file.path(config$data_dir, "11_protein_report_all_muscle.csv"),
          row.names = FALSE)

# Print summary
print(protein_summary)
# --- 9: Sample Assessments -----------------

invisible(lapply(names(dals), function(nm) {
  meta <- dals[[nm]]$metadata
  
  # Timepoint counts
  t1_count <- sum(meta$timepoint == "T1", na.rm = T)
  t2_count <- sum(meta$timepoint == "T2", na.rm = T)
  
  # Supplement counts
  CRE_count <- sum(meta$supp == "CRE", na.rm = T)
  PLA_count <- sum(meta$supp == "PLA", na.rm = T)
  
  # Sex counts
  m_count <- sum(meta$sex == "Male", na.rm = T)
  f_count <- sum(meta$sex == "Female", na.rm = T)
  
  # Supp_time counts
  CRE_T1_count <- sum(meta$supp_time == "CRE_T1", na.rm = T)
  CRE_T2_count <- sum(meta$supp_time == "CRE_T2", na.rm = T)
  PLA_T1_count <- sum(meta$supp_time == "PLA_T1", na.rm = T)
  PLA_T2_count <- sum(meta$supp_time == "PLA_T2", na.rm = T)
  
  # Sex_time counts
  m_T1_count <- sum(meta$sex_time == "Male_T1", na.rm = T)
  m_T2_count <- sum(meta$sex_time == "Male_T2", na.rm = T)
  f_T1_count <- sum(meta$sex_time == "Female_T1", na.rm = T)
  f_T2_count <- sum(meta$sex_time == "Female_T2", na.rm = T)
  
  # Paired vs unpaired subjects
  paired_subjects <- meta %>% 
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_cre <- meta %>%
    filter(supp == "CRE") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_pla <- meta %>%
    filter(supp == "PLA") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  unpaired_subjects <- meta %>%
    count(pid) %>%
    filter(n == 1) %>%
    nrow()
  
  # Controls (only relevant for SURV v CTL)
  n_ctl <- sum(meta$cancer == "CTL", na.rm = T)
  n_surv <- sum(meta$cancer == "SURV", na.rm = T)
  
  # T1 and T2 only
  t1_only <- t1_count - paired_subjects
  t2_only <- t2_count - paired_subjects
  
  report <- tibble(
    total_samples       = nrow(meta),
    n_CRE               = CRE_count,
    n_PLA               = PLA_count,
    n_male              = m_count,
    n_female            = f_count,
    n_male_T1           = m_T1_count,
    n_male_T2           = m_T2_count,
    n_female_T1         = f_T1_count,
    n_female_T2         = f_T2_count,
    n_CRE_T1            = CRE_T1_count,
    n_CRE_T2            = CRE_T2_count,
    n_PLA_T1            = PLA_T1_count,
    n_PLA_T2            = PLA_T2_count,
    n_T1                = t1_count,
    n_T2                = t2_count,
    n_paired            = paired_subjects,
    n_paired_CRE        = paired_cre,
    n_paired_PLA        = paired_pla,
    n_unpaired          = unpaired_subjects,
    n_T1_only           = t1_only,
    n_T2_only           = t2_only,
    n_SURV              = n_surv,
    n_CTL               = n_ctl
  )
  
  write.csv(report, file.path(report_dir, sprintf("12_sample_report_%s_muscle.csv", nm)))
}))

sample_summary <- map_dfr(names(dals), function(nm) {
  
  meta <- dals[[nm]]$metadata
  
  # --- Counts ---
  t1_count <- sum(meta$timepoint == "T1", na.rm = TRUE)
  t2_count <- sum(meta$timepoint == "T2", na.rm = TRUE)
  
  CRE_count <- sum(meta$supp == "CRE", na.rm = TRUE)
  PLA_count <- sum(meta$supp == "PLA", na.rm = TRUE)
  
  m_count <- sum(meta$sex == "Male", na.rm = TRUE)
  f_count <- sum(meta$sex == "Female", na.rm = TRUE)
  
  CRE_T1_count <- sum(meta$supp_time == "CRE_T1", na.rm = TRUE)
  CRE_T2_count <- sum(meta$supp_time == "CRE_T2", na.rm = TRUE)
  PLA_T1_count <- sum(meta$supp_time == "PLA_T1", na.rm = TRUE)
  PLA_T2_count <- sum(meta$supp_time == "PLA_T2", na.rm = TRUE)
  
  m_T1_count <- sum(meta$sex_time == "Male_T1", na.rm = TRUE)
  m_T2_count <- sum(meta$sex_time == "Male_T2", na.rm = TRUE)
  f_T1_count <- sum(meta$sex_time == "Female_T1", na.rm = TRUE)
  f_T2_count <- sum(meta$sex_time == "Female_T2", na.rm = TRUE)
  
  # --- Pairing ---
  paired_subjects <- meta %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_cre <- meta %>%
    filter(supp == "CRE") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  paired_pla <- meta %>%
    filter(supp == "PLA") %>%
    count(pid) %>%
    filter(n == 2) %>%
    nrow()
  
  unpaired_subjects <- meta %>%
    count(pid) %>%
    filter(n == 1) %>%
    nrow()
  
  # --- Groups ---
  n_ctl  <- sum(meta$cancer == "CTL",  na.rm = TRUE)
  n_surv <- sum(meta$cancer == "SURV", na.rm = TRUE)
  
  t1_only <- t1_count - paired_subjects
  t2_only <- t2_count - paired_subjects
  
  # --- Build row ---
  tibble(
    dataset = ifelse(str_detect(nm, "^SURV"), "SURV", "CRE"),
    filter  = str_extract(nm, "none|moderate|stringent"),
    ratio   = str_extract(nm, "BM[^_]+") |> str_remove("BM"),
    
    total_samples = nrow(meta),
    n_CRE = CRE_count,
    n_PLA = PLA_count,
    n_male = m_count,
    n_female = f_count,
    
    n_male_T1 = m_T1_count,
    n_male_T2 = m_T2_count,
    n_female_T1 = f_T1_count,
    n_female_T2 = f_T2_count,
    
    n_CRE_T1 = CRE_T1_count,
    n_CRE_T2 = CRE_T2_count,
    n_PLA_T1 = PLA_T1_count,
    n_PLA_T2 = PLA_T2_count,
    
    n_T1 = t1_count,
    n_T2 = t2_count,
    
    n_paired = paired_subjects,
    n_paired_CRE = paired_cre,
    n_paired_PLA = paired_pla,
    n_unpaired = unpaired_subjects,
    
    n_T1_only = t1_only,
    n_T2_only = t2_only,
    
    n_SURV = n_surv,
    n_CTL = n_ctl
  )
})

# Clean + order factors (optional but useful)
sample_summary <- sample_summary %>%
  mutate(
    filter = factor(filter, levels = c("none", "moderate", "stringent")),
    ratio  = factor(ratio, levels = c("none", "5", "10", "15", "20"))
  )

# Save combined file
write.csv(sample_summary,
          file.path(config$report_dir, "12_sample_report_all_muscle.csv"),
          row.names = FALSE)

# Print preview
print(sample_summary)
