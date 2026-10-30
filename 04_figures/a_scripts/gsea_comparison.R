# GSEA Pathway Comp

# --- 0: Setup ----------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, patchwork, cowplot, scales, enrichplot, grid,
               clusterProfiler, org.Hs.eg.db)

# --- Configuration -----------------------------------------------------------
config <- list(
  limma_muscle = "03_analysis/b_reports/muscle/all_data/limma",
  limma_plasma = "03_analysis/b_reports/plasma/all_data/limma",
  gsea_muscle  = "03_analysis/b_reports/muscle/SURVvCTL/gsea",
  gsea_plasma  = "03_analysis/b_reports/plasma/SURVvCTL/gsea",
  fig_dir      = "04_figures",
  data_dir     = "04_figures/data"
)
dir.create(config$fig_dir,  showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir, showWarnings = FALSE, recursive = TRUE)

# --- Color palette + theme ----------------------------------------------------
# This script needs two palettes (std = the shared pipeline palette, brt = a
# brighter variant for a specific overlay), so it nests them as colors$std$*
# and colors$brt$* instead of the flat colors$* used elsewhere.
source("04_figures/a_scripts/00_theme.R")

colors <- list(
  std = colors,
  brt = list(
    red        = "#E74C3C", red_light   = "#F1948A", red_dark   = "#922B21",
    blue       = "#2E86C1", blue_light  = "#85C1E9", blue_dark  = "#1A5276",
    green      = "#27AE60", green_light = "#82E0AA", green_dark = "#1D6A3A",
    pink       = "#E91E8C", pink_light  = "#F48FB1", pink_dark  = "#880E4F",
    purple     = "#8E44AD", purple_light= "#C39BD3", purple_dark= "#4A235A",
    yellow     = "#F1C40F", yellow_light= "#F9E79F", yellow_dark= "#B7950B",
    orange     = "#E67E22", orange_light= "#F0B27A", orange_dark= "#935116",
    aqua       = "#1ABC9C", aqua_light  = "#76D7C4", aqua_dark  = "#0E6655",
    brown      = "#A04000", brown_light = "#D7BEA8", brown_dark = "#6E2C00",
    red_pale   = "#FADBD8", blue_pale   = "#D6EAF8",
    gray_light = "#F5F5F5", gray        = "#9E9E9E", gray_dark  = "#4A4A4A",
    black      = "#000000", white       = "#FFFFFF"
  )
)

tissue_colors <- c(muscle = colors$std$pink_dark, plasma = colors$std$aqua_dark)

# --- Shared constants --------------------------------------------------------
adjp_thresh <- 0.05
pi_thresh   <- 0.05

comp_labels <- c(
  Baseline_SURVvCTL = "CRC_PRE\nvs\nCTL",
  Training_SURVvCTL = "CRC_POST\nvs\nCTL",
  Training_SURV     = "CRC_POST\nvs\nCRC_PRE"
)

# --- 1: Load Data ------------------------------------------------------------
load_gsea_rds <- function(dir) {
  files <- list.files(dir, pattern = "^gsea_results_.*\\.rds$", full.names = TRUE)
  nms   <- gsub("^gsea_results_|\\.rds$", "",
                list.files(dir, pattern = "^gsea_results_.*\\.rds$"))
  setNames(lapply(files, readRDS), nms)
}

gsea_muscle <- load_gsea_rds(config$gsea_muscle)
gsea_plasma <- load_gsea_rds(config$gsea_plasma)

limma_results <- list(
  muscle = list(
    Baseline_SURVvCTL = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle.xlsx"), sheet = 11),
    Training_SURVvCTL = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle.xlsx"), sheet = 12),
    Training_SURV     = read_excel(file.path(config$limma_muscle,
                                             "limma_results_all_muscle.xlsx"), sheet = 13)
  ),
  plasma = list(
    Baseline_SURVvCTL = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"), sheet = 11),
    Training_SURVvCTL = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"), sheet = 12),
    Training_SURV     = read_excel(file.path(config$limma_plasma,
                                             "limma_results_all_plasma.xlsx"), sheet = 13)
  )
)

# --- 2: Shared helper — gene symbol lookup -----------------------------------
id_to_gene <- setNames(
  limma_results$muscle$Baseline_SURVvCTL$gene,
  gsub("\\..*$", "", limma_results$muscle$Baseline_SURVvCTL$Ensembl)
)

to_symbols <- function(ids) {
  syms <- id_to_gene[ids]
  syms[is.na(syms)] <- ids[is.na(syms)]
  paste(syms, collapse = ", ")
}

# --- 3: Core functions -------------------------------------------------------
get_leading_edge <- function(gsea_result, comp_label, pathway) {
  df <- as.data.frame(gsea_result@result)
  hit <- df |>
    dplyr::filter(Description == pathway | ID == pathway) |>
    dplyr::slice(1)
  if (nrow(hit) == 0) {
    message("Pathway '", pathway, "' not found in ", comp_label)
    return(NULL)
  }
  proteins <- strsplit(hit$core_enrichment, "/")[[1]]
  tibble(
    comparison     = comp_label,
    pathway_id     = hit$ID,
    description    = hit$Description,
    NES            = hit$NES,
    p.adjust       = hit$p.adjust,
    n_leading_edge = length(proteins),
    protein        = list(proteins)
  )
}

get_da_stats <- function(ensembl_ids,
                         comparisons = c("Baseline_SURVvCTL", "Training_SURV")) {
  bind_rows(lapply(comparisons, function(comp) {
    limma_results$muscle[[comp]] |>
      dplyr::filter(gsub("\\..*$", "", Ensembl) %in% ensembl_ids) |>
      dplyr::select(gene, Ensembl, logFC, adj.P.Val, adj.Pi.Val) |>
      dplyr::mutate(comparison = comp)
  }))
}

# --- 4: Single pathway analysis ----------------------------------------------
pathway <- "actin cytoskeleton organization"

le_baseline <- get_leading_edge(gsea_muscle$Baseline_SURVvCTL$gsea_bp_clean,
                                "Baseline_SURVvCTL", pathway)
le_training <- get_leading_edge(gsea_muscle$Training_SURV$gsea_bp_clean,
                                "Training_SURV",     pathway)

proteins_base  <- gsub("\\..*$", "", le_baseline$protein[[1]])
proteins_train <- gsub("\\..*$", "", le_training$protein[[1]])

overlap    <- intersect(proteins_base,  proteins_train)
only_base  <- setdiff(proteins_base,   proteins_train)
only_train <- setdiff(proteins_train,  proteins_base)

cat("Shared:        ", length(overlap),    "\n")
cat("Baseline only: ", length(only_base),  "\n")
cat("Training only: ", length(only_train), "\n")
cat("Jaccard:       ",
    round(length(overlap) / length(union(proteins_base, proteins_train)), 3), "\n")

overlap_summary <- tibble(
  group    = c("Shared", "Baseline only", "Training only"),
  n        = c(length(overlap), length(only_base), length(only_train)),
  proteins = c(to_symbols(overlap), to_symbols(only_base), to_symbols(only_train))
)
print(overlap_summary)

all_stats <- bind_rows(
  get_da_stats(overlap)    |> dplyr::mutate(group = "Shared"),
  get_da_stats(only_base)  |> dplyr::mutate(group = "Baseline only"),
  get_da_stats(only_train) |> dplyr::mutate(group = "Training only")
)

# --- 5: Heatmap --------------------------------------------------------------
lfc_wide <- all_stats |>
  dplyr::select(gene, comparison, logFC, group) |>
  pivot_wider(names_from = comparison, values_from = logFC) |>
  dplyr::mutate(group = factor(group,
                               levels = c("Shared", "Baseline only", "Training only"))) |>
  arrange(group, desc(Baseline_SURVvCTL))

lfc_long <- lfc_wide |>
  pivot_longer(cols = c(Baseline_SURVvCTL, Training_SURV),
               names_to = "comparison", values_to = "logFC")

p_heatmap <- ggplot(lfc_long,
                    aes(x = comparison, y = reorder(gene, logFC), fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = round(logFC, 2)), size = 2.5, fontface = "bold") +
  scale_fill_gradientn(
    colors = c(colors$brt$blue, colors$brt$white, colors$brt$red),
    values = scales::rescale(c(
      min(lfc_long$logFC, na.rm = TRUE), 0, max(lfc_long$logFC, na.rm = TRUE)
    )),
    name = "log2FC"
  ) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y") +
  scale_x_discrete(labels = c(
    Baseline_SURVvCTL = "Baseline\nSURV vs CTL",
    Training_SURV     = "Training\nSURV"
  )) +
  labs(x = NULL, y = NULL, title = paste("Leading edge —", pathway)) +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(strip.text.y = element_text(angle = 0), legend.position = "right")

p_heatmap

# --- 6: Volcano --------------------------------------------------------------
p_volcano <- ggplot(all_stats,
                    aes(x = logFC, y = -log10(adj.P.Val),
                        color = comparison, shape = group)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0,            linetype = "dashed", color = "grey50") +
  ggrepel::geom_text_repel(aes(label = gene), size = 2.5, max.overlaps = 15) +
  scale_color_manual(values = c(
    Baseline_SURVvCTL = colors$std$blue,
    Training_SURV     = colors$std$pink
  )) +
  labs(x = "log2FC", y = "-log10(adj.p)",
       title = paste("Leading edge —", pathway),
       color = "Comparison", shape = "Group") +
  theme_cr(base_size = 10, title_hjust = 0)

p_volcano

# --- 7: ORA on leading edge subgroups ----------------------------------------
run_quick_ora <- function(ensembl_ids, label) {
  entrez <- bitr(ensembl_ids, fromType = "ENSEMBL",
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez) == 0) { message("No IDs mapped for ", label); return(NULL) }
  res <- enrichGO(gene = entrez$ENTREZID, OrgDb = org.Hs.eg.db,
                  ont = "BP", pAdjustMethod = "BH",
                  pvalueCutoff = 0.05, readable = TRUE)
  if (is.null(res) || nrow(as.data.frame(res)) == 0) {
    message("No enrichment for ", label); return(NULL)
  }
  as.data.frame(res) |> dplyr::mutate(group = label)
}

ora_base   <- run_quick_ora(only_base,  "Baseline only")
ora_train  <- run_quick_ora(only_train, "Training only")
ora_shared <- run_quick_ora(overlap,    "Shared")

ora_comb <- bind_rows(ora_base, ora_train, ora_shared) |>
  dplyr::select(group, Description, GeneRatio, p.adjust, geneID) |>
  arrange(group, p.adjust)
print(ora_comb)

# --- 8: Statistical tests ----------------------------------------------------
# Paired Wilcoxon on shared proteins
shared_wide <- all_stats |>
  dplyr::filter(group == "Shared") |>
  dplyr::select(gene, comparison, logFC) |>
  pivot_wider(names_from = comparison, values_from = logFC)

wilcox.test(shared_wide$Baseline_SURVvCTL, shared_wide$Training_SURV,
            paired = TRUE, alternative = "greater")

# NES comparison
pathway_stats <- bind_rows(
  as.data.frame(gsea_muscle$Baseline_SURVvCTL$gsea_bp_clean@result) |>
    dplyr::filter(Description == pathway) |>
    dplyr::mutate(comparison = "Baseline_SURVvCTL"),
  as.data.frame(gsea_muscle$Training_SURV$gsea_bp_clean@result) |>
    dplyr::filter(Description == pathway) |>
    dplyr::mutate(comparison = "Training_SURV")
) |>
  dplyr::select(comparison, NES, pvalue, p.adjust, setSize, leading_edge)
print(pathway_stats)

# Jaccard permutation test
observed_jaccard    <- length(overlap) / length(union(proteins_base, proteins_train))
all_pathway_prots   <- union(proteins_base, proteins_train)
set.seed(42)
perm_jaccards <- replicate(10000, {
  pb <- sample(all_pathway_prots, length(proteins_base))
  pt <- sample(all_pathway_prots, length(proteins_train))
  length(intersect(pb, pt)) / length(union(pb, pt))
})
p_perm <- mean(perm_jaccards <= observed_jaccard)
cat("Observed Jaccard: ", round(observed_jaccard, 3), "\n")
cat("Expected Jaccard: ", round(mean(perm_jaccards), 3), "\n")
cat("Permutation p:    ", p_perm, "\n")

# Fisher's exact test
all_detected <- nrow(limma_results$muscle$Baseline_SURVvCTL)
mat <- matrix(c(
  length(overlap),
  length(proteins_base)  - length(overlap),
  length(proteins_train) - length(overlap),
  all_detected - length(proteins_base) - length(proteins_train) + length(overlap)
), nrow = 2)
fisher.test(mat)

# --- 9: Sequenced analysis across all overlapping pathways -------------------
get_sig_pathways <- function(comp) {
  as.data.frame(gsea_muscle[[comp]]$gsea_bp_clean@result) |>
    dplyr::filter(p.adjust < adjp_thresh) |>
    dplyr::pull(ID)
}

overlapping_pathways <- intersect(
  get_sig_pathways("Baseline_SURVvCTL"),
  get_sig_pathways("Training_SURV")
)
cat("Overlapping significant pathways:", length(overlapping_pathways), "\n")

pathway_descriptions <- as.data.frame(gsea_muscle$Baseline_SURVvCTL$gsea_bp_clean@result) |>
  dplyr::select(ID, Description) |>
  distinct()

compare_leading_edge <- function(pathway_id) {
  desc <- pathway_descriptions$Description[pathway_descriptions$ID == pathway_id]
  if (length(desc) == 0) desc <- pathway_id
  
  le_base  <- get_leading_edge(gsea_muscle$Baseline_SURVvCTL$gsea_bp_clean,
                               "Baseline_SURVvCTL", pathway_id)
  le_train <- get_leading_edge(gsea_muscle$Training_SURV$gsea_bp_clean,
                               "Training_SURV",     pathway_id)
  
  if (is.null(le_base) | is.null(le_train)) {
    message("Skipping ", pathway_id)
    return(NULL)
  }
  
  pb <- gsub("\\..*$", "", le_base$protein[[1]])
  pt <- gsub("\\..*$", "", le_train$protein[[1]])
  
  ov  <- intersect(pb, pt)
  ob  <- setdiff(pb, pt)
  ot  <- setdiff(pt, pb)
  
  all_stats <- bind_rows(
    get_da_stats(ov, c("Baseline_SURVvCTL", "Training_SURV")) |>
      dplyr::mutate(group = "Shared"),
    get_da_stats(ob, c("Baseline_SURVvCTL", "Training_SURV")) |>
      dplyr::mutate(group = "Baseline only"),
    get_da_stats(ot, c("Baseline_SURVvCTL", "Training_SURV")) |>
      dplyr::mutate(group = "Training only")
  )
  
  if (nrow(all_stats) == 0) return(NULL)
  
  lfc_long <- all_stats |>
    dplyr::select(gene, comparison, logFC, group) |>
    dplyr::mutate(group = factor(group,
                                 levels = c("Shared", "Baseline only", "Training only")))
  
  list(
    pathway_id   = pathway_id,
    description  = desc,
    jaccard      = round(length(ov) / length(union(pb, pt)), 3),
    n_shared     = length(ov),
    n_base_only  = length(ob),
    n_train_only = length(ot),
    lfc_long     = lfc_long,
    stats        = all_stats
  )
}

plot_leading_edge_heatmap <- function(result) {
  if (is.null(result)) return(NULL)
  
  lfc_long <- result$lfc_long
  
  gene_order <- lfc_long |>
    dplyr::filter(comparison == "Baseline_SURVvCTL") |>
    arrange(group, desc(logFC)) |>
    dplyr::pull(gene)
  
  lfc_long <- lfc_long |>
    dplyr::mutate(
      gene       = factor(gene, levels = gene_order),
      comparison = factor(comparison,
                          levels = c("Baseline_SURVvCTL", "Training_SURV"),
                          labels = c("Baseline\nSURV vs CTL", "Training\nSURV"))
    )
  
  ggplot(lfc_long, aes(x = comparison, y = gene, fill = logFC)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = round(logFC, 2)), size = 2.5, fontface = "bold") +
    scale_fill_gradientn(
      colors = c(colors$brt$blue, colors$brt$white, colors$brt$red),
      values = scales::rescale(c(
        min(lfc_long$logFC, na.rm = TRUE), 0, max(lfc_long$logFC, na.rm = TRUE)
      )),
      name = "log2FC"
    ) +
    facet_grid(group ~ ., scales = "free_y", space = "free_y") +
    labs(
      x        = NULL,
      y        = NULL,
      title    = result$description,
      subtitle = sprintf(
        "Jaccard = %.3f | Shared: %d | Baseline only: %d | Training only: %d",
        result$jaccard, result$n_shared, result$n_base_only, result$n_train_only
      )
    ) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      strip.text.y    = element_text(angle = 0),
      legend.position = "right",
      plot.subtitle   = element_text(size = 8, face = "italic")
    )
}

# Run
pathway_results <- lapply(overlapping_pathways, compare_leading_edge)
names(pathway_results) <- overlapping_pathways
pathway_results <- Filter(Negate(is.null), pathway_results)
cat("Pathways successfully processed:", length(pathway_results), "\n")

# Save heatmaps
output_dir <- file.path(config$fig_dir, "leading_edge_heatmaps")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

lapply(names(pathway_results), function(pid) {
  result <- pathway_results[[pid]]
  p      <- plot_leading_edge_heatmap(result)
  if (is.null(p)) return(NULL)
  fname  <- substr(gsub("[^A-Za-z0-9_]", "_", result$description), 1, 60)
  ggsave(
    file.path(output_dir, paste0(fname, ".pdf")),
    p,
    width  = 5,
    height = max(4, nrow(result$lfc_long) / 2 * 0.2 + 2),
    units  = "in",
    limitsize = FALSE
  )
})
cat("Saved", length(pathway_results), "heatmaps to", output_dir, "\n")

# Summary table
pathway_summary <- bind_rows(lapply(pathway_results, function(r) {
  tibble(
    pathway_id   = r$pathway_id,
    description  = r$description,
    jaccard      = r$jaccard,
    n_shared     = r$n_shared,
    n_base_only  = r$n_base_only,
    n_train_only = r$n_train_only
  )
})) |> arrange(jaccard)

print(pathway_summary)

class_counts <- nes_df_three_way |>
  filter(padj_PRE < 0.05 | padj_POST < 0.05) |>
  count(normalization_class, go_set) |>
  filter(normalization_class != "Not dysregulated at baseline")

ggplot(class_counts, aes(x = reorder(normalization_class, n),
                         y = n, fill = normalization_class)) +
  geom_col(color = "black", linewidth = 0.3) +
  geom_text(aes(label = n, y = n + 0.3), 
            vjust = 0, hjust = 0, size = 3, fontface = "bold") +
  scale_fill_manual(values = c(
    "Normalized"            = colors$green,
    "Partial normalization"  = colors$green_light,
    "Exacerbated"            = colors$red,
    "Overcorrected"          = colors$orange,
    "Newly dysregulated"     = colors$purple,
    "Other"                  = colors$gray
  )) +
  facet_wrap(~ go_set, scales = "free_x") +
  coord_flip() +
  labs(x = NULL, y = "Number of pathways", title = "Pathway normalization summary") +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(legend.position = "none")

# --- REVERSAL 2 ---------------
nes_df_reversal <- nes_df_bp %>%
  mutate(
    # Magnitude of change
    delta_NES = NES_RET - NES_PRE,
    
    # Classify trajectory more precisely
    reversal_class = case_when(
      # Full reversal: significant PRE, flipped direction, significant RET
      padj_PRE < 0.05 & padj_RET < 0.05 & 
        sign(NES_PRE) != sign(NES_RET)                          ~ "Full reversal",
      
      # Partial reversal: significant PRE, flipped direction, NS in RET
      padj_PRE < 0.05 & padj_RET >= 0.05 & 
        sign(NES_PRE) != sign(NES_RET)                          ~ "Partial reversal",
      
      # Attenuation: same direction, significant PRE, reduced magnitude
      padj_PRE < 0.05 & 
        sign(NES_PRE) == sign(NES_RET) & 
        abs(NES_RET) < abs(NES_PRE)                             ~ "Attenuated",
      
      # Exacerbation: same direction, significant PRE, increased magnitude
      padj_PRE < 0.05 & 
        sign(NES_PRE) == sign(NES_RET) & 
        abs(NES_RET) > abs(NES_PRE)                             ~ "Exacerbated",
      
      # Training only: NS PRE, significant RET
      padj_PRE >= 0.05 & padj_RET < 0.05                       ~ "Training induced",
      
      TRUE                                                      ~ "Other"
    ),
    
    # Percent change in NES magnitude
    pct_change = ((abs(NES_RET) - abs(NES_PRE)) / abs(NES_PRE)) * 100
  )

# Join NES_RET onto the PRE vs POST data
nes_df_three_way <- nes_df_bp_post %>%
  left_join(
    nes_df_bp %>% dplyr::select(ID, NES_RET, padj_RET, trajectory),
    by = "ID"
  ) %>%
  mutate(
    # Direction of training effect
    training_direction = case_when(
      is.na(NES_RET)               ~ "Not tested",
      NES_RET > 0                  ~ "Up with training",
      NES_RET < 0                  ~ "Down with training"
    ),
    # Whether training moved toward or away from CTL
    # (PRE is the cancer phenotype — if PRE is up vs CTL and training is down, that's normalizing)
    training_effect = case_when(
      is.na(NES_RET)                                     ~ "Not tested",
      sign(NES_PRE) != sign(NES_RET)                     ~ "Normalizing",
      sign(NES_PRE) == sign(NES_RET) & 
        abs(NES_RET) > abs(NES_PRE) * 0.1                ~ "Exacerbating",
      TRUE                                               ~ "Neutral"
    )
  )

ggplot(nes_df_three_way,
       aes(x = NES_PRE, y = NES_POST, color = training_effect)) +
  geom_point(alpha = 0.8, size = 3.5) +
  scale_color_manual(values = c(
    "Normalizing"  = colors$aqua,
    "Exacerbating" = colors$red,
    "Not tested"   = colors$gray
  )) +
  theme_cr(base_size = 10, title_hjust = 0)

table(nes_df_three_way$training_effect, useNA = "always")
table(nes_df_three_way$go_set, useNA = "always")

nes_df_three_way <- nes_df_three_way |>
  mutate(
    normalization_class = case_when(
      # Not significant at baseline — training effect not relevant to attenuation
      is.na(padj_PRE) | padj_PRE >= 0.05                          ~ "Not sig. at BL",
      
      # Significant PRE, not significant POST — gap closed regardless of direction
      padj_PRE < 0.05 & (is.na(padj_POST) | padj_POST >= 0.05)   ~ "Attenuated",
      
      # Significant in both, same direction, smaller NES POST — partial attenuation
      padj_PRE < 0.05 & padj_POST < 0.05 &
        sign(NES_PRE) == sign(NES_POST) &
        abs(NES_POST) < abs(NES_PRE)                              ~ "Partial attenuation",
      
      # Significant in both, same direction, larger NES POST — gap widened
      padj_PRE < 0.05 & padj_POST < 0.05 &
        sign(NES_PRE) == sign(NES_POST) &
        abs(NES_POST) > abs(NES_PRE)                              ~ "Exacerbated",
      
      # Significant in both, opposite direction — reversal
      padj_PRE < 0.05 & padj_POST < 0.05 &
        sign(NES_PRE) != sign(NES_POST)                           ~ "Reversed",
      
      # Not significant PRE, significant POST — newly dysregulated
      padj_PRE >= 0.05 & padj_POST < 0.05                        ~ "Newly dysregulated",
      
      TRUE                                                        ~ "Other"
    )
  )

table(nes_df_three_way$normalization_class)

ggplot(nes_df_three_way |>
         filter(!is.na(NES_PRE), !is.na(NES_POST)),
       aes(x = NES_PRE, y = NES_POST, color = normalization_class)) +
  geom_point(
    aes(size = -log10(pmin(padj_PRE, padj_POST, na.rm = TRUE))),
    alpha = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = colors$gray) +
  geom_vline(xintercept = 0, linetype = "dashed", color = colors$gray) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = colors$gray) +
  ggrepel::geom_text_repel(
    data = \(d) filter(d, padj_PRE < 0.05 | padj_POST < 0.05,
                       normalization_class %in% c("Normalized",
                                                  "Partial normalization",
                                                  "Overcorrected")),
    aes(label = Description),
    size = 2.5, max.overlaps = 15
  ) +
  scale_color_manual(values = c(
    "Attenuated"                  = colors$green,
    "Partial attenuation"        = colors$green_light,
    "Exacerbated"                  = colors$red,
    "Reversed"                = colors$orange,
    "Newly dysregulated"           = colors$purple,
    "Not sig. at BL" = colors$gray,
    "Other"                        = colors$gray
  )) +
  scale_size_continuous(name = "-log10(best adj.p)", range = c(1, 5)) +
  labs(
    x     = "NES (CRC_PRE vs CTL)",
    y     = "NES (CRC_POST vs CTL)",
    color = NULL,
    title = "Training-induced normalization of pathway dysregulation"
  ) +
  theme_cr(base_size = 10, title_hjust = 0)

# Filter to pathways significant in at least one comparison
dumbbell_data <- nes_df_three_way |>
  filter(padj_PRE < 0.05 | padj_POST < 0.05) |>
  filter(!is.na(NES_PRE), !is.na(NES_POST)) |>
  arrange(normalization_class, NES_PRE) |>
  mutate(Description = factor(Description, levels = Description))

ggplot(dumbbell_data) +
  geom_segment(aes(x = NES_PRE, xend = NES_POST,
                   y = Description, yend = Description,
                   color = normalization_class),
               linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(x = NES_PRE,  y = Description), 
             shape = 21, fill = "white", size = 2.5, stroke = 0.8) +
  geom_point(aes(x = NES_POST, y = Description,
                 color = normalization_class), 
             size = 2.5) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c(
    "Normalized"                   = colors$green,
    "Partial normalization"         = colors$green_light,
    "Exacerbated"                   = colors$red,
    "Overcorrected"                 = colors$orange,
    "Newly dysregulated"            = colors$purple,
    "Not dysregulated at baseline"  = colors$gray,
    "Other"                         = colors$gray_dark
  )) +
  facet_grid(normalization_class ~ ., scales = "free_y", space = "free_y") +
  labs(x = "NES", y = NULL, color = NULL,
       title = "Pathway normalization after training") +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(strip.text.y   = element_text(angle = 0),
        legend.position = "none")

# Join in Training_SURV NES
heatmap_data <- nes_df_three_way |>
  filter(padj_PRE < 0.05 | padj_POST < 0.05) |>
  filter(!is.na(NES_PRE), !is.na(NES_POST)) |>
  dplyr::select(ID, Description, NES_PRE, NES_POST, NES_RET, normalization_class) |>
  pivot_longer(cols = c(NES_PRE, NES_POST, NES_RET),
               names_to = "comparison", values_to = "NES") |>
  mutate(
    comparison = factor(comparison,
                        levels = c("NES_PRE", "NES_POST", "NES_RET"),
                        labels = c("PRE vs CTL", "POST vs CTL", "POST vs PRE")),
    normalization_class = factor(normalization_class,
                                 levels = c("Normalized", "Partial normalization",
                                            "Exacerbated", "Overcorrected",
                                            "Newly dysregulated", "Other"))
  ) |>
  arrange(normalization_class, NES) |>
  mutate(Description = factor(Description, levels = unique(Description)))

ggplot(heatmap_data, aes(x = comparison, y = Description, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = c(colors$blue, colors$white, colors$red),
    values = scales::rescale(c(min(heatmap_data$NES, na.rm = TRUE),
                               0,
                               max(heatmap_data$NES, na.rm = TRUE))),
    na.value = "grey90",
    name = "NES"
  ) +
  facet_grid(normalization_class ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Pathway trajectory across comparisons") +
  theme_cr(base_size = 10, title_hjust = 0) +
  theme(
    strip.text.y    = element_text(angle = 0, size = 8),
    axis.text.y     = element_text(size = 7),
    legend.position = "right"
  )

# GSEA PATHWAY UPSETS

# Build pathway membership sets per database
make_pathway_upset <- function(gsea_obj, db_slot, db_label, adjp_thresh = 0.05) {
  
  comps <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")
  
  sig_sets <- lapply(comps, function(comp) {
    res <- gsea_obj[[comp]][[db_slot]]
    if (is.null(res)) return(character(0))
    as.data.frame(res) |>
      filter(p.adjust < adjp_thresh) |>
      pull(ID)
  })
  names(sig_sets) <- comps
  
  # Reuse your existing upset infrastructure
  comb_grid <- expand.grid(
    Baseline_SURVvCTL = c(TRUE, FALSE),
    Training_SURVvCTL = c(TRUE, FALSE),
    Training_SURV     = c(TRUE, FALSE)
  ) |> filter(rowSums(across(everything())) > 0)
  
  pmap_dfr(comb_grid, function(Baseline_SURVvCTL, Training_SURVvCTL, Training_SURV) {
    membership <- c(Baseline_SURVvCTL, Training_SURVvCTL, Training_SURV)
    in_sets    <- comps[membership]
    out_sets   <- comps[!membership]
    members    <- Reduce(intersect, sig_sets[in_sets])
    if (length(out_sets) > 0)
      members <- setdiff(members, unique(unlist(sig_sets[out_sets])))
    if (length(members) == 0) return(NULL)
    
    # NES direction per pathway per comparison
    dirs <- sapply(members, function(pid) {
      nes_vals <- sapply(in_sets, function(comp) {
        res <- gsea_obj[[comp]][[db_slot]]
        if (is.null(res)) return(NA)
        df <- as.data.frame(res)
        df$NES[df$ID == pid]
      })
      nes_vals <- nes_vals[!is.na(nes_vals)]
      if (all(nes_vals > 0)) "Up" else if (all(nes_vals < 0)) "Down" else "Mixed"
    })
    
    tibble(
      in_sets  = list(in_sets),
      degree   = length(in_sets),
      up       = sum(dirs == "Up"),
      down     = sum(dirs == "Down"),
      mixed    = sum(dirs == "Mixed"),
      total    = length(members),
      database = db_label
    )
  }) |>
    arrange(desc(total)) |>
    mutate(x = row_number())
}

# Run per database
upset_gobp     <- make_pathway_upset(gsea_muscle, "gsea_bp_simple", "GO:BP")
upset_hallmark <- make_pathway_upset(gsea_muscle, "gsea_hall",      "Hallmark")
upset_reactome <- make_pathway_upset(gsea_muscle, "gsea_reactome",  "Reactome")
upset_kegg     <- make_pathway_upset(gsea_muscle, "gsea_kegg",      "KEGG")


# NES HEATMAP
make_nes_heatmap <- function(gsea_obj, db_slot, adjp_thresh = 0.05) {
  
  comps <- c("Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV")
  
  all_sig <- unique(unlist(lapply(comps, function(comp) {
    res <- gsea_obj[[comp]][[db_slot]]
    if (is.null(res)) return(character(0))
    as.data.frame(res) |> filter(p.adjust < adjp_thresh) |> pull(ID)
  })))
  
  if (length(all_sig) == 0) stop("No significant pathways found in ", db_slot)
  
  nes_df <- bind_rows(lapply(comps, function(comp) {
    res <- gsea_obj[[comp]][[db_slot]]
    if (is.null(res)) return(NULL)
    as.data.frame(res) |>
      filter(ID %in% all_sig) |>
      dplyr::select(ID, Description, NES, p.adjust) |>
      mutate(comparison = comp, sig = p.adjust < adjp_thresh)
  }))
  
  sig_class <- nes_df |>
    group_by(ID, Description) |>
    summarise(
      sig_in = paste(sort(comparison[sig]), collapse = " + "),
      .groups = "drop"
    ) |>
    mutate(
      sig_group = case_when(
        sig_in == "Baseline_SURVvCTL"                                    ~ "Cancer only",
        sig_in == "Training_SURV"                                        ~ "Training only",
        sig_in == "Training_SURVvCTL"                                    ~ "POST vs CTL only",
        sig_in == "Baseline_SURVvCTL + Training_SURV"                    ~ "Cancer + Training",
        sig_in == "Baseline_SURVvCTL + Training_SURVvCTL"                ~ "Cancer + POST",
        sig_in == "Training_SURVvCTL + Training_SURV"                    ~ "POST + Training",
        sig_in == "Baseline_SURVvCTL + Training_SURVvCTL + Training_SURV"~ "All three",
        TRUE                                                              ~ "Other"
      ),
      sig_group = factor(sig_group, levels = c(
        "Cancer only", "Training only", "POST vs CTL only",
        "Cancer + Training", "Cancer + POST", "POST + Training", "All three"
      ))
    )
  
  plot_df <- nes_df |>
    left_join(sig_class, by = c("ID", "Description")) |>
    mutate(
      comparison = factor(comparison,
                          levels = comps,
                          labels = c("Cancer\n(PRE vs CTL)",
                                     "POST\nvs CTL",
                                     "Training\n(POST vs PRE)")),
      # Use NES for sig, 0 for NS so color scale still works
      NES_fill = ifelse(sig, NES, 0)
    )
  
  nes_range <- max(abs(plot_df$NES), na.rm = TRUE)
  
  pathway_order <- plot_df |>
    filter(comparison == "Cancer\n(PRE vs CTL)") |>
    arrange(sig_group, desc(NES)) |>
    pull(Description)
  
  plot_df <- plot_df |>
    mutate(Description = factor(Description, levels = unique(pathway_order)))
  
  cat("Sig rows:", sum(plot_df$sig), "\n")
  cat("NES range:", range(plot_df$NES[plot_df$sig], na.rm = TRUE), "\n")
  cat("nes_range value:", nes_range, "\n")
  
  ggplot(plot_df, aes(x = comparison, y = Description)) +
    # NS tiles first — grey background
    geom_tile(data = filter(plot_df, !sig),
              fill = "grey88", color = "white", linewidth = 0.2) +
    # Significant tiles on top
    geom_tile(data = filter(plot_df, sig),
              aes(fill = NES), color = "white", linewidth = 0.2) +
    scale_fill_gradientn(
      colors   = c(colors$brt$blue, colors$brt$white, colors$brt$red),
      values   = scales::rescale(c(-nes_range, 0, nes_range)),
      limits   = c(-nes_range, nes_range),
      name     = "NES"
    ) +
    facet_grid(sig_group ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_cr(base_size = 10, title_hjust = 0) +
    theme(
      strip.text.y    = element_text(angle = 0, size = 8),
      axis.text.y     = element_text(size = 7),
      legend.position = "right"
    )
}

p_nes_heatmap_gobp     <- make_nes_heatmap(gsea_muscle, "gsea_bp_simple")
p_nes_heatmap_hallmark <- make_nes_heatmap(gsea_muscle, "gsea_hall")
