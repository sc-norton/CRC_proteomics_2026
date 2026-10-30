# Muscle GSEA Quadrant Plot: PRE vs CTL (x) vs POST vs PRE (y)
setwd(rprojroot::find_rstudio_root_file())
# --- 0: Setup --------------------------
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, ggplot2, ggrepel, stringr)

# --- 1: Load & filter -------------------
xlsx_path <- file.path("03_analysis", "c_data", "blood_covar",
                       "gsea_consolidated", "gsea_pathway_data_supplement.xlsx")  

muscle <- read_excel(xlsx_path, sheet = "Muscle")

muscle_sub <- muscle %>%
  filter(database %in% c("GO Slim", "Hallmark"),
         comparison %in% c("PRE vs CTL", "POST vs PRE"))

# --- 2: Pivot to wide (one row per pathway) ---
muscle_sub <- muscle_sub %>%
  mutate(comparison_clean = str_replace_all(comparison, " ", "_"))  # "PRE vs CTL" -> "PRE_vs_CTL"

plot_df <- muscle_sub %>%
  select(database, ID, Description, comparison_clean, NES, p.adjust, significant) %>%
  pivot_wider(
    id_cols     = c(database, ID, Description),
    names_from  = comparison_clean,
    values_from = c(NES, p.adjust, significant)
  )

# check what got created before renaming
print(colnames(plot_df))

# --- 2: Pivot to wide (one row per pathway) ---
# "significant in at least one comparison" is evaluated across the two
# comparisons being plotted (PRE vs CTL, POST vs PRE) — adjust the filter
# below if you'd rather include POST vs CTL significance too.
plot_df <- plot_df %>%
  dplyr::rename(
    NES_PRE_CTL   = NES_PRE_vs_CTL,
    NES_POST_PRE  = NES_POST_vs_PRE,
    padj_PRE_CTL  = p.adjust_PRE_vs_CTL,
    padj_POST_PRE = p.adjust_POST_vs_PRE,
    sig_PRE_CTL   = significant_PRE_vs_CTL,
    sig_POST_PRE  = significant_POST_vs_PRE
  ) %>%
  dplyr::filter(sig_PRE_CTL | sig_POST_PRE) %>%
  dplyr::mutate(
    sig_status = dplyr::case_when(
      sig_PRE_CTL & sig_POST_PRE  ~ "Both",
      sig_PRE_CTL & !sig_POST_PRE ~ "PRE vs CTL only",
      !sig_PRE_CTL & sig_POST_PRE ~ "POST vs PRE only"
    ),
    quadrant = dplyr::case_when(
      NES_PRE_CTL > 0 & NES_POST_PRE > 0 ~ "Concordant up",
      NES_PRE_CTL < 0 & NES_POST_PRE < 0 ~ "Concordant down",
      TRUE ~ "Discordant"
    )
  )

cat(sprintf("Pathways plotted: %d (GO Slim: %d, Hallmark: %d)\n",
            nrow(plot_df),
            sum(plot_df$database == "GO Slim"),
            sum(plot_df$database == "Hallmark")))

# --- 3: Plot -----------------------------
p <- ggplot(plot_df, aes(x = NES_PRE_CTL, y = NES_POST_PRE)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(aes(shape = database, color = sig_status), size = 3, alpha = 0.85) +
  geom_text_repel(
    data = plot_df %>%
      filter(abs(NES_PRE_CTL) > 1.5 | abs(NES_POST_PRE) > 1.5),
    aes(label = str_wrap(Description, 25)),
    size = 2.8, max.overlaps = 20, segment.size = 0.3
  ) +
  scale_shape_manual(values = c("GO Slim" = 16, "Hallmark" = 17)) +
  scale_color_manual(values = c("Both" = "green4", "PRE vs CTL only" = "blue4", "POST vs PRE only" = "red4")) +
  labs(
    title    = "Muscle GSEA: PRE vs CTL vs. POST vs PRE",
    subtitle = "GO Slim & Hallmark pathways, significant (FDR < 0.05) in at least one comparison",
    x        = "NES (PRE vs CTL)",
    y        = "NES (POST vs PRE)",
    shape    = "Database",
    color    = "Significance"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

print(p)

# --- 4: Save ------------------------------
ggsave("muscle_GSEA_PREvCTL_vs_POSTvPRE.png", p, width = 9, height = 7, dpi = 300)
