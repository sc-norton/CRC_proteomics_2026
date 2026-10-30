# Protein reversal plot

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
# Load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2,
               ggrepel, ggvenn, patchwork, ggtext, cowplot, ggh4x, ggalluvial,
               scales, gtsummary, limma, Hmisc, nlme, ggpubr, tidytext, ggsankey,
               enrichplot)

# Configuration
config <- list(
  # Data files
  limma_muscle = "03_analysis/b_reports/muscle/all_data/limma",
  # Directories
  fig_dir = "04_figures/limma",
  data_dir = "04_figures/data"
)

dir.create(config$fig_dir, showWarnings = F, recursive = T)
dir.create(config$data_dir, showWarnings = F, recursive = T)

# Color palette + theme
source("04_figures/a_scripts/00_theme.R")

# Labels
# --- 0e: Comparison Labels ---------------------------------------
comp_labels <- c(
  "Baseline_SURVvCTL" = "SURV_T1 v CTL",
  "Training_SURVvCTL" = "SURV_T2 v CTL",
  "Training_SURV"     = "SURV_T2 v SURV_T1",
  "Baseline_SURVvCTL_str" = "SURV_T1 v CTL",
  "Training_SURVvCTL_str" = "SURV_T2 v CTL",
  "Training_SURV_str"     = "SURV_T2 v SURV_T1",
  "Baseline_CREvPLA"  = "CRE_T1 v PLA_T1",
  "Training_CREvPLA"  = "CRE_T2 v PLA_T2",
  "Training_CRE"      = "CRE_T2 v CRE_T1",
  "Training_PLA"      = "PLA_T2 v PLA_T1",
  "Interaction_supp"  = "Interaction\nCRE v PLA",
  "UP_both" = "UP in CRC at PRE & POST",
  "DOWN_both" = "DOWN in CRC at PRE & POST",
  "UP_PRE_only" = "UP in CRC PRE only",
  "UP_POST_only" = "UP in CRC POST only",
  "DOWN_PRE_only" = "DOWN in CRC PRE only",
  "DOWN_POST_only" = "DOWN in CRC POST only"
)

surv_labels <- c(
  "Baseline_SURVvCTL" = "SURV_T1 v CTL",
  "Training_SURVvCTL" = "SURV_T2 v CTL",
  "Training_SURV"     = "SURV_T2 v SURV_T1"
)

# --- 0f: Display Label Remapping ---------------------------------
remap_label <- function(x) {
  x %>%
    str_replace_all("SURV", "CRC") %>%
    str_replace_all("_T1",  "_PRE") %>%
    str_replace_all("_T2",  "_POST") %>%
    str_replace_all("T1",   "PRE") %>%
    str_replace_all("T2",   "POST")
}

comp_labels  <- setNames(remap_label(comp_labels),  names(comp_labels))
surv_labels  <- setNames(remap_label(surv_labels),  names(surv_labels))

# --- 1: Load Data ----------------------------
limma_results <- list(
  Baseline_SURVvCTL = read_excel(file.path(config$limma_muscle,
                                           "limma_results_all_muscle.xlsx"), 
                                 sheet = 11),
  Training_SURVvCTL = read_excel(file.path(config$limma_muscle,
                                           "limma_results_all_muscle.xlsx"),
                                 sheet = 12),
  Training_SURV     = read_excel(file.path(config$limma_muscle,
                                           "limma_results_all_muscle.xlsx"),
                                 sheet = 13)
)

# --- 2: Restructure for plotting -----------------
df_PRE <- limma_results$Baseline_SURVvCTL
df_RET <- limma_results$Training_SURV

df_PRE <- df_PRE %>%
  dplyr::select(uniprot_id, gene, logFC_CRC = logFC, adj.Pi.Val_CRC = adj.Pi.Val, adj.P.Val_CRC = adj.P.Val)

df_RET <- df_RET %>%
  dplyr::select(uniprot_id, gene, logFC_RET = logFC, adj.Pi.Val_RET = adj.Pi.Val, adj.P.Val_RET = adj.P.Val)

df_rev <- inner_join(df_PRE, df_RET, by = c("uniprot_id", "gene"))

# --- 3: Reversal Plot --------------------------
# --- 2. Significance & quadrant categories ---
sig_thresh <- 0.05

df_rev <- df_rev %>%
  mutate(
    sig = case_when(
      adj.Pi.Val_CRC < sig_thresh & adj.Pi.Val_RET < sig_thresh ~ "Both",
      adj.Pi.Val_CRC < sig_thresh                                    ~ "CRC only",
      adj.Pi.Val_RET < sig_thresh                                  ~ "RET only",
      TRUE                                                              ~ "Neither"
    ),
    quadrant = case_when(
      logFC_CRC > 0 & logFC_RET < 0 ~ "Reversed (↑CRC, ↓RET)",
      logFC_CRC < 0 & logFC_RET > 0 ~ "Reversed (↓CRC, ↑RET)",
      logFC_CRC > 0 & logFC_RET > 0 ~ "Exacerbated (↑Both)",
      logFC_CRC < 0 & logFC_RET < 0 ~ "Exacerbated (↓Both)",
      TRUE                                   ~ "On axis"
    ),
    # Combined factor for color — sig drives shade, quadrant drives hue
    color_group = paste(sig, quadrant, sep = " | ")
  )

reversal_genes <- df_rev %>%
  dplyr::filter(quadrant == "Reversed (↑CRC, ↓RET)" | quadrant == "Reversed (↓CRC, ↑RET)") 

ex_up <- df_rev %>% 
  dplyr::filter(quadrant == "Exacerbated (↑Both)")
ex_dn <- df_rev %>% 
  dplyr::filter(quadrant == "Exacerbated (↓Both)")
rev_dt <- df_rev %>% 
  dplyr::filter(quadrant == "Reversed (↑CRC, ↓RET)")
rev_ut <- df_rev %>% 
  dplyr::filter(quadrant == "Reversed (↓CRC, ↑RET)")

label_df <- df_rev %>%
  filter(sig != "Neither",
         abs(logFC_CRC) > 1 | abs(logFC_RET) > 1) %>%
  
  # avoid -Inf
  mutate(
    adj.Pi.Val_RET = pmax(adj.Pi.Val_RET, 1e-300),
    adj.Pi.Val_CRC = pmax(adj.Pi.Val_CRC, 1e-300)
  ) %>%
  
  # direction
  mutate(
    direction = ifelse(
      abs(logFC_CRC) >= abs(logFC_RET),
      ifelse(logFC_CRC > 0, "Up", "Down"),
      ifelse(logFC_RET > 0, "Up", "Down")
    )
  ) %>%
  
  # score
  mutate(
    score = pmax(-log10(adj.Pi.Val_RET),
                 -log10(adj.Pi.Val_CRC))
  ) %>%
  
  group_by(sig, direction) %>%
  slice_max(score, n = 5, with_ties = FALSE) %>%
  ungroup()

# --- 3. Color palette ---
# Hue = direction (reversal = blue/teal, exacerbation = red/orange, neither = grey)
# Saturation/shade = significance
pal <- c(
  # Reversed
  "Both"      = colors$green,
  "CRC only"      = colors$blue,
  "RET only" = colors$red,
  "Neither" = colors$gray
)

# --- 4. Plot ---
lim <- max(abs(c(df_rev$logFC_CRC, df_rev$logFC_RET)), na.rm = TRUE)

# Compute quadrant counts
quadrant_counts <- df_rev |>
  mutate(
    quadrant = case_when(
      logFC_CRC > 0 & logFC_RET > 0 ~ "PP",
      logFC_CRC < 0 & logFC_RET < 0 ~ "NN",
      logFC_CRC > 0 & logFC_RET < 0 ~ "PN",
      logFC_CRC < 0 & logFC_RET > 0 ~ "NP"
    )
  ) |>
  group_by(quadrant) |>
  summarise(
    n_total = n(),
    n_sig   = sum(sig != "Neither", na.rm = TRUE),
    .groups = "drop"
  )

# Build label data frame with corner positions
quad_labels <- tibble(
  quadrant = c("PP", "NN", "PN", "NP"),
  title    = c(
    "Exacerbated (Up)",
    "Exacerbated (Down)",
    "Reversed (Down RET)",
    "Reversed (Up RET)"
  ),
  x = c( lim,  -lim,  lim, -lim),
  y = c( lim,  -lim, -lim,  lim),
  hjust = c(1, 0, 1, 0),
  vjust = c(1, 0, 0, 1)
) |>
  left_join(quadrant_counts, by = "quadrant") |>
  mutate(
    label = sprintf("%s\nn = %d | sig = %d", title, n_total, n_sig)
  ) |>
  mutate(
    bg_color = case_when(
      quadrant %in% c("PP", "NN") ~ colors$blue_dark,
      quadrant %in% c("PN", "NP") ~ colors$red_dark
    )
  )

p1 <- ggplot(df_rev, aes(x = logFC_CRC, y = logFC_RET,
                    color = sig, fill = sig)) +
  # Quadrant shading — add before geom_point so points render on top
  annotate("rect", xmin = 0,  xmax = Inf, ymin = 0,  ymax = Inf,  fill = colors$blue, alpha = 0.2) +  # both positive
  annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,  fill = colors$blue, alpha = 0.2) +  # both negative
  annotate("rect", xmin = 0,  xmax = Inf, ymin = -Inf, ymax = 0,  fill = colors$red,   alpha = 0.2) +  # discordant
  annotate("rect", xmin = -Inf, xmax = 0, ymin = 0,  ymax = Inf,  fill = colors$red,   alpha = 0.2) +  # discordant
  # Reference lines
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  # Points — neither sig plotted first (bottom layer)
  geom_point(data = filter(df_rev, sig == "Neither"),
             size = 1.2, alpha = 0.8, shape = 16) +
  geom_point(data = filter(df_rev, sig != "Neither"),
             size = 2.2, alpha = 0.85, shape = 16) +
  scale_color_manual(values = pal, name = NULL) +
  scale_fill_manual(values  = pal, guide = "none") +
  labs(
    x = NULL,
    y = NULL,
    title = "CRC vs. RET Direction"
  ) +
  annotate("text",
           x = lim * -1, y = 0,
           label = "log[2]*FC~'(CRC_PRE vs CTL)'",
           hjust = 0, vjust = -0.3,
           size = 3, fontface = "bold", color = "grey30",
           parse = TRUE) +
  annotate("text",
           x = 0, y = lim * -1,
           label = "log[2]*FC~'(CRC_POST vs CRC_PRE)'",
           hjust = 0, vjust = -0.3,
           angle = 90,
           size = 3, fontface = "bold", color = "grey30",
           parse = TRUE) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) + 
  theme_cr()

p1

p2 <- p1 +
  geom_label(
    data        = quad_labels,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes = FALSE,
    fill        = quad_labels$bg_color,
    color       = colors$white,
    alpha       = 0.7,
    size        = 3,
    fontface    = "bold",
    lineheight  = 1.2,
    label.size  = NA,        # removes the border
    label.padding = unit(0.25, "lines")
  )
p2

p3 <- p2 +
  geom_label_repel(
    data = label_df,
    aes(label = gene, fill = sig),
    color = colors$white,
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.3,
    label.padding = unit(0.15, "lines"),
    label.size = NA,
    segment.color = "grey60",
    segment.size = 0.3,
    show.legend = FALSE
  )
p3

# --- SAVE PLOTS ----
# PDF
ggsave(file.path(config$fig_dir, "protein_reversal_plot.pdf"), p2, device = "pdf")

# PNG
ggsave(file.path(config$fig_dir, "protein_reversal_plot.png"), p2, device = "png")
