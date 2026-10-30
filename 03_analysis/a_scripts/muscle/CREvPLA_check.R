# =============================================================================
# CRE vs PLA Concordance Analysis
# =============================================================================

config <- list(
  gsea_muscle  = "03_analysis/b_reports/muscle/reformat/SURVvCTL/gsea",
  gsea_muscle_CRE = "03_analysis/b_reports/muscle/reformat/CREvPLA/gsea",
  limma_muscle = "03_analysis/b_reports/muscle/reformat/all_data/limma",
  fig_dir      = "04_figures/fig_05/new",
  data_dir     = "04_figures/data/fig_05/new"
)
# --- Load CRE and PLA training limma results ---------------------------------
# Adjust sheet numbers to match your actual excel file
limma_results_CRE <- list(
  Baseline_CREvPLA   = read_excel(file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 1),
  Training_CREvPLA   = read_excel(file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 2),
  Training_CRE       = read_excel(file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 3),
  Training_PLA       = read_excel(file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 4),
  Interaction_supp   = read_excel(file.path(config$limma_muscle, "limma_results_all_muscle.xlsx"), sheet = 5)
)

# --- Build concordance data frame --------------------------------------------
df_cre <- limma_results_CRE$Training_CRE |>
  dplyr::select(uniprot_id, gene,
                logFC_CRE     = logFC,
                Pi.Val_CRE    = Pi.Val,
                adj.P.Val_CRE = adj.P.Val)

df_pla <- limma_results_CRE$Training_PLA |>
  dplyr::select(uniprot_id, gene,
                logFC_PLA     = logFC,
                Pi.Val_PLA    = Pi.Val,
                adj.P.Val_PLA = adj.P.Val)

df_conc <- inner_join(df_cre, df_pla, by = c("uniprot_id", "gene")) |>
  # Also join the collapsed SURV training result for three-way comparison
  inner_join(
    limma_results$Training_SURV |>
      dplyr::select(uniprot_id, gene,
                    logFC_SURV     = logFC,
                    Pi.Val_SURV    = Pi.Val,
                    adj.P.Val_SURV = adj.P.Val),
    by = c("uniprot_id", "gene")
  ) |>
  mutate(
    # Direction agreement between CRE and PLA
    same_direction = sign(logFC_CRE) == sign(logFC_PLA),
    # Sig classification
    sig_CRE  = Pi.Val_CRE  < 0.05,
    sig_PLA  = Pi.Val_PLA  < 0.05,
    sig_SURV = Pi.Val_SURV < 0.05,
    sig_group = case_when(
      sig_CRE & sig_PLA  ~ "Both",
      sig_CRE & !sig_PLA ~ "CRE only",
      !sig_CRE & sig_PLA ~ "PLA only",
      TRUE               ~ "Neither"
    )
  )

# --- Summary stats -----------------------------------------------------------
# Overall direction concordance
conc_summary <- df_conc |>
  summarise(
    n_proteins        = n(),
    n_same_direction  = sum(same_direction, na.rm = TRUE),
    pct_same          = round(mean(same_direction, na.rm = TRUE) * 100, 1),
    r_cre_pla         = round(cor(logFC_CRE,  logFC_PLA,  use = "complete.obs"), 3),
    r_cre_surv        = round(cor(logFC_CRE,  logFC_SURV, use = "complete.obs"), 3),
    r_pla_surv        = round(cor(logFC_PLA,  logFC_SURV, use = "complete.obs"), 3)
  )

print(conc_summary)

# Direction concordance broken down by significance
conc_by_sig <- df_conc |>
  group_by(sig_group) |>
  summarise(
    n                = n(),
    pct_same_dir     = round(mean(same_direction, na.rm = TRUE) * 100, 1),
    mean_logFC_CRE   = round(mean(logFC_CRE,  na.rm = TRUE), 3),
    mean_logFC_PLA   = round(mean(logFC_PLA,  na.rm = TRUE), 3),
    mean_logFC_SURV  = round(mean(logFC_SURV, na.rm = TRUE), 3),
    r_cre_pla        = round(cor(logFC_CRE, logFC_PLA, use = "complete.obs"), 3),
    .groups = "drop"
  )

print(conc_by_sig)

# =============================================================================
# Plots
# =============================================================================

lim_conc <- max(abs(c(df_conc$logFC_CRE, df_conc$logFC_PLA)), na.rm = TRUE) * 1.05

# --- Plot A: CRE vs PLA logFC ------------------------------------------------
p_cre_pla <- ggplot(df_conc, aes(x = logFC_CRE, y = logFC_PLA, color = sig_group)) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin =    0, ymax =  Inf, fill = colors$blue, alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin = -Inf, ymax =    0, fill = colors$blue, alpha = 0.15) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin = -Inf, ymax =    0, fill = colors$red,  alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin =    0, ymax =  Inf, fill = colors$red,  alpha = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = colors$gray_dark,
              linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE,
              aes(x = logFC_CRE, y = logFC_PLA)) +
  geom_point(data = filter(df_conc, sig_group == "Neither"),
             size = 1.2, alpha = 0.5, shape = 16) +
  geom_point(data = filter(df_conc, sig_group != "Neither"),
             size = 2.2, alpha = 0.85, shape = 16) +
  geom_label_repel(
    data          = filter(df_conc, sig_group != "Neither"),
    aes(label     = gene, fill = sig_group),
    color         = colors$white,
    size          = 2.5,
    max.overlaps  = 15,
    box.padding   = 0.3,
    label.padding = unit(0.15, "lines"),
    label.size    = NA,
    segment.color = "grey60",
    segment.size  = 0.3,
    show.legend   = FALSE
  ) +
  annotate("text",
           x = -lim_conc * 0.95, y = lim_conc * 0.95,
           label = sprintf("r = %.3f", conc_summary$r_cre_pla),
           hjust = 0, vjust = 1, size = 4, fontface = "bold", color = colors$gray_dark) +
  scale_color_manual(values = c("Both"     = colors$green,
                                "CRE only" = colors$blue,
                                "PLA only" = colors$orange,
                                "Neither"  = colors$gray),
                     name = NULL) +
  scale_fill_manual(values  = c("Both"     = colors$green,
                                "CRE only" = colors$blue,
                                "PLA only" = colors$orange,
                                "Neither"  = colors$gray),
                    guide = "none") +
  coord_cartesian(xlim = c(-lim_conc, lim_conc), ylim = c(-lim_conc, lim_conc)) +
  labs(x = "logFC (CRE: POST vs PRE)",
       y = "logFC (PLA: POST vs PRE)",
       title = "CRE vs PLA Training Response") +
  theme_cr()

# --- Plot B: CRE vs collapsed SURV -------------------------------------------
lim_surv <- max(abs(c(df_conc$logFC_CRE, df_conc$logFC_SURV)), na.rm = TRUE) * 1.05

p_cre_surv <- ggplot(df_conc, aes(x = logFC_SURV, y = logFC_CRE, color = sig_group)) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin =    0, ymax =  Inf, fill = colors$blue, alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin = -Inf, ymax =    0, fill = colors$blue, alpha = 0.15) +
  annotate("rect", xmin =    0, xmax =  Inf, ymin = -Inf, ymax =    0, fill = colors$red,  alpha = 0.15) +
  annotate("rect", xmin = -Inf, xmax =    0, ymin =    0, ymax =  Inf, fill = colors$red,  alpha = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey50", linetype = "dashed") +
  geom_smooth(method = "lm", se = TRUE, color = colors$gray_dark,
              linewidth = 0.7, linetype = "dashed", inherit.aes = FALSE,
              aes(x = logFC_SURV, y = logFC_CRE)) +
  geom_point(data = filter(df_conc, sig_group == "Neither"),
             size = 1.2, alpha = 0.5, shape = 16) +
  geom_point(data = filter(df_conc, sig_group != "Neither"),
             size = 2.2, alpha = 0.85, shape = 16) +
  annotate("text",
           x = -lim_surv * 0.95, y = lim_surv * 0.95,
           label = sprintf("r = %.3f", conc_summary$r_cre_surv),
           hjust = 0, vjust = 1, size = 4, fontface = "bold", color = colors$gray_dark) +
  scale_color_manual(values = c("Both"     = colors$green,
                                "CRE only" = colors$blue,
                                "PLA only" = colors$orange,
                                "Neither"  = colors$gray),
                     name = NULL) +
  scale_fill_manual(values  = c("Both"     = colors$green,
                                "CRE only" = colors$blue,
                                "PLA only" = colors$orange,
                                "Neither"  = colors$gray),
                    guide = "none") +
  coord_cartesian(xlim = c(-lim_surv, lim_surv), ylim = c(-lim_surv, lim_surv)) +
  labs(x = "logFC (Collapsed SURV: POST vs PRE)",
       y = "logFC (CRE: POST vs PRE)",
       title = "CRE vs Collapsed SURV") +
  theme_cr()

# --- Plot C: PLA vs collapsed SURV -------------------------------------------
p_pla_surv <- p_cre_surv +
  aes(y = logFC_PLA) +
  annotate("text",
           x = -lim_surv * 0.95, y = lim_surv * 0.95,
           label = sprintf("r = %.3f", conc_summary$r_pla_surv),
           hjust = 0, vjust = 1, size = 4, fontface = "bold", color = colors$gray_dark) +
  labs(y = "logFC (PLA: POST vs PRE)",
       title = "PLA vs Collapsed SURV")

# --- Composite ---------------------------------------------------------------
fig_concordance <- plot_grid(
  p_cre_pla, p_cre_surv, p_pla_surv,
  labels = c("A", "B", "C"),
  ncol   = 3,
  align  = "h",
  axis   = "tb"
)

fig_concordance

save_plot(fig_concordance, "fig_supp_cre_pla_concordance", width = 18, height = 6)

# What fraction are in concordant vs discordant quadrants?
df_conc |>
  mutate(quadrant = case_when(
    logFC_CRE > 0 & logFC_PLA > 0 ~ "Both Up",
    logFC_CRE < 0 & logFC_PLA < 0 ~ "Both Down",
    logFC_CRE > 0 & logFC_PLA < 0 ~ "Discordant (CRE up)",
    logFC_CRE < 0 & logFC_PLA > 0 ~ "Discordant (PLA up)"
  )) |>
  count(quadrant) |>
  mutate(pct = round(n / sum(n) * 100, 1))

# Check if your key Training_SURV hits are concordant between CRE and PLA
key_hits <- df_conc |>
  filter(sig_SURV) |>
  mutate(
    concordant    = same_direction,
    logFC_diff    = abs(logFC_CRE - logFC_PLA),
    avg_logFC     = (logFC_CRE + logFC_PLA) / 2
  ) |>
  arrange(desc(logFC_diff)) |>
  dplyr::select(gene, logFC_CRE, logFC_PLA, logFC_SURV, concordant, logFC_diff)

print(key_hits, n = 30)

# Concordance rate among sig_SURV hits
key_hits |>
  summarise(
    n_sig          = n(),
    n_concordant   = sum(concordant),
    pct_concordant = round(mean(concordant) * 100, 1),
    # Is the collapsed logFC between the two subgroup estimates?
    n_surv_between = sum(
      (logFC_SURV >= pmin(logFC_CRE, logFC_PLA)) &
        (logFC_SURV <= pmax(logFC_CRE, logFC_PLA))
    ),
    # Mean absolute deviation of SURV from the CRE/PLA average
    mean_abs_dev   = round(mean(abs(logFC_SURV - (logFC_CRE + logFC_PLA) / 2)), 4),
    r_cre_pla_sig  = round(cor(logFC_CRE, logFC_PLA), 3)
  )

limma_results$Training_SURV |>
  filter(Pi.Val < 0.05) |>
  dplyr::select(gene, logFC, P.Value, adj.P.Val, Pi.Val, DE.pi) |>
  print(n = Inf)
