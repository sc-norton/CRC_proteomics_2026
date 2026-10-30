# ============================================================
# fig_06 — composite
# WGCNA co-expression modules and their eigengene response
# across the design (Proteins | Trajectory | Top ORA pathways)
#
# Skips the "Member response" (fry/fGSEA) column from the
# reference figure — not yet built (see conversation notes).
#
# Style helpers (scale_text, strip_for_composite, save_plot) come
# from fig_06_style.R, a self-contained file — not shared with
# fig_05, so FIG_WIDTH_MM/REFERENCE_WIDTH_MM here are independent
# of whatever fig_05 uses.
# ============================================================

pacman::p_load(tidyverse, ggplot2, patchwork, cowplot, ggtext, ggnewscale, writexl, rprojroot)

root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "04_figures", "a_scripts", "00_theme.R"))       # colors, theme_cr()
source(file.path(root, "04_figures", "a_scripts", "fig_06_style.R"))   # scale_text(), strip_for_composite(), save_plot()

fig_06_data <- readRDS(file.path(root, "04_figures", "data", "fig_06", "fig_06_data_plasma.rds"))
cfg <- fig_06_data$config
module_order <- fig_06_data$module_order  # size-descending, shared row order across all panels

group_labels <- c("CTL_T1" = "CTL", "SURV_T1" = "CRC", "SURV_T2" = "CRC+RET")
source_colors <- c("GO:BP" = colors$red_light, "GO:CC" = colors$yellow_light, "GO:MF" = colors$blue_light,
                   "Reactome" = colors$brown_light, "Hallmark" = colors$orange_light, "KEGG" = colors$purple_light)
# TODO: still hex fallbacks, not your real colors$ palette — swap in the
# actual named entries you want (e.g. colors$up_light / colors$down_light).

# ---- Panel widths (mm) for scale_text()'s panel_width_mm argument ----
# scale_text() needs the ACTUAL rendered width of the panel it's called in,
# not the total figure width — these follow the same relative widths passed
# to plot_layout() below, so keep the two in sync if you change either.
# FIG_WIDTH_MM comes from fig_06_style.R.
panel_rel_widths <- c(A = 1, B = 1.4, C = 2.2, D = 1.3)
panel_widths_mm <- FIG_WIDTH_MM * panel_rel_widths / sum(panel_rel_widths)

# Shared facet spacing constant, used by Panels B and C (both faceted).
# Panel A is a single non-faceted plot and doesn't use this.
ROW_SPACING <- unit(1, "lines")

# ---- Panel A: protein count per module ----
p_A_data <- fig_06_data$module_sizes %>%
  mutate(module = factor(module, levels = module_order),
         row_label = paste0(module, " \u00b7 ", n_proteins),
         y_dummy = 1)

p_A <- ggplot(p_A_data, aes(x = n_proteins, y = y_dummy)) +
  geom_col(aes(fill = as.character(module)), width = 0.6, orientation = "y") +
  geom_text(aes(x = 0, label = row_label), y = 1.4, hjust = 1,
            size = scale_text(7, panel_widths_mm["A"]), fontface = "bold") +
  # theme_cr()'s default panel.border draws a full box; with scale_x_reverse
  # its left edge lands at the bar TIPS, not the base. Blanked below and
  # replaced with a vline at x=0, which sits exactly at the bars' base
  # (the reversed axis's zero point) in every facet.
  geom_vline(xintercept = 0, color = "grey20", linewidth = 0.5) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0.15))) +
  scale_y_continuous(limits = c(0.5, 1.6), expand = c(0, 0)) +
  scale_fill_identity() +
  coord_cartesian(clip = "off") +
  # True facet_wrap -- same panel.spacing.y mechanism as B/C, for genuine
  # row alignment. Fixed (shared, default) scales rather than "free" so
  # every module's bar is measured against the same 0-max range, matching
  # the reference figure (smaller modules show visibly shorter bars, not
  # independently rescaled to fill their own row).
  # axes = "margins" (ggplot2 >= 3.5.0) suppresses the repeated x-axis on
  # every interior facet, leaving only the bottom-most one -- if your
  # ggplot2 version predates this, every facet will show its own x-axis
  # instead; upgrading ggplot2 is the clean fix, or tell me and I'll build
  # a manual per-facet axis-blanking workaround instead.
  facet_wrap(~ module, ncol = 1, axes = "margins") +
  labs(x = "count per module", y = NULL) +
  theme_cr() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.border = element_blank(),
        panel.spacing.y = ROW_SPACING)


p_A <- strip_for_composite(p_A) +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))

# ---- Panel B: module eigengene trajectory, one facet row per module ----
traj <- fig_06_data$trajectory_stats %>%
  mutate(
    module = factor(module, levels = module_order),
    group_label = factor(group_labels[.data[[cfg$group_col]]], levels = group_labels)
  )

# Display labels per contrast -- shown above each bracket alongside its
# p-value, since there are now 3 brackets and it's otherwise ambiguous
# which one represents which comparison. Matches the CRC/CRC+RET naming
# convention used consistently across the rest of this project's figures.
CONTRAST_BRACKET_LABELS <- c(
  Baseline_SURVvCTL = "CRC v CTL",
  Training_SURV     = "CRC+RET v CRC",
  Training_SURVvCTL = "CRC+RET v CTL"
)

# Per-module eigengene range, used to scale both the tick length and the
# two bracket-height tiers -- additive offsets relative to THIS module's
# own range (not a fixed/multiplicative value), same reasoning as the
# original single-bracket version: a fixed or multiplicative offset can
# push a bracket to the wrong side or off-panel for modules with very
# different eigengene scales.
bracket_range <- traj %>%
  group_by(module) %>%
  summarise(
    y_max = max(mean_eigengene + se_eigengene, na.rm = TRUE),
    y_min = min(mean_eigengene - se_eigengene, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(range = pmax(y_max - y_min, 1e-6))

TICK_FRAC  <- 0.04  # vertical end-tick length, fraction of module's own range
TIER1_FRAC <- 0.15  # short-span brackets (CTL-CRC, CRC-CRC+RET) height above y_max
TIER2_FRAC <- 0.42  # wide-span bracket (CTL-CRC+RET) height above y_max -- clears tier-1 labels

brackets <- fig_06_data$trajectory_brackets %>%
  left_join(bracket_range, by = "module") %>%
  filter(is.finite(range), sig) %>%  # only draw brackets for FDR-significant
  # contrasts -- a module can end up with
  # 0-3 brackets depending on how many of
  # its 3 pairwise contrasts clear q<0.05
  mutate(
    y_bracket     = y_max + if_else(tier == 1, TIER1_FRAC, TIER2_FRAC) * range,
    y_tick_bottom = y_bracket - TICK_FRAC * range,
    x_mid         = (x1 + x2) / 2,
    contrast_label = CONTRAST_BRACKET_LABELS[named_contrast],
    label         = sprintf("p=%.3f", p.value),
    fontface      = "bold"  # every bracket shown is significant by construction now
  )

p_B <- ggplot(traj, aes(x = group_label, y = mean_eigengene, group = module,
                        color = as.character(module))) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.4) +
  geom_line(linewidth = 0.9) +
  geom_errorbar(aes(ymin = mean_eigengene - se_eigengene, ymax = mean_eigengene + se_eigengene),
                width = 0.08) +
  geom_point(size = 1.6) +
  # bracket: horizontal span + a short vertical tick at each end
  geom_segment(data = brackets, aes(x = x1, xend = x2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, color = "grey20", linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = x1, xend = x1, y = y_bracket, yend = y_tick_bottom),
               inherit.aes = FALSE, color = "grey20", linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = x2, xend = x2, y = y_bracket, yend = y_tick_bottom),
               inherit.aes = FALSE, color = "grey20", linewidth = 0.4) +
  geom_text(data = brackets, aes(x = x_mid, y = y_bracket, label = label, fontface = fontface),
            inherit.aes = FALSE, vjust = -0.3, size = scale_text(3.2, panel_widths_mm["B"]),
            lineheight = 0.9) +
  scale_color_identity() +
  # extra top headroom for THREE stacked brackets (tier 1 + tier 2 + their
  # labels), well beyond what a single bracket needed
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.62))) +
  facet_wrap(~ module, ncol = 1, scales = "free_y") +
  coord_cartesian(clip = "off") +
  labs(x = NULL, y = "eigengene (a.u.)") +
  theme_cr() +
  theme(strip.text = element_blank(), panel.spacing.y = ROW_SPACING)

# strip_for_composite() resets plot.margin to a flat value -- apply any
# panel-specific margin AFTER calling it, or it gets silently overwritten
# (this is exactly what happened to Panel C's left margin before).
p_B <- strip_for_composite(p_B) +
  theme(plot.margin = margin(t = 2, r = 8, b = 2, l = 2, unit = "mm"))  # extra right
# margin so B
# doesn't crowd C

# ---- Panel C: top ORA pathways per module, one facet row per module ----
# Shorten specific pathway names here rather than wrapping every label --
# same pattern as ora_label_map in fig_05_panel_A_protein.R. Add an entry
# for the long name that's spilling out of its bar: left side is the exact
# Description string as it appears in fig_06_data$top_ora, right side is
# your shortened replacement. Tell me the actual long name and I'll fill
# this in for you instead of guessing.
pathway_label_overrides <- c(
  "Diseases Of Signal Transduction By Growth Factor Receptors And Second Messengers" 
  = "Growth Factor Receptor Signaling Diseases"
)

apply_pathway_labels <- function(description, overrides) {
  if_else(description %in% names(overrides), overrides[description], description)
}

ora_plot_data <- fig_06_data$top_ora %>%
  mutate(
    module = factor(module, levels = module_order),
    Description = apply_pathway_labels(Description, pathway_label_overrides)
  )

dark_modules <- c("black", "blue", "brown", "red")  # only these need white text for contrast
ora_plot_data <- ora_plot_data %>%
  mutate(text_color = if_else(as.character(module) %in% dark_modules, "white", "black"))

p_C <- ggplot(ora_plot_data, aes(x = neg_log_fdr, y = row_id)) +
  geom_col(aes(fill = as.character(module), alpha = bold_sig), width = 0.7) +
  # pathway name printed on the bar itself, right-justified inside it;
  # text color is fixed per module (dark modules get white, rest get black),
  # not tied to significance — the alpha shading already carries that signal
  geom_text(aes(x = neg_log_fdr / 2, label = Description, color = text_color), hjust = 0.5,
            size = scale_text(3, panel_widths_mm["C"]), fontface = "bold") +
  scale_fill_identity() +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  scale_color_identity() +
  ggnewscale::new_scale_color() +
  # database source label, left of the bar's start
  geom_text(aes(label = source, fontface = if_else(bold_sig, "bold", "plain"), color = source),
            x = 0, hjust = 1.05, size = scale_text(3.5, panel_widths_mm["C"])) +
  scale_color_manual(values = source_colors, guide = "none") +
  scale_y_discrete(expand = expansion(add = 0.5)) +  # a little breathing room
  # around wrapped 2-line labels
  facet_wrap(~ module, ncol = 1, scales = "free") +
  coord_cartesian(clip = "off") +
  labs(x = expression(-log[10]~FDR), y = NULL) +
  theme_cr() +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank(),
        panel.spacing.y = ROW_SPACING)

# Same margin-after-strip fix as Panel B -- this restores the left margin
# Panel C needs for the source labels (which sit left of x=0 via hjust=1.05),
# which strip_for_composite() would otherwise silently wipe out.
p_C <- strip_for_composite(p_C) +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 10, unit = "mm"))

# ---- Panel D: hub proteins per module, one facet row per module ----
# Bar length = kME with the protein's OWN module (0-1 scale); gene symbol
# printed on the bar, same on-bar/dark-module-text-contrast pattern as
# Panel C. dark_modules is already defined above (Panel C section) and
# reused here unchanged.
hub_plot_data <- fig_06_data$top_hubs %>%
  mutate(
    module = factor(module, levels = module_order),
    text_color = if_else(as.character(module) %in% dark_modules, "white", "black")
  )

p_D <- ggplot(hub_plot_data, aes(x = kME_own_module, y = row_id)) +
  geom_col(aes(fill = as.character(module)), width = 0.7) +
  geom_text(aes(x = kME_own_module / 2, label = display_name, color = text_color),
            hjust = 0.5, fontface = "italic",
            size = scale_text(3.5, panel_widths_mm["D"])) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  facet_wrap(~ module, ncol = 1, scales = "free") +
  coord_cartesian(clip = "off") +
  labs(x = "kME (own module)", y = NULL) +
  theme_cr() +
  theme(strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.spacing.y = ROW_SPACING)

p_D <- strip_for_composite(p_D) +
  theme(plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm"))

# ---- Stitch ----
composite <- p_A + p_B + p_C + p_D +
  plot_layout(widths = unname(panel_rel_widths))

title_row <- cowplot::ggdraw() +
  cowplot::draw_label("Plasma WGCNA co-expression modules and their eigengene response across the design",
                      fontface = "bold", size = scale_text(13, FIG_WIDTH_MM), hjust = 0, x = 0.01)

header_size <- scale_text(11, FIG_WIDTH_MM)  # same width arg for all -> same resulting size,
# unlike per-panel width which deliberately varies it
col_headers <- cowplot::plot_grid(
  cowplot::ggdraw() + cowplot::draw_label("A  Proteins", fontface = "bold", size = header_size, x = 0.05, hjust = 0),
  cowplot::ggdraw() + cowplot::draw_label("B  Module eigengene trajectory", fontface = "bold", size = header_size, x = 0.05, hjust = 0),
  cowplot::ggdraw() + cowplot::draw_label("C  Top ORA pathways", fontface = "bold", size = header_size, x = 0.05, hjust = 0),
  cowplot::ggdraw() + cowplot::draw_label("D  Hub proteins", fontface = "bold", size = header_size, x = 0.05, hjust = 0),
  nrow = 1, rel_widths = unname(panel_rel_widths)
)

final_fig <- cowplot::plot_grid(
  title_row, col_headers, composite,
  ncol = 1, rel_heights = c(0.05, 0.04, 1)
)

n_modules <- length(module_order)
save_plot(
  final_fig,
  file.path(cfg$out_root, paste0("fig_06_wgcna_modules_", cfg$tissue, ".pdf")),
  width = 380, height = max(120, n_modules * 45)  # mm; scale height to module count
)

write_xlsx(
  list(module_sizes = fig_06_data$module_sizes,
       trajectory_stats = traj,
       trajectory_brackets = brackets,
       top_ora = fig_06_data$top_ora,
       top_hubs = hub_plot_data),
  file.path(cfg$out_root, paste0("fig_06_wgcna_modules_", cfg$tissue, ".xlsx"))
)