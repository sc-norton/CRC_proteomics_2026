# =============================================================================
# fig_03_style.R
# Shared sizing/theme/helpers for the Figure 03 composite (volcano | Muscle
# union heatmap | Plasma union heatmap). Source this ONCE, before building
# any panel, from wherever fig_03's composite is assembled.
#
# Mirrors fig_05_style.R's sizing philosophy: every text element has a BASE
# size defined at a reference width (PANEL_MD), and scale_text() scales that
# base size by sqrt(actual_width / PANEL_MD) when a panel is embedded at
# some other width. Kept as a separate, self-contained script per panel
# (matching fig_05_style.R's own convention) rather than sharing one style
# file across figures, so a change tuned for one figure can't silently
# affect another.
# =============================================================================

library(ggplot2)
library(scales)
library(grid)

source("04_figures/a_scripts/00_theme.R")   # colors, theme_cr()

# --- Reference width + scaling -----------------------------------------------

PANEL_MD <- 178   # reference width (mm) all BASE_* constants below are tuned for -- matches fig_05_style.R's PANEL_MD so the two figures scale on the same baseline

scale_text <- function(base_size, panel_width_mm, ref_width = PANEL_MD) {
  base_size * sqrt(panel_width_mm / ref_width)
}

# --- Base text sizes at PANEL_MD ---------------------------------------------
# PT sizes -- used in theme(element_text(size = ...))
FIG_TITLE_SIZE    <- 10
FIG_SUBTITLE_SIZE <- 7
FIG_AXIS_TEXT     <- 7
FIG_LEGEND_TITLE  <- 8
FIG_LEGEND_TEXT   <- 7
FIG_STRIP_SIZE    <- 8
BASE_TAG          <- 11   # composite-level panel tags (A/B/C), used by the stitcher

# MM sizes -- used in geom_text()/annotate()/geom_label() size = ... (ggplot
# text geoms are in mm, roughly pt / 2.845)
BASE_GENE_LABEL <- 2.3   # volcano repelled gene-symbol labels (was hardcoded size = 3)
BASE_COUNT      <- 2.6   # volcano UP/DOWN facet counts (was hardcoded size = 3)
BASE_SIG_MARK   <- 2.6   # heatmap significance asterisk (was hardcoded size = 3.5)

# PT sizes for legend key/colorbar geometry -- used in guide_legend()/
# guide_colorbar() key/bar dimensions, scaled the same way as text.
BASE_LEGEND_KEY    <- 7    # discrete legend key width/height
BASE_COLORBAR_LEN  <- 45   # colorbar long-axis length
BASE_COLORBAR_WID  <- 6    # colorbar short-axis width

# --- Scaled theme, built on theme_cr() ---------------------------------------
fig_theme_scaled <- function(panel_width_mm, title_hjust = 0) {
  theme_cr(base_size = round(scale_text(FIG_AXIS_TEXT, panel_width_mm), 1),
           title_hjust = title_hjust) +
    theme(
      plot.title    = element_text(face = "bold", hjust = title_hjust,
                                   size = round(scale_text(FIG_TITLE_SIZE, panel_width_mm), 1)),
      plot.subtitle = element_text(face = "bold.italic", hjust = title_hjust,
                                   color = "grey30",
                                   size = round(scale_text(FIG_SUBTITLE_SIZE, panel_width_mm), 1)),
      legend.title  = element_text(face = "bold",
                                   size = round(scale_text(FIG_LEGEND_TITLE, panel_width_mm), 1)),
      legend.text   = element_text(size = round(scale_text(FIG_LEGEND_TEXT, panel_width_mm), 1)),
      strip.text    = element_text(face = "bold",
                                   size = round(scale_text(FIG_STRIP_SIZE, panel_width_mm), 1))
    )
}

# Strip title/subtitle/legend before embedding a standalone-authored plot
# into the composite -- the composite supplies its own tag+title band via
# draw_label(), and a shared legend via build_heatmap_legend() below, so
# panels being merged shouldn't carry their own copies of either.
#
# Uses & (not +), same reasoning as fig_05_style.R: & broadcasts into every
# nested plot of a composed patchwork object, not just the outer object's
# own theme/labs -- required here since a plot built with ggnewscale still
# behaves like a single ggplot for this purpose, but any future panel that
# IS a nested patchwork needs the broadcast to reach its own children too.
strip_for_composite <- function(p) {
  p & labs(title = NULL, subtitle = NULL) &
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm"))
}

# --- Shared heatmap legend, extracted from a real panel ----------------------
# Unlike fig_05_style.R's build_shared_legend() (which builds a SYNTHETIC
# plot from a named palette -- fine when every merged panel shares one
# simple categorical scale), our Muscle/Plasma legend combines a discrete
# manual scale (membership) with a continuous gradient scale (logFC) via
# ggnewscale::new_scale_fill(). That combination doesn't reduce to a single
# palette, so instead of synthesizing a reference plot, this pulls the
# legend straight off an already-built panel (before stripping it for
# embedding) via cowplot::get_legend(). Muscle and Plasma share an
# identical scale (same shared_range, same membership_colors), so either
# panel is a valid source.
build_heatmap_legend <- function(reference_plot, panel_width_mm = PANEL_MD) {
  p <- reference_plot +
    theme(legend.position = "bottom", legend.box = "horizontal", legend.direction = "vertical",
          legend.text  = element_text(size = round(scale_text(FIG_LEGEND_TEXT, panel_width_mm), 1)),
          legend.title = element_text(size = round(scale_text(FIG_LEGEND_TITLE, panel_width_mm), 1)))
  cowplot::get_legend(p)
}

# --- Save helper -- mm-native, cairo_pdf for Unicode glyphs -----------------
save_scaled <- function(p, dir, name, width_mm, height_mm) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".pdf")), p,
         width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  ggsave(file.path(dir, paste0(name, ".png")), p,
         width = width_mm, height = height_mm, units = "mm", dpi = 300)
}