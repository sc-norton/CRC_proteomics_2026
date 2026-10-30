# =============================================================================
# fig_05_style.R
# Shared sizing/theme/helpers for the Figure 05 composite (protein reversal |
# pathway reversal | fry | heatmap-sankey). Source this ONCE, before any
# panel script, from fig_05_stitch.R.
#
# Sizing philosophy: every text element has a BASE size defined at a
# reference width (PANEL_MD, a standard double-column-journal width in mm).
# When a panel is drawn at some other width -- because it's 2/3 of a 380mm
# composite row instead of a standalone 178mm figure -- scale_text() scales
# that base size by sqrt(actual_width / PANEL_MD), so panels embedded wider
# get modestly larger text and panels embedded narrower get modestly
# smaller text, instead of every panel guessing its own sizes independently.
# =============================================================================

library(ggplot2)
library(scales)
library(grid)

source("04_figures/a_scripts/00_theme.R")   # colors, theme_cr()

# --- Reference width + scaling -----------------------------------------------

PANEL_MD <- 178   # reference width (mm) all BASE_* constants below are tuned for

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
BASE_TAG          <- 11   # composite-level panel tags (A/B/C/D), used by the stitcher

# MM sizes -- used in geom_text()/annotate()/geom_label() size = ... (ggplot
# text geoms are in mm, roughly pt / 2.845)
BASE_QUADRANT <- 2.6   # quadrant corner labels (protein + pathway scatters)
BASE_GENE     <- 2.3   # repelled gene/pathway point labels
BASE_PATHWAY  <- 3.0   # ORA bar pathway-name labels
BASE_COUNT    <- 2.6   # ORA sig stars, UP/DOWN facet counts
BASE_STAT     <- 2.6   # in-panel stat annotations (rho, axis callouts)

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
# draw_label(), so panels shouldn't carry their own.
#
# Uses & (not +) so this reaches into every nested plot of a COMPOSED
# patchwork object (e.g. fig_compound = left_col | p_scatter | right_col,
# or the fry 5-row stack), not just the outer composite's own chrome. Using
# + here was the original bug that let Panel A's own bottom legend leak
# into the composite despite this function supposedly removing it -- + only
# touches the top-level patchwork object's own theme/labs, & broadcasts
# into every panel nested inside it.
strip_for_composite <- function(p) {
  p & labs(title = NULL, subtitle = NULL) &
    theme(legend.position = "none", plot.margin = margin(2, 2, 2, 2, "mm"))
}

# --- Shared reversal palette (protein scatter + pathway scatter + legend) ---
# Names here MUST exactly match the category strings in Panel A's `sig`
# column and Panel B's `sig_in` column (both use "CRC v CTL only" /
# "CRC+RET v CRC only", consistent with the Panel D sankey script's own
# CRC/CRC+RET terminology) -- scale_color_manual() matches by exact name,
# so any mismatch here silently drops colors for every point in both
# scatter plots and mislabels the shared legend built from these names.
pal_shared <- c(
  "Both"                = "#4E8B6F",
  "CRC v CTL only"       = "#2471A3",
  "CRC+RET v CRC only"   = "#C0392B",
  "Neither"              = "#9E9E9E"
)

# --- Quadrant-label helpers (protein + pathway reversal scatters) -----------
geom_quadrant_label <- function(df, size = BASE_QUADRANT) {
  geom_label(
    data          = df,
    aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
    inherit.aes   = FALSE,
    fill          = df$bg_color,
    color         = colors$white,
    alpha         = 0.7,
    size          = size,
    fontface      = "bold",
    lineheight    = 1.2,
    label.size    = NA,
    label.padding = unit(0.25, "lines"),
    show.legend   = FALSE
  )
}

build_quad_label_df <- function(counts, quad_ids, titles, label_fn,
                                blue_quads, red_quads) {
  tibble::tibble(
    q     = quad_ids,
    title = titles,
    x     = c(-Inf,  Inf,  Inf, -Inf),
    y     = c( Inf,  Inf, -Inf, -Inf),
    hjust = c(   0,    1,    1,    0),
    vjust = c(   1,    1,    0,    0)
  ) |>
    dplyr::left_join(counts, by = "q") |>
    dplyr::mutate(
      label    = label_fn(title, dplyr::pick(dplyr::everything())),
      bg_color = dplyr::if_else(q %in% blue_quads, colors$blue_dark, colors$red_dark)
    )
}

# --- Standalone single-row legend for the shared reversal palette -----------
# Decoupled from any panel's own data -- built from the palette names alone,
# so it never depends on nes_df or df_rev existing in a particular shape.
# breaks = names(pal) makes the legend order deterministic (matches pal's
# insertion order) rather than ggplot's default discrete-scale ordering.
build_shared_legend <- function(pal, panel_width_mm = PANEL_MD) {
  df <- data.frame(x = 1, y = 1, cat = factor(names(pal), levels = names(pal)))
  p <- ggplot(df, aes(x, y, color = cat)) +
    geom_point() +
    scale_color_manual(values = pal, breaks = names(pal), name = NULL) +
    guides(color = guide_legend(nrow = 1, override.aes = list(size = 4))) +
    theme_void() +
    theme(legend.position = "bottom",
          legend.text = element_text(size = round(scale_text(FIG_LEGEND_TEXT, panel_width_mm), 1)))
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