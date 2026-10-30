# =============================================================================
# 00_theme.R
# Shared color palette and ggplot theme for CRC_pro figures.
#
# Source this near the top of every figure script instead of redefining
# `colors` / `theme_cr()` locally:
#
#   source("04_figures/a_scripts/00_theme.R")
#
# Canonical defaults are base_size = 12, title_hjust = 0.5 (left over from
# the majority of existing figure scripts). Override either argument at a
# given call site (e.g. theme_cr(base_size = 9, title_hjust = 0)) for dense
# panels or left-aligned titles without touching this file.
#
# Figure-specific colors (e.g. a `tissue` or `trajectory` palette used by
# only one script) should stay local to that script, not live here.
#
# Each named hue has _light / _dark / _neon variants. _neon is a brighter,
# more saturated version of the same hue for use against dark backgrounds
# or when a color needs to visually pop (e.g. highlighting a single series).
# Note: brown_neon is a stretch — true "neon brown" isn't really a thing
# since neon implies high saturation/brightness that brown resists by
# definition. I used a vivid copper/amber as the closest approximation;
# swap it out if it doesn't read the way you want.
# =============================================================================
library(ggplot2)
colors <- list(
  red          = "#C0392B", red_light    = "#E6A19A", red_dark     = "#7B241C", red_neon    = "#FF3131",
  blue         = "#2471A3", blue_light   = "#A9CCE3", blue_dark    = "#154360", blue_neon   = "#1F51FF",
  green        = "#4E8B6F", green_light  = "#A9C9B8", green_dark   = "#2E5E4E", green_neon  = "#39FF14",
  pink         = "#D98C9C", pink_light   = "#F2C6CF", pink_dark    = "#A65C6B", pink_neon   = "#FF10F0",
  purple       = "#7B5EA7", purple_light = "#C5B4DD", purple_dark  = "#4A3A73", purple_neon = "#BC13FE",
  yellow       = "#D4A017", yellow_light = "#F1D77A", yellow_dark  = "#9A6F0A", yellow_neon = "#FFFF33",
  orange       = "#C4693A", orange_light = "#E5A07B", orange_dark  = "#8A4321", orange_neon = "#FF5F1F",
  aqua         = "#5DAFA3", aqua_light   = "#A9D6CF", aqua_dark    = "#3A7F76", aqua_neon   = "#0FF0FC",
  brown        = "#8C5A3C", brown_light  = "#C7A189", brown_dark   = "#5A3A28", brown_neon  = "#E8871E",
  red_pale     = "#FADBD8", blue_pale    = "#D6EAF8",
  gray_light   = "#F5F5F5", gray         = "#9E9E9E", gray_dark    = "#4A4A4A",
  black        = "#000000", white        = "#FFFFFF",
  auburn       = c(navy = "#0D234A", auburn = "#E96A2C")
)

theme_cr <- function(base_size = 12, title_hjust = 0.5) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(face = "bold"),
      axis.title      = element_text(face = "bold", size = base_size),
      axis.text       = element_text(face = "bold", size = base_size),
      strip.text      = element_text(face = "bold", size = base_size),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = base_size),
      legend.text     = element_text(size = base_size),
      panel.spacing   = unit(1.5, "lines"),
      panel.border    = element_blank(),
      plot.title      = element_text(face = "bold", size = base_size + 1, hjust = title_hjust)
    )
}