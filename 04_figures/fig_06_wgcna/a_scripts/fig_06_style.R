# ============================================================
# fig_06_style.R
# Composite-figure style helpers for fig_06 (WGCNA module summary)
#
# Self-contained version of the pattern used in fig_05_style.R —
# duplicated rather than shared because fig_06 has its own width
# and panel structure. If you'd rather unify these into one shared
# style file down the line, the functions below are drop-in
# compatible with the fig_05 call sites (same names/args).
# ============================================================

pacman::p_load(ggplot2, grid)

# ---- Reference geometry ----
# Text sizes are specified as "the size this would be at REFERENCE_WIDTH_MM"
# and scaled down/up for whatever panel actually renders at. sqrt scaling
# (rather than linear) keeps small panels from producing illegibly tiny text.
FIG_WIDTH_MM       <- 380   # total composite width
REFERENCE_WIDTH_MM <- 178   # single-panel reference width text sizes are specified against

# ---- Text scaling ----
# base_size: the font size you'd want at REFERENCE_WIDTH_MM
# panel_width_mm: actual rendered width of the panel this text lives in
scale_text <- function(base_size, panel_width_mm, reference_width_mm = REFERENCE_WIDTH_MM) {
  if (missing(panel_width_mm)) {
    stop("scale_text() requires panel_width_mm — pass the actual rendered width ",
         "of the panel this text is drawn in, not the full figure width.")
  }
  base_size * sqrt(panel_width_mm / reference_width_mm)
}

# ---- Strip redundant elements before compositing ----
# Removes title/subtitle/caption (rendered once by the stitcher instead) and
# collapses default plot margins so panels sit flush against each other.
# Use `&` rather than `+` when this is applied across a patchwork of >1 plot,
# so the theme change broadcasts to every patch rather than only the last one.
strip_for_composite <- function(p, keep_axis_titles = TRUE) {
  p <- p + theme(
    plot.title    = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption  = element_blank(),
    plot.margin   = margin(2, 2, 2, 2, unit = "mm")
  )
  if (!keep_axis_titles) {
    p <- p + theme(axis.title = element_blank())
  }
  p
}

# ---- Save composite to PDF ----
# width/height in mm, Cairo device for correct Unicode rendering
# (matches the cairo_pdf convention used elsewhere in the pipeline).
save_plot <- function(plot, path, width, height, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  cairo_pdf(path, width = width / 25.4, height = height / 25.4)
  print(plot)
  dev.off()
  invisible(path)
}
