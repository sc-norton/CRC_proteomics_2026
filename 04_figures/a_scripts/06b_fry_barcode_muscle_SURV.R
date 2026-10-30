## =============================================================================
## 06b_fry_barcode_muscle_SURV.R
##
## Barcode / running-ES plot for the signature_up / signature_down fry()
## reversal test, adapted from colleague's make_barcode() design.
##
## Proteins ranked by t-statistic in the TEST contrast (Training_SURV); tick
## marks show where each signature set's members fall in that ranking, with
## a running enrichment score curve and annotation of the fry() PValue,
## PValue.Mixed, Direction, and expected/consistent call.
##
## Depends on outputs of 06_fry_reversal_signature_muscle_SURV.R:
##   - fry_reversal (or its saved RDS)
##   - res_source / res_test result tables (re-loaded here for rank building)
##
## =============================================================================

pacman::p_load(ggplot2, dplyr, readxl, patchwork, scales)
source("04_figures/a_scripts/00_theme.R")

## -----------------------------------------------------------------------
## 0. Config -- reuse paths/labels from the analysis script
## -----------------------------------------------------------------------

config <- list(
  paths = list(
    limma_xlsx = "03_analysis/c_data/muscle/reformat/limma_results_all_muscle.xlsx",
    fry_rds    = "03_analysis/c_data/muscle/fry_reversal/fry_reversal_signature_muscle_SURV.RDS",
    report_dir = "03_analysis/b_reports/muscle/fry_reversal"
  ),
  sig_col    = "Pi.Val",
  sig_thresh = 0.05,
  source_contrast = "Baseline_SURVvCTL",
  test_contrast   = "Training_SURV",
  test_contrast_label = "POST vs PRE"   # matches your display-name convention
)

build_output_path <- function(dir, filename) file.path(dir, filename)  # >>> VERIFY <<< swap for your real helper if signature differs

## -----------------------------------------------------------------------
## 1. Reload result tables + fry results
## -----------------------------------------------------------------------

res_source <- read_excel(config$paths$limma_xlsx, sheet = config$source_contrast)  # >>> VERIFY <<< sheet name
res_test   <- read_excel(config$paths$limma_xlsx, sheet = config$test_contrast)     # >>> VERIFY <<< sheet name

fry_reversal <- readRDS(config$paths$fry_rds)

sig_source <- res_source %>%
  filter(.data[[config$sig_col]] < config$sig_thresh)

up_ids   <- sig_source$uniprot_id[sig_source$logFC > 0]
down_ids <- sig_source$uniprot_id[sig_source$logFC < 0]

## -----------------------------------------------------------------------
## 2. Rank table -- all proteins ranked by t-statistic in Training_SURV
## -----------------------------------------------------------------------

t_rank <- res_test %>%
  filter(!is.na(t)) %>%
  arrange(desc(t)) %>%
  mutate(
    rank    = row_number(),
    in_up   = uniprot_id %in% up_ids,
    in_down = uniprot_id %in% down_ids
  )

n_all <- nrow(t_rank)

## Running enrichment score, weighted by |t| for hits -- same construction
## as colleague's running_es(), applied here to signature_up / signature_down
## against the Training_SURV ranking.
running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w  <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

t_rank$es_up   <- running_es(t_rank$t, t_rank$in_up)
t_rank$es_down <- running_es(t_rank$t, t_rank$in_down)

## -----------------------------------------------------------------------
## 3. Colors -- signature_up / signature_down barcode colors, following the
##    same red = up / blue = down convention used elsewhere in the pipeline
## -----------------------------------------------------------------------

colors$up   <- colors$red
colors$down <- colors$blue

## -----------------------------------------------------------------------
## 4. Barcode panel builder
## -----------------------------------------------------------------------

make_barcode_panel <- function(set_name, in_col, es_col, color, fry_row) {

  marks  <- t_rank %>% filter(.data[[in_col]])
  is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
  line_color <- if (is_sig) color else scales::alpha(color, 0.45)

  label <- sprintf(
    "p = %.3f, p.Mixed = %.3f, n = %d%s",
    fry_row$PValue, fry_row$PValue.Mixed,
    fry_row$n, if (isTRUE(fry_row$consistent)) "" else ""
  )
  label_color <- if (isTRUE(fry_row$consistent)) "grey20" else "#DC2626"

  p_es <- ggplot(t_rank, aes(x = rank, y = .data[[es_col]])) +
    geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
    geom_line(color = line_color, linewidth = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
    annotate("text", x = Inf, y = Inf, label = label,
             hjust = 1.02, vjust = 1.5, size = 2.6, fontface = "bold", color = label_color) +
    labs(title = sprintf("%s (n = %d)", set_name, fry_row$n), x = NULL, y = "Running ES") +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    theme_cr() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
    geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_cr() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "grey97"))

  list(es = p_es, bc = p_bc)
}

fry_up_row   <- fry_reversal %>% filter(set == "signature_up")
fry_down_row <- fry_reversal %>% filter(set == "signature_down")

panel_up   <- make_barcode_panel("signature_up",   "in_up",   "es_up",   colors$up,   fry_up_row)
panel_down <- make_barcode_panel("signature_down", "in_down", "es_down", colors$down, fry_down_row)

## Shared rank axis at the bottom, showing the t-statistic distribution of
## the test contrast -- same role as colleague's p_t panel.
p_rank <- ggplot(t_rank, aes(x = rank, y = t)) +
  geom_area(fill = scales::alpha("grey40", 0.20), color = "grey40", linewidth = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  labs(x = sprintf("Rank (%s t-statistic, n = %d)", config$test_contrast_label, n_all), y = NULL) +
  scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
  theme_cr()

## -----------------------------------------------------------------------
## 5. Assemble compound figure
## -----------------------------------------------------------------------

fry_barcode_fig <- panel_up$es / panel_up$bc / panel_down$es / panel_down$bc / p_rank +
  plot_layout(heights = c(2, 0.25, 2, 0.25, 0.6)) +
  plot_annotation(
    title = "fry reversal test",
    subtitle = sprintf("Within-subject correlation = %.3f", fry_reversal$cor_within[1]),
    theme = theme(plot.title = element_text(face = "bold", size = 11),
                  plot.subtitle = element_text(color = "grey30", size = 9))
  )

## -----------------------------------------------------------------------
## 6. Save -- cairo_pdf per pipeline convention
## -----------------------------------------------------------------------

dir.create(config$paths$report_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  build_output_path(config$paths$report_dir, "fry_barcode_muscle_SURV.pdf"),
  fry_barcode_fig, width = 180, height = 160, units = "mm", device = cairo_pdf
)
ggsave(
  build_output_path(config$paths$report_dir, "fry_barcode_muscle_SURV.png"),
  fry_barcode_fig, width = 180, height = 160, units = "mm", dpi = 300
)

message(sprintf("Saved: %s", build_output_path(config$paths$report_dir, "fry_barcode_muscle_SURV.pdf")))
