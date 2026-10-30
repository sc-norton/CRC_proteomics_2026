# =============================================================================
# fig_05_panel_C_fry.R
# Panel C: fry rotation-test barcode/running-ES plot for the signature_up /
# signature_down reversal test. Depends on fig_05_style.R being sourced first.
# Also runnable standalone -- see bottom of file.
# =============================================================================

pacman::p_load(dplyr, ggplot2, patchwork, scales, readxl)

running_es <- function(t_vals, in_set) {
  n <- length(t_vals); n_h <- sum(in_set)
  if (n_h == 0) return(rep(0, n))
  hit_w <- ifelse(in_set, abs(t_vals), 0)
  miss_w <- 1 / (n - n_h)
  cumsum(ifelse(in_set, hit_w / sum(hit_w), -miss_w))
}

#' Build Panel C: fry reversal barcode plot
#'
#' @param limma_baseline Raw Baseline_SURVvCTL sheet (needs uniprot, logFC, Pi.Val)
#' @param limma_training_surv Raw Training_SURV sheet (needs uniprot, t, Pi.Val)
#' @param fry_reversal Data frame from fry_reversal_signature RDS
#' @param panel_width_mm Rendered width of Panel C in the final composite (mm)
#' @param sig_col / sig_thresh Significance column + threshold used to define signature sets
#' @param test_contrast_label Display label for the x-axis (e.g. "POST vs PRE")
build_panel_C <- function(limma_baseline, limma_training_surv, fry_reversal,
                          panel_width_mm = PANEL_MD,
                          sig_col = "Pi.Val", sig_thresh = 0.05,
                          test_contrast_label = "POST vs PRE") {

  sig_source <- limma_baseline %>% filter(.data[[sig_col]] < sig_thresh)
  up_ids   <- sig_source$uniprot[sig_source$logFC > 0]
  down_ids <- sig_source$uniprot[sig_source$logFC < 0]

  t_rank <- limma_training_surv %>%
    filter(!is.na(t)) %>%
    arrange(desc(t)) %>%
    mutate(rank = row_number(), in_up = uniprot %in% up_ids, in_down = uniprot %in% down_ids)

  n_all <- nrow(t_rank)
  t_rank$es_up   <- running_es(t_rank$t, t_rank$in_up)
  t_rank$es_down <- running_es(t_rank$t, t_rank$in_down)

  col_up   <- colors$red
  col_down <- colors$blue

  make_barcode_panel <- function(set_name, in_col, es_col, color, fry_row) {
    marks  <- t_rank %>% filter(.data[[in_col]])
    is_sig <- !is.na(fry_row$PValue) && fry_row$PValue < 0.05
    line_color <- if (is_sig) color else scales::alpha(color, 0.45)

    label <- sprintf("p = %.3f, p.Mixed = %.3f, n = %d", fry_row$PValue, fry_row$PValue.Mixed, fry_row$n)
    label_color <- if (isTRUE(fry_row$consistent)) "grey20" else "#DC2626"

    th <- fig_theme_scaled(panel_width_mm)

    p_es <- ggplot(t_rank, aes(x = rank, y = .data[[es_col]])) +
      geom_area(fill = scales::alpha(line_color, 0.15), color = NA) +
      geom_line(color = line_color, linewidth = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.3) +
      annotate("text", x = Inf, y = Inf, label = label, hjust = 1.02, vjust = 1.5,
               size = scale_text(BASE_COUNT, panel_width_mm), fontface = "bold", color = label_color) +
      labs(title = sprintf("%s (n = %d)", set_name, fry_row$n), x = NULL, y = "Running ES") +
      scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
      th + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
                 plot.title = element_text(color = color, face = "bold"))

    p_bc <- ggplot(marks, aes(x = rank, xend = rank, y = 0, yend = 1)) +
      geom_segment(color = line_color, linewidth = 0.3, alpha = 0.7) +
      scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      th +
      theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank(), panel.background = element_rect(fill = "grey97"))

    list(es = p_es, bc = p_bc)
  }

  fry_up_row   <- fry_reversal %>% filter(set == "signature_up")
  fry_down_row <- fry_reversal %>% filter(set == "signature_down")

  panel_up   <- make_barcode_panel("Significant UP",   "in_up",   "es_up",   col_up,   fry_up_row)
  panel_down <- make_barcode_panel("Significant DOWN", "in_down", "es_down", col_down, fry_down_row)

  p_rank <- ggplot(t_rank, aes(x = rank, y = t)) +
    geom_area(fill = scales::alpha("grey40", 0.20), color = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
    labs(x = sprintf("Rank (%s t-statistic, n = %d)", test_contrast_label, n_all), y = NULL) +
    scale_x_continuous(limits = c(1, n_all), expand = c(0.005, 0)) +
    fig_theme_scaled(panel_width_mm)

  fig <- panel_up$es / panel_up$bc / panel_down$es / panel_down$bc / p_rank +
    plot_layout(heights = c(2, 0.25, 2, 0.25, 0.6))

  ttl <- "fry: Reversal Rotation Test"
  sub <- sprintf("Within-subject correlation = %.3f | n = %d proteins", fry_reversal$cor_within[1], n_all)

  list(
    plot           = fig + plot_annotation(title = ttl, subtitle = sub,
                        theme = theme(plot.title = element_text(face = "bold", size = scale_text(FIG_TITLE_SIZE, panel_width_mm)),
                                     plot.subtitle = element_text(color = "grey30", size = scale_text(FIG_SUBTITLE_SIZE, panel_width_mm)))),
    plot_composite = strip_for_composite(fig),   # was: plot_composite = fig -- fig has no title/
                                                  # subtitle of its own yet at this point, but still
                                                  # needs legend.position="none" broadcast via &
                                                  # into all 5 nested rows for composite embedding
    title = ttl, subtitle = sub
  )
}

# --- Standalone run --------------------------------------------------------
if (sys.nframe() == 0 || (exists("RUN_STANDALONE") && isTRUE(RUN_STANDALONE))) {
  source("04_figures/a_scripts/fig_05_style.R")

  paths <- list(
    limma_xlsx = "03_analysis/c_data/muscle/limma_results_all_muscle_cv.xlsx",
    fry_rds    = "03_analysis/c_data/muscle/fry_reversal/fry_reversal_signature_muscle.RDS",
    report_dir = "03_analysis/b_reports/muscle/fry_reversal"
  )
  limma_baseline      <- read_excel(paths$limma_xlsx, sheet = "Cancer")
  limma_training_surv <- read_excel(paths$limma_xlsx, sheet = "Training")
  fry_reversal         <- readRDS(paths$fry_rds)

  result <- build_panel_C(limma_baseline, limma_training_surv, fry_reversal, panel_width_mm = 180)
  save_scaled(result$plot, paths$report_dir, "fry_barcode_muscle_SURV", width_mm = 180, height_mm = 160)
}
