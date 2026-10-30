# =============================================================================
# GSEA / ORA Dotplot PDF Report Generator
# Standalone script — loads saved .rds results and produces a single PDF
# with dotplots for every database and every comparison.
# =============================================================================

# --- 0: Setup ----------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()


if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  ggplot2,          # plotting
  patchwork,        # layout (available if needed for future multi-panel pages)
  grid,             # title pages
  clusterProfiler,  # dotplot() method for gseaResult / enrichResult
  enrichplot,       # dotplot() extensions
  org.Hs.eg.db,     # OrgDb (needed if dotplot calls readable = TRUE internally)
  dplyr
)

# --- 1: Configuration --------------------------------------------------------

config <- list(
  report_dir = "03_analysis/b_reports/plasma/reformat"
)

# Folder each comparison lives under (mirrors folder_map in the analysis script)
folder_map <- c(
  Baseline_SURVvCTL = "SURVvCTL",
  Training_SURVvCTL = "SURVvCTL",
  Training_SURV     = "SURVvCTL",
  UP_PRE_only       = "SURVvCTL",
  UP_POST_only      = "SURVvCTL",
  UP_both           = "SURVvCTL",
  DOWN_PRE_only     = "SURVvCTL",
  DOWN_POST_only    = "SURVvCTL",
  DOWN_both         = "SURVvCTL",
  REVERSAL          = "SURVvCTL",
  EXACERBATION      = "SURVvCTL",
  Baseline_CREvPLA  = "CREvPLA",
  Training_CREvPLA  = "CREvPLA",
  Training_CRE      = "CREvPLA",
  Training_PLA      = "CREvPLA",
  Interaction_supp  = "CREvPLA",
  Baseline_FvM      = "FvM",
  Training_FvM      = "FvM",
  Training_F        = "FvM",
  Training_M        = "FvM",
  Interaction_sex   = "FvM"
)

# GSEA comparisons  →  filename pattern: gsea_results_15_<label>.rds
gsea_labels <- c(
  "Baseline_SURVvCTL", "Training_SURVvCTL", "Training_SURV",
  "Baseline_CREvPLA",  "Training_CREvPLA",  "Training_CRE",
  "Training_PLA",      "Interaction_supp"
)

# ORA-only comparisons  →  filename pattern: gsea_results_<label>.rds
ora_labels <- c(
  "UP_PRE_only", "UP_POST_only", "UP_both",
  "DOWN_PRE_only", "DOWN_POST_only", "DOWN_both",
  "REVERSAL", "EXACERBATION"
)

# Build full paths
gsea_rds_paths <- setNames(
  vapply(gsea_labels, function(lbl) {
    file.path(config$report_dir, folder_map[[lbl]], "gsea",
              sprintf("gsea_results_%s.rds", lbl))
  }, character(1)),
  gsea_labels
)

ora_rds_paths <- setNames(
  vapply(ora_labels, function(lbl) {
    file.path(config$report_dir, "SURVvCTL", "gsea",
              sprintf("gsea_results_%s.rds", lbl))
  }, character(1)),
  ora_labels
)

all_rds_paths <- c(gsea_rds_paths, ora_rds_paths)

# Slots present in GSEA results and their display labels
gsea_slots <- c(
  gsea_hall      = "Hallmark (GSEA)",
  gsea_reactome  = "Reactome (GSEA)",
  gsea_bp_simple = "GO Biological Process — simplified (GSEA)",
  gsea_bp_clean  = "GO Biological Process — full (GSEA)",
  gsea_cc_simple = "GO Cellular Component — simplified (GSEA)",
  gsea_cc_clean  = "GO Cellular Component — full (GSEA)",
  gsea_mf_simple = "GO Molecular Function — simplified (GSEA)",
  gsea_mf_clean  = "GO Molecular Function — full (GSEA)",
  gsea_bp_slim   = "GO Slim (GSEA)",
  gsea_kegg      = "KEGG (GSEA)"
)

# Slots present in both GSEA and ORA-only results
ora_slots <- c(
  ora_bp_up   = "GO BP — UP-regulated (ORA)",
  ora_cc_up   = "GO CC — UP-regulated (ORA)",
  ora_mf_up   = "GO MF — UP-regulated (ORA)",
  ora_bp_down = "GO BP — DOWN-regulated (ORA)",
  ora_cc_down = "GO CC — DOWN-regulated (ORA)",
  ora_mf_down = "GO MF — DOWN-regulated (ORA)"
)

all_slots <- c(gsea_slots, ora_slots)

# --- 2: Helper Functions -----------------------------------------------------

#' Draw a plain title page for a comparison
draw_title_page <- function(label) {
  grid.newpage()
  grid.text(
    label = sprintf("Comparison: %s", label),
    x     = 0.5,
    y     = 0.5,
    gp    = gpar(fontsize = 24, fontface = "bold", col = "grey20")
  )
}

#' Return a ggplot dotplot, or NULL if nothing significant / object absent
safe_dotplot <- function(obj, title, showCategory = 20) {
  if (is.null(obj)) return(NULL)
  
  # Must be a recognised result class
  if (!inherits(obj, c("gseaResult", "enrichResult"))) return(NULL)
  
  res <- obj@result
  if (is.null(res) || nrow(res) == 0) return(NULL)
  
  # Only plot if at least one term passes FDR threshold
  if (!any(res$p.adjust < 0.05, na.rm = TRUE)) return(NULL)
  
  tryCatch({
    dotplot(obj, showCategory = showCategory, title = title) +
      theme_bw(base_size = 9) +
      theme(
        plot.title   = element_text(size = 10, face = "bold", hjust = 0.5),
        axis.text.y  = element_text(size = 7),
        axis.text.x  = element_text(size = 7),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 8)
      )
  }, error = function(e) {
    message(sprintf("   dotplot() failed for '%s': %s", title, e$message))
    NULL
  })
}

#' Draw a blank placeholder page with a message
draw_blank_page <- function(msg) {
  p <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg,
             size = 5, colour = "grey60", hjust = 0.5, vjust = 0.5) +
    theme_void()
  print(p)
}

# --- 3: Main Report Function -------------------------------------------------

#' Generate a single PDF with dotplots for every comparison × database
#'
#' @param rds_paths   Named character vector: names = comparison labels,
#'                    values = paths to .rds files produced by the analysis script.
#' @param output_pdf  Full path for the output PDF.
#' @param showCategory Number of top terms to display per dotplot (default 20).
#' @param page_width  PDF page width in inches (default 11, landscape).
#' @param page_height PDF page height in inches (default 8.5, landscape).
#' @param skip_missing If TRUE, silently skip comparisons whose .rds file does
#'                     not exist (default TRUE).

generate_gsea_pdf_report <- function(
    rds_paths,
    output_pdf,
    showCategory  = 20,
    page_width    = 11,
    page_height   = 8.5,
    skip_missing  = TRUE
) {
  cat(sprintf("\n>> Starting PDF report: %s\n", output_pdf))
  
  # Ensure output directory exists
  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)
  
  pdf(output_pdf, width = page_width, height = page_height)
  on.exit(dev.off(), add = TRUE)
  
  n_comparisons <- length(rds_paths)
  
  for (i in seq_along(rds_paths)) {
    label <- names(rds_paths)[i]
    path  <- rds_paths[[i]]
    
    cat(sprintf("   [%d/%d] %s\n", i, n_comparisons, label))
    
    # --- Handle missing file ---
    if (!file.exists(path)) {
      if (skip_missing) {
        message(sprintf("   Skipping '%s': file not found (%s)", label, path))
        next
      } else {
        stop(sprintf("File not found for '%s': %s", label, path))
      }
    }
    
    res <- tryCatch(readRDS(path), error = function(e) {
      message(sprintf("   Failed to load '%s': %s", label, e$message))
      NULL
    })
    if (is.null(res)) next
    
    # --- Section title page ---
    draw_title_page(label)
    
    # --- One page per database slot ---
    for (slot_name in names(all_slots)) {
      db_label <- all_slots[[slot_name]]
      obj      <- res[[slot_name]]
      
      if (is.null(obj)) next   # slot absent in this result type (e.g. ORA-only)
      
      full_title <- sprintf("%s  |  %s", label, db_label)
      p <- safe_dotplot(obj, title = full_title, showCategory = showCategory)
      
      if (is.null(p)) {
        draw_blank_page(sprintf("No significant terms (p.adjust \u2265 0.05)\n%s", full_title))
      } else {
        print(p)
      }
    }
  }
  
  cat(sprintf("\n>> Report saved: %s\n", output_pdf))
  invisible(output_pdf)
}

# --- 4: Run the Report -------------------------------------------------------

generate_gsea_pdf_report(
  rds_paths    = all_rds_paths,
  output_pdf   = file.path(config$report_dir, "GSEA_dotplot_report.pdf"),
  showCategory = 20
)
