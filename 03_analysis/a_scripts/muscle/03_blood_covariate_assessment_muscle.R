# ============================================================
# Blood covariate choice — Muscle
# CR Proteomics
#
# Compares candidate blood covariates for the limma DA model
# (~ 0 + cancer_time + covariate, blocking on pid) on criteria that do NOT
# depend on how many DAPs each model produces:
#
#   A. Validity   : does the covariate track RBC content (RBC panel score)
#                   WITHIN each group, and how collinear is it with group?
#   B. Fit        : how much protein variance does it explain (residual SD,
#                   share of proteins with a covariate effect, AIC)?
#   C. Controls   : does it strip contrast signal from RBC panel proteins
#                   (positive controls) while leaving non-blood proteins'
#                   t-statistics largely unchanged (negative controls)?
#   D. Residuals  : after adjustment, how much do protein residuals still
#                   track RBC content? Lower = better adjustment.
#
# RBC content is measured from a fixed list of RBC-specific proteins
# (config$rbc_genes) in the normalized data itself -- no dependence on
# WGCNA. Run this BEFORE removing those proteins upstream: they are both
# the RBC score and the positive controls here.
# ============================================================

# --- 0: Setup --------------------------
setwd(rprojroot::find_rstudio_root_file())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, tibble, ggplot2, patchwork, limma, writexl)

source(file.path("04_figures", "a_scripts", "00_theme.R"))  # colors, theme_cr()

# --- 1: Config -------------------------
config <- list(
  norm_RDS     = "01_normalization/c_data/muscle/02_normalized_data_muscle.RDS",
  report_dir   = "03_analysis/b_reports/muscle/covariate_choice",
  group_col    = "cancer_time",
  group_levels = c("CTL_T1", "CRC_T1", "CRC_T2"),
  subject_col  = "pid",
  # RBC-specific proteins: membrane skeleton / band 3 complex, RBC-specific
  # enzymes, heme synthesis/catabolism. Deliberately excludes proteins with
  # substantial muscle expression (e.g. CAT, CA2, PRDX2, SLC2A1) so the
  # score reflects RBC content rather than muscle biology. Edit as needed.
  rbc_genes = c("SPTA1", "SPTB", "ANK1", "EPB41", "EPB42", "SLC4A1", "ADD2",
                "MPP1", "GYPA", "RHAG", "BPGM", "CA1", "ALAD", "HMBS",
                "BLVRB", "SNCA"),
  min_rbc_present = 5,
  # proteins whose within-group |r| with the RBC score is below this count
  # as "non-blood" negative controls
  nonblood_r = 0.2
)
dir.create(config$report_dir, recursive = TRUE, showWarnings = FALSE)

# --- 2: Load ---------------------------
CRC  <- readRDS(config$norm_RDS)
expr <- as.matrix(CRC$data)
meta <- as_tibble(CRC$metadata)
stopifnot(all(colnames(expr) == meta$sample_id))
stopifnot(nrow(CRC$annotation) == nrow(expr))

# Candidate covariates. log(B_M_ratio) is included because a skewed ratio
# gives a few high-blood samples a lot of leverage on a linear covariate.
meta <- meta %>% dplyr::mutate(log_B_M_ratio = log(B_M_ratio))
candidates <- c("blood_pct", "B_M_ratio", "log_B_M_ratio")

missing <- setdiff(c(candidates, config$group_col, config$subject_col), colnames(meta))
if (length(missing) > 0) stop("Missing metadata column(s): ", paste(missing, collapse = ", "))
na_cov <- vapply(meta[candidates], function(x) sum(is.na(x)), integer(1))
if (any(na_cov > 0)) stop("NA in candidate covariate(s): ",
                          paste(names(na_cov)[na_cov > 0], collapse = ", "))

grp <- factor(meta[[config$group_col]], levels = config$group_levels)
center_by_group <- function(x, g) x - ave(x, g, FUN = function(v) mean(v, na.rm = TRUE))

# --- 2b: RBC panel score ---------------
# z-score each panel protein across samples, then average, so the most
# abundant panel proteins don't dominate the score.
rbc_rows <- CRC$annotation$gene %in% config$rbc_genes
rbc_members <- rownames(expr)[rbc_rows]

cat(sprintf("RBC panel: %d of %d genes present: %s\n",
            sum(rbc_rows), length(config$rbc_genes),
            paste(sort(CRC$annotation$gene[rbc_rows]), collapse = ", ")))
if (sum(rbc_rows) < config$min_rbc_present) {
  stop("Fewer than ", config$min_rbc_present, " RBC panel proteins in the data -- ",
       "were they already filtered out? This script needs them present.")
}

rbc_z <- t(scale(t(expr[rbc_rows, , drop = FALSE])))
meta$rbc_index <- colMeans(rbc_z, na.rm = TRUE)

# Negative controls: non-panel proteins that don't track the RBC score
# within groups
expr_c <- t(apply(expr, 1, center_by_group, g = grp))
rbc_c  <- center_by_group(meta$rbc_index, grp)
r_rbc  <- as.vector(stats::cor(t(expr_c), rbc_c, method = "spearman",
                               use = "pairwise.complete.obs"))
nonblood <- rownames(expr)[!rbc_rows & !is.na(r_rbc) & abs(r_rbc) < config$nonblood_r]

cat(sprintf("%d RBC panel proteins (positive controls); %d non-blood proteins (negative controls)\n",
            length(rbc_members), length(nonblood)))

# --- 3: Fit one limma model per candidate --

fit_model <- function(cov) {
  design <- if (is.null(cov)) {
    stats::model.matrix(~ 0 + grp)
  } else {
    stats::model.matrix(~ 0 + grp + meta[[cov]])
  }
  colnames(design) <- c(levels(grp), if (!is.null(cov)) "cov")
  
  cm <- limma::makeContrasts(
    Cancer   = CRC_T1 - CTL_T1,
    Residual = CRC_T2 - CTL_T1,
    Training = CRC_T2 - CRC_T1,
    levels   = design
  )
  rho <- limma::duplicateCorrelation(expr, design,
                                     block = meta[[config$subject_col]])$consensus.correlation
  fit <- limma::lmFit(expr, design, block = meta[[config$subject_col]], correlation = rho)
  
  list(cov    = cov,
       design = design,
       rho    = rho,
       fit_e  = limma::eBayes(fit),
       fit_c  = limma::eBayes(limma::contrasts.fit(fit, cm)),
       resid  = stats::residuals(fit, expr))
}

models <- purrr::map(c(list(none = NULL), as.list(purrr::set_names(candidates))), fit_model)

# --- 4A: Validity ----------------------

collinearity <- purrr::map_dfr(candidates, function(cv) {
  tibble(covariate = cv,
         r2_with_group = summary(stats::lm(meta[[cv]] ~ grp))$r.squared)
})

rbc_agreement <- purrr::map_dfr(candidates, function(cv) {
  meta %>%
    dplyr::mutate(group = grp) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      r_with_rbc_score = stats::cor(.data[[cv]], rbc_index, method = "spearman",
                                    use = "complete.obs"),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(covariate = cv, .before = 1)
})

# --- 4B: Fit ---------------------------

sigma_none <- models$none$fit_e$sigma

fit_summary <- purrr::imap_dfr(models, function(m, nm) {
  out <- tibble(model = nm,
                within_subject_rho = m$rho,
                median_sigma = stats::median(m$fit_e$sigma, na.rm = TRUE),
                median_sigma_ratio_vs_none = stats::median(m$fit_e$sigma / sigma_none, na.rm = TRUE))
  if (!is.null(m$cov)) {
    tt <- limma::topTable(m$fit_e, coef = "cov", number = Inf, sort.by = "none")
    pv <- tt$P.Value[!is.na(tt$P.Value)]
    out <- out %>% dplyr::mutate(
      n_cov_fdr05                   = sum(tt$adj.P.Val < 0.05, na.rm = TRUE),
      prop_proteins_with_cov_effect = 1 - limma::propTrueNull(pv)
    )
  }
  out
})

# Per-protein AIC preference between designs. selectModel() doesn't take
# the subject block, and needs complete rows, so this is descriptive:
# complete-case proteins only, no repeated-measures correction.
expr_cc <- expr[rowSums(is.na(expr)) == 0, , drop = FALSE]
aic_sel <- limma::selectModel(expr_cc, purrr::map(models, "design"), criterion = "aic")
aic_pref <- tibble(model = names(models)) %>%
  dplyr::mutate(n_proteins_preferring = as.integer(table(factor(aic_sel$pref, levels = model))),
                prop_preferring = n_proteins_preferring / nrow(expr_cc))

# --- 4C: Positive / negative controls --

t_long <- purrr::imap_dfr(models, function(m, nm) {
  purrr::map_dfr(colnames(m$fit_c$t), function(cn) {
    tibble(model = nm, contrast = cn,
           protein = rownames(m$fit_c$t),
           t = m$fit_c$t[, cn],
           logFC = m$fit_c$coefficients[, cn])
  })
}) %>%
  dplyr::mutate(protein_class = dplyr::case_when(
    protein %in% rbc_members ~ "RBC panel",
    protein %in% nonblood    ~ "non-blood",
    TRUE                     ~ "other"
  ))

ref_t <- t_long %>%
  dplyr::filter(model == "none") %>%
  dplyr::select(contrast, protein, t_none = t, logFC_none = logFC)

control_summary <- t_long %>%
  dplyr::left_join(ref_t, by = c("contrast", "protein")) %>%
  dplyr::group_by(model, contrast) %>%
  dplyr::summarise(
    # positive controls: lower = more blood signal removed
    rbc_median_abs_t           = stats::median(abs(t[protein_class == "RBC panel"]), na.rm = TRUE),
    # negative controls: closer to 1 / 0 = less collateral change
    nonblood_r_t_vs_none       = stats::cor(t[protein_class == "non-blood"],
                                            t_none[protein_class == "non-blood"],
                                            use = "complete.obs"),
    nonblood_median_abs_dlogFC = stats::median(abs(logFC - logFC_none)[protein_class == "non-blood"],
                                               na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(model = factor(model, levels = names(models))) %>%
  dplyr::arrange(contrast, model)

# --- 4D: Residual blood tracking -------
# Residuals already have group effects (and the covariate) removed, so
# correlating them with the group-centered RBC score asks: how much
# within-group blood signal is LEFT after this adjustment? Read the
# "all other" rows -- panel proteins track the score by construction.
resid_long <- purrr::imap_dfr(models, function(m, nm) {
  r <- as.vector(stats::cor(t(m$resid), rbc_c, method = "spearman",
                            use = "pairwise.complete.obs"))
  tibble(model = nm, protein = rownames(m$resid), r_resid_rbc = r)
}) %>%
  dplyr::mutate(protein_class = dplyr::if_else(protein %in% rbc_members, "RBC panel", "all other"),
                model = factor(model, levels = names(models)))

resid_summary <- resid_long %>%
  dplyr::group_by(model, protein_class) %>%
  dplyr::summarise(
    median_abs_r         = stats::median(abs(r_resid_rbc), na.rm = TRUE),
    prop_abs_r_over_0.35 = mean(abs(r_resid_rbc) > 0.35, na.rm = TRUE),
    .groups = "drop"
  )

# --- 5: Plots --------------------------

p_cov_p <- purrr::imap_dfr(models[candidates], function(m, nm) {
  tibble(covariate = nm,
         P.Value = limma::topTable(m$fit_e, coef = "cov", number = Inf)$P.Value)
}) %>%
  ggplot(aes(P.Value)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = colors$primary, color = "white") +
  facet_wrap(~ covariate, nrow = 1) +
  labs(x = "Covariate coefficient p-value", y = "Proteins",
       title = "Covariate effect across the proteome",
       subtitle = "Taller spike near 0 = covariate explains variance in more proteins") +
  theme_cr()

p_resid <- ggplot(resid_long, aes(r_resid_rbc, color = model)) +
  geom_vline(xintercept = 0, color = "grey60") +
  geom_density(linewidth = 0.6) +
  facet_wrap(~ protein_class, scales = "free_y") +
  labs(x = "Spearman r: residual vs. RBC score (within group)", y = "Density",
       color = "Model",
       title = "Blood signal left after adjustment",
       subtitle = "Narrower around 0 = better adjustment") +
  theme_cr()

p_controls <- t_long %>%
  dplyr::filter(model != "none", protein_class != "other") %>%
  dplyr::left_join(ref_t, by = c("contrast", "protein")) %>%
  ggplot(aes(t_none, t, color = protein_class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 0.6, alpha = 0.5) +
  facet_grid(model ~ contrast) +
  scale_color_manual(values = c("RBC panel" = colors$red, "non-blood" = "grey40")) +
  labs(x = "t (no covariate)", y = "t (adjusted)", color = NULL,
       title = "Contrast t before vs. after adjustment",
       subtitle = "Good covariate: red (RBC) points pulled toward 0; grey on the diagonal") +
  theme_cr()

cairo_pdf(file.path(config$report_dir, "covariate_choice_muscle.pdf"),
          width = 12, height = 18)
print(p_cov_p / p_resid / p_controls + plot_layout(heights = c(1, 1, 2.2)))
dev.off()

# --- 6: Export -------------------------

write_xlsx(
  list(A_collinearity     = collinearity,
       A_rbc_agreement    = rbc_agreement,
       B_fit_summary      = fit_summary,
       B_aic_preference   = aic_pref,
       C_controls         = control_summary,
       D_residual_summary = resid_summary),
  file.path(config$report_dir, "covariate_choice_muscle.xlsx")
)

cat("\n>> A. Agreement with RBC score (within group):\n"); print(rbc_agreement)
cat("\n>> A. Collinearity with group:\n");                print(collinearity)
cat("\n>> B. Fit:\n");                                    print(fit_summary)
cat("\n>> B. AIC preference:\n");                         print(aic_pref)
cat("\n>> C. Controls:\n");                               print(control_summary)
cat("\n>> D. Residual blood tracking:\n");                print(resid_summary)

