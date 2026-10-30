# ============================================================
# WGCNA Module Analysis — Muscle
# CR Proteomics — full dataset (CTL_T1 + CRC_T1 + CRC_T2)
#
# Network construction -> module eigengenes -> ME mapping across
# groups -> ME vs contrasts (6c lmer, 6d fry, 6e MM-vs-GS) ->
# ME-phenotype correlation -> ORA on module proteins
#
# To run on plasma: duplicate this file as 05_wgcna_plasma.R,
# set config$tissue <- "plasma" and point dalist_path at the
# plasma normalized DAList. Everything else is parameterized.
# ============================================================

# ---- 0. Setup ----
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  WGCNA, tidyverse, ggplot2, patchwork, cowplot, ggtext, ggnewscale,
  writexl, rprojroot, lmerTest, emmeans
)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install("limma")  # 6d/6e fry + MM-vs-GS
library(clusterProfiler)
library(org.Hs.eg.db)

options(stringsAsFactors = FALSE)
WGCNA::enableWGCNAThreads()

root <- rprojroot::find_rstudio_root_file()
source(file.path(root, "04_figures", "a_scripts", "00_theme.R"))  # colors, theme_cr()

# ---- 1. Config ----
config <- list(
  tissue            = "muscle",
  dalist_path       = file.path(root, "02_imputation", "c_data", "muscle",
                                "muscle_imputed_dalist.RDS"),
  group_col         = "cancer_time",   # metadata col: CTL_T1 / CRC_T1 / CRC_T2
  group_levels      = c("CTL_T1", "CRC_T1", "CRC_T2"),
  subject_col       = "pid",           # subject ID, for repeated-measures blocking
  # TODO: character vector of metadata columns to correlate with MEs.
  # NULL = auto-detect all numeric columns in metadata (broad, unfiltered pass)
  pheno_cols        = NULL,
  # Base variable names -- your metadata already has separate pre_<name>/
  # post_<name> columns (e.g. pre_fat_kg, post_fat_kg), repeated identically
  # across a subject's CRC_T1 and CRC_T2 rows (confirmed from your actual
  # dal$metadata: CR006_T1 and CR006_T2 both show pre_fat_kg = 26.4). Delta
  # is computed directly as post_<name> - pre_<name> after deduplicating to
  # one row per subject (see section 7b) -- no cancer_time pivot needed for
  # these, unlike the module eigengene delta below, which DOES need one
  # since ME values genuinely differ per sample row. CTL subjects never
  # contribute here regardless (no PRE/POST pair to begin with).
  delta_pheno_cols  = c("fat_kg", "fat_pct", "ALM_kg", "ALM_index", "BMI",
                        "leg_ext_lbs", "chest_press_lbs", "sts_max_pwr",
                        "grip_lbs"),  # example -- edit to the variables you actually want
  power_range       = c(1:10, seq(12, 30, 2)),
  network_type      = "signed hybrid",
  tom_type          = "signed",
  min_module_size   = 30,
  merge_cut_height  = 0.25,
  deep_split        = 2,
  # modulePreservation() permutation count (section 5c). WGCNA's own
  # recommended default is 200 for results you'd actually report, but that
  # can be slow with ~2000 proteins -- drop to ~50 for a fast first look,
  # then bump back to 200 (or higher) for the number that goes in the
  # manuscript. More permutations = more stable/precise Zsummary estimate,
  # not a different test.
  n_permutations    = 200,
  ora_pvalue_cutoff = 0.05,
  ora_qvalue_cutoff = 0.2,
  fdr_method        = "BH",
  out_root          = file.path(root, "05_wgcna")
)

build_output_path <- function(...) {
  p <- file.path(config$out_root, ...)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- 2. Load data ----
dal <- readRDS(config$dalist_path)

# proteins x samples -> samples x proteins for WGCNA
# TODO: confirm the correct DAList accessor for your build of proteoDA
# (e.g. dal$data vs dal$e_data) — swap in whichever holds normalized abundance
expr_mat <- t(dal$data)
if ("sample_id" %in% colnames(dal$metadata)) {
  meta <- as_tibble(dal$metadata)
} else {
  meta <- dal$metadata %>% as_tibble(rownames = "sample_id")
}

stopifnot(all(rownames(expr_mat) == meta$sample_id))

# ---- 3. QC: good samples/genes + outlier check ----
gsg <- goodSamplesGenes(expr_mat, verbose = 3)
if (!gsg$allOK) {
  expr_mat <- expr_mat[gsg$goodSamples, gsg$goodGenes]
  meta <- meta %>% filter(sample_id %in% rownames(expr_mat))
}

sample_tree <- hclust(dist(expr_mat), method = "average")

pdf_path <- build_output_path("00_qc", paste0("sample_clustering_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 10, height = 6)
plot(sample_tree, main = paste("Sample clustering —", config$tissue),
     sub = "", xlab = "", cex = 0.7)
dev.off()

# Inspect the plot; if there's a clear outlier branch, cut it manually:
# cut_h <- 120
# keep <- cutreeStatic(sample_tree, cutHeight = cut_h, minSize = 10) == 1
# expr_mat <- expr_mat[keep, ]; meta <- meta %>% filter(sample_id %in% rownames(expr_mat))

# ---- 4. Soft-thresholding power ----
sft <- pickSoftThreshold(
  expr_mat, powerVector = config$power_range,
  networkType = config$network_type, verbose = 5
)

sft_df <- sft$fitIndices %>%
  dplyr::mutate(scale_free_r2 = -sign(slope) * SFT.R.sq)

p_sft <- ggplot(sft_df, aes(Power, scale_free_r2, label = Power)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50") +
  geom_text(color = colors$primary) +
  labs(x = "Soft threshold power", y = expression(Signed~R^2),
       title = paste("Scale-free topology fit —", config$tissue)) +
  theme_cr()

p_conn <- ggplot(sft_df, aes(Power, mean.k., label = Power)) +
  geom_text(color = colors$primary) +
  labs(x = "Soft threshold power", y = "Mean connectivity",
       title = paste("Mean connectivity —", config$tissue)) +
  theme_cr()

pdf_path <- build_output_path("00_qc", paste0("soft_threshold_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 10, height = 5)
print(p_sft + p_conn)
dev.off()

# lowest power reaching R^2 >= 0.8; falls back to max tested power otherwise
soft_power <- sft_df$Power[which(sft_df$scale_free_r2 >= 0.8)][1]
if (is.na(soft_power)) soft_power <- max(config$power_range)
message("Selected soft power: ", soft_power)

# ---- 5. Network construction ----
# clusterProfiler/org.Hs.eg.db (loaded above, for the ORA step later in this
# script) pull in Bioconductor's BiocGenerics, which redefines cor() as an
# S4 generic. blockwiseModules() internally calls cor() expecting WGCNA's
# own C-implemented version -- the only one that accepts the extra
# weights.x/weights.y/cosine arguments it passes -- so once BiocGenerics'
# cor sits ahead of WGCNA's in the search path, that internal call breaks
# with "unused arguments". Re-point cor at WGCNA's version for exactly this
# call, then restore it immediately after so nothing later in the script is
# silently affected by a lingering global override.
cor <- WGCNA::cor

net <- blockwiseModules(
  expr_mat,
  power             = soft_power,
  networkType       = config$network_type,
  TOMType           = config$tom_type,
  minModuleSize     = config$min_module_size,
  mergeCutHeight    = config$merge_cut_height,
  deepSplit         = config$deep_split,
  numericLabels     = TRUE,
  saveTOMs          = TRUE,
  saveTOMFileBase   = build_output_path("01_network", paste0("TOM_", config$tissue)),
  pamRespectsDendro = FALSE,
  verbose           = 3
)

cor <- stats::cor  # restore immediately -- see note above

module_labels <- net$colors
module_colors <- labels2colors(module_labels)
names(module_colors) <- colnames(expr_mat)

MEs <- net$MEs %>% orderMEs()
# net$MEs column names are "ME0", "ME1", ... because blockwiseModules() was
# run with numericLabels = TRUE. module_colors (used for module_membership)
# converts those same integer labels via labels2colors() — MEs needs the
# identical conversion, or its column names and module_membership$module
# silently refer to two different label vocabularies for the same modules.
me_numeric_labels <- as.numeric(gsub("^ME", "", colnames(MEs)))
colnames(MEs) <- labels2colors(me_numeric_labels)

pdf_path <- build_output_path("01_network", paste0("dendrogram_modules_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 10, height = 6)
plotDendroAndColors(
  net$dendrograms[[1]], module_colors[net$blockGenes[[1]]],
  groupLabels = "Module", dendroLabels = FALSE,
  addGuide = TRUE, hang = 0.03
)
dev.off()

module_membership <- tibble(
  protein = colnames(expr_mat),
  module  = unname(module_colors)  # module_colors carries protein-ID names from
  # labels2colors() above; strip them here so
  # this column is a plain character vector
)

write_xlsx(
  list(module_membership = module_membership),
  build_output_path("01_network", paste0("module_membership_", config$tissue, ".xlsx"))
)

# ---- 5b. Intramodular connectivity (kME) and hub proteins ----
# kME = correlation of each protein's own expression with each module
# eigengene. A protein's kME with its OWN assigned module is the standard
# WGCNA measure of how "hub-like" it is -- high kME means its expression
# pattern closely tracks the module's overall trend. Hub proteins = top N
# per module ranked by kME_own_module (raw, not absolute value -- module
# assignment already aligns direction for the large majority of member
# proteins under a signed/signed-hybrid network; a NEGATIVE kME_own_module
# would mean a protein assigned to a module while trending opposite its
# eigengene, which is unusual and worth a manual look if it shows up among
# your "top" hits rather than silently ranking it as a strong hub).
# NOT using WGCNA::signedKME() here: it assumes MEs columns start with the
# literal "ME" prefix (WGCNA's standard convention, e.g. "MEred") and
# strips what it assumes are those first 2 characters before prepending
# outputColumnName. Our MEs columns were deliberately renamed to bare color
# names earlier in this pipeline (fixing the numeric-vs-color-label
# mismatch), so signedKME() instead strips the first 2 LETTERS of each
# color name (e.g. "red" -> "d", "black" -> "ack") -- silently mangled
# column names, not an error. kME is mathematically just the correlation
# of each protein's expression with each module eigengene ("signed" only
# means sign-preserving vs. an absolute-value variant elsewhere in WGCNA),
# so computing it directly sidesteps signedKME()'s naming assumption
# entirely rather than working around it.
kME <- as.data.frame(stats::cor(expr_mat, MEs, use = "pairwise.complete.obs"))
rownames(kME) <- colnames(expr_mat)
colnames(kME) <- colnames(MEs)

HUB_TOP_N <- 10  # top N proteins per module by kME with their own module

hub_proteins <- module_membership %>%
  filter(module != "grey") %>%
  rowwise() %>%
  dplyr::mutate(kME_own_module = kME[protein, module]) %>%
  ungroup() %>%
  group_by(module) %>%
  slice_max(kME_own_module, n = HUB_TOP_N, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(module, desc(kME_own_module))

write_xlsx(
  list(hub_proteins = hub_proteins),
  build_output_path("01_network", paste0("hub_proteins_", config$tissue, ".xlsx"))
)

# ---- 5c. Module preservation: whole-dataset modules in CRC-only subset ----
# Network construction (section 5) uses the FULL dataset (CTL + CRC_T1 +
# CRC_T2) rather than a CRC-only subset, primarily because the CRC-only
# subset is only ~12-13 independent SUBJECTS (CRC_T1/CRC_T2 are repeated
# measures on the same people, not independent samples -- WGCNA's network
# estimation has no equivalent of limma's duplicateCorrelation() blocking,
# so those 25 rows behave statistically closer to 12-13 than 25), below
# WGCNA's rough floor of ~15-20 for stable module detection. The whole
# dataset (~22 independent subjects: 10 CTL + 12 CRC) is better-supported.
#
# But the whole-dataset network was built across a real structural
# confound: CTL used Bergstrom needle biopsy, CRC used microneedle --
# module topology could partly reflect that difference rather than pure
# biology. Rather than assume either way, this step tests it directly:
# do the whole-dataset modules hold up (same proteins staying correlated
# with each other) when evaluated ONLY on the CRC-only subset, where that
# confound doesn't exist?
#
# Interpretation (Langfelder & Horvath, 2011 preservation statistics):
#   Zsummary.pres > 10 : strong evidence of preservation
#   Zsummary.pres 2-10 : weak-to-moderate evidence
#   Zsummary.pres < 2  : no evidence of preservation -- this module's
#     structure may be substantially shaped by the CTL-vs-CRC split
#     (including the biopsy-method confound) rather than biology that
#     holds specifically within the CRC-only subset. Treat outcome-variable
#     correlations for a non-preserved module with real caution.
crc_sample_ids <- meta %>%
  dplyr::filter(.data[[config$group_col]] %in% c("CRC_T1", "CRC_T2")) %>%
  dplyr::pull(sample_id)

expr_mat_crc_only <- expr_mat[rownames(expr_mat) %in% crc_sample_ids, , drop = FALSE]

message(sprintf(
  "Module preservation: testing %d whole-dataset modules (excluding grey) against a CRC-only subset (n = %d samples, %d proteins). nPermutations = %d -- this can take a while.",
  length(unique(module_colors)) - 1L,
  nrow(expr_mat_crc_only), ncol(expr_mat_crc_only), config$n_permutations
))

multi_expr <- list(
  whole_dataset = list(data = as.data.frame(expr_mat)),
  crc_only      = list(data = as.data.frame(expr_mat_crc_only))
)
multi_color <- list(whole_dataset = module_colors)

# Same cor() masking issue as section 5 (BiocGenerics S4 generic from
# clusterProfiler/org.Hs.eg.db) -- modulePreservation() also calls cor()
# internally and needs WGCNA's version specifically.
cor <- WGCNA::cor

set.seed(26)
mod_preservation <- modulePreservation(
  multi_expr,
  multi_color,
  referenceNetworks = 1,
  nPermutations     = config$n_permutations,
  randomSeed        = 26,
  quickCor          = 0,
  verbose           = 3
)

cor <- stats::cor  # restore immediately -- see note above

preservation_Z <- mod_preservation$preservation$Z$ref.whole_dataset$inColumnsAlsoPresentIn.crc_only %>%
  tibble::as_tibble(rownames = "module") %>%
  dplyr::select(module, moduleSize, Zsummary.pres)

# medianRank.pres lives in a SEPARATE component from Zsummary.pres -- WGCNA
# splits Z-transformed statistics ($preservation$Z$...) from raw/observed
# statistics ($preservation$observed$...) into two different tables, not
# one combined data frame.
preservation_observed <- mod_preservation$preservation$observed$ref.whole_dataset$inColumnsAlsoPresentIn.crc_only %>%
  tibble::as_tibble(rownames = "module") %>%
  dplyr::select(module, medianRank.pres)

preservation_stats <- preservation_Z %>%
  dplyr::left_join(preservation_observed, by = "module") %>%
  # "gold" is WGCNA's internal random-reference module (used to calibrate
  # the null distribution), not a real co-expression module -- exclude it
  # alongside "grey" (unassigned proteins), same convention as ORA/hub
  # proteins elsewhere in this script.
  dplyr::filter(!module %in% c("gold", "grey")) %>%
  dplyr::arrange(dplyr::desc(Zsummary.pres))

message("\nModule preservation summary (Zsummary.pres > 10 = strong, 2-10 = weak/moderate, <2 = none):")
print(preservation_stats)

p_preservation <- ggplot(preservation_stats,
                         aes(x = reorder(module, Zsummary.pres), y = Zsummary.pres, fill = module)) +
  geom_col() +
  geom_hline(yintercept = c(2, 10), linetype = "dashed", color = "grey40") +
  scale_fill_identity() +
  coord_flip() +
  labs(x = NULL, y = "Zsummary.pres",
       title = paste("Module preservation in CRC-only subset —", config$tissue),
       subtitle = "Dashed lines: Zsummary.pres = 2 (weak) and 10 (strong) preservation thresholds") +
  theme_cr()

pdf_path <- build_output_path("01_network", paste0("module_preservation_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 8, height = max(4, nrow(preservation_stats) * 0.5))
print(p_preservation)
dev.off()

write_xlsx(
  list(preservation_stats = preservation_stats),
  build_output_path("01_network", paste0("module_preservation_", config$tissue, ".xlsx"))
)

# ---- 6. Module eigengenes across groups ----
me_long <- MEs %>%
  as_tibble(rownames = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  pivot_longer(cols = all_of(colnames(MEs)), names_to = "module", values_to = "eigengene") %>%
  dplyr::mutate(!!config$group_col := factor(.data[[config$group_col]], levels = config$group_levels))

p_me_group <- ggplot(me_long, aes(x = .data[[config$group_col]], y = eigengene, fill = module)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_fill_identity() +
  facet_wrap(~ module, scales = "free_y") +
  labs(x = NULL, y = "Module eigengene",
       title = paste("Module eigengenes across groups —", config$tissue)) +
  theme_cr() +
  theme(legend.position = "none")

pdf_path <- build_output_path("02_eigengenes", paste0("me_by_group_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 12, height = 8)
print(p_me_group)
dev.off()

# ---- 6b. ME ~ group, blocking on subject (repeated-measures aware) ----
# CRC_T1/CRC_T2 are paired within subject; CTL_T1 is baseline-only.
# lmer with a subject random intercept respects that pairing, unlike a
# naive one-way ANOVA on the pooled groups.
test_me_group <- function(module_name) {
  df <- me_long %>% filter(module == module_name)
  fm <- reformulate(c(config$group_col, paste0("(1|", config$subject_col, ")")),
                    response = "eigengene")
  fit <- lmerTest::lmer(fm, data = df)
  anova_tab <- anova(fit) %>% as_tibble(rownames = "term") %>% dplyr::mutate(module = module_name)
  pairwise <- emmeans::emmeans(fit, specs = config$group_col) %>%
    pairs() %>% as_tibble() %>% dplyr::mutate(module = module_name)
  list(anova = anova_tab, pairwise = pairwise)
}

me_group_tests <- map(colnames(MEs), test_me_group)
me_group_anova <- map_dfr(me_group_tests, "anova") %>%
  dplyr::mutate(p_adj = p.adjust(`Pr(>F)`, method = config$fdr_method))
me_group_pairwise <- map_dfr(me_group_tests, "pairwise")
# p.value here is ALREADY Tukey HSD-adjusted -- emmeans::pairs() defaults
# to adjust = "tukey" (no adjust= argument was passed above), correcting
# for the 3-way pairwise family WITHIN each module. Applying a second,
# unrelated BH correction on top of an already-adjusted p-value would
# stack two different correction philosophies and be needlessly
# conservative. Per project decision, Tukey alone is the correction used
# for significance calls on these pairwise contrasts -- p_adj is kept as
# a passthrough (not a second adjustment) so downstream code that reads
# $p_adj (mc_long, trajectory_brackets, Panel B's bold/sig flag) doesn't
# need to change what column it looks at.
me_group_pairwise <- me_group_pairwise %>%
  dplyr::mutate(p_adj = p.value)

write_xlsx(
  list(anova = me_group_anova, pairwise_contrasts = me_group_pairwise),
  build_output_path("02_eigengenes", paste0("me_group_stats_", config$tissue, ".xlsx"))
)

# ---- 6c. Module eigengene vs. named contrast (repeated-measures aware) ----
# "Correlating modules with contrasts" isn't a Pearson-style correlation --
# contrasts are group comparisons, not continuous per-sample traits. This
# reuses me_group_pairwise (the lmer/subject-blocked emmeans contrasts from
# 6b) rather than computing a separate point-biserial correlation, which
# would reintroduce the same independence problem already flagged for the
# phenotype heatmap in step 7.
#
# SIGN CONVENTION: emmeans::pairs() reports "earlier factor level minus
# later level" (e.g. "CTL_T1 - CRC_T1"), the OPPOSITE of this pipeline's
# limma contrasts (Baseline_CRCvCTL = CRC_T1 - CTL_T1, confirmed from
# 06_fry_reversal_signature_muscle_CRC.R's makeContrasts() call). Flipping
# sign here so a positive estimate means the same direction as a positive
# logFC in your limma output -- get this wrong and every module's direction
# in this heatmap reads backwards relative to your actual DE results.
contrast_definition <- list(
  Cancer = c(from = "CTL_T1",  to = "CRC_T1"),
  Residual = c(from = "CTL_T1",  to = "CRC_T2"),
  Training = c(from = "CRC_T1", to = "CRC_T2")
)

contrast_display_labels <- c(
  Cancer = "CRC vs CTL",
  Residual = "CRC+RET vs CTL",
  Training = "CRC+RET vs CRC"
)

# Matches one emmeans contrast string (e.g. "CTL_T1 - CRC_T1") to one of
# the 3 named limma contrasts regardless of which order emmeans printed the
# two levels in, and returns the sign flip needed to reorient it to match
# limma's "to - from" convention.
map_emmeans_contrast <- function(contrast_str) {
  parts <- trimws(strsplit(contrast_str, "-")[[1]])
  match_idx <- purrr::detect_index(contrast_definition, function(def) setequal(parts, unname(def)))
  if (match_idx == 0) return(tibble(named_contrast = NA_character_, sign_flip = NA_real_))
  def <- contrast_definition[[match_idx]]
  # emmeans prints "from - to" when its first term equals our `from` --
  # negate to reorient as "to - from", matching limma's convention.
  sign_flip <- if (parts[1] == unname(def["from"])) -1 else 1
  tibble(named_contrast = names(contrast_definition)[match_idx], sign_flip = sign_flip)
}

me_contrast <- me_group_pairwise %>%
  rowwise() %>%
  dplyr::mutate(map_emmeans_contrast(contrast)) %>%
  ungroup() %>%
  filter(!is.na(named_contrast)) %>%
  dplyr::mutate(
    estimate_oriented = estimate * sign_flip,
    t_oriented        = t.ratio * sign_flip,
    # r-equivalent effect size from the subject-blocked model:
    # r = t / sqrt(t^2 + df). Correlation-scale (bounded -1..1) like a
    # classic WGCNA module-trait cell, but still respects CRC_T1/T2
    # pairing -- unlike a naive point-biserial r on a 0/1 group code.
    # Approximate, since df comes from the mixed model (Satterthwaite/KR).
    r_equiv           = t_oriented / sqrt(t_oriented^2 + df),
    named_contrast    = factor(named_contrast, levels = names(contrast_definition))
  )

n_expected <- length(unique(me_group_pairwise$module)) * length(contrast_definition)
message(sprintf(
  "Module x contrast: mapped %d of %d pairwise rows to the 3 named contrasts (expected %d = %d modules x %d contrasts)",
  nrow(me_contrast), nrow(me_group_pairwise), n_expected,
  length(unique(me_group_pairwise$module)), length(contrast_definition)
))
if (nrow(me_contrast) != n_expected) {
  warning("Not every module x contrast combination mapped cleanly -- check emmeans contrast string format against contrast_definition above before trusting this heatmap.")
}

mc_long <- me_contrast %>%
  dplyr::mutate(
    label = sprintf("%.2f\n(%.3f)", estimate_oriented, p.value),
    sig   = p_adj < 0.05
  )

p_mc <- ggplot(mc_long, aes(x = named_contrast, y = module, fill = estimate_oriented)) +
  geom_tile(color = "white") +
  geom_tile(data = filter(mc_long, sig), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red, midpoint = 0) +
  scale_x_discrete(labels = contrast_display_labels[levels(mc_long$named_contrast)]) +
  labs(x = NULL, y = NULL, fill = "\u0394 ME",
       title = paste("Module eigengene vs. contrast \u2014", config$tissue),
       subtitle = "lmer estimate (p); subject-blocked; black outline = Tukey p < 0.05") +
  theme_cr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf_path <- build_output_path("02_eigengenes", paste0("me_contrast_heatmap_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 8, height = max(6, ncol(MEs) * 0.4))
print(p_mc)
dev.off()

# Same heatmap on the r-equivalent scale -- reads like a classic WGCNA
# module-trait heatmap (fixed -1..1 limits), so cells are comparable
# across modules regardless of each eigengene's raw spread.
p_mc_r <- ggplot(mc_long, aes(x = named_contrast, y = module, fill = r_equiv)) +
  geom_tile(color = "white") +
  geom_tile(data = filter(mc_long, sig), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = sprintf("%.2f\n(%.3f)", r_equiv, p.value)), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red,
                       midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = contrast_display_labels[levels(mc_long$named_contrast)]) +
  labs(x = NULL, y = NULL, fill = "r (equiv.)",
       title = paste("Module eigengene vs. contrast —", config$tissue),
       subtitle = "r from lmer t (p); subject-blocked; black outline = Tukey p < 0.05") +
  theme_cr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf_path <- build_output_path("02_eigengenes", paste0("me_contrast_heatmap_r_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 8, height = max(6, ncol(MEs) * 0.4))
print(p_mc_r)
dev.off()

write_xlsx(
  list(module_contrast = mc_long %>% dplyr::select(-sign_flip)),
  build_output_path("02_eigengenes", paste0("me_contrast_", config$tissue, ".xlsx"))
)

# ---- 6d. Module protein-set tests on the named contrasts (limma fry) ----
# Complements 6c at the PROTEIN level: 6c asks whether the module's single
# summary number (the eigengene) shifts in each contrast; fry asks whether
# the module's member proteins, using every protein's own limma statistic,
# shift together in that contrast. Agreement in direction + significance
# between the two is the robustness check worth reporting.
#
# DESIGN: fit here, from this script's own expr_mat/meta, as
#   ~ 0 + cancer_time, blocking on pid via duplicateCorrelation()
# -- the limma analogue of 6b's ME ~ cancer_time + (1|pid). Building it
# locally (rather than loading the 06_fry script's objects) guarantees the
# exact same samples, protein set (post-goodSamplesGenes), and group coding
# as the network. If your main DE model includes covariates, add them to
# `limma_design` below so protein-level estimates match your primary DE.
#
# SIGN: contrasts are built as "to - from" from contrast_definition, the
# same orientation as estimate_oriented in 6c, so fry's Up/Down and 6c's
# sign are directly comparable (checked explicitly below).
#
# fry's own FDR is BH across modules WITHIN each contrast (sets tested
# together in one call), matching the project's BH default.
stopifnot(all(rownames(expr_mat) == meta$sample_id))

limma_expr  <- t(expr_mat)  # proteins x samples, same QC'd protein set as the network
limma_group <- factor(meta[[config$group_col]], levels = config$group_levels)
limma_block <- meta[[config$subject_col]]

limma_design <- stats::model.matrix(~ 0 + limma_group)
colnames(limma_design) <- levels(limma_group)

contrast_strings <- vapply(contrast_definition,
                           function(d) paste0(d[["to"]], " - ", d[["from"]]),
                           character(1))
contrast_matrix <- limma::makeContrasts(contrasts = contrast_strings, levels = limma_design)
colnames(contrast_matrix) <- names(contrast_definition)

limma_corfit <- limma::duplicateCorrelation(limma_expr, limma_design, block = limma_block)
limma_rho    <- limma_corfit$consensus.correlation
message(sprintf("6d: duplicateCorrelation consensus within-subject correlation = %.3f", limma_rho))

limma_fit <- limma::lmFit(limma_expr, limma_design, block = limma_block, correlation = limma_rho) %>%
  limma::contrasts.fit(contrast_matrix) %>%
  limma::eBayes()

module_sets <- split(module_membership$protein, module_membership$module)
module_sets <- module_sets[names(module_sets) != "grey"]
module_idx  <- limma::ids2indices(module_sets, rownames(limma_expr))

module_fry <- purrr::map_dfr(names(contrast_definition), function(cn) {
  res <- limma::fry(
    limma_expr, index = module_idx, design = limma_design,
    contrast = contrast_matrix[, cn],
    block = limma_block, correlation = limma_rho,
    sort = "none"
  )
  if (!"FDR" %in% colnames(res)) res$FDR <- res$PValue  # fry omits FDR when only 1 set
  tibble::as_tibble(res, rownames = "module") %>%
    dplyr::transmute(
      module,
      named_contrast = cn,
      n_proteins     = NGenes,
      direction      = Direction,
      p              = PValue,
      fdr            = FDR,
      p_mixed        = PValue.Mixed,
      fdr_mixed      = FDR.Mixed
    )
}) %>%
  # attach 6c's eigengene estimate for the same module x contrast
  dplyr::left_join(
    me_contrast %>%
      dplyr::transmute(module, named_contrast = as.character(named_contrast),
                       me_estimate_oriented = estimate_oriented, me_p = p.value),
    by = c("module", "named_contrast")
  ) %>%
  dplyr::mutate(
    direction_agrees_with_me = sign(me_estimate_oriented) == dplyr::if_else(direction == "Up", 1, -1),
    signed_neg_log_fdr       = dplyr::if_else(direction == "Up", 1, -1) * -log10(fdr),
    sig                      = fdr < 0.05,
    named_contrast           = factor(named_contrast, levels = names(contrast_definition))
  )

message(sprintf(
  "6d: fry direction agrees with eigengene direction (6c) in %d of %d module x contrast cells; %d cells FDR < 0.05.",
  sum(module_fry$direction_agrees_with_me, na.rm = TRUE), nrow(module_fry), sum(module_fry$sig)
))

p_fry <- ggplot(module_fry, aes(x = named_contrast, y = module, fill = signed_neg_log_fdr)) +
  geom_tile(color = "white") +
  geom_tile(data = dplyr::filter(module_fry, sig), fill = NA, color = "black", linewidth = 1) +
  geom_text(aes(label = sprintf("%s\n(%.3f)", direction, fdr)), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red, midpoint = 0) +
  scale_x_discrete(labels = contrast_display_labels[levels(module_fry$named_contrast)]) +
  labs(x = NULL, y = NULL, fill = "signed\n-log10 FDR",
       title = paste("Module member response (fry) —", config$tissue),
       subtitle = "Direction (FDR); subject-blocked limma; black outline = FDR < 0.05") +
  theme_cr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf_path <- build_output_path("02_eigengenes", paste0("module_fry_heatmap_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = 8, height = max(6, length(module_sets) * 0.4))
print(p_fry)
dev.off()

# Protein-level DE from the same fit -- needed for 6e, exported for
# reference. These are the WGCNA-matched-design estimates; they will only
# match your primary DE tables if that model has no extra covariates.
limma_de <- purrr::map_dfr(names(contrast_definition), function(cn) {
  limma::topTable(limma_fit, coef = cn, number = Inf, sort.by = "none") %>%
    tibble::as_tibble(rownames = "protein") %>%
    dplyr::mutate(named_contrast = cn)
})

write_xlsx(
  c(list(module_fry = module_fry),
    split(limma_de, limma_de$named_contrast)),
  build_output_path("02_eigengenes", paste0("module_fry_", config$tissue, ".xlsx"))
)

# ---- 6e. Module membership (kME) vs. protein significance (GS) ----
# Classic WGCNA "MM vs GS": within each module, correlate each member
# protein's kME with its OWN module against its moderated t for a contrast
# (signed GS -- t rather than logFC so per-protein variance is accounted
# for). A strong positive r means the proteins most central to the module
# are also the ones most changed in that contrast -> the module behaves as
# a coherent unit, not a few DE proteins that happened to land in it.
#
# CAVEAT: cor.test() treats member proteins as independent, which they are
# not (they were grouped BECAUSE they co-vary), so these p-values are
# anti-conservative. Treat r as descriptive; 6d (fry) is the inferential
# test of module-level response.
kME_mat <- as.matrix(kME)

mm_gs_points <- module_membership %>%
  dplyr::filter(module != "grey") %>%
  dplyr::mutate(kME_own_module = kME_mat[cbind(match(protein, rownames(kME_mat)),
                                               match(module, colnames(kME_mat)))]) %>%
  dplyr::inner_join(limma_de %>% dplyr::select(protein, named_contrast, logFC, gs_t = t),
                    by = "protein", relationship = "many-to-many") %>%
  dplyr::mutate(named_contrast = factor(named_contrast, levels = names(contrast_definition)))

mm_gs <- mm_gs_points %>%
  dplyr::group_by(module, named_contrast) %>%
  dplyr::summarise(
    n = dplyr::n(),
    r = if (n >= 5) stats::cor(kME_own_module, gs_t) else NA_real_,
    p = if (n >= 5) stats::cor.test(kME_own_module, gs_t)$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  dplyr::mutate(fdr = p.adjust(p, method = config$fdr_method))

for (cn in names(contrast_definition)) {
  pts  <- dplyr::filter(mm_gs_points, named_contrast == cn)
  labs_df <- dplyr::filter(mm_gs, named_contrast == cn) %>%
    dplyr::mutate(label = sprintf("r = %.2f\np = %.3g", r, p))
  
  p_mmgs <- ggplot(pts, aes(kME_own_module, gs_t)) +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_point(aes(color = module), size = 0.9, alpha = 0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "grey20", linewidth = 0.5) +
    geom_text(data = labs_df, aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 2.8, inherit.aes = FALSE) +
    scale_color_identity() +
    facet_wrap(~ module, scales = "free") +
    labs(x = "kME (own module)", y = "moderated t (GS)",
         title = paste0("Module membership vs. protein significance — ",
                        contrast_display_labels[[cn]], " — ", config$tissue),
         subtitle = "Per-module Pearson r across member proteins; descriptive (proteins are not independent)") +
    theme_cr()
  
  pdf_path <- build_output_path("02_eigengenes", "mm_vs_gs",
                                paste0("mm_vs_gs_", cn, "_", config$tissue, ".pdf"))
  cairo_pdf(pdf_path, width = 12, height = 9)
  print(p_mmgs)
  dev.off()
}

write_xlsx(
  list(mm_gs_summary = mm_gs, mm_gs_points = mm_gs_points),
  build_output_path("02_eigengenes", paste0("mm_vs_gs_", config$tissue, ".xlsx"))
)

# ---- 7. Module eigengene – phenotype correlation ----
if (is.null(config$pheno_cols)) {
  pheno_cols <- meta %>% dplyr::select(where(is.numeric)) %>% colnames()
} else {
  pheno_cols <- config$pheno_cols
}

pheno_mat <- meta %>%
  filter(sample_id %in% rownames(MEs)) %>%
  arrange(match(sample_id, rownames(MEs))) %>%
  dplyr::select(all_of(pheno_cols)) %>%
  as.data.frame()
rownames(pheno_mat) <- rownames(MEs)

# NOTE: WGCNA::cor / corPvalueStudent treat samples as independent. With
# paired CRC_T1/CRC_T2 in the sample set, p-values here are anti-conservative
# for that subset — treat as descriptive/screening, not confirmatory. For a
# confirmatory version, refit trait~ME per module with the same lmer/subject
# blocking pattern used in 6b.
module_trait_cor  <- WGCNA::cor(MEs, pheno_mat, use = "pairwise.complete.obs")

# nSamples must be computed PER TRAIT COLUMN, not once for the whole matrix.
# Phenotype columns with missingness (e.g. years-since-chemo, which is
# structurally NA for CTL) get their correlation computed on a smaller
# pairwise-complete subset via `use = "pairwise.complete.obs"` above, but a
# single nrow(MEs) passed to corPvalueStudent() would apply the FULL sample
# size to every column's p-value regardless -- understating p for any
# column with missingness (same r looks "more significant" at a larger n
# than it actually was computed on).
n_per_trait <- sapply(pheno_mat, function(x) sum(!is.na(x)))
module_trait_pval <- module_trait_cor  # same dimensions, overwritten column-by-column below
for (tr in colnames(module_trait_cor)) {
  module_trait_pval[, tr] <- WGCNA::corPvalueStudent(module_trait_cor[, tr], nSamples = n_per_trait[tr])
}
module_trait_fdr  <- matrix(
  p.adjust(module_trait_pval, method = config$fdr_method),
  nrow = nrow(module_trait_pval), dimnames = dimnames(module_trait_pval)
)

mt_long <- as_tibble(module_trait_cor, rownames = "module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "r") %>%
  left_join(
    as_tibble(module_trait_pval, rownames = "module") %>%
      pivot_longer(-module, names_to = "trait", values_to = "p"),
    by = c("module", "trait")
  ) %>%
  dplyr::mutate(label = sprintf("%.2f\n(%.3f)", r, p))

p_mt <- ggplot(mt_long, aes(trait, module, fill = r)) +
  geom_tile(color = "white") +
  geom_tile(data = filter(mt_long, p < 0.05), fill = NA,
            color = "black", linewidth = 1) +
  geom_text(aes(label = label), size = 2.5) +
  scale_fill_gradient2(low = colors$blue, mid = "white", high = colors$red,
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "r",
       title = paste("Module–trait correlations —", config$tissue),
       subtitle = "r (p); pairwise.complete — see caveat on repeated measures in script") +
  theme_cr() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf_path <- build_output_path("03_module_trait", paste0("module_trait_heatmap_", config$tissue, ".pdf"))
cairo_pdf(pdf_path, width = max(8, length(pheno_cols) * 0.6), height = max(6, ncol(MEs) * 0.4))
print(p_mt)
dev.off()

write_xlsx(
  list(correlation = as_tibble(module_trait_cor, rownames = "module"),
       pvalue      = as_tibble(module_trait_pval, rownames = "module"),
       fdr         = as_tibble(module_trait_fdr, rownames = "module")),
  build_output_path("03_module_trait", paste0("module_trait_stats_", config$tissue, ".xlsx"))
)

# ---- 7b. Module eigengene delta vs. phenotype delta (POST - PRE) ----
# Unlike section 7's raw-value correlation, this correlates CHANGE SCORES:
# for each subject who has BOTH a CRC_T1 and CRC_T2 sample, delta_ME =
# ME(CRC_T2) - ME(CRC_T1) (computed via cancer_time pivot, since MEs
# genuinely differ per sample row). delta_pheno for each variable in
# config$delta_pheno_cols is post_<name> - pre_<name>, taken directly from
# your metadata's existing pre_/post_ columns (which are subject-level
# constants repeated across both of a subject's rows, NOT per-sample
# values) -- no cancer_time pivot needed on that side.
#
# This is a genuinely different question from section 7: "does the SIZE OF
# THE CHANGE in this module track the size of the change in this
# phenotype variable" rather than "does the module's raw level track the
# phenotype's raw level." Mixing a per-subject delta trait against
# per-sample (CTL+PRE+POST all mixed) raw MEs wouldn't align cleanly --
# both sides need to be the same per-subject delta quantity.
#
# STRUCTURALLY CRC-ONLY: only subjects with both a PRE and POST sample
# contribute a delta at all. CTL subjects (single timepoint) can never
# appear here, same restriction already established for years-since-chemo.
delta_modules <- setdiff(colnames(MEs), "grey")
delta_trait_long <- NULL  # stays NULL if delta_pheno_cols has no valid matches this run

build_subject_delta <- function(long_df, value_col, id_col_name) {
  # long_df must have: config$subject_col, config$group_col, and value_col
  wide <- long_df %>%
    dplyr::filter(.data[[config$group_col]] %in% c("CRC_T1", "CRC_T2")) %>%
    dplyr::select(subject = dplyr::all_of(config$subject_col),
                  time = dplyr::all_of(config$group_col),
                  value = dplyr::all_of(value_col)) %>%
    dplyr::distinct(subject, time, .keep_all = TRUE) %>%
    tidyr::pivot_wider(names_from = time, values_from = value)
  if (!all(c("CRC_T1", "CRC_T2") %in% names(wide))) {
    return(tibble::tibble(subject = character(0), delta = numeric(0)) %>%
             stats::setNames(c("subject", id_col_name)))
  }
  tibble::tibble(subject = wide$subject, delta = wide$CRC_T2 - wide$CRC_T1) %>%
    stats::setNames(c("subject", id_col_name))
}

# Module eigengene deltas -- one column per module, joined on subject
me_delta_list <- lapply(delta_modules, function(mod) {
  build_subject_delta(me_long %>% dplyr::filter(module == mod), "eigengene", mod)
})
me_delta <- Reduce(function(a, b) dplyr::full_join(a, b, by = "subject"), me_delta_list)

# Phenotype deltas: metadata already carries pre_<name>/post_<name> as
# subject-level constants repeated across BOTH a subject's CRC_T1 and
# CRC_T2 rows -- delta = post_<name> - pre_<name> directly, no pivot on
# cancer_time needed (unlike the module eigengene delta above). Dedup to
# one row per subject FIRST, or a subject with 2 rows would get counted
# twice in the correlation below.
meta_one_row_per_subject <- meta %>%
  dplyr::distinct(.data[[config$subject_col]], .keep_all = TRUE)

pheno_delta_list <- lapply(config$delta_pheno_cols, function(v) {
  pre_col  <- paste0("pre_", v)
  post_col <- paste0("post_", v)
  if (!pre_col %in% colnames(meta) || !post_col %in% colnames(meta)) {
    warning(sprintf("delta_pheno_cols entry '%s' -- couldn't find both '%s' and '%s' in metadata. Skipping.",
                    v, pre_col, post_col))
    return(NULL)
  }
  meta_one_row_per_subject %>%
    dplyr::transmute(
      subject = .data[[config$subject_col]],
      !!v := .data[[post_col]] - .data[[pre_col]]
    )
})
pheno_delta_list <- purrr::compact(pheno_delta_list)

if (length(pheno_delta_list) == 0) {
  warning("No valid delta_pheno_cols found in metadata -- section 7b skipped. ",
          "Update config$delta_pheno_cols with real column names and re-run from here.")
} else {
  pheno_delta <- Reduce(function(a, b) dplyr::full_join(a, b, by = "subject"), pheno_delta_list)
  
  delta_merged <- dplyr::inner_join(me_delta, pheno_delta, by = "subject")
  
  delta_pheno_names <- setdiff(names(pheno_delta), "subject")
  
  delta_me_mat    <- delta_merged %>% dplyr::select(dplyr::all_of(delta_modules)) %>% as.data.frame()
  delta_pheno_mat <- delta_merged %>% dplyr::select(dplyr::all_of(delta_pheno_names)) %>% as.data.frame()
  rownames(delta_me_mat) <- delta_merged$subject
  rownames(delta_pheno_mat) <- delta_merged$subject
  
  delta_trait_cor <- WGCNA::cor(delta_me_mat, delta_pheno_mat, use = "pairwise.complete.obs")
  
  # Per-column nSamples (same fix as section 7's module_trait_pval) -- a
  # delta variable with more missingness than another must not borrow the
  # full-matrix sample size for its p-value.
  n_per_delta_trait <- sapply(delta_pheno_mat, function(x) sum(!is.na(x)))
  delta_trait_pval <- delta_trait_cor
  for (tr in colnames(delta_trait_cor)) {
    delta_trait_pval[, tr] <- WGCNA::corPvalueStudent(delta_trait_cor[, tr], nSamples = n_per_delta_trait[tr])
  }
  delta_trait_fdr <- matrix(
    p.adjust(delta_trait_pval, method = config$fdr_method),
    nrow = nrow(delta_trait_pval), dimnames = dimnames(delta_trait_pval)
  )
  
  delta_trait_long <- tibble::as_tibble(delta_trait_cor, rownames = "module") %>%
    tidyr::pivot_longer(-module, names_to = "variable", values_to = "r") %>%
    dplyr::left_join(
      tibble::as_tibble(delta_trait_pval, rownames = "module") %>%
        tidyr::pivot_longer(-module, names_to = "variable", values_to = "p"),
      by = c("module", "variable")
    ) %>%
    dplyr::left_join(
      tibble::as_tibble(delta_trait_fdr, rownames = "module") %>%
        tidyr::pivot_longer(-module, names_to = "variable", values_to = "fdr"),
      by = c("module", "variable")
    ) %>%
    dplyr::mutate(sig = fdr < 0.05, n = n_per_delta_trait[variable])
  
  message(sprintf(
    "Module eigengene delta vs. phenotype delta: %d subjects with a CRC_T1+CRC_T2 pair contributed at least one delta variable.",
    dplyr::n_distinct(delta_merged$subject)
  ))
  
  write_xlsx(
    list(delta_trait_long = delta_trait_long),
    build_output_path("03_module_trait", paste0("module_delta_trait_", config$tissue, ".xlsx"))
  )
}

# ---- 8. ORA on module proteins ----
all_proteins <- colnames(expr_mat)

# TODO: confirm fromType matches your protein ID space (SYMBOL/UNIPROT/ENSEMBL)
id_map <- bitr(all_proteins, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db) %>%
  distinct(UNIPROT, .keep_all = TRUE)

universe_entrez <- id_map$ENTREZID

run_module_ora <- function(module_name) {
  mod_proteins <- module_membership %>% filter(module == module_name) %>% pull(protein)
  mod_entrez <- id_map %>% filter(UNIPROT %in% mod_proteins) %>% pull(ENTREZID)
  if (length(mod_entrez) < 5) return(NULL)
  
  enrichGO(
    gene          = mod_entrez,
    universe      = universe_entrez,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    pAdjustMethod = config$fdr_method,
    pvalueCutoff  = config$ora_pvalue_cutoff,
    qvalueCutoff  = config$ora_qvalue_cutoff,
    readable      = TRUE
  )
}

# "grey" = unassigned proteins in WGCNA convention — excluded from ORA
module_names <- unique(module_membership$module) %>% setdiff("grey")
ora_results <- set_names(map(module_names, run_module_ora), module_names) %>% compact()

ora_table <- imap_dfr(ora_results, ~ as_tibble(.x@result) %>% dplyr::mutate(module = .y))

write_xlsx(
  c(list(all_modules = ora_table),
    set_names(map(ora_results, ~ as_tibble(.x@result)), names(ora_results))),
  build_output_path("04_ora", paste0("module_ora_", config$tissue, ".xlsx"))
)

walk(names(ora_results), function(m) {
  res <- ora_results[[m]]
  if (nrow(res@result) == 0) return(invisible(NULL))
  p <- dotplot(res, showCategory = 15) +
    labs(title = paste0(config$tissue, " — ", m, " module ORA (GO:BP)")) +
    theme_cr()
  pdf_path <- build_output_path("04_ora", "dotplots", paste0("ora_dotplot_", m, "_", config$tissue, ".pdf"))
  cairo_pdf(pdf_path, width = 8, height = 6)
  print(p)
  dev.off()
})

# ---- 9. Bundle for reporting ----
report_data <- list(
  config             = config,
  soft_power         = soft_power,
  net                = net,
  module_membership  = module_membership,
  kME                = kME,
  hub_proteins       = hub_proteins,
  preservation_stats = preservation_stats,
  mod_preservation   = mod_preservation,
  MEs                = MEs,
  me_long            = me_long,
  me_group_anova     = me_group_anova,
  me_group_pairwise  = me_group_pairwise,
  mc_long                  = mc_long,
  contrast_display_labels  = contrast_display_labels,
  contrast_definition      = contrast_definition,
  module_fry         = module_fry,
  limma_de           = limma_de,
  limma_rho          = limma_rho,
  mm_gs              = mm_gs,
  mm_gs_points       = mm_gs_points,
  module_trait_cor   = module_trait_cor,
  module_trait_pval  = module_trait_pval,
  module_trait_fdr   = module_trait_fdr,
  delta_trait_long   = delta_trait_long,
  ora_results        = ora_results,
  ora_table          = ora_table
)
saveRDS(report_data, build_output_path("report_data_muscle.rds"))