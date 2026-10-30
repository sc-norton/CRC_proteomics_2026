# === Contamination Assessment ===============
# --- 0: Setup ---------------------
# load packages
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, stringr, purrr, ggplot2,
               ggrepel, patchwork, org.Hs.eg.db, AnnotationDbi, tibble, ggpubr,
               UniprotR)

# set working directory
setwd(rprojroot::find_rstudio_root_file())

# Configuration
config <- list(
  # Raw data files - muscle
  raw_muscle = "00_input/muscle/CRm_raw.csv",
  meta_muscle = "00_input/muscle/CRm_meta.csv",
  hpa_muscle = "00_input/muscle/HPA_muscle.csv",
  go_map_muscle = "00_input/muscle/CRm_raw_map.csv",
  # Normalized data files - muscle
  norm_non = "01_normalization/muscle/",
  norm_mod = "01_normalization/muscle/",
  norm_str = "01_normalization/muscle/",
  imp_non = "02_imputation/muscle/",
  imp_mod = "02_imputation/muscle/",
  imp_str = "02_imputation/muscle/",
  # Directories
)

# --- 1: Load & Restructure Raw Data -------------------------------
# Load data
raw <- read.csv(file.path(base_dir, "01_input", "muscle", "CRm_raw.csv"))
meta <- read.csv(file.path(base_dir, "01_input", "muscle", "CRm_meta.csv"))

# Rename annotation column names
annot_cols <- c("uniprot_id", "protein", "gene", "description", "n_seq", "sum_spectra")
colnames(raw)[1:6] <- annot_cols

# Make Uniprot_ID the row name (or gene)


# Split input data into protein intensity data and annotation data
intensity <- raw[ , setdiff(names(raw), annot_cols)] 
annotation <- raw[ ,annot_cols] 

cat(sprintf("\n Raw data: %d proteins x %d samples\n", nrow(raw), ncol(intensity)))
cat(sprintf("\n Meta data: %d samples\n", nrow(meta)))

# Align row names of metadata with column names of data
cat(" Aligning data and metadata...\n")
rownames(meta) <- meta[ ,1]

meta <- meta[order(meta$sample_id), ]

cat(" Data aligned to metadata.\n")

# Validate alignment
cat(" Validating alignment...\n")

data_samples <- colnames(intensity)
meta_samples <- meta$sample_id

stopifnot(
  "Sample mismatch between data and metadata" =
    setequal(data_samples, meta_samples)
)

# Reorder intensity columns to match metadata row order
intensity <- intensity[, meta_samples]

cat("   Validation passed: all", length(meta_samples), "samples aligned.\n")

# --- 1a: Load Normalized and Imputed Data ---------------------------
norm_non <- read.csv(file.path(base_dir, "02_normalization", "c_data", "muscle_stringent",
                  "01_normalized_data_SURV_none_muscle.csv"))
norm <- read.csv(file.path(base_dir, "02_normalization", "c_data", "muscle",
                            "01_normalized_data_SURV_muscle.csv"))
norm_str <- read.csv(file.path(base_dir, "02_normalization", "c_data", "muscle_stringent",
                              "ratio_15", "01_normalized_data_SURV_stringent_muscle.csv"))
imp <- read.csv(file.path(base_dir, "03_imputation", "c_data", "muscle", "01_imputed_SURV.csv"))

# --- 2a: Identify Muscle & Blood Proteins ----------------------------------
cat("\n>> 2 - Contaminant Assessment\n")

canonical_muscle <- c(
  # Contractile
  "MYH1", "MYH2", "MYH4", "MYH7", "ACTA1", "TPM1", "TPM2", "TPM3",
  "TNNT1", "TNNT3", "TNNC1", "TNNC2", "TNNI1", "TNNI2",
  # Sarcomeric/structural
  "TTN", "NEB", "MYOM1", "MYOM2", "FLNC",
  # Metabolic
  "CKM", "CKMT2", "PYGM", "PFKM", "ENO3", "ALDOA",
  # Mitochondrial
  "COX4I1", "SDHA", "UQCRC1", "ATP5F1A",
  # Already in your top detected
  "MB", "CA3", "FABP3"
)

immunoglobulins <- annotation[grepl("^IGH|^IGK|^IGL", annotation$gene), ] %>%
  pull(gene)

keratins <- annotation[grepl("^KRT", annotation$gene), ] %>%
  pull(gene)

serpins <- annotation[grepl("^SERPIN", annotation$gene), ] %>%
  pull(gene)

hemoglobins <- annotation[grepl("^HB", annotation$gene), ] %>%
  pull(gene)

blood_blacklist_stringent <- c(
  # Albumin & transport
  "ALB", "TTR", "APOA1", "APOA2", "APOA4", "APOB", "APOC1", "APOC2", "APOC3", 
  "APOD", "APOE", "APOH", "APOL1", "APOL2", "APOC4", "APOM", "CLU",
  # Coagulation
  "FGA", "FGB", "FGG", "F2", "F5", "F7", "F10", "F13A1", "F13B",
  "SERPINC1", "SERPIND1", "PROS1", "PROC", "VWF",
  # Complement
  "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "C4BPA", "C4BPB",
  "C5", "C6", "C7", "C8A", "C8B", "C8G", "C9",
  # Hemoglobin
  "HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBG2",
  # Immunoglobulins (Fc regions / constant chains)
  "IGHA1", "IGHA2", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "IGHM",
  "IGKC", "IGLC1", "IGLC2",
  # Acute phase / other plasma
  "CRP", "SAA1", "SAA2", "ORM1", "ORM2", "HP", "HPX", "TF",
  "A1BG", "A2M", "GC", "PLG", "KNG1", "ITIH1", "ITIH2", "ITIH3", "ITIH4",
  # Platelet markers 
  "PPBP", "PF4", "THBS1", "SELP",
  # Neutrophil/immune - NEW ADDITIONS***
  "S100A8", "S100A9", "S100A6", "LYZ", "DEFA3", "JCHAIN",
  # Complement regulators - NEW ADDITIONS***
  "CFH", "CFB", "CFI", "SERPING1", "SERPINA1", "SERPINF2", "SERPINA3", 
  "SERPINA6", "SERPIND1", "SERPINF1", "SERPINH1", "SERPINB1", "SERPINB6",
  # Platelet/exosome - NEW ADDITIONS***
  "FLOT1", "FLOT2", "CD9", "CD81", "ANXA5",
  # Keratins
  "KRT1", "KRT2", "KRT8", "KRT9", "KRT10", "KRT19", "KRT6B", "KRT74",
  immunoglobulins,
  keratins,
  serpins,
  hemoglobins
)

blood_hiab_genes <- c("ALB", "HBA1", "HBB")

# Check which are actually detected in your dataset
detected_canonical <- intersect(canonical_muscle, raw$gene)
cat("Detected canonical muscle proteins:", length(detected_canonical), "of", length(canonical_muscle), "\n")

# Map blacklist gene symbols to uniprot IDs using your raw annotation
detected_blacklist <- raw %>%
  filter(gene %in% blood_blacklist_stringent) %>%
  pull(gene)

# --- 2b: Assess Protein Intensities -----------------------------
# Identify sample columns
sample_cols <- setdiff(colnames(raw), annot_cols)

# Total intensity per sample (calculate separately first)
total_intensity <- raw %>%
  dplyr::select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

alb_intensity <- raw %>%
  filter(gene == "ALB") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

hba1_intensity <- raw %>%
  filter(gene == "HBA1") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

hbb_intensity <- raw %>%
  filter(gene == "HBB") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

mb_intensity <- raw %>%
  filter(gene == "MB") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

ckm_intensity <- raw %>%
  filter(gene == "CKM") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh7_intensity <- raw %>%
  filter(gene == "MYH7") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh1_intensity <- raw %>%
  filter(gene == "MYH1") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh2_intensity <- raw %>%
  filter(gene == "MYH2") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

int_index <- data.frame(
  sample_id = names(ckm_intensity),
  total_intensity = total_intensity,
  alb_intensity = alb_intensity,
  hba1_intensity = hba1_intensity,
  hbb_intensity = hbb_intensity,
  mb_intensity = mb_intensity,
  ckm_intensity = ckm_intensity,
  myh7_intensity = myh7_intensity,
  myh1_intensity = myh1_intensity,
  myh2_intensity = myh2_intensity
)

int_index <- int_index %>%
  mutate(
  alb_pct = alb_intensity / total_intensity,
  hba1_pct = hba1_intensity / total_intensity,
  hbb_pct = hbb_intensity / total_intensity,
  mb_pct = mb_intensity / total_intensity,
  ckm_pct = ckm_intensity / total_intensity,
  myh7_pct = myh7_intensity / total_intensity,
  myh1_pct = myh1_intensity / total_intensity,
  myh2_pct = myh2_intensity / total_intensity
)

int_index <- int_index %>%
  mutate(group = case_when(
    grepl("_T1$", sample_id) ~ "SURV T1",
    grepl("_T2$", sample_id) ~ "SURV T2",
    grepl("^PPS", sample_id) ~ "CTL"
  ))

int_index_long <- int_index %>%
  dplyr::select(sample_id, group,
         alb_pct, hba1_pct, hbb_pct, mb_pct, ckm_pct, myh7_pct, myh1_pct, myh2_pct) %>%
  pivot_longer(
    cols = -c(sample_id, group),
    names_to = "pro_pct",
    values_to = "value"
  )

int_index_long <- int_index_long %>%
  mutate(protein = case_when(
    grepl("hba1_pct", pro_pct) ~ "HBA1",
    grepl("hbb_pct", pro_pct) ~ "HBB",
    grepl("alb_pct", pro_pct) ~ "ALB",
    grepl("mb_pct", pro_pct) ~ "MB",
    grepl("ckm_pct", pro_pct) ~ "CKM",
    grepl("myh7_pct", pro_pct) ~ "MYH7",
    grepl("myh1_pct", pro_pct) ~ "MYH1",
    grepl("myh2_pct", pro_pct) ~ "MYH2"
  ))

# Plot intensities
ggplot(int_index_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~protein, scales = "free_y") +
  theme_classic() +
  labs(x = NULL, y = "Protein % of total intensity") +
  scale_fill_manual(values = c(
    "SURV T1" = "#1b9e77",
    "SURV T2" = "#d95f02",
    "CTL" = "#7570b3"
  )) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

int_index_long2 <- int_index_long %>%
  filter(protein %in% c("ALB", "HBA1", "HBB"))

summary(lm(value ~ group, data = int_index_long2))

# --- 2c: Assess Contaminant Intensities --------------------------

# Select only raw for retained participants after normalization:
#raw2 <- raw[,-c(13,15)]

# Contamination intensity per sample (only hi abundance)
contam_intensity_hiab <- raw %>%
  filter(gene %in% blood_hiab_genes) %>%
  dplyr::select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)


# Contamination intsensity per sample (whole blacklist)
contam_intensity_all <- raw %>%
  filter(gene %in% blood_blacklist_stringent) %>%
  dplyr::select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

# Combine into index df
contam_index <- data.frame(
  sample_id      = names(contam_intensity_hiab),
  contam_intensity_hiab   = contam_intensity_hiab,
  contam_intensity_all = contam_intensity_all,
  total_intensity = total_intensity[names(contam_intensity_hiab)]
) %>%
  mutate(contam_pct_hiab = contam_intensity_hiab / total_intensity,
         contam_pct_all  = contam_intensity_all / total_intensity)

rownames(contam_index) <- NULL

# Join to metadata
meta <- meta %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_intensity_hiab, contam_pct_hiab, 
                            contam_intensity_all, contam_pct_all), 
            by = "sample_id")

# Quick check — compare between groups
meta %>%
  group_by(cancer_time) %>%
  summarise(
    mean_contam_pct_hiab  = mean(contam_pct_hiab,  na.rm = TRUE),
    median_contam_pct_hiab = median(contam_pct_hiab, na.rm = TRUE),
    sd_contam_pct_hiab    = sd(contam_pct_hiab,    na.rm = TRUE),
    mean_contam_pct_all  = mean(contam_pct_all,  na.rm = TRUE),
    median_contam_pct_all = median(contam_pct_all, na.rm = TRUE),
    sd_contam_pct_all    = sd(contam_pct_all,    na.rm = TRUE)
  )

sample_cols <- colnames(norm_str[, 7:41])

# --- 2d: Create Dataframes of Muscle Proteins & Blood Contaminants ----------
canon_df_norm <- norm_str %>%
  filter(gene %in% canonical_muscle) %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
  group_by(sample_id) %>%
  summarise(
    mean_canonical  = mean(intensity, na.rm = TRUE),
    total_canonical = sum(intensity,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    norm_str %>%
      dplyr::select(all_of(sample_cols)) %>%
      pivot_longer(everything(), names_to = "sample_id", values_to = "intensity") %>%
      group_by(sample_id) %>%
      summarise(total_norm = sum(intensity, na.rm = TRUE), .groups = "drop"),
    by = "sample_id"
  ) %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_intensity_hiab, contam_pct_hiab, 
                            contam_intensity_all, contam_pct_all, total_intensity), 
            by = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(canonical_pct = total_canonical / total_norm)

canon_df_none <- norm_non %>%
  filter(gene %in% canonical_muscle) %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
  group_by(sample_id) %>%
  summarise(
    mean_canonical  = mean(intensity, na.rm = TRUE),
    total_canonical = sum(intensity,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    norm_non %>%
      dplyr::select(all_of(sample_cols)) %>%
      pivot_longer(everything(), names_to = "sample_id", values_to = "intensity") %>%
      group_by(sample_id) %>%
      summarise(total_norm = sum(intensity, na.rm = TRUE), .groups = "drop"),
    by = "sample_id"
  ) %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_intensity_hiab, contam_pct_hiab, 
                            contam_intensity_all, contam_pct_all, total_intensity), 
            by = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(canonical_pct = total_canonical / total_norm)

canon_df_str <- norm_str %>%
  filter(gene %in% canonical_muscle) %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
  group_by(sample_id) %>%
  summarise(
    mean_canonical  = mean(intensity, na.rm = TRUE),
    total_canonical = sum(intensity,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    norm_str %>%
      dplyr::select(all_of(sample_cols)) %>%
      pivot_longer(everything(), names_to = "sample_id", values_to = "intensity") %>%
      group_by(sample_id) %>%
      summarise(total_norm = sum(intensity, na.rm = TRUE), .groups = "drop"),
    by = "sample_id"
  ) %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_intensity_hiab, contam_pct_hiab, 
                            contam_intensity_all, contam_pct_all, total_intensity), 
            by = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(canonical_pct = total_canonical / total_norm)

canon_df_imp <- imp %>%
  filter(gene %in% canonical_muscle) %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
  group_by(sample_id) %>%
  summarise(
    mean_canonical  = mean(intensity, na.rm = TRUE),
    total_canonical = sum(intensity,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    imp %>%
      dplyr::select(all_of(sample_cols)) %>%
      pivot_longer(everything(), names_to = "sample_id", values_to = "intensity") %>%
      group_by(sample_id) %>%
      summarise(total_norm = sum(intensity, na.rm = TRUE), .groups = "drop"),
    by = "sample_id"
  ) %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_intensity_hiab, contam_pct_hiab, 
                            contam_intensity_all, contam_pct_all, total_intensity), 
            by = "sample_id") %>%
  left_join(meta, by = "sample_id") %>%
  mutate(canonical_pct = total_canonical / total_norm)

# --- 2e: Plot Correlations of Contaminant with Muscle --------------------
# These plots will identify whether there is a "suppression" effect due to
# blood protein contamination.

# Normalized data - original analysis
p_norm <- ggplot(canon_df_norm,
                        aes(x = contam_pct_hiab.x, y = canonical_pct, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Standard Blood Filter - Original Analysis",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Proportion (% of normalized intensity)"
  ) +
  theme_minimal()

# Normalized data - no blood filter
p_none <- ggplot(canon_df_none,
                 aes(x = contam_pct_hiab.x, y = canonical_pct, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "No Blood Filter",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Proportion (% of normalized intensity)"
  ) +
  theme_minimal()

# Normalized data - stringent blood filter
p_str <- ggplot(canon_df_str,
                 aes(x = contam_pct_hiab.x, y = canonical_pct, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c(colors$green, colors$blue, colors$red), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Stringent Blood Filter",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Proportion (% of normalized intensity)"
  ) +
  theme_minimal()

# Normalized data - imputed dataset, original analysis
p_imp <- ggplot(canon_df_imp,
                 aes(x = contam_pct_hiab.x, y = canonical_pct, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Standard Blood Filter - Imputed Data",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Proportion (% of normalized intensity)"
  ) +
  theme_minimal()

# --- Plot against normalized intensity instead of proportion
# Normalized data - original analysis
p_norm_int <- ggplot(canon_df_norm,
                 aes(x = contam_pct_hiab.x, y = mean_canonical, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Standard Blood Filter - Original Analysis",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Intensity (mean normalized)"
  ) +
  theme_minimal()

# Normalized data - no blood filter
p_none_int <- ggplot(canon_df_none,
                 aes(x = contam_pct_hiab.x, y = mean_canonical, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "No Blood Filter",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Intensity (mean normalized)"
  ) +
  theme_minimal()

# Normalized data - stringent blood filter
p_str_int <- ggplot(canon_df_str,
                aes(x = contam_pct_hiab.x, y = mean_canonical, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Stringent Blood Filter",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Intensity (mean normalized)"
  ) +
  theme_minimal()

# Normalized data - imputed dataset, original analysis
p_imp_int <- ggplot(canon_df_imp,
                aes(x = contam_pct_hiab.x, y = mean_canonical, color = cancer_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "grey2",
              fill = "grey", linewidth = 0.8) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  scale_color_manual(values = c("darkgreen", "blue3", "red4"), name = "Group") +
  stat_cor(method = "spearman",
           label.x.npc = 0.05, label.y.npc = 0.95) +
  labs(
    title    = "Contamination vs Canonical Muscle Protein Abundance",
    subtitle = "Standard Blood Filter - Imputed Data",
    x        = "Contamination Proportion (% of raw intensity)",
    y        = "Canonical Muscle Protein Intensity (mean normalized)"
  ) +
  theme_minimal()

# --- 2x: Assess cellular compartments ----------------------------
raw_map <- read.csv(file.path(base_dir, "01_input", "muscle", "CRm_raw_map.csv"))
head(raw_map)

raw_map <- raw_map %>%
  # Split by ";" into a list column
  mutate(comp_list = str_split(Gene.Ontology..cellular.component., ";")) %>%
  # Remove GO IDs and trim whitespace
  mutate(comp_list = lapply(comp_list, function(x) str_trim(str_remove(x, "\\s*\\[GO:\\d+\\]"))))

raw_map <- raw_map %>%
  mutate(
    comp1 = sapply(comp_list, function(x) if(length(x) >= 1) x[1] else NA),
    comp2 = sapply(comp_list, function(x) if(length(x) >= 2) x[2] else NA),
    comp3 = sapply(comp_list, function(x) if(length(x) >= 3) x[3] else NA)
  )

raw <- raw %>%
  left_join(raw_map %>% 
              dplyr::select(From, comp1, comp2, comp3), by = c("uniprot_id" = "From"))

norm <- norm %>%
  left_join(raw_map %>% 
              dplyr::select(From, comp1, comp2, comp3), by = c("uniprot_id" = "From"))

norm_non <- norm_non %>%
  left_join(raw_map %>% 
              dplyr::select(From, comp1, comp2, comp3), by = c("uniprot_id" = "From"))

norm_str <- norm_str %>%
  left_join(raw_map %>% 
              dplyr::select(From, comp1, comp2, comp3), by = c("uniprot_id" = "From"))

# Create a map for like terms
comp_map <- tibble(
  pattern = c("acrosomal", 
              "actin", 
              "axon", 
              "autolysosome",
              "autophagosome", 
              "azurophil",
              "blood microparticle",
              "brush border",
              "calcium",
              "cardiac",
              "chromosome",
              "Golgi",
              "clathrin",
              "collagen",
              "COP",
              "Cul",
              "cytoplasm",
              "cytosol",
              "endosome",
              "endoplasmic reticulum",
              "extracellular",
              "lysosomal",
              "mitochondrial",
              "nuclear",
              "sarcoplasmic reticulum"),
  category = c("acrosomal", 
               "actin", 
               "axon", 
               "autolysosome",
               "autophagosome", 
               "azurophil",
               "blood microparticle",
               "brush border",
               "calcium",
               "cardiac",
               "chromosome",
               "Golgi",
               "clathrin",
               "collagen",
               "COP",
               "Cul",
               "cytoplasm",
               "cytosol",
               "endosome",
               "endoplasmic reticulum",
               "extracellular",
               "lysosomal",
               "mitochondrial",
               "nuclear",
               "sarcoplasmic reticulum")
)

condense_compartment <- function(x) {
  dplyr::case_when(
    # Cytoplasm / Cytosol
    grepl("cytoplasm|cytosol", x, ignore.case = TRUE) ~ "Cytoplasm/Cytosol",
    
    # Mitochondria
    grepl("mitochondr|MICOS complex|m-AAA|endopeptidase Clp", x, ignore.case = TRUE) ~ "Mitochondria",
    
    # Skeletal muscle-specific
    grepl("sarcomere|sarcoplasmic reticulum|contractile muscle fiber|cardiac myofibril|M band|A band|costamere|longitudinal sarcoplasmic|cardiac Troponin|ankyrin-1|dystrophin", x, ignore.case = TRUE) ~ "Sarcoplasmic/Myofibrillar",  # + dystrophin
    
    # Actin / Cytoskeleton
    grepl("actin|cytoskeleton|actomyosin|intermediate filament|microtubule", x, ignore.case = TRUE) ~ "Cytoskeleton/Contractile",  # + intermediate filament, microtubule
    
    # Blood microparticle
    grepl("blood microparticle", x, ignore.case = TRUE) ~ "Blood Microparticle",
    
    # Extracellular / Secreted
    grepl("extracellular|basement membrane|collagen|elastic fiber|discoidal|extracellular region|extracellular space", x, ignore.case = TRUE) ~ "Extracellular/Secreted",  # + elastic fiber
    
    # Membrane / Cell surface
    grepl("membrane|cell surface|caveola|plasma membrane|cell cortex|cell periphery|cell leading edge", x, ignore.case = TRUE) ~ "Membrane/Cell Surface",
    
    # Nucleus / Chromatin
    grepl("nucleus|chromatin|chromosome|nucleoplasm|Cajal|fibrillar|Barr body|CENP|nuclear speck|cleavage body|chromatoid body", x, ignore.case = TRUE) ~ "Nucleus/Chromatin",  # + nuclear speck, cleavage body, chromatoid body
    
    # ER / Golgi / Vesicle
    grepl("endoplasmic reticulum|COPII|COPI|Golgi|endomembrane|annulate|endosome|clathrin|amphisome|autophagosome|autolysosome|lysosom|endocytic|early endosome", x, ignore.case = TRUE) ~ "ER/Golgi/Vesicle",  # + COPI, Golgi apparatus, endosome, clathrin, amphisome, autophagosome variants
    
    # Junctions / Adhesion
    grepl("junction|adherens|anchoring|focal adhesion|tight junction|cell-cell|fascia adherens|contact zone", x, ignore.case = TRUE) ~ "Cell Junction/Adhesion",  # + fascia adherens, contact zone
    
    # Granules / Lysosomes
    grepl("granule|aggresome|vesicle|lamellar body|stress granule", x, ignore.case = TRUE) ~ "Granule/Lysosome/Vesicle",  # + stress granule; lysosome/autophagosome moved to ER/Golgi/Vesicle above
    
    # Cilia / Centrosome / Neuron
    grepl("cili|centrosome|centriole|axoneme|basal body|axolemma|axon|dendrite|synapse|GABA|cerebellar|postsynaptic|cleavage furrow", x, ignore.case = TRUE) ~ "Cilia/Centrosome/Neuron",  # + postsynaptic, cleavage furrow
    
    # Named complexes
    grepl("complex|spliceosome|signalosome|ubiquitin|ligase|chaperonin|ATPase|synthase|synthetase|isocitrate|flavoprotein|glycine cleavage|calcineurin|calpain|CD40|angiogenin|AnxA2|MIB|CRD|iron-sulfur|aminoacyl|calmodulin|kinase|Arp2|beta-catenin|Cul[0-9]|aryl hydrocarbon|AP-[0-9]|BRISC|BRCA|branched-chain|acetolactate|chloride channel|commitment|core-binding|Cdc73|DNA replication|snoRNP|bub[0-9]|CAAX|B-WICH|Atg|adaptor", x, ignore.case = TRUE) ~ "Protein Complex",  # + AP-1/2, BRISC, BRCA1-A, branched-chain dehydrogenase, chloride channel, and other new named complexes
    
    # Acrosome
    grepl("acrosom", x, ignore.case = TRUE) ~ "Acrosome",
    
    # Empty
    x == "" ~ "Unknown",
    
    TRUE ~ "Other"
  )
}

norm_str$comp1_grouped <- condense_compartment(norm_str$comp1)
raw$comp1_grouped <- condense_compartment(raw$comp1)
norm$comp1_grouped <- condense_compartment(norm$comp1)
norm_non$comp1_grouped <- condense_compartment(norm_non$comp1)
# Apply the map to comp columns
# Function to map a term to a broad category
map_category <- function(term, comp_map) {
  if(is.na(term)) return(NA)
  matched <- comp_map$category[sapply(comp_map$pattern, function(p) str_detect(term, p))]
  if(length(matched) > 0) return(matched[1]) 
  else return(term)  # keep original term if no match
}

raw <- raw %>%
  mutate(
    comp1_cat = sapply(comp1, map_category, comp_map = comp_map),
    comp2_cat = sapply(comp2, map_category, comp_map = comp_map),
    comp3_cat = sapply(comp3, map_category, comp_map = comp_map)
  )

norm <- norm %>%
  mutate(
    comp1_cat = sapply(comp1, map_category, comp_map = comp_map),
    comp2_cat = sapply(comp2, map_category, comp_map = comp_map),
    comp3_cat = sapply(comp3, map_category, comp_map = comp_map)
  )

norm_non <- norm_non %>%
  mutate(
    comp1_cat = sapply(comp1, map_category, comp_map = comp_map),
    comp2_cat = sapply(comp2, map_category, comp_map = comp_map),
    comp3_cat = sapply(comp3, map_category, comp_map = comp_map)
  )

norm_str <- norm_str %>%
  mutate(
    comp1_cat = sapply(comp1, map_category, comp_map = comp_map),
    comp2_cat = sapply(comp2, map_category, comp_map = comp_map),
    comp3_cat = sapply(comp3, map_category, comp_map = comp_map)
  )

# raw
raw_long <- raw %>%
  pivot_longer(
    cols = starts_with("CR") | starts_with("PPS"),  # all intensity columns
    names_to = "sample_id",
    values_to = "intensity"
  )

raw_long <- raw_long %>%
  left_join(meta, by = "sample_id")  # metadata contains sample → group/tissue info

# standard (moderate)
norm_long <- norm %>%
  pivot_longer(
    cols = starts_with("CR") | starts_with("PPS"),  # all intensity columns
    names_to = "sample_id",
    values_to = "intensity"
  )

norm_long <- norm_long %>%
  left_join(meta, by = "sample_id")

# no filter
norm_non_long <- norm_non %>%
  pivot_longer(
    cols = starts_with("CR") | starts_with("PPS"),  # all intensity columns
    names_to = "sample_id",
    values_to = "intensity"
  )

norm_non_long <- norm_non_long %>%
  left_join(meta, by = "sample_id")

# stringent
norm_str_long <- norm_str %>%
  pivot_longer(
    cols = starts_with("CR") | starts_with("PPS"),  # all intensity columns
    names_to = "sample_id",
    values_to = "intensity"
  )

norm_str_long <- norm_str_long %>%
  left_join(meta, by = "sample_id")

# Create summaries
global_summary <- raw_long %>%
  group_by(comp1_cat) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  mutate(prop = total_intensity / sum(total_intensity))

gp_sum_raw <- raw_long %>%
  group_by(cancer_time, comp1_grouped) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  group_by(cancer_time) %>%
  mutate(prop = total_intensity / sum(total_intensity))

gp_sum_norm_str <- norm_str_long %>%
  group_by(cancer_time, comp1_grouped) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  group_by(cancer_time) %>%
  mutate(prop = total_intensity / sum(total_intensity))

sample_summary <- raw_long %>%
  group_by(sample_id, comp1_cat) %>%
  summarise(total_intensity = sum(intensity, na.rm = TRUE)) %>%
  group_by(sample_id) %>%
  mutate(prop = total_intensity / sum(total_intensity))

highlight_term <- "blood microparticle"  # the category you want to highlight

group_summary_blood <- group_summary %>%
  mutate(color_group = ifelse(comp1_cat == highlight_term, highlight_term, "Other"))

sample_summary_blood <- sample_summary %>%
  mutate(color_group = ifelse(comp1_cat == highlight_term, highlight_term, "Other"))


ggplot(gp_sum_raw, aes(x = cancer_time, y = prop, fill = comp1_grouped)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Uniprot GO-CC Proportions by Group",
       subtitle = "Raw Data",
       x = "Group",
       y = "Proportion of Total Raw Intensity")

ggplot(gp_sum_raw, aes(x = cancer_time, y = total_intensity, fill = comp1_grouped)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Uniprot GO-CC Proportions by Group",
       subtitle = "Raw Data",
       x = "Group",
       y = "Total Raw Intensity")

ggplot(gp_sum_norm_str, aes(x = cancer_time, y = prop, fill = comp1_grouped)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Uniprot GO-CC Proportions by Group",
       subtitle = "Normalized Data - Stringent Filter",
       x = "Group",
       y = "Proportion of Total Normalized Intensity")

ggplot(gp_sum_norm_str, aes(x = cancer_time, y = total_intensity, fill = comp1_grouped)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Uniprot GO-CC Proportions by Group",
       subtitle = "Normalized Data - Stringent Filter",
       x = "Group",
       y = "Total Normalized Intensity")

ggplot(sample_summary_blood, aes(x = sample_id, y = prop, fill = color_group)) +
  geom_col() +
  scale_fill_manual(values = c("blood microparticle" = colors$red, 
                               Other = colors$gray)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  ) +
  labs(fill = "Cellular Compartment")

# --- 2y: Plot correlations of compartments with contam --------------------

# 1. Filter to mitochondria group
mito_long <- norm_str_long %>%
  filter(comp1_grouped == "Mitochondria")

# Summarise to one mean intensity per sample
mito_sample <- mito_long %>%
  group_by(sample_id) %>%
  summarise(
    mean_intensity    = mean(intensity, na.rm = TRUE),
    contam_pct_hiab   = dplyr::first(contam_pct_hiab),
    cancer_time       = dplyr::first(cancer_time),
    .groups = "drop"
  )

# 2. Correlation stat for annotation
cor_res <- cor.test(mito_sample$mean_intensity, mito_sample$contam_pct_hiab,
                    method = "spearman", exact = FALSE)

cor_label <- sprintf("r = %.2f, p = %.3f", cor_res$estimate, cor_res$p.value)

# 3. Plot
ggplot(mito_sample, aes(x = contam_pct_hiab, y = mean_intensity, colour = cancer_time)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = TRUE, aes(fill = cancer_time), alpha = 0.15) +
  scale_colour_manual(values = c(colors$gray, colors$red, colors$blue)) +
  scale_fill_manual(values = c(colors$gray, colors$red, colors$blue)) +
  annotate("text", x = Inf, y = Inf, label = cor_label,
           hjust = 1.1, vjust = 1.8, size = 3.5,
           colour = "grey30") +
  labs(x = "Contamination index (contam_pct_hiab)",
       y = "Mean mitochondrial protein intensity",
       colour = NULL, fill = NULL,
       title = "Mitochondrial intensity ~ contamination") +
  theme_cr()

# --- 2z: Correlation plots by compartment --------------------
# --- 1. Helper: summarise mean intensity per sample per compartment ----------
make_comp_summary <- function(long_df) {
  long_df %>%
    filter(!is.na(comp1_grouped), comp1_grouped != "Unknown", comp1_grouped != "Other") %>%
    group_by(sample_id, comp1_grouped) %>%
    summarise(
      mean_intensity = mean(intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(
      long_df %>%
        distinct(sample_id, contam_pct_hiab, cancer_time),
      by = "sample_id"
    )
}

comp_norm     <- make_comp_summary(norm_long)
comp_none     <- make_comp_summary(norm_non_long)
comp_str      <- make_comp_summary(norm_str_long)
comp_raw      <- make_comp_summary(raw_long)

# --- 2. Helper: make all compartment scatter plots for one pipeline ----------
make_comp_plots <- function(comp_df, subtitle) {
  compartments <- unique(comp_df$comp1_grouped)
  
  purrr::map(compartments, function(comp) {
    dat <- filter(comp_df, comp1_grouped == comp)
    
    ggplot(dat, aes(x = contam_pct_hiab, y = mean_intensity, color = cancer_time)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, colour = "grey30",
                  fill = "grey70", linewidth = 0.8) +
      geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
      scale_colour_manual(values = c(colors$gray, colors$red, colors$blue), name = NULL) +
      ggpubr::stat_cor(method = "spearman",
                       label.x.npc = 0.02, label.y.npc = 0.98, size = 3) +
      labs(
        title    = subtitle,
        x        = "Contamination index (% of raw intensity)",
        y        = "Mean normalized intensity"
      ) +
      theme_cr()
  }) %>%
    purrr::set_names(compartments)
}

plots_norm <- make_comp_plots(comp_norm, "Standard filter")
plots_none <- make_comp_plots(comp_none, "No filter")
plots_str  <- make_comp_plots(comp_str,  "Stringent filter")
plots_raw  <- make_comp_plots(comp_raw,  "Raw data")

# --- 3. Access individual plots ----------------------------------------------
# e.g. plots_str[["Mitochondria"]]
# e.g. plots_norm[["Cytoskeleton/Contractile"]]

# --- 4. Optional: patchwork panel for one compartment across all pipelines ---
compare_pipelines <- function(comp_name) {
  (plots_norm[[comp_name]] | plots_none[[comp_name]]) /
    (plots_str[[comp_name]]  | plots_raw[[comp_name]]) +
    patchwork::plot_annotation(
      title = comp_name,
      subtitle = "Contamination index vs mean intensity across pipeline versions"
    ) 
}

# e.g.
compare_pipelines("Mitochondria")
compare_pipelines("Blood Microparticle")
compare_pipelines("Extracellular/Secreted")
compare_pipelines("Skeletal Muscle")

# --- 2zz: Comparison by proportion correlations ----------------------
# --- 1. Helper: summarise total intensity per sample first -------------------
make_comp_summary_pct <- function(long_df) {
  
  # total normalized intensity per sample (denominator)
  sample_totals <- long_df %>%
    group_by(sample_id) %>%
    summarise(total_intensity = sum(intensity, na.rm = TRUE), .groups = "drop")
  
  # sum intensity per sample per compartment, then divide by total
  long_df %>%
    filter(!is.na(comp1_grouped), comp1_grouped != "Unknown", comp1_grouped != "Other") %>%
    group_by(sample_id, comp1_grouped) %>%
    summarise(
      comp_intensity = sum(intensity, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(sample_totals, by = "sample_id") %>%
    mutate(comp_pct = comp_intensity / total_intensity) %>%
    left_join(
      long_df %>% distinct(sample_id, contam_pct_hiab, cancer_time),
      by = "sample_id"
    )
}

comp_norm <- make_comp_summary_pct(norm_long)
comp_none <- make_comp_summary_pct(norm_non_long)
comp_str  <- make_comp_summary_pct(norm_str_long)
comp_raw  <- make_comp_summary_pct(raw_long)

# --- 2. Helper: make all compartment scatter plots for one pipeline ----------
make_comp_plots <- function(comp_df, subtitle, col_ctl, col_surv_bl, col_surv_t) {
  compartments <- unique(comp_df$comp1_grouped)
  
  purrr::map(compartments, function(comp) {
    dat <- filter(comp_df, comp1_grouped == comp)
    
    ggplot(dat, aes(x = contam_pct_hiab, y = comp_pct, color = cancer_time)) +
      geom_point(size = 2.5, alpha = 0.8) +
      geom_smooth(method = "lm", se = TRUE, colour = "grey30",
                  fill = "grey70", linewidth = 0.8) +
      geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
      scale_colour_manual(
        values = c(col_ctl, col_surv_bl, col_surv_t),
        name   = NULL
      ) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 0.1)) +
      ggpubr::stat_cor(method = "spearman",
                       label.x.npc = 0.05, label.y.npc = 0.95, size = 3) +
      labs(
        title    = comp,
        subtitle = subtitle,
        x        = "Contamination index (%)",
        y        = "Compartment proportion (%)"
      ) +
      theme_cr()
  }) %>%
    purrr::set_names(compartments)
}

# --- 3. Build plot lists — swap in your actual colours ----------------------
plots_norm <- make_comp_plots(comp_norm, "Standard filter",   colors$gray, colors$blue, colors$red)
plots_none <- make_comp_plots(comp_none, "No filter",         colors$gray, colors$blue, colors$red)
plots_str  <- make_comp_plots(comp_str,  "Stringent filter",  colors$gray, colors$blue, colors$red)
plots_raw  <- make_comp_plots(comp_raw,  "Raw (unfiltered)",  colors$gray, colors$blue, colors$red)

# --- 4. 2x2 panel for one compartment across all four pipelines -------------
compare_pipelines <- function(comp_name) {
  (plots_norm[[comp_name]] | plots_none[[comp_name]]) /
    (plots_str[[comp_name]]  | plots_raw[[comp_name]]) +
    patchwork::plot_annotation(
      title    = comp_name,
      subtitle = "Contamination index vs compartment proportion across pipeline versions"
    )
}

# e.g.
compare_pipelines("Mitochondria")
compare_pipelines("Skeletal Muscle")
compare_pipelines("Blood Microparticle")
compare_pipelines("Extracellular/Secreted")
# --- 2f: Additional Correlations and plots ---------------
canonical_cors <- raw %>%
  filter(gene %in% detected_canonical) %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "intensity") %>%
  left_join(contam_index %>% 
              dplyr::select(sample_id, contam_pct), by = "sample_id") %>%
  group_by(gene) %>%
  summarise(
    cor     = cor(intensity, contam_pct,
                  use = "pairwise.complete.obs", method = "spearman"),
    p_value = cor.test(intensity, contam_pct, method = "spearman")$p.value,
    .groups = "drop"
  ) %>%
  arrange(cor)

p_suppression_lollipop <- ggplot(canonical_cors,
                                 aes(x = cor, y = reorder(gene, cor),
                                     color = cor < 0)) +
  geom_segment(aes(x = 0, xend = cor,
                   y = reorder(gene, cor), yend = reorder(gene, cor)),
               linewidth = 0.5) +
  geom_point(size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "red"),
                     labels = c("TRUE" = "Suppressed", "FALSE" = "Enhanced"),
                     name   = "Direction") +
  labs(
    title = "Canonical Muscle Proteins vs Contamination Index",
    x     = "Spearman r",
    y     = NULL
  ) +
  theme_minimal()

# --- 3: Linear Modeling Assessment ------------------------------
lm(canonical_pct ~ contam_pct_hiab.x + cancer_time, data = canon_df_norm) %>%
  summary()

# --- 4: PCA Assessment --------------------------------------
# --- 4a: PCA - Imputed dataset -------------------------------
# Create matrix and run PCA
imp_mat <- imp[, 6:44]
rownames(imp_mat) <- imp$gene
pca_imp     <- prcomp(t(imp_mat), center = TRUE, scale. = TRUE)
pc_scores   <- pca_imp$x[, 1:3]
mahal_dist  <- mahalanobis(pc_scores, colMeans(pc_scores), cov(pc_scores))

# Extract variance
var_exp_imp <- round(summary(pca_imp)$importance[2, 1:2] * 100, 1)

# Create dataframe for plotting
pc_df_imp <- as.data.frame(pca_imp$x[, 1:2])
pc_df_imp$sample_id <- rownames(pc_df)

pc_df_imp <- pc_df_imp %>%
  mutate(
    group = case_when(
      grepl("^PPS", sample_id) ~ "CTL",
      grepl("^CR.*_T1$", sample_id) ~ "SURV_T1",
      grepl("^CR.*_T2$", sample_id) ~ "SURV_T2"
    )
  )
pc_df_imp$group <- factor(pc_df_imp$group,
                      levels = c("CTL", "SURV_T1", "SURV_T2"))

# --- 4b: PCA - stringent dataset (normalized) -------------------------
# Create matrix
data_for_pca_mat <- norm_str %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  mutate(across(-gene, ~ as.numeric(.))) %>%
  column_to_rownames("gene")  # attach gene names before converting to matrix

# Remove rows with all NA or zero variance
data_for_pca_mat <- data_for_pca_mat[complete.cases(data_for_pca_mat), ]
data_for_pca_mat <- data_for_pca_mat[apply(data_for_pca_mat, 1, var, na.rm = TRUE) > 0, ]

# Run PCA on transposed matrix (samples as rows, proteins as columns)
pca_str <- prcomp(t(data_for_pca_mat), scale. = TRUE, center = TRUE)

# Extract explained variance 
var_exp_str <- round(summary(pca_str)$importance[2, 1:2] * 100, 1)

# Create dataframe for plotting
pc_df_str <- as.data.frame(pca_str$x[, 1:2])
pc_df_str$sample_id <- rownames(pc_df_str)

pc_df_str <- pc_df_str %>%
  mutate(
    group = case_when(
      grepl("^PPS", sample_id) ~ "CTL",
      grepl("^CR.*_T1$", sample_id) ~ "SURV_T1",
      grepl("^CR.*_T2$", sample_id) ~ "SURV_T2"
    )
  )
pc_df_str$group <- factor(pc_df_str$group,
                          levels = c("CTL", "SURV_T1", "SURV_T2"))


# --- 4c: PCA - original analysis (normalized) ---------------------------
# Create matrix
data_for_pca_norm <- norm %>%
  dplyr::select(gene, all_of(sample_cols)) %>%
  mutate(across(-gene, ~ as.numeric(.))) %>%
  column_to_rownames("gene")  # attach gene names before converting to matrix

# Remove rows with all NA or zero variance
data_for_pca_norm <- data_for_pca_norm[complete.cases(data_for_pca_norm), ]
data_for_pca_norm <- data_for_pca_norm[apply(data_for_pca_norm, 1, var, na.rm = TRUE) > 0, ]

# Run PCA on transposed matrix (samples as rows, proteins as columns)
pca_norm <- prcomp(t(data_for_pca_norm), scale. = TRUE, center = TRUE)

# Extract variance
var_exp_norm <- round(summary(pca_norm)$importance[2, 1:2] * 100, 1)

pc_df_norm <- as.data.frame(pca_norm$x[, 1:2])
pc_df_norm$sample_id <- rownames(pc_df_norm)

pc_df_norm <- pc_df_norm %>%
  mutate(
    group = case_when(
      grepl("^PPS", sample_id) ~ "CTL",
      grepl("^CR.*_T1$", sample_id) ~ "SURV_T1",
      grepl("^CR.*_T2$", sample_id) ~ "SURV_T2"
    )
  )
pc_df_norm$group <- factor(pc_df_norm$group,
                          levels = c("CTL", "SURV_T1", "SURV_T2"))

# --- 4d: PCA Plotting ----------------------------------------------
# Imputed
p_pca_imp <- ggplot(pc_df_imp, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "CTL" = "gray40",
    "SURV_T1" = "#377EB8",
    "SURV_T2" = "#E41A1C"
  )) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp_imp[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp_imp[2]),
       title = "PCA - Standard Blood Filter",
       subtitle = "Imputed Data",
       color = "Group") +
  stat_ellipse(
    geom  = "polygon",
    alpha = 0.1,
    level = 0.95,
    type  = "t"
  ) +
  theme_minimal(base_size = 11)

# Stringent
p_pca_str <- ggplot(pc_df_str, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "CTL" = "gray40",
    "SURV_T1" = "#377EB8",
    "SURV_T2" = "#E41A1C"
  )) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp_str[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp_str[2]),
       title = "PCA - Stringent Blood Filter",
       subtitle = "Normalized Data",
       color = "Group") +
  stat_ellipse(
    geom  = "polygon",
    alpha = 0.1,
    level = 0.95,
    type  = "t"
  ) +
  theme_minimal(base_size = 11)

# Original
p_pca_norm <- ggplot(pc_df_norm, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  scale_color_manual(values = c(
    "CTL" = "gray40",
    "SURV_T1" = "#377EB8",
    "SURV_T2" = "#E41A1C"
  )) +
  geom_text_repel(aes(label = sample_id), size = 2.5, max.overlaps = 10) +
  labs(x = sprintf("PC1 (%.1f%%)", var_exp_norm[1]),
       y = sprintf("PC2 (%.1f%%)", var_exp_norm[2]),
       title = "PCA - Standard Blood Filter",
       subtitle = "Normalized Data - Original Analysis",
       color = "Group") +
  stat_ellipse(
    geom  = "polygon",
    alpha = 0.1,
    level = 0.95,
    type  = "t"
  ) +
  theme_minimal(base_size = 11)

# --- 5: Barplot of Contamination by Group ---------------------------
# Sort data hi-low
cont_df <- meta %>%
  arrange(desc(contam_pct_hiab)) %>%
  mutate(sample_id = factor(sample_id, levels = sample_id))

# plot
p_bar_contam <- ggplot(cont_df, 
             aes(x = sample_id, y = contam_pct_hiab, fill = cancer_time)) +
  geom_col(color = "black") +
  scale_fill_manual(values = c(
    "CTL_T1" = "white",
    "SURV_T1" = "grey",
    "SURV_T2" = "grey2"
  )) +
  # geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", linewidth = 1) +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# --- 6: PCA Correlation with Contamination Index -------------------------
# Check PCA corr with contamination index
# Stringent filtering
pca_cor_str <- data.frame(
  sample_id = rownames(pca_str$x),
  PC1 = pca_str$x[,1],
  PC2 = pca_str$x[,2]
) %>%
  left_join(meta, by = "sample_id")

cor.test(pca_cor_str$PC1, pca_cor_str$contam_pct_hiab)

# Standard - original analysis
pca_cor_norm <- data.frame(
  sample_id = rownames(pca_norm$x),
  PC1 = pca_norm$x[,1],
  PC2 = pca_norm$x[,2]
) %>%
  left_join(meta, by = "sample_id")

cor.test(pca_cor_norm$PC1, pca_cor_norm$contam_pct_hiab)

# --- 7: Assess PCA_1 top contributors ---------------------------------
# Extract loadings
loadings <- pca_str$rotation
head(rownames(loadings))  # should show gene names 

# Top PC1 contributors
pc1_loadings <- loadings[, "PC1"]
top_pc1_df <- data.frame(
  gene    = names(sort(abs(pc1_loadings), decreasing = TRUE)[1:20]),
  loading = pc1_loadings[names(sort(abs(pc1_loadings), decreasing = TRUE)[1:20])]
) %>%
  arrange(desc(abs(loading)))



# --- 8: NORMALIZED INTENSITIES ----------------------------
sample_cols <- colnames(norm[, 7:45])

# Total intensity per sample (calculate separately first)
total_intensity_n <- norm %>%
  dplyr::select(all_of(sample_cols)) %>%
  colSums(na.rm = TRUE)

mb_intensity_n <- norm %>%
  filter(gene == "MB") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

ckm_intensity_n <- norm %>%
  filter(gene == "CKM") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh7_intensity_n <- norm %>%
  filter(gene == "MYH7") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh1_intensity_n <- norm %>%
  filter(gene == "MYH1") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

myh2_intensity_n <- norm %>%
  filter(gene == "MYH2") %>%
  dplyr::select(all_of(sample_cols)) %>%
  unlist(use.names = TRUE)

int_index_n <- data.frame(
  sample_id = names(ckm_intensity_n),
  total_intensity_n = total_intensity_n,
  mb_intensity_n = mb_intensity_n,
  ckm_intensity_n = ckm_intensity_n,
  myh7_intensity_n = myh7_intensity_n,
  myh1_intensity_n = myh1_intensity_n,
  myh2_intensity_n = myh2_intensity_n
)

int_index_n <- int_index_n %>%
  mutate(
    mb_pct = mb_intensity_n / total_intensity_n,
    ckm_pct = ckm_intensity_n / total_intensity_n,
    myh7_pct = myh7_intensity_n / total_intensity_n,
    myh1_pct = myh1_intensity_n / total_intensity_n,
    myh2_pct = myh2_intensity_n / total_intensity_n
  )

int_index_n <- int_index_n %>%
  mutate(group = case_when(
    grepl("_T1$", sample_id) ~ "SURV T1",
    grepl("_T2$", sample_id) ~ "SURV T2",
    grepl("^PPS", sample_id) ~ "CTL"
  ))

int_index_n_long <- int_index_n %>%
  dplyr::select(sample_id, group,
                mb_pct, ckm_pct, myh7_pct, myh1_pct, myh2_pct) %>%
  pivot_longer(
    cols = -c(sample_id, group),
    names_to = "pro_pct",
    values_to = "value"
  )

int_index_n_long <- int_index_n_long %>%
  mutate(protein = case_when(
    grepl("mb_pct", pro_pct) ~ "MB",
    grepl("ckm_pct", pro_pct) ~ "CKM",
    grepl("myh7_pct", pro_pct) ~ "MYH7",
    grepl("myh1_pct", pro_pct) ~ "MYH1",
    grepl("myh2_pct", pro_pct) ~ "MYH2"
  ))

ggplot(int_index_n_long, aes(x = group, y = value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
  facet_wrap(~protein, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Normalized Data",
    subtitle = "Standard Blood Filter",
    x = NULL, y = "Protein % of total intensity") +
  scale_fill_manual(values = c(
    "SURV T1" = "#1b9e77",
    "SURV T2" = "#d95f02",
    "CTL" = "#7570b3"
  )) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  )

# --- Correlation Plot (Kruger et al.) -----------------------------------
raw2 <- raw
raw2$gene <- make.unique(raw2$gene)
rownames(raw2) <- raw2$gene

cor_mat <- raw2[ , setdiff(names(raw2), annot_cols)]
dim(cor_mat)
alb <- cor_mat["ALB", ]
length(alb)
mb <- cor_mat["MB", ]

ALB_cor <- apply(cor_mat, 1, function(x) {
  cor(x, alb, use = "pairwise.complete.obs")
})

MB_cor <- apply(cor_mat, 1, function(x) {
  cor(x, mb, use = "pairwise.complete.obs")
})

alb <- as.numeric(cor_mat["ALB", ])
mb  <- as.numeric(cor_mat["MB", ])

# Ensure ordering matches
raw2 <- raw2 %>%
  mutate(
    ALB_cor = ALB_cor[gene],
    MB_cor  = MB_cor[gene]
  )

# Plot

top_proteins <- raw2 %>% 
  arrange(desc(ALB_cor)) %>% 
  slice_head(n = 10) %>% 
  pull(gene)

raw2 <- raw2 %>%
  mutate(
    cor_group = case_when(
      ALB_cor > 0 & MB_cor < 0 ~ "ALB+ / MB-",
      ALB_cor < 0 & MB_cor > 0 ~ "MB+ / ALB-",
      TRUE                     ~ "Same direction"  # both positive or both negative
    )
  )

ggplot(raw2, aes(x = ALB_cor, y = MB_cor, color = cor_group)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c(
    "ALB+ / MB-" = "red",
    "MB+ / ALB-" = "blue",
    "Same direction" = "gray"
  )) +
  theme_minimal() +
  labs(
    x = "Correlation with Albumin (ALB)",
    y = "Correlation with Myoglobin (MB)",
    color = "Correlation type"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")

muscle_gene_list <- raw2 %>%
  filter(ALB_cor < 0 & MB_cor > 0) %>%
  pull(gene)

