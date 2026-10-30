# R21 full cohort -- mixed-model ANOVA on phenotypic outcomes (standalone)
# Model per outcome: value ~ group * time + (1 | subject), Kenward-Roger df.
# All available observations are used (non-completers contribute their PRE values).
#
# Outputs (config$out_dir):
#   fig_prepost_by_group.pdf/.png -- model means +/- 95% CI by group, individuals joined
#   fig_prepost_pooled.pdf/.png   -- pooled (CRE + PLA) means +/- 95% CI, Time p bracket
#   R21_full_anova.xlsx           -- main table, model means, contrasts, F tests, data

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readr, dplyr, tidyr, purrr, ggplot2, lme4, lmerTest, pbkrtest, emmeans, writexl)

config <- list(
  data_file = "00_input/R21/data_R21_full.csv",   # <-- CHECK path
  out_dir   = "04_figures/fig_01_pheno_full",
  group_map = c(`1` = "PLA", `0` = "CRE")         # <-- CHECK: R21 group code -> label
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

group_colors <- c(PLA = "#1F4E79", CRE = "#E69F00")
pooled_color <- "grey45"

# R21 column stem, display label, unit conversion
# (R21 units: weight lb; lbm/alm/fat g -> converted to kg)
vars <- tribble(
  ~source,         ~label,                         ~mult,
  "weight",        "Weight (kg)",                  0.45359237,
  "BMI",           "BMI",                          1,
  "lbm",           "Lean Body Mass (kg)",          1 / 1000,
  "alm",           "Appendicular Lean Mass (kg)",  1 / 1000,
  "alm_height",    "ALM Index",                    1,
  "fat",           "Fat Mass (kg)",                1 / 1000,
  "bf",            "Fat (%)",                      1,
  "sts_max",       "Sit-to-Stand Max Power",       1,
  "grip",          "Grip Strength (lbs)",          1,
  "chest_press",   "Chest Press (lbs)",            1,
  "leg_extension", "Leg Extension (lbs)",          1
)

# --- 1: Load and reshape ---------------------------------------------------------
raw <- read_csv(config$data_file, show_col_types = FALSE)

# One row per subject x outcome
wide <- map_dfr(seq_len(nrow(vars)), function(i) {
  tibble(
    subject  = raw$id,
    group    = unname(config$group_map[as.character(raw$group)]),
    variable = vars$label[i],
    PRE      = as.numeric(raw[[paste0("pre_",  vars$source[i])]]) * vars$mult[i],
    POST     = as.numeric(raw[[paste0("post_", vars$source[i])]]) * vars$mult[i]
  )
}) %>%
  mutate(group    = factor(group, levels = c("PLA", "CRE")),
         variable = factor(variable, levels = vars$label),
         CHG      = POST - PRE,
         PCT      = if_else(PRE > 0, 100 * CHG / PRE, NA_real_))

# One row per subject x outcome x timepoint (what the model is fit to)
long <- wide %>%
  pivot_longer(c(PRE, POST), names_to = "time", values_to = "value") %>%
  filter(!is.na(value)) %>%
  mutate(time = factor(time, levels = c("PRE", "POST")))

cat("Participants by group:\n"); print(table(distinct(wide, subject, group)$group))

# --- 2: Mixed-model ANOVA per outcome --------------------------------------------
emm_options(msg.interaction = FALSE)
KR <- "kenward-roger"

tidy_emm <- function(x) {
  as.data.frame(summary(x, infer = TRUE)) %>% rename(any_of(c(estimate = "emmean")))
}
std <- function(x, label) {
  tidy_emm(x) %>% transmute(effect = label, estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value)
}

run_one <- function(v) {
  d <- filter(long, variable == v)
  m <- lmer(value ~ group * time + (1 | subject), data = d)
  
  emm_gt <- emmeans(m, ~ group * time, lmer.df = KR)
  emm_t  <- emmeans(m, ~ time,         lmer.df = KR)
  emm_g  <- emmeans(m, ~ group,        lmer.df = KR)
  
  by_grp <- tidy_emm(contrast(emm_gt, "revpairwise", by = "group", adjust = "none"))
  a      <- anova(m, ddf = "Kenward-Roger")
  
  list(
    means      = tidy_emm(emm_gt),
    time_means = tidy_emm(emm_t),
    contrasts  = bind_rows(
      by_grp %>% transmute(effect = paste("Change", group), estimate, SE, df,
                           lower.CL, upper.CL, t.ratio, p.value),
      std(contrast(emm_t,  "revpairwise", adjust = "none"),            "Pooled change"),
      std(contrast(emm_gt, interaction = "revpairwise", adjust = "none"), "Diff in change (CRE - PLA)"),
      std(contrast(emm_g,  "revpairwise", adjust = "none"),            "Group main (CRE - PLA)"),
      tidy_emm(contrast(emm_gt, "revpairwise", by = "time", adjust = "none")) %>%
        filter(time == "PRE") %>%
        transmute(effect = "Baseline diff (CRE - PLA)", estimate, SE, df,
                  lower.CL, upper.CL, t.ratio, p.value)
    ),
    ftests = tibble(term = rownames(a), F_value = a$`F value`, df_num = a$NumDF,
                    df_den = a$DenDF, p.value = a$`Pr(>F)`, singular_fit = isSingular(m))
  )
}

results <- map(set_names(vars$label), run_one)

stack <- function(part) {
  imap_dfr(results, ~ as_tibble(.x[[part]]) %>% mutate(variable = .y, .before = 1)) %>%
    mutate(variable = factor(variable, levels = vars$label))
}
model_means  <- stack("means")
pooled_means <- stack("time_means")
contrasts    <- stack("contrasts")
ftests       <- stack("ftests")

# --- 3: Main table ---------------------------------------------------------------
fmt_p    <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
fmt_plab <- function(p) ifelse(!is.na(p) & p < 0.001, "p < 0.001", paste0("p = ", fmt_p(p)))
fnum     <- function(x, d) sprintf("%.*f", d, x)

digits <- wide %>% group_by(variable) %>%
  summarise(d = if_else(max(abs(PRE), na.rm = TRUE) >= 20, 1L, 2L), .groups = "drop")

n_tbl <- wide %>% group_by(variable) %>%
  summarise(
    n        = paste0(sum(group == "CRE" & (!is.na(PRE) | !is.na(POST))), " / ",
                      sum(group == "PLA" & (!is.na(PRE) | !is.na(POST)))),
    n_paired = paste0(sum(group == "CRE" & !is.na(CHG)), " / ",
                      sum(group == "PLA" & !is.na(CHG))),
    .groups = "drop")

obs_tbl <- long %>%
  group_by(variable, group, time) %>%
  summarise(m = mean(value), s = sd(value), .groups = "drop") %>%
  left_join(digits, by = "variable") %>%
  mutate(col = paste(group, time), txt = paste0(fnum(m, d), " \u00b1 ", fnum(s, d))) %>%
  dplyr::select(variable, col, txt) %>%
  pivot_wider(names_from = col, values_from = txt)

pct_tbl <- bind_rows(wide %>% mutate(group = as.character(group)),
                     wide %>% mutate(group = "Pooled")) %>%
  filter(!is.na(PCT)) %>%
  group_by(variable, group) %>%
  summarise(txt = sprintf("%.1f \u00b1 %.1f", mean(PCT), sd(PCT)), .groups = "drop") %>%
  mutate(col = paste("pct", group)) %>%
  dplyr::select(variable, col, txt) %>%
  pivot_wider(names_from = col, values_from = txt)

ct_tbl <- contrasts %>%
  left_join(digits, by = "variable") %>%
  mutate(ci = paste0(fnum(estimate, d), " [", fnum(lower.CL, d), ", ", fnum(upper.CL, d), "]"),
         p  = fmt_p(p.value)) %>%
  dplyr::select(variable, effect, ci, p) %>%
  pivot_wider(names_from = effect, values_from = c(ci, p), names_sep = "|")

Main_Table <- tibble(variable = factor(vars$label, levels = vars$label)) %>%
  left_join(n_tbl,   by = "variable") %>%
  left_join(obs_tbl, by = "variable") %>%
  left_join(pct_tbl, by = "variable") %>%
  left_join(ct_tbl,  by = "variable") %>%
  transmute(
    Outcome                                         = as.character(variable),
    `n (CRE / PLA)`                                 = n,
    `n paired (CRE / PLA)`                          = n_paired,
    `CRE PRE (mean +/- SD)`                         = `CRE PRE`,
    `CRE POST (mean +/- SD)`                        = `CRE POST`,
    `PLA PRE (mean +/- SD)`                         = `PLA PRE`,
    `PLA POST (mean +/- SD)`                        = `PLA POST`,
    `Baseline diff, CRE - PLA [95% CI]`             = `ci|Baseline diff (CRE - PLA)`,
    `Change CRE [95% CI]`                           = `ci|Change CRE`,
    `Change PLA [95% CI]`                           = `ci|Change PLA`,
    `% change CRE (mean +/- SD)`                    = `pct CRE`,
    `% change PLA (mean +/- SD)`                    = `pct PLA`,
    `% change Pooled (mean +/- SD)`                 = `pct Pooled`,
    `Diff in change, CRE - PLA (Group x Time) [95% CI]` = `ci|Diff in change (CRE - PLA)`,
    `Group x Time p`                                = `p|Diff in change (CRE - PLA)`,
    `Pooled change [95% CI]`                        = `ci|Pooled change`,
    `Time p`                                        = `p|Pooled change`,
    `Group main, CRE - PLA [95% CI]`                = `ci|Group main (CRE - PLA)`,
    `Group p`                                       = `p|Group main (CRE - PLA)`
  )

print(as.data.frame(Main_Table))

# --- 4: Figures --------------------------------------------------------------------
fig_theme <- theme_bw(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey92", color = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        plot.title.position = "plot")

n_cols      <- ceiling(sqrt(nrow(vars)))
plot_width  <- max(8, 3.2 * n_cols)
plot_height <- max(5, 3.0 * ceiling(nrow(vars) / n_cols))

save_fig <- function(p, name) {
  ggsave(file.path(config$out_dir, paste0(name, ".pdf")), p,
         width = plot_width, height = plot_height, device = cairo_pdf)
  ggsave(file.path(config$out_dir, paste0(name, ".png")), p,
         width = plot_width, height = plot_height, dpi = 300, bg = "white")
}

set.seed(1)
subj_jit <- distinct(long, subject) %>% mutate(jit = runif(n(), -0.12, 0.12))

p_lab <- ftests %>%
  dplyr::select(variable, term, p.value) %>%
  pivot_wider(names_from = term, values_from = p.value)

# 4a: PRE/POST by group
pos <- tibble(group = factor(c("PLA", "PLA", "CRE", "CRE"), levels = c("PLA", "CRE")),
              time  = factor(c("PRE", "POST", "PRE", "POST"), levels = c("PRE", "POST")),
              x     = c(1, 1.8, 2.8, 3.6))

bars_g <- model_means %>% left_join(pos, by = c("group", "time"))
pts_g  <- long %>% left_join(pos, by = c("group", "time")) %>%
  left_join(subj_jit, by = "subject") %>% mutate(x = x + jit)
ann_g  <- p_lab %>%
  mutate(label = paste0("Group ", fmt_plab(group), "\nTime ", fmt_plab(time),
                        "\nGroup\u00d7Time ", fmt_plab(`group:time`)))

p_group <- ggplot(bars_g, aes(x = x, y = estimate)) +
  geom_col(aes(fill = group, alpha = time), width = 0.7, color = "black") +
  geom_line(data = pts_g, aes(y = value, group = subject), color = "grey30", alpha = 0.35, linewidth = 0.3) +
  geom_point(data = pts_g, aes(y = value), shape = 21, size = 1.4, fill = "white", color = "grey30") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.5) +
  geom_text(data = ann_g, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
            hjust = -0.05, vjust = 1.1, size = 2.3, lineheight = 0.95) +
  scale_fill_manual(values = group_colors, name = "Group") +
  scale_alpha_manual(values = c(PRE = 0.4, POST = 1), name = "Time") +
  scale_x_continuous(breaks = pos$x, labels = paste0(pos$time, "\n", pos$group)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.40))) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "R21 full cohort: CRE vs PLA (mixed-model ANOVA)",
       subtitle = "Bars = model-estimated means \u00b1 95% CI; points/lines = individual participants",
       x = NULL, y = NULL) +
  fig_theme

# 4b: Pooled PRE/POST with Time p bracket
bars_p <- pooled_means %>% mutate(x = if_else(time == "PRE", 1, 2))
pts_p  <- long %>% left_join(subj_jit, by = "subject") %>%
  mutate(x = if_else(time == "PRE", 1, 2) + jit)
brk <- bind_rows(transmute(bars_p, variable, y = upper.CL), transmute(pts_p, variable, y = value)) %>%
  group_by(variable) %>% summarise(y_top = max(y), .groups = "drop") %>%
  left_join(p_lab, by = "variable") %>%
  mutate(y_tick = y_top * 1.03, y_br = y_top * 1.06, y_lab = y_top * 1.08,
         label = fmt_plab(time),
         face  = if_else(time < 0.05, "bold", "plain"))

p_pooled <- ggplot(bars_p, aes(x = x, y = estimate)) +
  geom_col(aes(alpha = time), fill = pooled_color, width = 0.6, color = "black") +
  geom_line(data = pts_p, aes(y = value, group = subject), color = "grey30", alpha = 0.35, linewidth = 0.3) +
  geom_point(data = pts_p, aes(y = value), shape = 21, size = 1.4, fill = "white", color = "grey30") +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.15, linewidth = 0.5) +
  geom_segment(data = brk, aes(x = 1, xend = 2, y = y_br,   yend = y_br),   inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = brk, aes(x = 1, xend = 1, y = y_tick, yend = y_br),   inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = brk, aes(x = 2, xend = 2, y = y_tick, yend = y_br),   inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = brk, aes(x = 1.5, y = y_lab, label = label, fontface = face),
            inherit.aes = FALSE, vjust = 0, size = 2.5) +
  scale_alpha_manual(values = c(PRE = 0.4, POST = 1), guide = "none") +
  scale_x_continuous(breaks = c(1, 2), labels = c("PRE", "POST")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  facet_wrap(~ variable, scales = "free_y") +
  labs(title = "R21 full cohort: pooled CRE + PLA (mixed-model ANOVA)",
       subtitle = "Bars = model-estimated means \u00b1 95% CI, averaged over groups; points/lines = individuals; p = Time effect",
       x = NULL, y = NULL) +
  fig_theme

save_fig(p_group,  "fig_prepost_by_group")
save_fig(p_pooled, "fig_prepost_pooled")

# --- 5: Excel ------------------------------------------------------------------------
notes <- tibble(Note = c(
  "Model: value ~ group * time + (1 | subject); Kenward-Roger df; all available observations.",
  "Change and Pooled change are POST - PRE. Pooled change averages CRE and PLA with equal weight (Time effect).",
  "Diff in change (Group x Time) is the supplement effect. Group main averages over PRE and POST and is not a supplement effect.",
  "% change = mean +/- SD of individual 100 x (POST - PRE) / PRE, participants with both timepoints. Descriptive only.",
  "R21 units converted: weight lb -> kg; lean, appendicular lean and fat mass g -> kg.",
  "p-values uncorrected (pre-specified outcome battery)."
))

write_xlsx(
  list(
    Main_Table   = Main_Table,
    Notes        = notes,
    Model_Means  = model_means,
    Pooled_Means = pooled_means,
    Contrasts    = contrasts,
    F_Tests      = ftests,
    Data         = wide
  ),
  path = file.path(config$out_dir, "R21_full_anova.xlsx")
)

cat("\nSaved to:", config$out_dir, "\n")
