# ==============================
# Diversity indices & H′ boxplot
# ==============================

# Clean workspace (optional)
remove(list = ls())
# If you use renv, uncomment:
# renv::restore()

library(tidyverse)
library(janitor)
library(readr)
library(vegan)
library(ggpubr)

# -------- Helper to prepare one period --------
prep_period <- function(url, years_keep, period_label) {
  readr::read_csv(url, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      fish_species = stringr::str_squish(tolower(fish_species)), # standardize names
      abundance = 1
    ) %>%
    filter(
      !is.na(location_id), !is.na(sampling_year), !is.na(fish_species),
      sampling_year %in% years_keep,
      location_id %in% paste0("M", 1:9)  # keep M1–M9 (adjust if needed)
    ) %>%
    # counts per site × year × species
    group_by(location_id, sampling_year, fish_species) %>%
    summarise(n = sum(abundance), .groups = "drop") %>%
    # wide species matrix per site-year
    pivot_wider(names_from = fish_species, values_from = n, values_fill = 0) %>%
    mutate(period = period_label) %>%
    relocate(period, location_id, sampling_year)
}

# -------- Load both periods --------
url_past <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv"
url_curr <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv"

past_wide  <- prep_period(url_past, years_keep = c(2013, 2014, 2016), period_label = "Past (2013–2016)")
curr_wide  <- prep_period(url_curr, years_keep = c(2021, 2022),        period_label = "Current (2021–2022)")

# species columns (exclude id/meta)
species_cols <- function(df) setdiff(names(df), c("period","location_id","sampling_year"))

all_species <- sort(unique(c(species_cols(past_wide), species_cols(curr_wide))))

# helper: add any missing species columns as zeros and order columns
add_missing_species <- function(df, all_species) {
  miss <- setdiff(all_species, names(df))
  if (length(miss) > 0) {
    df[miss] <- 0
  }
  df %>% select(period, location_id, sampling_year, all_of(all_species))
}

past_wide <- add_missing_species(past_wide, all_species)
curr_wide <- add_missing_species(curr_wide, all_species)

# combine rows
div_mat <- dplyr::bind_rows(past_wide, curr_wide)

# -------- Compute diversity indices per site-year --------
species_mat <- div_mat %>%
  select(all_of(all_species)) %>%
  as.data.frame()

# IMPORTANT: replace NAs with zeros for diversity calcs
species_mat[is.na(species_mat)] <- 0

H  <- vegan::diversity(species_mat, index = "shannon")  # ln base
S  <- vegan::specnumber(species_mat)
J  <- ifelse(S > 0, H / log(S), NA_real_)
EN <- exp(H)

div_indices <- div_mat %>%
  select(period, location_id, sampling_year) %>%
  mutate(
    shannon_H   = H,
    richness_S  = S,
    pielou_J    = J,
    effective_N = EN
  )

# --- Clean for plotting/tests: keep finite H and both periods present
div_indices_clean <- div_indices %>%
  filter(is.finite(shannon_H)) %>%
  mutate(period = factor(period)) %>%
  droplevels()

# ----- Focus on sites M4, M7, M9 and tidy factors -----
sites_of_interest <- c("M4","M7","M9")

div_m <- div_indices_clean %>%
  filter(location_id %in% sites_of_interest) %>%
  mutate(
    period = factor(period, levels = c("Past (2013–2016)", "Current (2021–2022)")),
    location_id = factor(location_id, levels = sites_of_interest)
  )

# Quick sanity check: counts per site × period
print(div_m %>% count(location_id, period))

# ----- Per-site Wilcoxon (Past vs Current) with Bonferroni adjust -----
library(rstatix)

pvals <- div_m %>%
  group_by(location_id) %>%
  wilcox_test(shannon_H ~ period) %>%
  adjust_pvalue(method = "bonferroni") %>%
  mutate(p.label = ifelse(p.adj < .001, "***",
                          ifelse(p.adj < .01, "**",
                                 ifelse(p.adj < .05, "*", "ns"))),
         group1 = "Past (2013–2016)",
         group2 = "Current (2021–2022)") %>%
  ungroup()

# y-positions for p-value labels, per facet
ypos <- div_m %>%
  group_by(location_id) %>%
  summarise(y.position = max(shannon_H, na.rm = TRUE) * 1.05, .groups = "drop")

pvals <- pvals %>%
  left_join(ypos, by = "location_id")

# ----- Boxplot (year-level H′ as observations), faceted by site -----
# Ensure factors are ordered as requested
div_m <- div_m %>%
  mutate(
    location_id = factor(location_id, levels = c("M4","M7","M9")),
    period = factor(period, levels = c("Past (2013–2016)", "Current (2021–2022)"))
  )

# Dodging helpers for neat separation of the two periods at each site
pd  <- position_dodge(width = 0.65)
pjd <- position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65)

ggplot(div_m, aes(x = location_id, y = shannon_H, fill = period)) +
  geom_boxplot(position = pd, width = 0.55, alpha = 0.9,
               outlier.shape = NA, color = "black") +
  geom_jitter(position = pjd, size = 2, alpha = 0.6) +
  labs(
    title = "Shannon–Wiener Diversity (H′)",
    x = "Sampling Site",
    y = "H′ (Shannon–Wiener)",
    fill = "Sampling Period"
  ) +
  # optional custom colors; remove this block if you prefer defaults
  scale_fill_manual(values = c("Past (2013–2016)" = "#A6CEE3",
                               "Current (2021–2022)" = "#1F78B4")) +
  scale_x_discrete(drop = FALSE, limits = c("M4","M7","M9"),
                   labels = c("M4","M7","M9")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )
###################################################################
# ==============================
# Diversity indices & H′ boxplot
# ==============================

remove(list = ls())

library(tidyverse)
library(janitor)
library(readr)
library(vegan)
library(ggpubr)
library(rstatix)

# -------- Helper to prepare one period --------
prep_period <- function(url, years_keep, period_label) {
  readr::read_csv(url, show_col_types = FALSE) %>%
    clean_names() %>%
    mutate(
      fish_species = stringr::str_squish(tolower(fish_species)),
      abundance = 1
    ) %>%
    filter(
      !is.na(location_id), !is.na(sampling_year), !is.na(fish_species),
      sampling_year %in% years_keep,
      location_id %in% paste0("M", 1:9)
    ) %>%
    group_by(location_id, sampling_year, fish_species) %>%
    summarise(n = sum(abundance), .groups = "drop") %>%
    pivot_wider(names_from = fish_species, values_from = n, values_fill = 0) %>%
    mutate(period = period_label) %>%
    relocate(period, location_id, sampling_year)
}

# -------- Load both periods --------
url_past <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv"
url_curr <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv"

past_wide  <- prep_period(url_past, years_keep = c(2013, 2014, 2016), period_label = "Past (2013–2016)")
curr_wide  <- prep_period(url_curr, years_keep = c(2021, 2022),        period_label = "Current (2021–2022)")

# Species columns and combine
species_cols <- function(df) setdiff(names(df), c("period","location_id","sampling_year"))
all_species <- sort(unique(c(species_cols(past_wide), species_cols(curr_wide))))

add_missing_species <- function(df, all_species) {
  miss <- setdiff(all_species, names(df))
  if (length(miss) > 0) df[miss] <- 0
  df %>% select(period, location_id, sampling_year, all_of(all_species))
}

past_wide <- add_missing_species(past_wide, all_species)
curr_wide <- add_missing_species(curr_wide, all_species)
div_mat   <- bind_rows(past_wide, curr_wide)

# -------- Compute diversity indices per site-year --------
species_mat <- div_mat %>% select(all_of(all_species)) %>% as.data.frame()
species_mat[is.na(species_mat)] <- 0

H  <- vegan::diversity(species_mat, index = "shannon")
S  <- vegan::specnumber(species_mat)
J  <- ifelse(S > 0, H / log(S), NA_real_)
EN <- exp(H)

div_indices <- div_mat %>%
  select(period, location_id, sampling_year) %>%
  mutate(
    shannon_H   = H,
    richness_S  = S,
    pielou_J    = J,
    effective_N = EN
  ) %>%
  filter(is.finite(shannon_H)) %>%
  mutate(period = factor(period)) %>%
  droplevels()

# ----- Focus on sites M4, M7, M9 and tidy factors -----
sites_of_interest <- c("M4","M7","M9")
div_m <- div_indices %>%
  filter(location_id %in% sites_of_interest) %>%
  mutate(
    period = factor(period, levels = c("Past (2013–2016)", "Current (2021–2022)")),
    location_id = factor(location_id, levels = sites_of_interest)
  )

# -------- Global significance test (Past vs Current) --------
# Non-parametric Wilcoxon rank-sum (Mann–Whitney)
global_test <- wilcox.test(shannon_H ~ period, data = div_m, exact = FALSE)
p_glob <- global_test$p.value
test_label <- sprintf("Wilcoxon rank-sum: p = %.3f", p_glob)

# For positioning the annotation
y_top <- max(div_m$shannon_H, na.rm = TRUE)
# -------- Global significance test (Past vs Current) --------
# Non-parametric Wilcoxon rank-sum (Mann–Whitney)

library(rstatix)

# Tidy Wilcoxon test between periods for Shannon H'
global_test <- div_m %>%
  wilcox_test(shannon_H ~ period, detailed = TRUE)

global_test
# This gives a tibble with:
# .y., group1, group2, n1, n2, statistic, p, method, alternative

# Extract statistic and p-value
W_glob <- global_test$statistic[1]
p_glob <- global_test$p[1]

# Optional: round nicely for reporting
W_glob_round <- round(W_glob, 2)
p_glob_round <- signif(p_glob, 3)

# Label to use in the plot (or in text)
test_label <- sprintf("Wilcoxon rank-sum: W = %.2f, p = %.3f",
                      W_glob_round, p_glob_round)

# For positioning the annotation in the plot
y_top <- max(div_m$shannon_H, na.rm = TRUE)

# You can check the test results in the console:
print(global_test)
print(test_label)

# -------- Plot --------
pd  <- position_dodge(width = 0.65)
pjd <- position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65)

ggplot(div_m, aes(x = location_id, y = shannon_H, fill = period)) +
  geom_boxplot(position = pd, width = 0.55, alpha = 0.9,
               outlier.shape = NA, color = "black") +
  geom_jitter(position = pjd, size = 2, alpha = 0.6) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "H′ (Shannon–Wiener)",
    fill = "Sampling Period"
  ) +
  scale_fill_manual(values = c("Past (2013–2016)" = "#A6CEE3",
                               "Current (2021–2022)" = "#1F78B4")) +
  expand_limits(y = y_top * 1.18) +
  # Top-left annotation with test name and p-value
  annotate("text",
           x = 0.55, y = y_top * 1.15,
           hjust = 0, vjust = 1,
           label = test_label, size = 4.2) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )
###########################################################################
# ==========================================
# Combined barplots (Past & Current together)
# Kruskal–Wallis + Dunn(BH) letters per period
# ==========================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(FSA)            # dunnTest
library(multcompView)   # multcompLetters

sites_focus   <- c("M4","M7","M9")
period_levels <- c("Past (2013–2016)", "Current (2021–2022)")

# Keep only the sites and periods of interest
div_m2 <- div_m %>%
  filter(location_id %in% sites_focus,
         period %in% period_levels) %>%
  mutate(
    location_id = factor(location_id, levels = sites_focus),
    period      = factor(period, levels = period_levels)
  )

# ---------- helper: Dunn letters among sites for ONE period ----------
letters_for_period <- function(dat_period, sites_order) {
  dat_period <- dat_period %>% filter(location_id %in% sites_order)
  if (length(unique(dat_period$location_id)) < 2) {
    return(list(
      letters_df = data.frame(location_id = sites_order, letter = NA_character_),
      kw_p = NA_real_
    ))
  }
  # Global KW
  kw <- kruskal.test(shannon_H ~ location_id, data = dat_period)
  # Pairwise Dunn (BH)
  dunn <- FSA::dunnTest(shannon_H ~ location_id, data = dat_period, method = "bh")
  pw <- dunn$res %>% select(Comparison, P.adj)
  
  present <- sort(unique(dat_period$location_id))
  pmat <- matrix(1, nrow = length(present), ncol = length(present),
                 dimnames = list(present, present))
  if (nrow(pw) > 0) {
    pairs <- do.call(rbind, strsplit(pw$Comparison, " - "))
    colnames(pairs) <- c("g1","g2")
    p_tab <- data.frame(g1 = pairs[,1], g2 = pairs[,2], p = pw$P.adj,
                        stringsAsFactors = FALSE)
    for (i in seq_len(nrow(p_tab))) {
      if (p_tab$g1[i] %in% present && p_tab$g2[i] %in% present) {
        pmat[p_tab$g1[i], p_tab$g2[i]] <- p_tab$p[i]
        pmat[p_tab$g2[i], p_tab$g1[i]] <- p_tab$p[i]
      }
    }
  }
  letters_obj <- multcompView::multcompLetters(pmat < 0.05)
  letters_df  <- data.frame(location_id = names(letters_obj$Letters),
                            letter      = unname(letters_obj$Letters),
                            stringsAsFactors = FALSE)
  letters_df <- data.frame(location_id = sites_order, stringsAsFactors = FALSE) %>%
    left_join(letters_df, by = "location_id")
  list(letters_df = letters_df, kw_p = kw$p.value)
}

# ---------- letters & summaries per period ----------
# Past
res_past <- letters_for_period(filter(div_m2, period == period_levels[1]), sites_focus)
# Current
res_curr <- letters_for_period(filter(div_m2, period == period_levels[2]), sites_focus)

letters_all <- bind_rows(
  res_past$letters_df %>% mutate(period = period_levels[1]),
  res_curr$letters_df %>% mutate(period = period_levels[2])
) %>%
  mutate(
    period      = factor(period, levels = period_levels),
    location_id = factor(location_id, levels = sites_focus)
  )

# Summaries for plotting (mean ± SE of H′ per site × period)
summ_H <- div_m2 %>%
  group_by(period, location_id) %>%
  summarise(
    n_years = dplyr::n(),
    mean_H  = mean(shannon_H, na.rm = TRUE),
    sd_H    = sd(shannon_H,   na.rm = TRUE),
    se_H    = sd_H / sqrt(n_years),
    .groups = "drop"
  ) %>%
  left_join(letters_all, by = c("period","location_id")) %>%
  group_by(period) %>%
  mutate(
    y_lab = mean_H + se_H + 0.05 * max(mean_H + se_H, na.rm = TRUE),
    kw_label = dplyr::case_when(
      period == period_levels[1] ~ paste0("Kruskal–Wallis p = ",
                                          formatC(res_past$kw_p, format = "f", digits = 3)),
      period == period_levels[2] ~ paste0("Kruskal–Wallis p = ",
                                          formatC(res_curr$kw_p, format = "f", digits = 3)),
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup()

# ---------- single combined figure (two panels) ----------
ymax <- summ_H %>%
  group_by(period) %>%
  summarise(ym = max(y_lab, na.rm = TRUE), .groups = "drop") %>%
  pull(ym) %>% max(na.rm = TRUE)

ggplot(summ_H, aes(x = location_id, y = mean_H)) +
  geom_col(width = 0.65, color = "black", fill = "grey70") +
  geom_errorbar(aes(ymin = pmax(mean_H - se_H, 0), ymax = mean_H + se_H),
                width = 0.22) +
  geom_text(aes(y = y_lab, label = letter),
            vjust = 0, size = 6, fontface = "bold", na.rm = TRUE) +
  facet_wrap(~ period, ncol = 2) +
  scale_y_continuous(limits = c(0, ymax * 1.05),
                     breaks = scales::pretty_breaks(n = 6)) +
  labs(
    title = "",
    x = "Site",
    y = "Shannon–Wiener (H′; mean ± SE)"
  ) +
  # KW p-value annotation per panel
  geom_text(
    data = summ_H %>% group_by(period) %>% slice(1),
    aes(x = 0.6, y = max(summ_H$y_lab, na.rm = TRUE) * 1.02, label = kw_label),
    inherit.aes = FALSE, hjust = 0, vjust = 0, size = 4.2
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    strip.text  = element_text(face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 12)
  )
###########################################################################
# ==========================================================
# Fish diversity (Current only: 2021–2022)
# Shannon–Wiener H′ per site-year + boxplot
# Global test among sites: Friedman test (year as block)
# (Optional) Pairwise paired Wilcoxon among selected sites
# ==========================================================

rm(list = ls())

library(tidyverse)
library(janitor)
library(readr)
library(vegan)
library(rstatix)
library(ggpubr)

# ----------------------------
# 1) Load + prepare current data (2021–2022)
# ----------------------------
url_curr <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv"
years_keep <- c(2021, 2022)

curr_wide <- read_csv(url_curr, show_col_types = FALSE) %>%
  clean_names() %>%
  mutate(
    fish_species = str_squish(tolower(fish_species)),
    abundance = 1
  ) %>%
  filter(
    !is.na(location_id),
    !is.na(sampling_year),
    !is.na(fish_species),
    sampling_year %in% years_keep,
    location_id %in% paste0("M", 1:9)
  ) %>%
  group_by(location_id, sampling_year, fish_species) %>%
  summarise(n = sum(abundance), .groups = "drop") %>%
  pivot_wider(names_from = fish_species, values_from = n, values_fill = 0) %>%
  relocate(location_id, sampling_year)

# ----------------------------
# 2) Compute Shannon H′ per site-year
# ----------------------------
species_cols <- setdiff(names(curr_wide), c("location_id", "sampling_year"))
species_mat  <- curr_wide %>% select(all_of(species_cols)) %>% as.data.frame()
species_mat[is.na(species_mat)] <- 0

H <- vegan::diversity(species_mat, index = "shannon")

div_curr <- curr_wide %>%
  select(location_id, sampling_year) %>%
  mutate(
    shannon_H = H,
    sampling_year = factor(sampling_year),
    location_id = factor(location_id, levels = paste0("M", 1:9))
  ) %>%
  filter(is.finite(shannon_H))

# ----------------------------
# 3) Diagnostics: confirm you truly have 2 obs per site (2021 & 2022)
# ----------------------------
cat("\nCounts per site-year:\n")
print(div_curr %>% count(location_id, sampling_year) %>% arrange(location_id, sampling_year), n = 100)

cat("\nCounts per site (should be 2 each):\n")
print(div_curr %>% count(location_id) %>% arrange(n, location_id), n = 100)

# ----------------------------
# 4) Keep only sites with BOTH years (balanced design)
#    This prevents 'not enough observations' errors
# ----------------------------
div_bal <- div_curr %>%
  group_by(location_id) %>%
  filter(n_distinct(sampling_year) == 2) %>%
  ungroup()

# If any sites were dropped, print which ones
dropped <- setdiff(levels(div_curr$location_id), unique(as.character(div_bal$location_id)))
if (length(dropped) > 0) {
  message("\nDropped sites (missing one of the years 2021/2022 or non-finite H′): ",
          paste(dropped, collapse = ", "))
}

# Reset factor levels to the remaining sites only (helps plotting)
div_bal <- div_bal %>%
  mutate(location_id = factor(as.character(location_id),
                              levels = sort(unique(as.character(location_id)))))

# ----------------------------
# 5) Global test among sites accounting for year: Friedman test
#    (Nonparametric repeated-measures; year is the block)
# ----------------------------
fried_test <- div_bal %>%
  friedman_test(shannon_H ~ location_id | sampling_year)

cat("\nFriedman test:\n")
print(fried_test)

p_fried <- fried_test$p[1]

# ----------------------------
# 6) OPTIONAL: Pairwise site comparisons (paired)
#    WARNING: For many sites this is lots of comparisons; plot will get messy.
#    Best practice: restrict to a small set of sites.
# ----------------------------

# Toggle this to TRUE if you want pairwise tests on the plot
ADD_PAIRWISE <- FALSE

# If you want pairwise comparisons, choose a manageable subset of sites:
sites_of_interest <- c("M4", "M7", "M9")  # edit as needed

pairwise_tbl <- NULL
div_pair <- NULL

if (ADD_PAIRWISE) {
  div_pair <- div_bal %>%
    filter(as.character(location_id) %in% sites_of_interest) %>%
    mutate(location_id = factor(as.character(location_id), levels = sites_of_interest))
  
  # Ensure still balanced within chosen sites
  div_pair <- div_pair %>%
    group_by(location_id) %>%
    filter(n_distinct(sampling_year) == 2) %>%
    ungroup()
  
  # Paired pairwise Wilcoxon among sites (paired by year)
  pairwise_tbl <- div_pair %>%
    pairwise_wilcox_test(
      shannon_H ~ location_id,
      paired = TRUE,
      p.adjust.method = "BH"
    ) %>%
    add_xy_position(x = "location_id")
  
  cat("\nPairwise paired Wilcoxon (BH-adjusted):\n")
  print(pairwise_tbl)
  
  # Keep only significant (optional)
  pairwise_tbl <- pairwise_tbl %>% filter(p.adj <= 0.05)
}

# ----------------------------
# 7) Plot: boxplot + points + global Friedman p-value
# ----------------------------
y_top <- max(div_bal$shannon_H, na.rm = TRUE)

p <- ggplot(div_bal, aes(x = location_id, y = shannon_H)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(aes(shape = sampling_year), width = 0.08, size = 2.4, alpha = 0.75) +
  labs(
    title = "Fish diversity (Shannon–Wiener H′), 2021–2022",
    x = "Sampling Site",
    y = "H′ (Shannon–Wiener)",
    shape = "Year"
  ) +
  annotate(
    "text",
    x = 1, y = y_top * 1.13,
    hjust = 0,
    label = paste0("Friedman test : p = ", signif(p_fried, 3)),
    size = 4.2
  ) +
  expand_limits(y = y_top * 1.18) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.position = "right"
  )

# Add pairwise brackets only if enabled and available
if (ADD_PAIRWISE && !is.null(pairwise_tbl) && nrow(pairwise_tbl) > 0) {
  p <- p + stat_pvalue_manual(
    pairwise_tbl,
    label = "p.adj.signif",
    tip.length = 0.01,
    hide.ns = TRUE
  )
}

print(p)

# ----------------------------
# 8) Save figure (optional)
# ----------------------------
# ggsave("fish_shannon_2021_2022_sites.png", p, width = 10, height = 5.5, dpi = 300)

