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

print(p_box_sites)

# (Optional) Also print the exact p-values table
pvals %>% select(location_id, p, p.adj, p.label) %>% arrange(location_id) %>% print()
