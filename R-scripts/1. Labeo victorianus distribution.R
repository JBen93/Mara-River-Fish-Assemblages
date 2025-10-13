#Labeo victorianus distribution analysis in the Mara River
# clear everything in memory (of R)
remove(list=ls())

renv::restore()
# =========================================================
# Labeo victorianus (M1–M9, 2021–2022)
# Distribution and abundance analysis
# =========================================================

# ---- Load libraries ----
libs <- c("tidyverse", "janitor", "MASS", "emmeans", "performance", "DHARMa", "rstatix")
to_install <- libs[!libs %in% installed.packages()[,"Package"]]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(libs, require, character.only = TRUE)

# ---- Read data (Google Sheet) ----
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE)
  
# Check column names
names(dat1)
# should now include: location_name, location_id, sampling_year, sampling_month, fish_species, fish_weight, total_length, standard_length

# ---- Filter for Labeo victorianus (M1–M9, years 2021–2022) ----
dat2 <- dat1 %>%
  mutate(
    site = as.character(location_ID),
    year = as.integer(sampling_year),
    species = tolower(trimws(fish_species))
  ) %>%
  filter(
    site %in% paste0("M", 1:9),
    year %in% c(2021, 2022),
    species == "labeo victorianus"
  )
dat2

# ---- Count individuals per Site × Year ----
dat_summary <- dat2 %>%
  count(site, year, name = "abundance") %>%
  complete(site = paste0("M", 1:9), year = c(2021, 2022), fill = list(abundance = 0)) %>%
  mutate(
    site = factor(site, levels = paste0("M", 1:9)),
    year = factor(year)
  )

# ---- Descriptive summary ----
dat_summary %>%
  group_by(site, year) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop")

## Non-parametric test of abundance by year
wilcox.test(abundance ~ year, data = dat_summary)

# =========================================================
# 1. Plot: Observed counts
# =========================================================

# ---- Labeo victorianus mean abundance per site (2021–2022) ----

library(tidyverse)
library(readr)

# Load data from Google Sheets
currentfish <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# Filter for Labeo victorianus, sites M1–M9, years 2021–2022
dat_lv <- currentfish %>%
  filter(fish_species == "Labeo victorianus",
         location_ID %in% paste0("M", 1:9),
         sampling_year %in% c(2021, 2022))

# Summarize abundance (count per site per year)
dat_summary <- dat_lv %>%
  group_by(location_ID, sampling_year) %>%
  summarise(abundance = n(), .groups = "drop")

# Compute the mean abundance across both years for each site
dat_mean <- dat_summary %>%
  group_by(location_ID) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop")

# ---- Create bar plot of mean abundance (grey color, no labels, Y-axis up to 50) ----
ggplot(dat_mean, aes(x = location_ID, y = mean_abundance, fill = "grey")) +
  geom_col(width = 0.7, color = "black") +
  labs(
    title = "Mean abundance of Labeo victorianus across sites (2021–2022)",
    x = "Sampling site",
    y = "Mean abundance (number of individuals)"
  ) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
  theme_minimal(base_size = 13) +
  scale_fill_manual(values = c("grey" = "grey70"), guide = "none") +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    axis.ticks.length = unit(0.25, "cm")
  )
####################################################################################
#Barplot for only three sites sampled in both periods, (M4,M7,M9)
##################################################################################
# ---- Mean abundance of Labeo victorianus (Sites M4, M7, M9; 2021–2022) ----
remove(list=ls())
library(tidyverse)
library(readr)

# Load data from Google Sheets
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# Filter for Labeo victorianus at sites M4, M7, and M9 for years 2021–2022
dat_lv <- dat1 %>%
  filter(
    fish_species == "Labeo victorianus",
    location_ID %in% c("M4", "M7", "M9"),
    sampling_year %in% c(2021, 2022)
  )

# Summarize abundance (count per site per year)
dat_summary <- dat_lv %>%
  group_by(location_ID, sampling_year) %>%
  summarise(abundance = n(), .groups = "drop")

# Compute mean abundance across both years for each site
dat_mean <- dat_summary %>%
  group_by(location_ID) %>%
  summarise(mean_abundance = mean(abundance), .groups = "drop")

# ---- Create bar plot (grey fill, Y-axis up to 50) ----
ggplot(dat_mean, aes(x = location_ID, y = mean_abundance, fill = "grey")) +
  geom_col(width = 0.7, color = "black") +
  labs(
    title = expression(paste("Mean Abundance of ", italic("Labeo victorianus"), 
                             " at Sites M4, M7, and M9 (2021–2022)")),
    x = "Sampling Site",
    y = "Mean Abundance (Number of Individuals)"
  ) +
  scale_y_continuous(limits = c(0, 30), expand = c(0, 0)) +
  scale_fill_manual(values = c("grey" = "grey70"), guide = "none") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
####################################################################################
#Relative abundance of the Labeo victorianus species at each site 
##################################################################################
# ---- Relative Abundance of Labeo victorianus (Sites M4, M7, M9; 2021–2022) ----
remove(list = ls())
library(tidyverse)
library(readr)

# Load data from Google Sheets
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# Filter data for sites M4, M7, M9 and years 2021–2022
dat_filtered <- dat1 %>%
  filter(
    location_ID %in% c("M4", "M7", "M9"),
    sampling_year %in% c(2021, 2022)
  )

# ---- Compute relative abundance of Labeo victorianus per site ----
dat_rel <- dat_filtered %>%
  group_by(location_ID, fish_species) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(location_ID) %>%
  mutate(
    total_count = sum(count),
    relative_abundance = (count / total_count) * 100
  ) %>%
  ungroup() %>%
  filter(fish_species == "Labeo victorianus")

# ---- Create bar plot (relative abundance %) ----
ggplot(dat_rel, aes(x = location_ID, y = relative_abundance, fill = "grey")) +
  geom_col(width = 0.7, color = "black") +
  labs(
    title = expression(paste("Relative Abundance of ", italic("Labeo victorianus"), 
                             "  (2021–2022)")),
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  scale_fill_manual(values = c("grey" = "grey70"), guide = "none") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
####################################################################################
#Relative abundance of the Labeo victorianus for all the sites in 2021-2022
##################################################################################
# ---- Relative Abundance of Labeo victorianus (Sites M4, M7, M9; 2021–2022) ----
remove(list = ls())
library(tidyverse)
library(readr)

# Load data from Google Sheets
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# Filter data for all the sites in 2021–2022
dat_filtered <- dat1 %>%
  filter(
    location_ID %in% c("M2","M3","M4","M5","M6","M7","M8","M9"),
    sampling_year %in% c(2021, 2022)
  )

# ---- Compute relative abundance of Labeo victorianus per site ----
dat_rel <- dat_filtered %>%
  group_by(location_ID, fish_species) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(location_ID) %>%
  mutate(
    total_count = sum(count),
    relative_abundance = (count / total_count) * 100
  ) %>%
  ungroup() %>%
  filter(fish_species == "Labeo victorianus")

# ---- Create bar plot (relative abundance %) ----
ggplot(dat_rel, aes(x = location_ID, y = relative_abundance, fill = "grey")) +
  geom_col(width = 0.7, color = "black") +
  labs(
    title = expression(paste("Relative Abundance of ", italic("Labeo victorianus"), 
                             "  (2021–2022)")),
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  scale_y_continuous(limits = c(0, 75), expand = c(0, 0)) +
  scale_fill_manual(values = c("grey" = "grey70"), guide = "none") +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 1, hjust = 0.5),
    axis.ticks.length = unit(0.25, "cm")
  )
####################################################################################
#Species accumulation curve#################################################################
# ---- Species Accumulation Curve (Sites M2, M3, M4,M5,M6, M7, M9; 2021–2022) ----
# ---- Species Accumulation Curve (Sites M2, M3, M4, M5, M6, M7, M9; 2021–2022) ----
remove(list = ls())

library(tidyverse)
library(readr)
library(vegan)

# Load data from Google Sheets
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# -----------------------------
# CONFIG: which sites to include
# -----------------------------
sites_keep <- c("M2","M3","M4","M5","M6","M7","M9")  # (M8 omitted per your list)
years_keep <- c(2021, 2022)

# Clean & filter (each row = one fish)
df <- dat1 %>%
  mutate(
    location_ID = as.character(location_ID),
    fish_species = stringr::str_squish(fish_species)
  ) %>%
  filter(location_ID %in% sites_keep,
         sampling_year %in% years_keep)

# Build site x species abundance matrix (counts of individuals)
#   - Rows: sites
#   - Cols: species
comm_mat <- df %>%
  count(location_ID, fish_species, name = "abund") %>%
  tidyr::pivot_wider(names_from = fish_species, values_from = abund, values_fill = 0) %>%
  tibble::column_to_rownames("location_ID") %>%
  as.data.frame()

# ---- Species Accumulation (randomized order, 1000 perms) ----
set.seed(123)  # for reproducibility
spec_acc <- vegan::specaccum(comm_mat, method = "random", permutations = 1000)

# Prepare data for ggplot (95% CI ≈ mean ± 1.96 * sd)
acc_df <- tibble(
  sites = spec_acc$sites,
  richness = spec_acc$richness,
  sd = spec_acc$sd
) %>%
  mutate(
    lower = pmax(richness - 1.96 * sd, 0),
    upper = richness + 1.96 * sd
  )

# ---- Plot: Species Accumulation Curve (All sites combined across 2021–2022) ----
ggplot(acc_df, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, fill = "grey70") +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(
    title = "Species Accumulation Curve",
    x = "Number of Sites",
    y = "Cumulative Species Richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
####################################################################################
#Species cumulative curve for only 3 sites M4,M7 and M9 during 2021-2022
#################################################################
remove(list = ls())

library(tidyverse)
library(readr)
library(vegan)

# Load data
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# -----------------------------
# Include only M4, M7, M9 (years 2021–2022)
# -----------------------------
sites_keep <- c("M4","M7","M9")
years_keep <- c(2021, 2022)

# Clean + filter
df <- dat1 %>%
  mutate(
    location_ID  = as.character(location_ID),
    fish_species = stringr::str_squish(fish_species)
  ) %>%
  filter(location_ID %in% sites_keep,
         sampling_year %in% years_keep,
         !is.na(fish_species) & fish_species != "")

# Build site × species matrix (counts)
comm_mat <- df %>%
  count(location_ID, fish_species, name = "abund") %>%
  tidyr::pivot_wider(names_from = fish_species, values_from = abund, values_fill = 0) %>%
  # ensure rows appear in the desired site order
  mutate(location_ID = factor(location_ID, levels = sites_keep)) %>%
  arrange(location_ID) %>%
  tibble::column_to_rownames("location_ID") %>%
  as.data.frame()

# (Optional) Use incidence (presence/absence) instead of counts:
# comm_mat[comm_mat > 0] <- 1

# ---- Species Accumulation (randomized order, 1000 perms) ----
set.seed(123)
spec_acc <- vegan::specaccum(comm_mat, method = "random", permutations = 1000)

# Prep for plotting (95% CI ≈ mean ± 1.96*sd)
acc_df <- tibble(
  sites    = spec_acc$sites,
  richness = spec_acc$richness,
  sd       = spec_acc$sd
) %>%
  mutate(
    lower = pmax(richness - 1.96 * sd, 0),
    upper = richness + 1.96 * sd
  )

# ---- Plot ----
ggplot(acc_df, aes(x = sites, y = richness)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.25) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:length(sites_keep), limits = c(1, length(sites_keep))) +
  labs(
    title = "Species Accumulation Curve",
    x     = "Number of Sites",
    y     = "Cumulative Species Richness"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5)
  )


##################################################################################
# Fish Biomass for all the sites in 2021 and 2022
##################################################################################
# ---- Fish Biomass (M2–M9; 2021–2022) + Length–Weight Regression Plot ----
remove(list = ls())

library(tidyverse)
library(readr)
library(janitor)

# Load data
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>% 
  janitor::clean_names()   # standardize: fish_weight, total_length, standard_length, location_id, sampling_year, fish_species, ...

# Filter sites and years
df <- dat1 %>%
  filter(location_id %in% paste0("M", 2:9),
         sampling_year %in% c(2021, 2022)) %>%
  mutate(
    # pick a length column (prefer total_length; fall back to standard_length if missing)
    length_mm = standard_length,
    weight_g  = fish_weight
  ) %>%
  # keep valid numeric observations
  filter(!is.na(weight_g), !is.na(length_mm), weight_g > 0, length_mm > 0)

# -----------------------------
# 1) Biomass (sum of weights) per site
# -----------------------------
biomass_site <- df %>%
  group_by(location_id) %>%
  summarise(biomass_g = sum(weight_g, na.rm = TRUE), .groups = "drop") %>%
  arrange(location_id)

print(biomass_site)

# (Optional) quick bar plot of biomass by site
ggplot(biomass_site, aes(x = location_id, y = biomass_g)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  labs(title = "Fish Biomass by Site (M2–M9; 2021–2022)",
       x = "Site", y = "Biomass (g)") +
  theme_minimal(base_size = 12)

# -----------------------------
# 2) Length–Weight regression (log–log)
#    log10(weight_g) ~ log10(length_mm)
# -----------------------------
df_log <- df %>%
  mutate(
    log_w = log10(weight_g),
    log_l = log10(length_mm)
  )

lm_log <- lm(log_w ~ log_l, data = df_log)
summary_lm <- summary(lm_log)
r2_log <- summary_lm$r.squared
p_log  <- coef(summary_lm)[2, 4]

# Plot: log–log regression with R² on the panel
ggplot(df_log, aes(x = log_l, y = log_w)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  annotate("text",
           x = min(df_log$log_l, na.rm = TRUE),
           y = max(df_log$log_w, na.rm = TRUE),
           hjust = 0, vjust = 1,
           label = paste0("R² = ", round(r2_log, 3),
                          "\nP-value = ", format.pval(p_log, digits = 3, eps = 1e-3))) +
  labs(title = "Length–Weight Relationship (log–log)",
       x = "log10(Length, mm)", y = "log10(Weight, g)") +
  theme_minimal(base_size = 12)

############################################################
# ---- Biomass by species-year, summed to site totals (M2–M9; 2021–2022) ----
remove(list = ls())

library(tidyverse)
library(readr)
library(janitor)
library(stringr)

# Load data
dat1 <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names()   # location_id, sampling_year, fish_species, fish_weight, total_length, standard_length, ...

# Filter + choose length column; keep valid rows
df <- dat1 %>%
  filter(location_id %in% paste0("M", 2:9),
         sampling_year %in% c(2021, 2022)) %>%
  mutate(
    length_mm = dplyr::coalesce(total_length, standard_length),
    weight_g  = fish_weight
  ) %>%
  filter(!is.na(fish_species), fish_species != "",
         !is.na(length_mm), !is.na(weight_g),
         length_mm > 0, weight_g > 0)

# --------------------------------------------------------------
# Biomass per species-year (your definition)
#   biomass_species_year = n_fish(species, site, year) * mean_weight(species, site, year)
# --------------------------------------------------------------
biomass_spp_year <- df %>%
  group_by(location_id, fish_species, sampling_year) %>%
  summarise(
    n_fish        = dplyr::n(),
    mean_weight_g = mean(weight_g, na.rm = TRUE),
    biomass_g     = n_fish * mean_weight_g,
    .groups = "drop"
  )

# Collapse across years to species totals per site
biomass_spp_site <- biomass_spp_year %>%
  group_by(location_id, fish_species) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop")

# Site-level biomass = sum across species (both years combined)
biomass_site <- biomass_spp_site %>%
  group_by(location_id) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(location_id = factor(location_id, levels = paste0("M", 2:9))) %>%
  arrange(location_id)

print(biomass_site)

# ---- Regression: Biomass vs Site Order (M2 -> M9) ----
site_stats2 <- biomass_site %>%
  mutate(
    site_order = as.numeric(stringr::str_remove(as.character(location_id), "^M")),
    location_id = factor(location_id, levels = paste0("M", 2:9))
  ) %>%
  arrange(site_order)

# Fit linear model: biomass ~ site_order
m_site <- lm(biomass_g ~ site_order, data = site_stats2)
sm_site <- summary(m_site)
r2_site <- sm_site$r.squared
p_site  <- coef(sm_site)[2, 4]

# Plot regression with R² and P, x-axis labeled as M2–M9
ggplot(site_stats2, aes(x = site_order, y = biomass_g)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_continuous(
    breaks = site_stats2$site_order,
    labels = site_stats2$location_id
  ) +
  annotate(
    "text",
    x = min(site_stats2$site_order, na.rm = TRUE),
    y = max(site_stats2$biomass_g, na.rm = TRUE),
    hjust = 0, vjust = 1,
    label = paste0("R² = ", round(r2_site, 3),
                   "\nP = ", format.pval(p_site, digits = 3, eps = 1e-3))
  ) +
  labs(
    title = "Biomass vs Site",
    x = "Site",
    y = "Biomass (g) = n × mean weight"
  ) +
  theme_minimal(base_size = 12)

# ---- Optional: also show a bar plot for quick comparison ----
ggplot(biomass_site, aes(x = location_id, y = biomass_g)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  labs(
    title = "Fish Biomass by Site (M2–M9; 2021–2022)",
    x = "Site",
    y = "Biomass (g) = n × mean weight"
  ) +
  theme_minimal(base_size = 12)

# -----------------------------
