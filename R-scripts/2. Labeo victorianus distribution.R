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




