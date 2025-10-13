####################################################################################
#Species accumulation curve#################################################################
# ---- Species Accumulation Curve (Sites M2, M3, M4,M5,M6, M7, M9; 2021–2022) ----
# ---- Species Accumulation Curve (Sites M2, M3, M4, M5, M6, M7, M9; 2021–2022) ----
# clear everything in memory (of R)
remove(list=ls())
# load the renv package
renv::restore()
library(tidyverse)
library(readr)
library(vegan)

#data URL source if you need to inspect for the whole dataset
#browseURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pubhtml")

# Load data from Google Sheets
currentspec <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
)

# -----------------------------
# CONFIG: which sites to include
# -----------------------------
sites_keep <- c("M2","M3","M4","M5","M6","M7","M9")  # (M8 omitted per your list)
years_keep <- c(2021, 2022)

# Clean & filter (each row = one fish)
currentspec2<- currentspec%>%
  mutate(
    location_ID = as.character(location_ID),
    fish_species = stringr::str_squish(fish_species)
  ) %>%
  filter(location_ID %in% sites_keep,
         sampling_year %in% years_keep)

# Build site x species abundance matrix (counts of individuals)
#   - Rows: sites
#   - Cols: species
comm_mat <- currentspec2 %>%
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




