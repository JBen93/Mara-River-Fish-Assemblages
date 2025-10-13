#Fish Biomass in the Mara River, Kenya, during past sampling (2013,2014, 2016) and current sampling (2021-2022)
# clear everything in memory (of R)
remove(list=ls())
# load the renv package
renv::restore()
#load libararies 
library(tidyverse) # for dplyr, ggplot2, tidyr, etc.
library(vegan) # for read_csv
library(janitor) # for clean_names
library(stringr) # for str_remove

# --------------------------------------------------------------
# Past fish data (2013, 2014, 2016)
#data URL source if you need to inspect for the whole dataset
#browseURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pubhtml")

# ---- Load past data ----
pastfish <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names()   # -> location_id, sampling_year, fish_species, fish_weight, total_length, standard_length, ...

# ---- Filter: sites & years ----
df <- pastfish %>%
  filter(location_id %in% paste0("M", 2:9),
         sampling_year %in% c(2013, 2014, 2016)) %>%
  mutate(
    # length is optional for biomass; keep if present but DON'T filter by it
    length_mm = dplyr::coalesce(standard_length),
    weight_g  = fish_weight
  ) %>%
  # Only require valid weights (biomass uses weight)
  filter(!is.na(fish_species), fish_species != "",
         !is.na(weight_g), weight_g > 0)

# ---- Biomass per species × site × year: biomass = n * mean(weight) ----
biomass_spp_year <- df %>%
  group_by(location_id, fish_species, sampling_year) %>%
  summarise(
    n_fish        = dplyr::n(),
    mean_weight_g = mean(weight_g, na.rm = TRUE),
    biomass_g     = n_fish * mean_weight_g,
    .groups = "drop"
  )

# ---- Collapse across years to species totals per site ----
biomass_spp_site <- biomass_spp_year %>%
  group_by(location_id, fish_species) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop")

# ---- Site-level biomass (sum across species, all past years combined) ----
biomass_site <- biomass_spp_site %>%
  group_by(location_id) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(location_id = factor(location_id, levels = paste0("M", 2:9))) %>%
  arrange(location_id)

print(biomass_site)

# ---- Build numeric site index for regression (M2 -> 2, ..., M9 -> 9) ----
site_stats2 <- biomass_site %>%
  mutate(
    site_order = as.numeric(str_remove(as.character(location_id), "^M")),
    location_id = factor(location_id, levels = paste0("M", 2:9))
  ) %>%
  arrange(site_order)

# ---- Linear model: biomass ~ site_order ----
m_site   <- lm(biomass_g ~ site_order, data = site_stats2)
sm_site  <- summary(m_site)
r2_site  <- sm_site$r.squared
p_site   <- coef(sm_site)[2, 4]

# ---- Plot: regression (biomass vs site order) ----
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
    title = "2013–2016",
    x = "Sampling Site",
    y = expression("Biomass (g) = n × mean weight")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

# ---- Optional: bar plot of site biomass ----
ggplot(biomass_site, aes(x = location_id, y = biomass_g)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  labs(
    title = "2013–2016",
    x = "Site",
    y = expression("Biomass (g) = n × mean weight")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
# --------------------------------------------------------------
# --------------------------------------------------------------
# Current fish data (2021-2022)
# clear everything in memory (of R)
remove(list=ls())
# load the renv package
renv::restore()
#load libararies 
library(tidyverse) # for dplyr, ggplot2, tidyr, etc.
library(readr) # for read_csv
library(janitor) # for clean_names
library(stringr) # for str_remove
# Load data
currentfish<- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>% 
  janitor::clean_names()   # standardize: fish_weight, total_length, standard_length, location_id, sampling_year, fish_species, ...


# Filter + choose length column; keep valid rows
df1 <- currentfish %>%
  filter(location_id %in% paste0("M", 2:9),
         sampling_year %in% c(2021, 2022)) %>%
  mutate(
    length_mm = dplyr::coalesce(standard_length),
    weight_g  = fish_weight
  ) %>%
  filter(!is.na(fish_species), fish_species != "",
         !is.na(length_mm), !is.na(weight_g),
         length_mm > 0, weight_g > 0)

# --------------------------------------------------------------
# Biomass per species, site, year 
# biomass_species_year = n_fish(species, site, year) * mean_weight(species, site, year)
# --------------------------------------------------------------
biomass_spp_year <- df1 %>%
  group_by(location_id, fish_species, sampling_year) %>%
  summarise(
    n_fish        = dplyr::n(),
    mean_weight_g = mean(weight_g, na.rm = TRUE),
    biomass_g     = n_fish * mean_weight_g,
    .groups = "drop"
  )

# Collapse across years to species totals per site by summing the biomass within the 2 years (2021 + 2022)
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
    title = "2021–2022",
    x = "Sampling Site",
    y = expression("Biomass (g) = n × mean weight")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 14,
      hjust = 0.5          # Center the title
    ),
    axis.text.x = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12)
  )

# ---- Optional: also show a bar plot for quick comparison ----
ggplot(biomass_site, aes(x = location_id, y = biomass_g)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  labs(
    title = "2021–2022",
    x = "Site",
    y = expression("Biomass (g) = n × mean weight")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      size = 14,
      hjust = 0.5,     # Centers the title horizontally
      vjust = 1        # Adjusts vertical placement closer to top
    ),
    plot.margin = margin(10, 10, 10, 10)  # Add balanced margins around plot
  )


# -----------------------------
