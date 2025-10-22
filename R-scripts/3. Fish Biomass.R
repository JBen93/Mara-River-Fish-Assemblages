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
pastbiomass_spp_year <- df %>%
  group_by(location_id, fish_species, sampling_year) %>%
  summarise(
    n_fish        = dplyr::n(),
    mean_weight_g = mean(weight_g, na.rm = TRUE),
    biomass_g     = n_fish * mean_weight_g,
    .groups = "drop"
  )

# ---- Collapse across years to species totals per site ----
pastbiomass_spp_site <- pastbiomass_spp_year %>%
  group_by(location_id, fish_species) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop")

# ---- Site-level biomass (sum across species, all past years combined) ----
pastbiomass_site <- pastbiomass_spp_site %>%
  group_by(location_id) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(location_id = factor(location_id, levels = paste0("M", 2:9))) %>%
  arrange(location_id)

print(pastbiomass_site)

# ---- Build numeric site index for regression (M2 -> 2, ..., M9 -> 9) ----
site_stats2 <- pastbiomass_site %>%
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

# --- Compute biomass per site × year (replicates) ---
site_year_biomass <- pastbiomass_spp_year %>%
  group_by(location_id, sampling_year) %>%                 # replicate = year within site
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop")

# --- Summarize to mean ± SE across years for each site ---
biomass_summary_se <- site_year_biomass %>%
  group_by(location_id) %>%
  summarise(
    n_years      = dplyr::n(),
    mean_biomass = mean(biomass_g, na.rm = TRUE),
    sd_biomass   = sd(biomass_g,   na.rm = TRUE),
    se_biomass   = sd_biomass / sqrt(n_years),
    .groups = "drop"
  ) %>%
  mutate(
    location_id = factor(location_id, levels = paste0("M", 2:9)),
    se_biomass  = ifelse(is.na(se_biomass), 0, se_biomass) # if only 1 year for a site
  )

print(biomass_summary_se)

# --- Bar plot: mean ± SE ---
pd <- position_dodge(width = 0.7)

ggplot(biomass_summary_se, aes(x = location_id, y = mean_biomass)) +
  geom_col(width = 0.65, color = "black", fill = "grey70") +
  geom_errorbar(
    aes(ymin = mean_biomass - se_biomass, ymax = mean_biomass + se_biomass),
    width = 0.22, position = pd
  ) +
  labs(
    title = "Fish Biomass (2013, 2014,2016)",
    x = "Site",
    y = "Biomass (g; mean ± SE)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0.5),
        axis.text.x = element_text(size = 12),
        axis.title = element_text(size = 13))
  


# --- Sum Biomass - total Yield: total biomass per site ---(this is less common due to variation brought by sampling years)

ggplot(pastbiomass_site, aes(x = location_id, y = biomass_g)) +
  geom_col(fill = "grey70", color = "black", width = 0.7) +
  labs(
    title = "2013–2016",
    x = "Site",
    y = expression("Biomass (g) = n × mean weight")
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Current fish data (2021-2022)
# clear everything in memory (of R)
# load the renv package
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
  filter(location_id %in% paste0("M", 1:9),
         sampling_year %in% c(2021, 2022)) %>%
  mutate(
    length_mm = dplyr::coalesce(total_length),
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
currentbiomass_site <- biomass_spp_site %>%
  group_by(location_id) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  mutate(location_id = factor(location_id, levels = paste0("M", 2:9))) %>%
  arrange(location_id)

print(currentbiomass_site)

# ---- Regression: Biomass vs Site Order (M2 -> M9) ----
site_stats2 <- currentbiomass_site %>%
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
#-- --------------------------------------------------------------
# --- 1) Biomass per site × year (replicates for SE) ---
site_year_biomass <- biomass_spp_year %>%
  group_by(location_id, sampling_year) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop")

# --- 2) Mean ± SE across years (2021 & 2022) for each site ---
current_biomass_mean_se <- site_year_biomass %>%
  group_by(location_id) %>%
  summarise(
    n_years      = dplyr::n(),                          # should be 2 if both years present
    mean_biomass = mean(biomass_g, na.rm = TRUE),
    sd_biomass   = sd(biomass_g,   na.rm = TRUE),
    se_biomass   = sd_biomass / sqrt(n_years),
    .groups = "drop"
  ) %>%
  mutate(
    # keep your site ordering M2..M9
    location_id = factor(location_id, levels = paste0("M", 2:9)),
    # if a site only has 1 year, SE would be NA; set to 0 so the plot still draws
    se_biomass  = dplyr::coalesce(se_biomass, 0)
  ) %>%
  arrange(location_id)

print(current_biomass_mean_se)

# --- 3) Barplot: Mean biomass ± SE (2021–2022) ---
plot_df <- current_biomass_mean_se %>%
  dplyr::mutate(
    mean_biomass = as.numeric(mean_biomass),
    se_biomass   = tidyr::replace_na(se_biomass, 0)
  ) %>%
  tidyr::drop_na(mean_biomass)

ggplot(plot_df, aes(x = location_id, y = mean_biomass)) +
  geom_col(width = 0.65, color = "black", fill = "grey70") +
  geom_errorbar(aes(ymin = pmax(mean_biomass - se_biomass, 0),
                    ymax = mean_biomass + se_biomass),
                width = 0.22) +
  scale_y_continuous(breaks = seq(0, 14000, by = 2000)) +
  coord_cartesian(ylim = c(0, 14000)) +   # clamps display without dropping data
  labs(
    title = "Fish Biomass (2021–2022)",
    x = "Site",
    y = "Biomass (g; mean ± SE)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.text  = element_text(size = 12),
    axis.title = element_text(size = 13)
  )

# --- 4) Regression on MEAN biomass (not totals) across sites ---
site_stats_mean <- current_biomass_mean_se %>%
  mutate(site_order = as.numeric(stringr::str_remove(as.character(location_id), "^M"))) %>%
  arrange(site_order)

m_site_mean  <- lm(mean_biomass ~ site_order, data = site_stats_mean)
sm_site_mean <- summary(m_site_mean)
r2_site_mean <- sm_site_mean$r.squared
p_site_mean  <- coef(sm_site_mean)[2, 4]

# Plot regression using mean biomass per site
ggplot(site_stats_mean, aes(x = site_order, y = mean_biomass)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  scale_x_continuous(
    breaks = site_stats_mean$site_order,
    labels = site_stats_mean$location_id
  ) +
  annotate(
    "text",
    x = min(site_stats_mean$site_order, na.rm = TRUE),
    y = max(site_stats_mean$mean_biomass, na.rm = TRUE),
    hjust = 0, vjust = 1,
    label = paste0("R² = ", round(r2_site_mean, 3),
                   "\nP = ", format.pval(p_site_mean, digits = 3, eps = 1e-3))
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Mean Biomass (g; 2021–2022)"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

# ---- Sum Biomass - total Yield: total biomass per site ---(this is less common due to variation brought by sampling years) ----
ggplot(currentbiomass_site, aes(x = location_id, y = biomass_g)) +
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
# --------------------------------------------------------------
# Compare past (2013–2016) vs current (2021–2022) biomass at M4, M7, M9
# --------------------------------------------------------------
#make an object for the current biomass data
curr_biomass_spp_year <- biomass_spp_year
library(dplyr)
library(ggplot2)

sites_focus <- c("M4","M7","M9")

# --- Build site×year totals for each period (years = replicates) ---
site_year_past <- pastbiomass_spp_year %>%
  group_by(location_id, sampling_year) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  filter(location_id %in% sites_focus) %>%
  mutate(period = "Past (2013–2016)")

site_year_curr <- curr_biomass_spp_year %>%
  group_by(location_id, sampling_year) %>%
  summarise(biomass_g = sum(biomass_g, na.rm = TRUE), .groups = "drop") %>%
  filter(location_id %in% sites_focus) %>%
  mutate(period = "Current (2021–2022)")

site_year_all <- bind_rows(site_year_past, site_year_curr)

# --- Mean ± SE per site × period ---
summary_site_period <- site_year_all %>%
  group_by(location_id, period) %>%
  summarise(
    n_years      = dplyr::n(),
    mean_biomass = mean(biomass_g, na.rm = TRUE),
    sd_biomass   = sd(biomass_g,   na.rm = TRUE),
    se_biomass   = sd_biomass / sqrt(n_years),
    .groups = "drop"
  ) %>%
  mutate(
    location_id = factor(location_id, levels = sites_focus),
    period      = factor(period, levels = c("Past (2013–2016)", "Current (2021–2022)")),
    se_biomass  = dplyr::coalesce(se_biomass, 0)
  )

# --- Grouped barplot: mean biomass ± SE, fill by period ---
pd <- position_dodge(width = 0.7)

ggplot(summary_site_period, aes(x = location_id, y = mean_biomass, fill = period)) +
  geom_col(position = pd, width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = pmax(mean_biomass - se_biomass, 0),
                    ymax = mean_biomass + se_biomass),
                width = 0.22, position = pd) +
  scale_fill_manual(values = c("Past (2013–2016)" = "#A6CEE3",
                               "Current (2021–2022)" = "#1F78B4")) +
  scale_y_continuous(limits = c(0, 14000), breaks = seq(0, 14000, by = 2000)) +
  labs(
    title = "Mean Fish Biomass",
    x = "Site",
    y = "Biomass (g; mean ± SE)",
    fill = "Sampling Period"
  ) +
  # ---- Add p-value annotation ----
annotate("text",
         x = 2,                     # horizontal position (M7 = middle site)
         y = 13500,               # vertical position (top of y-axis)
         label = "(Welch t-test: p = 0.1444)",
         size = 5,
         fontface = "italic") +
  # ---- Styling ----
theme_minimal(base_size = 13) +
  theme(
    plot.title  = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(size = 12),
    axis.title  = element_text(size = 12)
  )


# --------------------------------------------------------------
# --- Combine past and current datasets for comparison ---

# Add a "period" column to each dataset
past_biomass_spp_year <- pastbiomass_spp_year %>%
  mutate(period = "Past (2013–2016)")

current_biomass_spp_year <- biomass_spp_year %>%
  mutate(period = "Current (2021–2022)")

# Combine the two into one data frame
combined_df <- bind_rows(past_biomass_spp_year, current_biomass_spp_year)

t.test(mean_weight_g ~ period, data = combined_df)


