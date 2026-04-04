# ================================
# Relative abundance of macroinvertebrate families
# Sites: M4–M9; Years: 2021–2022
# Legend shows only families that hit ≥10% at any site
# ================================

# Clear environment
remove(list = ls())
renv::restore()
# Libraries
library(tidyverse)
library(readr)

# ---- Load data (use exact column headers from the sheet) ----
macros <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Filter to sites M4–M9 and years 2021–2022 ----
sites_focus <- paste0("M", 1:9)

macro_filt <- macros %>%
  filter(
    Location_ID %in% sites_focus,
    year %in% c(2021, 2022),
    !is.na(Family), Family != ""
  ) %>%
  mutate(
    Location_ID = factor(Location_ID, levels = sites_focus),
    Count = as.numeric(Count)  # ensure numeric
  )

# ---- Totals per site (denominator) ----
site_totals <- macro_filt %>%
  group_by(Location_ID) %>%
  summarise(Total_N = sum(Count, na.rm = TRUE), .groups = "drop")

# ---- Counts per family × site (numerator) ----
family_counts <- macro_filt %>%
  group_by(Location_ID, Family) %>%
  summarise(N_Family = sum(Count, na.rm = TRUE), .groups = "drop")

# ---- Relative abundance (%) per family per site ----
rel_abund <- family_counts %>%
  left_join(site_totals, by = "Location_ID") %>%
  mutate(Rel_Percent = if_else(Total_N > 0, 100 * N_Family / Total_N, 0))

# ---- Order families by overall abundance (for legend order) ----
family_order_all <- rel_abund %>%
  group_by(Family) %>%
  summarise(total_rel = sum(Rel_Percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_rel)) %>%
  pull(Family)

rel_abund <- rel_abund %>%
  mutate(Family = factor(Family, levels = family_order_all))

# ---- Legend: only families that ever reach ≥10% at any site ----
legend_families <- rel_abund %>%
  group_by(Family) %>%
  summarise(max_rel = max(Rel_Percent, na.rm = TRUE), .groups = "drop") %>%
  filter(max_rel >= 10) %>%
  pull(Family)

# ---- Plot ----
ggplot(rel_abund, aes(x = Location_ID, y = Rel_Percent, fill = Family)) +
  geom_col(position = "stack", color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 100)) +  # clip without dropping rows
  scale_fill_viridis_d(option = "turbo", breaks = legend_families) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.text.x   = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  )

# ================================
# OPTIONAL: lump very rare families into "Other (<10%)" per site
# (uncomment this block if you want simpler bars)
rel_abund_lumped <- rel_abund %>%
  group_by(Location_ID) %>%
 mutate(Family2 = if_else(Rel_Percent < 10, "Other (<10%)", as.character(Family))) %>%
group_by(Location_ID, Family2) %>%
summarise(Rel_Percent = sum(Rel_Percent), .groups = "drop")

ggplot(rel_abund_lumped, aes(Location_ID, Rel_Percent, fill = Family2)) +
 geom_col(color = "black", linewidth = 0.25) +
scale_y_continuous(labels = scales::percent_format(scale = 1)) +
coord_cartesian(ylim = c(0, 100)) +
scale_fill_viridis_d(option = "turbo") +
labs(title = "",
    x = "Site", y = "Relative Abundance (%)", fill = "Family") +
theme_minimal(base_size = 13) +
 theme(plot.title = element_text(face = "bold", hjust = 0.5))
# ================================
# ================================
# Relative abundance of macroinvertebrate families
# Sites: M4–M9; Years: 2021–2022
# Legend shows only families that hit ≥10% at any site
# ================================

# Clear environment
remove(list = ls())
renv::restore()
# Libraries
library(tidyverse)
library(readr)

# ---- Load data (use exact column headers from the sheet) ----
macros <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vR9TMKMzDZtRRS5WAsC1N-8lcQyAB7FM5IInNfD7kDp-AtWM1tG57aLG2Hgq3RVrRFNE8VQq8mrqbhl/pub?gid=1254679428&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Filter to sites M4–M9 and years 2021–2022 ----
sites_focus <- paste0("M", 4:9)

macro_filt <- macros %>%
  filter(
    Location_ID %in% sites_focus,
    year %in% c(2021, 2022),
    !is.na(Family), Family != ""
  ) %>%
  mutate(
    Location_ID = factor(Location_ID, levels = sites_focus),
    Count = as.numeric(Count)  # ensure numeric
  )

# ---- Totals per site (denominator) ----
site_totals <- macro_filt %>%
  group_by(Location_ID) %>%
  summarise(Total_N = sum(Count, na.rm = TRUE), .groups = "drop")

# ---- Counts per family × site (numerator) ----
family_counts <- macro_filt %>%
  group_by(Location_ID, Family) %>%
  summarise(N_Family = sum(Count, na.rm = TRUE), .groups = "drop")

# ---- Relative abundance (%) per family per site ----
rel_abund <- family_counts %>%
  left_join(site_totals, by = "Location_ID") %>%
  mutate(Rel_Percent = if_else(Total_N > 0, 100 * N_Family / Total_N, 0))

# ---- Order families by overall abundance (for legend order) ----
family_order_all <- rel_abund %>%
  group_by(Family) %>%
  summarise(total_rel = sum(Rel_Percent, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_rel)) %>%
  pull(Family)

rel_abund <- rel_abund %>%
  mutate(Family = factor(Family, levels = family_order_all))

# ---- Legend: only families that ever reach ≥10% at any site ----
legend_families <- rel_abund %>%
  group_by(Family) %>%
  summarise(max_rel = max(Rel_Percent, na.rm = TRUE), .groups = "drop") %>%
  filter(max_rel >= 10) %>%
  pull(Family)

# ---- Plot ----
ggplot(rel_abund, aes(x = Location_ID, y = Rel_Percent, fill = Family)) +
  geom_col(position = "stack", color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 100)) +  # clip without dropping rows
  scale_fill_viridis_d(option = "turbo", breaks = legend_families) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)",
    fill = "Family"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.text.x   = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 10)
  )

# ================================
# OPTIONAL: lump very rare families into "Other (<10%)" per site
# (uncomment this block if you want simpler bars)
rel_abund_lumped <- rel_abund %>%
  group_by(Location_ID) %>%
  mutate(Family2 = if_else(Rel_Percent < 10, "Other (<10%)", as.character(Family))) %>%
  group_by(Location_ID, Family2) %>%
  summarise(Rel_Percent = sum(Rel_Percent), .groups = "drop")

ggplot(rel_abund_lumped, aes(Location_ID, Rel_Percent, fill = Family2)) +
  geom_col(color = "black", linewidth = 0.25) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  coord_cartesian(ylim = c(0, 100)) +
  scale_fill_viridis_d(option = "turbo") +
  labs(title = "",
       x = "Site", y = "Relative Abundance (%)", fill = "Family") +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
# ================================
