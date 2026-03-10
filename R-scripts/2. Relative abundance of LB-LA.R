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
    title = "Labeo victorianus 2021–2022)",
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
#Relative abundance of the Labeo victorianus and Labeobarbus altinialis fish for all the sites in 2021-2022
##################################################################################

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
# Relative abundance of Labeo victorianus & Labeobarbus altianalis (M2–M9; 2021–2022)
remove(list = ls())

library(tidyverse)
library(readr)
library(janitor)
library(stringr)

# --- Load & clean data ---
dat <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    fish_species = str_squish(str_to_lower(fish_species))
  ) %>%
  filter(
    location_id %in% paste0("M", 4:9),
    sampling_year %in% c(2021, 2022),
    !is.na(fish_species), fish_species != ""
  )

# --- Canonicalize species names ---
dat <- dat %>%
  mutate(
    species_canon = case_when(
      str_detect(fish_species, "^labeo\\s+victorianus$") ~ "Labeo victorianus",
      str_detect(fish_species, "^labeobarbus\\s+alti?analis$") ~ "Labeobarbus altianalis",
      TRUE ~ NA_character_
    )
  )

# --- Count total catch per site (all species) ---
site_totals <- dat %>%
  count(location_id, name = "n_total")

# --- Count focal species per site ---
focal_counts <- dat %>%
  filter(!is.na(species_canon)) %>%
  count(location_id, species_canon, name = "n_focal") %>%
  complete(
    location_id = factor(paste0("M", 4:9), levels = paste0("M", 4:9)),
    species_canon = factor(c("Labeo victorianus", "Labeobarbus altianalis"),
                           levels = c("Labeo victorianus", "Labeobarbus altianalis")),
    fill = list(n_focal = 0)
  )

# --- Relative abundance (%) ---
rel_abund <- focal_counts %>%
  left_join(site_totals, by = "location_id") %>%
  mutate(
    rel_percent = if_else(n_total > 0, 100 * n_focal / n_total, 0)
  )

# --- Custom colours ---
pal_species <- c(
  "Labeo victorianus"      = "#1F78B4",  # deep blue
  "Labeobarbus altianalis" = "#E41A1C"   # vivid red
)

# --- Order legend by higher overall abundance ---
species_order <- rel_abund %>%
  group_by(species_canon) %>%
  summarise(mean_rel = mean(rel_percent, na.rm = TRUE)) %>%
  arrange(desc(mean_rel)) %>%
  pull(species_canon)

rel_abund <- rel_abund %>%
  mutate(species_canon = factor(species_canon, levels = species_order))
# -----------------------------
# Test whether relative abundance differs between the two species
# (paired across sites)
# -----------------------------
library(tidyr)
library(dplyr)
library(rstatix)
library(coin)
# Wide format: one row per site, one column per species
rel_wide <- rel_abund %>%
  select(location_id, species_canon, rel_percent) %>%
  pivot_wider(names_from = species_canon, values_from = rel_percent)

# Quick look
print(rel_wide)

# Paired Wilcoxon signed-rank test (recommended)
wilcox_res <- wilcox.test(
  rel_wide$`Labeo victorianus`,
  rel_wide$`Labeobarbus altianalis`,
  paired = TRUE,
  exact = FALSE
)

# Effect size (rank-biserial) for paired Wilcoxon
eff_res <- rel_abund %>%
  wilcox_effsize(rel_percent ~ species_canon, paired = TRUE)

cat("\nPaired Wilcoxon signed-rank test (across sites)\n")
cat("V =", wilcox_res$statistic, " p =", wilcox_res$p.value, "\n")
print(eff_res)

# Optional: paired t-test (only if you want a parametric comparison too)
ttest_res <- t.test(
  rel_wide$`Labeo victorianus`,
  rel_wide$`Labeobarbus altianalis`,
  paired = TRUE
)

cat("\nPaired t-test (across sites) [optional]\n")
print(ttest_res)


# --- Plot: stacked relative abundance ---
ggplot(rel_abund, aes(x = location_id, y = rel_percent, fill = species_canon)) +
  geom_col(width = 0.7, color = "black") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(
    name = "Species",
    values = pal_species,
    breaks = species_order
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x   = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )
# --- Custom colours (same as before) ---
pal_species <- c(
  "Labeo victorianus"      = "#1F78B4",  # deep blue
  "Labeobarbus altianalis" = "#E41A1C"   # vivid red
)

# --- Order legend by higher overall mean relative abundance ---
species_order <- rel_abund %>%
  dplyr::group_by(species_canon) %>%
  dplyr::summarise(mean_rel = mean(rel_percent, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mean_rel)) %>%
  dplyr::pull(species_canon)

rel_abund_plot <- rel_abund %>%
  dplyr::mutate(species_canon = factor(species_canon, levels = species_order))

# --- Grouped bar plot (two bars per site) with Y-axis 0–75% ---
ggplot(rel_abund_plot,
       aes(x = location_id, y = rel_percent, fill = species_canon)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "black") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 75),              # Start at 0%, end at 75%
    expand = expansion(mult = c(0, 0))
  ) +
  scale_fill_manual(
    name   = "Species",
    values = c(
      "Labeo victorianus"      = "#1F78B4",  # deep blue
      "Labeobarbus altianalis" = "#E41A1C"   # vivid red
    ),
    breaks = species_order
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.text.x   = element_text(size = 11),
    axis.title.y  = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
    )


#######add the statistical annotation to the plot #########
# ---- Build Wilcoxon label (make sure wilcox_res exists) ----
p_label <- paste0(
  "Paired Wilcoxon\n",
  "V = ", wilcox_res$statistic, "\n",
  "p = ", formatC(wilcox_res$p.value, format = "f", digits = 3)
)

label_size <- 4  # <-- adjust label text size here

# ---- Choose a safe top-right position in data coordinates ----
x_pos <- length(unique(rel_abund_plot$location_id)) + 0.5   # a bit to the right of last site
y_pos <- 75                                                # top of your y-limit

ggplot(rel_abund_plot,
       aes(x = location_id, y = rel_percent, fill = species_canon)) +
  geom_col(position = position_dodge(width = 0.72),
           width = 0.62, color = "black") +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 75),
    expand = expansion(mult = c(0, 0))
  ) +
  # give a little room on the right for the annotation
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(5.5, 40, 5.5, 5.5)) +
  scale_fill_manual(
    name   = "Species",
    values = c(
      "Labeo victorianus"      = "#1F78B4",
      "Labeobarbus altianalis" = "#E41A1C"
    ),
    breaks = species_order
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    axis.text.x   = element_text(size = 11),
    axis.title.y  = element_text(size = 12),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11)
  ) +
  # ---- Add annotation (top-right) ----
annotate(
  "text",
  x = Inf, y = Inf,
  label = "Paired Wilcoxon: V = 1, p = 0.059",
  hjust = 1.05, vjust = 1.2,
  size = 3.5
)



################################################################################
remove(list = ls())

library(tidyverse)
library(readr)
library(janitor)
library(stringr)
library(forcats)
library(scales)

# --- Load & clean ---
dat <- readr::read_csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    # standardize species names (lowercase, trim)
    fish_species = str_squish(str_to_lower(fish_species)),
    # keep only M2–M9 and years 2021–2022
    location_id  = factor(location_id, levels = paste0("M", 2:9))
  ) %>%
  filter(
    !is.na(location_id), location_id %in% levels(location_id),
    sampling_year %in% c(2021, 2022),
    !is.na(fish_species), fish_species != ""
  )

# Optional: nice printing names (Title Case) for legend/plot
nice_species <- function(x) {
  x %>%
    str_replace_all("_", " ") %>%
    str_replace_all("\\s+", " ") %>%
    str_trim() %>%
    str_to_sentence()
}

# Canonicalize focal species display names exactly as requested
to_display <- function(sp) {
  sp_lc <- str_squish(str_to_lower(sp))
  case_when(
    str_detect(sp_lc, "^labeo\\s+victorianus$") ~ "Labeo victorianus",
    str_detect(sp_lc, "^labeobarbus\\s+alti?analis$") ~ "Labeobarbus altianalis",
    TRUE ~ nice_species(sp_lc)
  )
}

dat <- dat %>%
  mutate(species_label = to_display(fish_species))

# ----- Count per site × species, then convert to relative abundance (%) -----
site_species_counts <- dat %>%
  count(location_id, species_label, name = "n_ind")

site_totals <- site_species_counts %>%
  group_by(location_id) %>%
  summarise(n_total = sum(n_ind), .groups = "drop")

rel_abund_all <- site_species_counts %>%
  left_join(site_totals, by = "location_id") %>%
  mutate(rel_percent = if_else(n_total > 0, 100 * n_ind / n_total, 0)) %>%
  ungroup()

# ----- Order species by overall abundance (helps legend readability) -----
species_order <- rel_abund_all %>%
  group_by(species_label) %>%
  summarise(overall = sum(rel_percent), .groups = "drop") %>%
  arrange(desc(overall)) %>%
  pull(species_label)

rel_abund_all <- rel_abund_all %>%
  mutate(species_label = factor(species_label, levels = species_order))

# ----- Build a color palette for ALL species, locking colors for the two focal spp. -----
all_species <- levels(rel_abund_all$species_label)

# start with a generic distinct palette (avoid blue hues since you asked earlier)
base_cols <- hue_pal(h = c(10, 350), c = 100, l = 60)(length(all_species))
pal <- setNames(base_cols, all_species)

# overwrite the two focal species with your fixed colors
pal["Labeo victorianus"]      <- "#1F78B4"  # deep blue
pal["Labeobarbus altianalis"] <- "#E41A1C"  #vivid red

# If any of those two labels are absent in the data, ignore the warning
pal <- pal[!is.na(names(pal))]

# ----- Plot: Stacked bars (relative abundance %) for all species across M2–M9 -----
ggplot(rel_abund_all, aes(x = location_id, y = rel_percent, fill = species_label)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.25) +
  scale_y_continuous(
    labels = percent_format(scale = 1),
    limits = c(0, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = pal, name = "Species") +
  labs(
    title = "Relative Abundance (2021–2022)",
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x   = element_text(size = 11),
    legend.title  = element_text(face = "bold"),
    legend.key.height = unit(0.45, "cm")
  )
##################################################################################
############################################################
# Relative abundance (%) of Labeo victorianus vs Labeobarbus altianalis
# Sites M4–M9, Years 2021–2022
# Mean ± SE across years (SE = 0 if only 1 year available)
# Grouped bar plot with black SE error bars + paired Wilcoxon across sites
############################################################

remove(list = ls())

# ---- Packages ----
libs <- c("tidyverse", "readr", "janitor", "stringr", "rstatix", "scales")
to_install <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(libs, library, character.only = TRUE))

# ---- Data link (Google Sheet CSV) ----
data_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv"

# ---- Load data ----
dat <- readr::read_csv(data_url, show_col_types = FALSE) %>%
  janitor::clean_names()

# ---- Canonicalize/clean key fields ----
dat <- dat %>%
  mutate(
    location_id   = as.character(location_id),
    sampling_year = suppressWarnings(as.integer(sampling_year)),
    fish_species  = stringr::str_squish(stringr::str_to_lower(fish_species))
  )

# ---- Keep only sites + years of interest ----
sites_keep <- paste0("M", 4:9)
years_keep <- c(2021, 2022)

dat_sub <- dat %>%
  filter(
    location_id %in% sites_keep,
    sampling_year %in% years_keep,
    !is.na(fish_species), fish_species != ""
  )

# ---- Map species names to the two focal taxa ----
dat_sub <- dat_sub %>%
  mutate(
    species_canon = case_when(
      str_detect(fish_species, "^labeo\\s+victorianus$") ~ "Labeo victorianus",
      str_detect(fish_species, "^labeobarbus\\s+alti?analis$") ~ "Labeobarbus altianalis",
      TRUE ~ NA_character_
    ),
    location_id = factor(location_id, levels = sites_keep),
    sampling_year = factor(sampling_year)
  )

# ============================================================
# STEP 1: total catch per site-year (ALL species)
# ============================================================
tot_site_year <- dat_sub %>%
  count(location_id, sampling_year, name = "n_total") %>%
  tidyr::complete(
    location_id = factor(sites_keep, levels = sites_keep),
    sampling_year = factor(years_keep),
    fill = list(n_total = 0)
  )

# ============================================================
# STEP 2: focal counts per site-year (two focal species only)
# ============================================================
focal_site_year <- dat_sub %>%
  filter(!is.na(species_canon)) %>%
  count(location_id, sampling_year, fish_species = species_canon, name = "n_focal") %>%
  tidyr::complete(
    location_id = factor(sites_keep, levels = sites_keep),
    sampling_year = factor(years_keep),
    fish_species = factor(
      c("Labeo victorianus", "Labeobarbus altianalis"),
      levels = c("Labeobarbus altianalis", "Labeo victorianus")
    ),
    fill = list(n_focal = 0)
  )

# ============================================================
# STEP 3: relative abundance (%) per site-year-species
# ============================================================
rel_site_year <- focal_site_year %>%
  left_join(tot_site_year, by = c("location_id", "sampling_year")) %>%
  mutate(
    n_total = replace_na(n_total, 0),
    rel_percent = if_else(n_total > 0, 100 * n_focal / n_total, 0)
  )

# ============================================================
# STEP 4: mean ± SE across years for each site × species
#   - IMPORTANT FIX: if only 1 year exists, set SE = 0
# ============================================================
rel_summary <- rel_site_year %>%
  group_by(location_id, fish_species) %>%
  summarise(
    n_rep    = sum(!is.na(rel_percent)),
    mean_rel = mean(rel_percent, na.rm = TRUE),
    se_rel   = if_else(n_rep > 1, sd(rel_percent, na.rm = TRUE) / sqrt(n_rep), 0),
    .groups = "drop"
  )

# Optional: check why an error bar might be missing (M4 example)
rel_site_year %>%
  filter(location_id == "M4", fish_species == "Labeobarbus altianalis") %>%
  arrange(sampling_year) %>%
  print(n = Inf)

# ============================================================
# STEP 5: Paired Wilcoxon across sites
#   - use site-level mean relative abundance (paired by site)
# ============================================================
rel_wide <- rel_summary %>%
  select(location_id, fish_species, mean_rel) %>%
  tidyr::pivot_wider(names_from = fish_species, values_from = mean_rel)

wilcox_res <- wilcox.test(
  rel_wide$`Labeo victorianus`,
  rel_wide$`Labeobarbus altianalis`,
  paired = TRUE,
  exact = FALSE
)

p_label <- paste0(
  "Paired Wilcoxon: V = ", unname(wilcox_res$statistic),
  ", p = ", formatC(wilcox_res$p.value, format = "f", digits = 3)
)

# ============================================================
# STEP 6: Plot (grouped bars) with black SE error bars
#   - Labeobarbus altianalis = red
#   - Labeo victorianus      = blue
# ============================================================
pal_species <- c(
  "Labeobarbus altianalis" = "#E41A1C",
  "Labeo victorianus"      = "#1F78B4"
)

rel_summary <- rel_summary %>%
  mutate(
    fish_species = factor(fish_species, levels = c("Labeobarbus altianalis", "Labeo victorianus"))
  )

p_rel <- ggplot(rel_summary, aes(x = location_id, y = mean_rel, fill = fish_species)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = pmax(mean_rel - se_rel, 0), ymax = mean_rel + se_rel),
    position = position_dodge(width = 0.72),
    width = 0.18,
    linewidth = 0.6,
    color = "black"
  ) +
  scale_fill_manual(name = "Species", values = pal_species) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 75),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Relative Abundance (%)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = p_label,
    hjust = 1.05, vjust = 1.2,
    size = 3.8
  )

print(p_rel)

# ============================================================
# OPTIONAL: Save publication-quality figure (JPEG)
# ============================================================
dir.create("Figures", showWarnings = FALSE)

ggsave(
  filename = "Figures/RelativeAbundance_LV_vs_LA_M4-M9_2021-2022.jpg",
  plot = p_rel,
  width = 8,
  height = 6,
  dpi = 300
)
##############################################################################
############################################################
# Relative abundance (%) of Labeo victorianus vs Labeobarbus altianalis
# Sites M4–M9, Years 2021–2022
# Mean ± SE across years (SE = SD/sqrt(n_years_with_data); SE = 0 if n=1)
# Grouped bar plot with black SE error bars + paired Wilcoxon across sites
############################################################

remove(list = ls())

# ---- Packages ----
libs <- c("tidyverse", "readr", "janitor", "stringr", "rstatix", "scales")
to_install <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)
invisible(lapply(libs, library, character.only = TRUE))

# ---- Data link (Google Sheet CSV) ----
data_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv"

# ---- Load + clean ----
dat <- readr::read_csv(data_url, show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    location_id   = as.character(location_id),
    sampling_year = suppressWarnings(as.integer(sampling_year)),
    fish_species  = stringr::str_squish(stringr::str_to_lower(fish_species))
  ) %>%
  filter(!is.na(location_id), !is.na(sampling_year), !is.na(fish_species), fish_species != "")

# ---- Keep only sites + years of interest ----
sites_keep <- paste0("M", 4:9)
years_keep <- c(2021, 2022)

dat_sub <- dat %>%
  filter(
    location_id %in% sites_keep,
    sampling_year %in% years_keep
  ) %>%
  mutate(
    location_id = factor(location_id, levels = sites_keep),
    sampling_year = factor(sampling_year, levels = years_keep)
  )

# ---- Canonicalize species names to 2 focal taxa ----
dat_sub <- dat_sub %>%
  mutate(
    species_canon = case_when(
      str_detect(fish_species, "^labeo\\s+victorianus$") ~ "Labeo victorianus",
      str_detect(fish_species, "^labeobarbus\\s+alti?analis$") ~ "Labeobarbus altianalis",
      TRUE ~ NA_character_
    )
  )

# ============================================================
# STEP 1: total catch per site-year (ALL species)
# ============================================================
tot_site_year <- dat_sub %>%
  count(location_id, sampling_year, name = "n_total") %>%
  tidyr::complete(
    location_id   = factor(sites_keep, levels = sites_keep),
    sampling_year = factor(years_keep, levels = years_keep),
    fill = list(n_total = 0)
  )

# ============================================================
# STEP 2: focal counts per site-year (ONLY the two species)
# ============================================================
focal_site_year <- dat_sub %>%
  filter(!is.na(species_canon)) %>%
  count(location_id, sampling_year, fish_species = species_canon, name = "n_focal") %>%
  tidyr::complete(
    location_id   = factor(sites_keep, levels = sites_keep),
    sampling_year = factor(years_keep, levels = years_keep),
    fish_species  = factor(c("Labeobarbus altianalis", "Labeo victorianus"),
                           levels = c("Labeobarbus altianalis", "Labeo victorianus")),
    fill = list(n_focal = 0)
  )

# ============================================================
# STEP 3: relative abundance (%) per site-year-species
# ============================================================
rel_site_year <- focal_site_year %>%
  left_join(tot_site_year, by = c("location_id", "sampling_year")) %>%
  mutate(
    n_total = replace_na(n_total, 0),
    rel_percent = if_else(n_total > 0, 100 * n_focal / n_total, 0)
  )

# ============================================================
# STEP 4: mean ± SE across years for each site × species
#   SE is computed across years (2021, 2022) within a site
# ============================================================
rel_summary <- rel_site_year %>%
  group_by(location_id, fish_species) %>%
  summarise(
    n_years  = sum(!is.na(rel_percent)),                   # will be 2 after complete()
    mean_rel = mean(rel_percent, na.rm = TRUE),
    sd_rel   = sd(rel_percent, na.rm = TRUE),
    se_rel   = if_else(n_years > 1, sd_rel / sqrt(n_years), 0),
    .groups = "drop"
  )

# (Optional) diagnose “missing” error bars: show the year-level values
# print(rel_site_year %>% filter(location_id == "M4", fish_species == "Labeobarbus altianalis"))

# ============================================================
# STEP 5: paired Wilcoxon across sites using site means (paired by site)
# ============================================================
rel_wide <- rel_summary %>%
  select(location_id, fish_species, mean_rel) %>%
  tidyr::pivot_wider(names_from = fish_species, values_from = mean_rel)

wilcox_res <- wilcox.test(
  rel_wide$`Labeo victorianus`,
  rel_wide$`Labeobarbus altianalis`,
  paired = TRUE,
  exact = FALSE
)

p_label <- paste0(
  "Paired Wilcoxon: V = ", unname(wilcox_res$statistic),
  ", p = ", formatC(wilcox_res$p.value, format = "f", digits = 3)
)

# ============================================================
# STEP 6: Plot (grouped bars) with black SE error bars
# ============================================================
pal_species <- c(
  "Labeobarbus altianalis" = "#E41A1C",  # red
  "Labeo victorianus"      = "#1F78B4"   # blue
)

rel_summary <- rel_summary %>%
  mutate(fish_species = factor(fish_species,
                               levels = c("Labeobarbus altianalis", "Labeo victorianus")))

p_rel <- ggplot(rel_summary, aes(x = location_id, y = mean_rel, fill = fish_species)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    color = "black"
  ) +
  geom_errorbar(
    aes(ymin = pmax(mean_rel - se_rel, 0), ymax = mean_rel + se_rel),
    position = position_dodge(width = 0.72),
    width = 0.18,
    linewidth = 0.6,
    color = "black"
  ) +
  scale_fill_manual(name = "Species", values = pal_species) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    limits = c(0, 80),          # <-- changed from 75 to 80
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "",
    x = "Sampling Site",
    y = "Mean relative abundance (%) 2021–2022 "
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    plot.margin = margin(5.5, 18, 5.5, 5.5)
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = p_label,
    hjust = 1.05, vjust = 1.2,
    size = 3.8
  )

print(p_rel)

# ============================================================
# STEP 7: Save publication-quality JPEG
# ============================================================
dir.create("Figures", showWarnings = FALSE)

ggsave(
  filename = "Figures/RelativeAbundance_LV_vs_LA_M4-M9_2021-2022.jpg",
  plot = p_rel,
  width = 8,
  height = 6,
  dpi = 300
)
###########################################
