# =========================
# Fish gut proportions plot
# =========================
remove(list = ls())

# Setup renv (optional)
renv::restore()

# ---- Packages ----
library(tidyverse)  # includes dplyr, ggplot2, stringr, etc.
library(readr)
library(janitor)
library(scales)

# ---- Data source (CSV via readr) ----
csv_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=2033839056&single=true&output=csv"

# ---- Allowed food codes + nice labels ----
allowed_codes <- c("FH","NAB","FA","Diatom","AD","LD","Grass","Insect","VM","IO")
food_labels <- c(
  FH       = "Fungal hyphae",
  NAB      = "Non-algal biofilm",
  FA       = "Filamentous algae",
  Diatom   = "Diatoms",
  AD       = "Amorphous detritus",
  LD       = "Leaf detritus",
  Grass    = "Grass",
  Insect   = "Insect",
  IO       = "Inorganic",
  VM       = "Vertebrate matter"
)

# ---- Load & basic clean ----
raw <- readr::read_csv(csv_url, show_col_types = FALSE) |>
  clean_names()   # -> site_name, site_code, fish_species, fish_replicant, food_item, picture, food_type, area_mm2, picture_number

# Safety check for required columns
stopifnot(all(c("site_code","fish_species","food_type","area_mm2") %in% names(raw)))

# ---- Food_type cleaning ----
dat <- raw |>
  mutate(
    food_type = stringr::str_squish(stringr::str_trim(food_type)),
    food_type = stringr::str_replace(food_type, "#\\d+$", ""),
    food_type_lower = tolower(food_type),
    food_type = dplyr::recode(
      food_type_lower,
      "fh" = "FH",
      "nab" = "NAB",
      "fa" = "FA",
      "diatom" = "Diatom",
      "diatoms" = "Diatom",
      "ad" = "AD",
      "ld" = "LD",
      "grass" = "Grass",
      "insect" = "Insect",
      "insects" = "Insect",
      "io" = "IO",
      "vm" = "VM",
      .default = food_type
    )
  ) |>
  select(-food_type_lower)

# ---- Filter allowed food codes ----
dat <- dat |>
  filter(food_type %in% allowed_codes)

# ---- Identify two dominant species by total gut area ----
top_species <- dat |>
  group_by(fish_species) |>
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total_area)) |>
  slice_head(n = 2) |>
  pull(fish_species)

message("Selected dominant species: ", paste(top_species, collapse = " | "))

# If you want to fix them manually:
# top_species <- c("Labeobarbus altianalis", "Labeo victorianus")

dat2 <- dat |>
  filter(fish_species %in% top_species)

# ---- Aggregate & compute proportions ----
agg <- dat2 |>
  group_by(site_code, fish_species, food_type) |>
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop")

props <- agg |>
  group_by(site_code, fish_species) |>
  mutate(
    denom = sum(total_area, na.rm = TRUE),
    proportion = ifelse(denom > 0, 100 * total_area / denom, 0)
  ) |>
  ungroup() |>
  mutate(
    food_type_label = recode(food_type, !!!food_labels, .default = food_type),
    food_type_label = factor(food_type_label, levels = unname(food_labels))
  )

# ---- Optional: site order (alphabetical by code) ----
site_order <- raw |> distinct(site_code) |> arrange(site_code) |> pull(site_code)
props <- props |> mutate(site_code = factor(site_code, levels = site_order))

# ---- Plot using site_code ----
p <- ggplot(props, aes(x = site_code, y = proportion, fill = food_type_label)) +
  geom_col() +
  facet_wrap(~ fish_species, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  labs(
    title = "Gut content composition (% area) by species and site",
    x = "Site code",
    y = "Proportion of gut content (%)",
    fill = "Food type"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "right",
    strip.text = element_text(face = "bold")
  )

print(p)

#########################################################################################
# =========================
# Fish gut proportions plot
# Distinct colors + ordered legend + fixed facet order
# =========================
remove(list = ls())

# ---- Packages ----
library(tidyverse)
library(readr)
library(janitor)
library(scales)
library(RColorBrewer)

# ---- Data source ----
csv_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=2033839056&single=true&output=csv"

# ---- Allowed food codes + labels ----
allowed_codes <- c("FH","NAB","FA","Diatom","AD","LD","Grass","Insect","VM","IO")
food_labels <- c(
  FH       = "Fungal hyphae",
  NAB      = "Non-algal biofilm",
  FA       = "Filamentous algae",
  Diatom   = "Diatom",
  AD       = "Amorphous detritus",
  LD       = "Leaf detritus",
  Grass    = "Grass",
  Insect   = "Insect",
  IO       = "Inorganic",
  VM       = "Vertebrate matter"
)

# ---- Load & clean ----
raw <- read_csv(csv_url, show_col_types = FALSE) |>
  clean_names()

stopifnot(all(c("site_code","fish_species","food_type","area_mm2") %in% names(raw)))

# ---- Clean food_type ----
dat <- raw |>
  mutate(
    food_type = str_squish(str_trim(food_type)),
    food_type = str_replace(food_type, "#\\d+$", ""),
    food_type_lower = tolower(food_type),
    food_type = recode(
      food_type_lower,
      "fh" = "FH",
      "nab" = "NAB",
      "fa" = "FA",
      "diatom" = "Diatom",
      "diatoms" = "Diatom",
      "ad" = "AD",
      "ld" = "LD",
      "grass" = "Grass",
      "insect" = "Insect",
      "insects" = "Insect",
      "io" = "IO",
      "vm" = "VM",
      .default = food_type
    )
  ) |>
  select(-food_type_lower)

# ---- Filter valid food codes ----
dat <- dat |> filter(food_type %in% allowed_codes)

# ---- Pick top 2 dominant species ----
top_species <- dat |>
  group_by(fish_species) |>
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(total_area)) |>
  slice_head(n = 2) |>
  pull(fish_species)

message("Selected dominant species: ", paste(top_species, collapse = " | "))

dat2 <- dat |> filter(fish_species %in% top_species)

# ---- Aggregate and calculate proportions ----
agg <- dat2 |>
  group_by(site_code, fish_species, food_type) |>
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop")

props <- agg |>
  group_by(site_code, fish_species) |>
  mutate(
    denom = sum(total_area, na.rm = TRUE),
    proportion = ifelse(denom > 0, 100 * total_area / denom, 0)
  ) |>
  ungroup() |>
  mutate(food_type_label = recode(food_type, !!!food_labels, .default = food_type))

# ---- Order legend by overall mean proportion ----
legend_order <- props |>
  group_by(food_type_label) |>
  summarise(mean_prop = mean(proportion, na.rm = TRUE)) |>
  arrange(desc(mean_prop)) |>
  pull(food_type_label)

props <- props |>
  mutate(food_type_label = factor(food_type_label, levels = legend_order))

# ---- Set site order ----
site_order <- raw |> distinct(site_code) |> arrange(site_code) |> pull(site_code)
props <- props |> mutate(site_code = factor(site_code, levels = site_order))

# ---- Fix facet order (species order) ----
props <- props |> 
  mutate(fish_species = factor(fish_species, 
                               levels = c("Labeobarbus altianalis", "Labeo victorianus")))

# ---- Define a vivid, colorblind-friendly palette ----
palette_colors <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#a65628", "#8dd3c7"
)

# ---- Plot ----
p <- ggplot(props, aes(x = site_code, y = proportion, fill = food_type_label)) +
  geom_col(color = "black", width = 0.8) +
  facet_wrap(~ fish_species, nrow = 1, scales = "free_x") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0.05))
  ) +
  scale_fill_manual(
    values = palette_colors,
    guide = guide_legend(ncol = 1)
  ) +
  labs(
    title = "",
    x = "Site",
    y = "Proportion of gut content (%)",
    fill = "Food type"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 11),
    legend.text = element_text(size = 11),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 12)
  )

print(p)

#############################################################
###########################################################

# Clear workspace
remove(list = ls())

# If you're using renv for project reproducibility:
# renv::restore()


###########################################################
# 1. LOAD PACKAGES
###########################################################

library(tidyverse)   # dplyr, ggplot2, tidyr, etc.
library(readr)       # read_csv
library(janitor)     # clean_names
library(vegan)       # vegdist, adonis2, betadisper, metaMDS
library(rstatix)     # anova_test, kruskal_test, etc.


###########################################################
# 2. LOAD AND CLEAN RAW GUT CONTENT DATA
###########################################################

# ---- Data source (CSV via readr) ----
csv_url <- "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=2033839056&single=true&output=csv"

# ---- Allowed food codes ----
allowed_codes <- c("FH", "NAB", "FA", "Diatom", "AD", "LD",
                   "Grass", "Insect", "VM", "IO")

# ---- Load and basic clean ----
raw <- read_csv(csv_url, show_col_types = FALSE) %>%
  clean_names()   # -> site_code, fish_species, fish_replicant, food_type, area_mm2, etc.

# Safety check for required columns
stopifnot(all(c("site_code", "fish_species", "fish_replicant",
                "food_type", "area_mm2") %in% names(raw)))

# ---- Food_type cleaning and standardisation ----
dat <- raw %>%
  mutate(
    food_type = stringr::str_squish(stringr::str_trim(food_type)),
    # remove #number suffix (e.g. "NAB#1" -> "NAB")
    food_type = stringr::str_replace(food_type, "#\\d+$", ""),
    food_type_lower = tolower(food_type),
    food_type = dplyr::recode(
      food_type_lower,
      "fh"       = "FH",
      "nab"      = "NAB",
      "fa"       = "FA",
      "diatom"   = "Diatom",
      "diatoms"  = "Diatom",
      "ad"       = "AD",
      "ld"       = "LD",
      "grass"    = "Grass",
      "insect"   = "Insect",
      "insects"  = "Insect",
      "io"       = "IO",
      "vm"       = "VM",
      .default   = food_type  # fallback to cleaned original
    )
  ) %>%
  select(-food_type_lower) %>%
  # keep only desired food categories
  filter(food_type %in% allowed_codes)

# OPTIONAL: restrict to the two dominant species only
top_species <- dat %>%
  group_by(fish_species) %>%
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_area)) %>%
  slice_head(n = 2) %>%
  pull(fish_species)

message("Selected dominant species: ", paste(top_species, collapse = " | "))

dat <- dat %>%
  filter(fish_species %in% top_species)


###########################################################
# 3. COMPUTE REPLICATE-LEVEL GUT COMPOSITION (% AREA)
###########################################################

# For each fish replicate, compute total area per food_type,
# then convert to percentage of total gut area for that fish.

rep_prop <- dat %>%
  group_by(site_code, fish_species, fish_replicant, food_type) %>%
  summarise(total_area = sum(area_mm2, na.rm = TRUE), .groups = "drop") %>%
  group_by(site_code, fish_species, fish_replicant) %>%
  mutate(
    denom = sum(total_area, na.rm = TRUE),
    proportion = 100 * total_area / denom
  ) %>%
  ungroup()

# Inspect the replicate-level proportions
print(rep_prop)


###########################################################
# 4. PREPARE WIDE MATRIX FOR MULTIVARIATE ANALYSES
###########################################################

# Rows = individual fish (site × species × replicate),
# Columns = food types, Values = % area

rep_wide <- rep_prop %>%
  select(site_code, fish_species, fish_replicant, food_type, proportion) %>%
  tidyr::pivot_wider(
    names_from  = food_type,
    values_from = proportion,
    values_fill = 0
  )

# Metadata and matrix for vegan functions
meta <- rep_wide %>%
  select(site_code, fish_species, fish_replicant)

diet_matrix <- rep_wide %>%
  select(-site_code, -fish_species, -fish_replicant)

# Quick check
head(rep_wide)


###########################################################
# 5. PERMANOVA: EFFECT OF SPECIES AND SITE ON GUT COMPOSITION
###########################################################

# Compute Bray–Curtis distance on diet composition
dist_bc <- vegdist(diet_matrix, method = "bray")

# PERMANOVA model:
#   - main effects: fish_species, site_code
#   - interaction: fish_species:site_code
adonis_result <- adonis2(
  dist_bc ~ fish_species * site_code,
  data        = meta,
  permutations = 999
)

print(adonis_result)

adonis2(dist_bc ~ fish_species + site_code + fish_species:site_code,
        data = meta, by = "margin")


###########################################################
# 6. CHECK HOMOGENEITY OF MULTIVARIATE DISPERSION (BETADISPER)
###########################################################

# This is analogous to checking homogeneity of variances in ANOVA

# (a) Dispersion among species
disp_species <- betadisper(dist_bc, group = meta$fish_species)
anova(disp_species)

# (b) Dispersion among sites
disp_site <- betadisper(dist_bc, group = meta$site_code)
anova(disp_site)


###########################################################
# 8. UNIVARIATE TESTS (ANOVA / KRUSKAL–WALLIS) FOR FOOD TYPES
###########################################################

# Here we work on rep_prop (replicate-level % area).
# Example 1: test if % Amorphous detritus (AD) differs between species
rep_prop %>%
  filter(food_type == "AD") %>%
  anova_test(proportion ~ fish_species)

# Example 2: test if % Amorphous detritus (AD) differs among sites
rep_prop %>%
  filter(food_type == "AD") %>%
  kruskal_test(proportion ~ site_code)

# Example 3: loop over food types and test species differences (ANOVA)
# NOTE: with small n you may prefer non-parametric tests.
species_anova_results <- rep_prop %>%
  group_by(food_type) %>%
  anova_test(proportion ~ fish_species)

print(species_anova_results)

# Example 4: loop over food types and test site differences (Kruskal–Wallis)
site_kw_results <- rep_prop %>%
  group_by(food_type) %>%
  kruskal_test(proportion ~ site_code)

print(site_kw_results)

