# Compare the historical and current fish community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.
# =====================================================================
# PERMANOVA + NMDS for changes in fish community
# =====================================================================

# =====================================================================
# PERMANOVA + NMDS for changes in fish community
# =====================================================================
remove(list = ls())
renv::restore()  # uncomment if you use renv

library(tidyverse)
library(readr)
library(janitor)
library(stringr)
library(vegan)
library(ggpubr)

# -----------------------------
# 1) Load LONG data and prepare
# -----------------------------
fish_long <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    fish_species   = str_squish(fish_species),
    # be robust if your file sometimes calls it `reach`
    river_reach    = coalesce(river_reach),
    abundance      = 1
  ) %>%
  # keep essential fields
  filter(!is.na(location_id), !is.na(sampling_year), !is.na(sampling_month), !is.na(fish_species))

# -----------------------------
# 2) Summarize and LONG -> WIDE
# -----------------------------
# Define sampling unit = site × year × month
fish_sum <- fish_long %>%
  group_by(location_id, river_reach, sampling_year, sampling_month, fish_species) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop")

fish_wide <- fish_sum %>%
  pivot_wider(
    names_from  = fish_species,
    values_from = total_abundance,
    values_fill = 0
  ) %>%
  # Create a unit_id that will be the row id
  mutate(unit_id = paste(location_id, sampling_year, sampling_month, sep = "_"))

# -----------------------------
# 3) Community matrix & metadata
# -----------------------------
# Community matrix (species only), rows = unit_id
comm <- fish_wide %>%
  select(-location_id, -river_reach, -sampling_year, -sampling_month, -unit_id)

rownames(comm) <- fish_wide$unit_id

# Drop rows with all zeros (if any)
#comm <- comm[rowSums(comm) > 0, , drop = FALSE]

# Matching metadata in the same row order as `comm`
meta <- fish_wide %>%
  select(location_id, river_reach, sampling_year, sampling_month, unit_id) %>%
  filter(unit_id %in% rownames(comm)) %>%
  arrange(match(unit_id, rownames(comm))) %>%
  mutate(
    river_reach    = factor(replace_na(river_reach, "Unknown")),
    sampling_year  = factor(sampling_year),
    sampling_month = factor(sampling_month)
  )

stopifnot(identical(meta$unit_id, rownames(comm)))

# -----------------------------
# 4) Distance & dispersion test
# -----------------------------
set.seed(123)
dist_bray <- vegdist(comm, method = "bray")

# Check homogeneity of dispersion (by sampling_year here; test others as needed)
bd_year   <- betadisper(dist_bray, meta$sampling_year)
disp_test <- permutest(bd_year, permutations = 999)
print(disp_test)  # if significant, interpret PERMANOVA with caution

# ---------------
# 5) PERMANOVA
# ---------------
# IMPORTANT: use `data = meta` (same rows as `dist_bray`) and refer to meta columns in the formula
perm <- adonis2(
  dist_bray ~ sampling_year + sampling_month + river_reach,
  data = meta,
  permutations = 999,
  by = "margin",           # marginal (type II) tests
  method = "bray"
)
print(perm)

# -------------------------------------
# 6) NMDS (optional) — visualize groups
# -------------------------------------
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100, autotransform = FALSE)
cat("NMDS stress:", nmds$stress, "\n")

nmds_scores <- scores(nmds, display = "sites") %>% as.data.frame()
nmds_scores$unit_id <- rownames(nmds_scores)
nmds_df <- nmds_scores %>% left_join(meta, by = "unit_id")

p_nmds <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = sampling_year, shape = river_reach)) +
  geom_point(size = 2.6, alpha = 0.9) +
  stat_ellipse(aes(group = sampling_year), linetype = "dashed", linewidth = 0.5, alpha = 0.6) +
  labs(title = paste0("NMDS (Bray–Curtis) — Stress = ", round(nmds$stress, 3)),
       color = "Year", shape = "Reach") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5))
print(p_nmds)

