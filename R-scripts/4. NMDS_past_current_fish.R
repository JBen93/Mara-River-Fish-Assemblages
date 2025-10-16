# Compare the historical and current fish community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.
# =====================================================================
# Fish Community Analysis — NMDS, PERMANOVA, PCoA, and DCA
# =====================================================================

# Clean workspace
remove(list = ls())

# Restore environment (if using renv)
renv::restore()

# Load libraries
library(tidyverse)
library(vegan)
library(janitor)
library(readr)
library(ggpubr)
library(ggplot2)

# ---------------------------------------------------------------------
# 1. Load and Prepare Data
# ---------------------------------------------------------------------
fish_long <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv",
  show_col_types = FALSE
) %>%
  clean_names() %>%
  mutate(
    fish_species = str_squish(fish_species),
    abundance = 1
  ) %>%
  filter(!is.na(location_id), !is.na(sampling_year), !is.na(fish_species))

# ---------------------------------------------------------------------
# 2. Summarize and Convert from LONG → WIDE
# ---------------------------------------------------------------------
fish_sum <- fish_long %>%
  group_by(location_id, river_reach, sampling_year, sampling_month, fish_species) %>%
  summarise(total_abundance = sum(abundance), .groups = "drop")

fish_wide <- fish_sum %>%
  pivot_wider(
    names_from  = fish_species,
    values_from = total_abundance,
    values_fill = 0
  ) %>%
  mutate(unit_id = paste(location_id, sampling_year, sep = "_"))

# ---------------------------------------------------------------------
# 3. Community Matrix and Metadata
# ---------------------------------------------------------------------
comm <- fish_wide %>%
  select(-location_id, -river_reach, -sampling_year, -sampling_month, -unit_id) %>%
  as.data.frame()

meta <- fish_wide %>%
  select(unit_id, location_id, river_reach, sampling_year, sampling_month) %>%
  mutate(
    sampling_year  = as.factor(sampling_year),
    river_reach    = as.factor(replace_na(river_reach, "Unknown")),
    sampling_month = as.factor(sampling_month)
  )

rownames(comm) <- meta$unit_id
comm <- comm[rowSums(comm) > 0, , drop = FALSE]

meta_aligned <- meta %>%
  filter(unit_id %in% rownames(comm)) %>%
  arrange(match(unit_id, rownames(comm)))

stopifnot(identical(meta_aligned$unit_id, rownames(comm)))

# ---------------------------------------------------------------------
# 4. Distance Matrix + PERMANOVA
# ---------------------------------------------------------------------
set.seed(123)
dist_bray <- vegdist(comm, method = "bray")

# Test for homogeneity of dispersion (important before PERMANOVA)
bd_year <- betadisper(dist_bray, meta_aligned$sampling_year)
disp_test <- permutest(bd_year, permutations = 999)
print(disp_test)

# PERMANOVA test
perm <- adonis2(
  dist_bray ~ sampling_year + sampling_month + river_reach,
  data = meta_aligned,
  permutations = 999,
  by = "margin"
)
print(perm)

# ---------------------------------------------------------------------
# 5. NMDS Ordination
# ---------------------------------------------------------------------
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

scores_nmds <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("unit_id")

nmds_df <- inner_join(scores_nmds, meta_aligned, by = "unit_id")

# Plot NMDS
ggplot(nmds_df, aes(NMDS1, NMDS2, color = sampling_year, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = sampling_year),
    linetype = "dashed",
    linewidth = 0.6,
    alpha = 0.5,
    na.rm = TRUE,
    show.legend = FALSE
  ) +
  labs(
    title = paste0("NMDS (Bray–Curtis) — Stress = ", round(nmds$stress, 3)),
    color = "Year",
    shape = "Reach"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# ---------------------------------------------------------------------
# 6. PCoA Ordination
# ---------------------------------------------------------------------
pcoa_res <- cmdscale(dist_bray, eig = TRUE, k = 2)
pcoa_df <- data.frame(
  Axis1 = pcoa_res$points[, 1],
  Axis2 = pcoa_res$points[, 2],
  meta_aligned
)

ggplot(pcoa_df, aes(Axis1, Axis2, color = sampling_year, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "PCoA (Bray–Curtis Distance)",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

# ---------------------------------------------------------------------
# 7. DCA (Detrended Correspondence Analysis)
# ---------------------------------------------------------------------
dca_res <- decorana(comm)

dca_scores <- scores(dca_res, display = "sites") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("unit_id")

dca_df <- inner_join(dca_scores, meta_aligned, by = "unit_id")

ggplot(dca_df, aes(DCA1, DCA2, color = sampling_year, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    title = "DCA — Temporal Variation in Fish Assemblages",
    x = "Axis 1",
    y = "Axis 2"
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))
