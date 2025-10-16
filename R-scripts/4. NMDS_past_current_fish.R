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
  dist_bray ~ sampling_year + river_reach,
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
# 6. PCoA Ordination (with ellipses)
# ---------------------------------------------------------------------

# Compute Bray–Curtis distance and PCoA
dist_bray <- vegan::vegdist(comm, method = "bray")
pcoa_res <- cmdscale(dist_bray, eig = TRUE, k = 2)

# Combine ordination axes with metadata
pcoa_df <- data.frame(
  Axis1 = pcoa_res$points[, 1],
  Axis2 = pcoa_res$points[, 2],
  meta_aligned
)

# Plot PCoA with ellipses grouped by sampling_year
library(ggplot2)

ggplot(pcoa_df, aes(Axis1, Axis2, color = sampling_year, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(aes(group = sampling_year),
               linetype = "dashed", linewidth = 0.6, alpha = 0.6) +
  labs(
    title = "PCoA (Bray–Curtis Distance)",
    x = paste0("Axis 1 (", round(pcoa_res$eig[1] / sum(pcoa_res$eig) * 100, 1), "%)"),
    y = paste0("Axis 2 (", round(pcoa_res$eig[2] / sum(pcoa_res$eig) * 100, 1), "%)"),
    color = "Sampling Year",
    shape = "River Reach"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )

# ---------------------------------------------------------------------
# 7. DCA Ordination — Sites (black) and Species (red)
# ---------------------------------------------------------------------

# Run DCA
dca_res <- vegan::decorana(comm)

# Extract site and species scores
site_scores <- as.data.frame(scores(dca_res, display = "sites")) %>%
  tibble::rownames_to_column("unit_id")

species_scores <- as.data.frame(scores(dca_res, display = "species")) %>%
  tibble::rownames_to_column("fish_species")

# Combine site scores with metadata
dca_df <- inner_join(site_scores, meta_aligned, by = "unit_id")

# --- compute padded limits from both sites and species ---
x_all <- c(dca_df$DCA1, species_scores$DCA1)
y_all <- c(dca_df$DCA2, species_scores$DCA2)
pad_x <- diff(range(x_all)) * 0.15  # 15% padding
pad_y <- diff(range(y_all)) * 0.15

library(ggplot2)
# If you have ggrepel installed, uncomment the next line and the two geom_text_repel layers
# library(ggrepel)

ggplot() +
  # Sites
  geom_point(
    data = dca_df,
    aes(DCA1, DCA2, shape = river_reach, color = sampling_year),
    size = 3.5, alpha = 0.9
  ) +
  # Site labels (black)
  # geom_text_repel(
  #   data = dca_df,
  #   aes(DCA1, DCA2, label = location_id),
  #   color = "black", size = 3, max.overlaps = Inf
  # ) +
  geom_text(
    data = dca_df,
    aes(DCA1, DCA2, label = location_id),
    color = "black", size = 3.5, vjust = -0.8
  ) +
  
  # Species labels (red, italic)
  # geom_text_repel(
  #   data = species_scores,
  #   aes(DCA1, DCA2, label = fish_species),
  #   color = "red", fontface = "italic", size = 3, max.overlaps = Inf
  # ) +
  geom_text(
    data = species_scores,
    aes(DCA1, DCA2, label = fish_species),
    color = "red", fontface = "italic", size = 3.5
  ) +
  
  # Axes padding
  scale_x_continuous(limits = c(min(x_all) - pad_x, max(x_all) + pad_x), expand = expansion(0)) +
  scale_y_continuous(limits = c(min(y_all) - pad_y, max(y_all) + pad_y), expand = expansion(0)) +
  
  labs(
    title = "DCA",
    x = "DCA Axis 1",
    y = "DCA Axis 2",
    color = "Year",
    shape = "River Reach"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(10, 20, 10, 20),  # extra room around the panel
    panel.spacing = unit(6, "pt")
  ) +
  coord_cartesian(clip = "off")  # lets labels extend beyond panel if needed
#---------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------
# 2021-2022 Fish Community Analysis — NMDS, PERMANOVA, PCoA, and DCA
# =====================================================================

# Clean workspace
remove(list = ls())

# Load libraries
library(tidyverse)
library(vegan)
library(janitor)
library(readr)
library(ggpubr)
library(ggplot2)

#load data
# Load and filter data for sites M1–M9 and years 2021–2022
fish_long <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names() %>%
  mutate(
    fish_species = stringr::str_squish(fish_species),
    abundance = 1
  ) %>%
  filter(
    location_id %in% paste0("M", 1:9),         # keep only M1 to M9
    sampling_year %in% c(2021, 2022),          # keep only years 2021 and 2022
    !is.na(fish_species)                       # ensure species info is not missing
  )

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
  dist_bray ~ sampling_year + river_reach,
  data = meta_aligned,
  permutations = 999,
  by = "margin"
)
print(perm)

# ---------------------------------------------------------------------
# 5. NMDS Ordination
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# 5. NMDS Ordination — Ellipses by River Reach
# ---------------------------------------------------------------------
set.seed(123)
nmds <- metaMDS(comm, distance = "bray", k = 2, trymax = 100)

# Extract NMDS site scores and join metadata
scores_nmds <- scores(nmds, display = "sites") %>%
  as.data.frame() %>%
  tibble::rownames_to_column("unit_id")

nmds_df <- inner_join(scores_nmds, meta_aligned, by = "unit_id") %>%
  mutate(
    river_reach = factor(river_reach, levels = c("Upstream", "Midstream", "Downstream"))
  )

# NMDS plot with ellipses around river reach groups
ggplot(nmds_df, aes(NMDS1, NMDS2, color = river_reach, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = river_reach),   # 👈 ellipse by reach
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.6,
    na.rm = TRUE
  ) +
  labs(
    title = paste0("NMDS (Bray–Curtis) — Stress = ", round(nmds$stress, 3)),
    color = "River Reach",
    shape = "River Reach"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right"
  )
#more tight ellipses
ggplot(nmds_df, aes(NMDS1, NMDS2, color = river_reach, shape = river_reach)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(
    aes(group = river_reach),
    type = "norm",      # use covariance ellipse
    level = 0.68,       # ~1 SD; try 0.50–0.75 if you want smaller/larger
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.6,
    na.rm = TRUE
  ) +
  labs(
    title = paste0("NMDS (Bray–Curtis) — Stress = ", round(nmds$stress, 3)),
    color = "River Reach",
    shape = "River Reach"
  ) +
  coord_equal() +        # prevents ellipse distortion
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5))

###########################################################################

# ---------------------------------------------------------------------
# 7. DCA Ordination — Sites (black) and Species (red)
# ---------------------------------------------------------------------

# Run DCA
dca_res <- vegan::decorana(comm)

# Extract site and species scores
site_scores <- as.data.frame(scores(dca_res, display = "sites")) %>%
  tibble::rownames_to_column("unit_id")

species_scores <- as.data.frame(scores(dca_res, display = "species")) %>%
  tibble::rownames_to_column("fish_species")

# Combine site scores with metadata
dca_df <- inner_join(site_scores, meta_aligned, by = "unit_id")

# --- compute padded limits from both sites and species ---
x_all <- c(dca_df$DCA1, species_scores$DCA1)
y_all <- c(dca_df$DCA2, species_scores$DCA2)
pad_x <- diff(range(x_all)) * 0.15  # 15% padding
pad_y <- diff(range(y_all)) * 0.15

library(ggplot2)
# If you have ggrepel installed, uncomment the next line and the two geom_text_repel layers
# library(ggrepel)

ggplot() +
  # Sites
  geom_point(
    data = dca_df,
    aes(DCA1, DCA2, shape = river_reach, color = sampling_year),
    size = 3.5, alpha = 0.9
  ) +
  # Site labels (black)
  # geom_text_repel(
  #   data = dca_df,
  #   aes(DCA1, DCA2, label = location_id),
  #   color = "black", size = 3, max.overlaps = Inf
  # ) +
  geom_text(
    data = dca_df,
    aes(DCA1, DCA2, label = location_id),
    color = "black", size = 3.5, vjust = -0.8
  ) +
  
  # Species labels (red, italic)
  # geom_text_repel(
  #   data = species_scores,
  #   aes(DCA1, DCA2, label = fish_species),
  #   color = "red", fontface = "italic", size = 3, max.overlaps = Inf
  # ) +
  geom_text(
    data = species_scores,
    aes(DCA1, DCA2, label = fish_species),
    color = "red", fontface = "italic", size = 3.5
  ) +
  
  # Axes padding
  scale_x_continuous(limits = c(min(x_all) - pad_x, max(x_all) + pad_x), expand = expansion(0)) +
  scale_y_continuous(limits = c(min(y_all) - pad_y, max(y_all) + pad_y), expand = expansion(0)) +
  
  labs(
    title = "DCA",
    x = "DCA Axis 1",
    y = "DCA Axis 2",
    color = "Year",
    shape = "River Reach"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title   = element_text(hjust = 0.5, size = 10),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 10),
    plot.margin  = margin(10, 20, 10, 20),  # extra room around the panel
    panel.spacing = unit(6, "pt")
  ) +
  coord_cartesian(clip = "off")  # lets labels extend beyond panel if needed
