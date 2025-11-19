# ===============================================
# Fish isotopes using SIBER (non-Bayesian)
# ===============================================

# (optional) clear & restore
remove(list = ls())

#setup renv
renv::restore()

# ---- Packages ----
library(tidyverse)
library(readr)
library(SIBER)

# ---- Load data ----
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code   %in% target_sites)

# ---- Robust column chooser (handles exact header + variants) ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) {
    stop("None of these columns were found: ", paste(candidates, collapse = " | "))
  }
  hit[[1]]
}

c13_name <- pick_first_col(df, c(
  "d13C (permil, vs VPDB)",
  "Normalized d13C",
  "d13C (‰, vs VPDB)",
  "d13C", "d13C_corrected", "C13"
))
n15_name <- pick_first_col(df, c(
  "d15N (permil, vs AIR)",
  "d15N (‰, vs AIR)",
  "d15N", "N15"
))

# ---- Construct analysis columns & basic QC ----
df <- df %>%
  mutate(
    d13C_use = .data[[c13_name]],
    d15N_use = .data[[n15_name]],
    trophic_group = case_when(
      Fish_species == "Labeobarbus altianalis" ~ "Omnivore–benthivore",
      Fish_species == "Labeo victorianus"      ~ "Detritivore–herbivore",
      TRUE ~ "Other"
    )
  ) %>%
  drop_na(d13C_use, d15N_use)

# ---- Ensure minimum sample size per species × site (SIBER requirement) ----
grp_sizes <- df %>% count(Site_code, Fish_species, name = "n")
df_ok <- df %>%
  inner_join(grp_sizes %>% filter(n >= 3),
             by = c("Site_code","Fish_species"))

if (nrow(df_ok) == 0) stop("No groups with n >= 3 after filtering.")

# ---- Keys to map species/site to SIBER integers ----
species_key <- df_ok %>%
  distinct(Fish_species) %>% arrange(Fish_species) %>%
  mutate(group_id = row_number())

site_key <- df_ok %>%
  distinct(Site_code) %>% arrange(Site_code) %>%
  mutate(comm_id = row_number())

df_id <- df_ok %>%
  left_join(species_key, by = "Fish_species") %>%
  left_join(site_key,   by = "Site_code")

# ---- Build SIBER object ----
siber_df <- df_id %>%
  transmute(
    iso1      = d13C_use,             # δ13C
    iso2      = d15N_use,             # δ15N
    group     = as.integer(group_id), # species id
    community = as.integer(comm_id)   # site id
  ) %>%
  as.data.frame()

siber_obj <- createSiberObject(siber_df)

# ---- SEAc (non-Bayesian standard ellipse areas) ----
SEAc <- siberEllipses(siber_obj)

cat("\nSEAc (rows = communities/sites per site_key order; cols = groups/species per species_key order):\n")
print(SEAc)
cat("\nSpecies key:\n"); print(species_key)
cat("\nSite key:\n");    print(site_key)

# ---- Plot (ggplot): points + 40% normal ellipses, facets by site ----
sp_cols <- c("Labeo victorianus" = "#1F78B4",      # blue
             "Labeobarbus altianalis" = "#E41A1C") # red

df_ok <- df_ok %>% mutate(Site_code = factor(Site_code, levels = paste0("M", 4:9)))

p <- ggplot(df_ok, aes(x = d13C_use, y = d15N_use, color = Fish_species)) +
  geom_point(size = 2.2, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.40, linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~ Site_code, nrow = 2) +
  scale_color_manual(values = sp_cols, name = "Species") +
  labs(
    title = expression(paste(delta^13, "C vs ", delta^15, "N by site (M4–M9)")),
    x = expression(paste(delta^13, "C (‰)")),
    y = expression(paste(delta^15, "N (‰)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.background= element_rect(fill = "grey92"),
    strip.text      = element_text(face = "bold")
  )

print(p)

# ---- Optional summaries ----
means_by_trophic <- df_ok %>%
  group_by(trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups   = "drop"
  )
cat("\nMeans by trophic group (n >= 3):\n"); print(means_by_trophic)

means_by_trophic_site <- df_ok %>%
  group_by(Site_code, trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups   = "drop"
  )
cat("\nMeans by trophic group and site (n >= 3):\n"); print(means_by_trophic_site)
# ---- End of script ----
#######################################################################
# ===============================================
# Fish isotopes using SIBER (non-Bayesian) — robust & ordered legend
# ===============================================

# (optional) clear env
remove(list = ls())
# renv::restore()

# ---- Packages ----
library(tidyverse)
library(readr)
library(SIBER)

# ---- Load data ----
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code   %in% target_sites)

# ---- Robust column chooser (handles exact header + variants) ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (!length(hit)) stop("None of these columns were found: ", paste(candidates, collapse = " | "))
  hit[[1]]
}

c13_name <- pick_first_col(df, c(
  "d13C (permil, vs VPDB)",
  "Normalized d13C",
  "d13C (‰, vs VPDB)",
  "d13C"
))
n15_name <- pick_first_col(df, c(
  "d15N (permil, vs AIR)",
  "d15N (‰, vs AIR)",
  "d15N"
))

# ---- Analysis columns & QC ----
df <- df %>%
  mutate(
    d13C_use = .data[[c13_name]],
    d15N_use = .data[[n15_name]],
    trophic_group = case_when(
      Fish_species == "Labeobarbus altianalis" ~ "Omnivore–benthivore",
      Fish_species == "Labeo victorianus"      ~ "Detritivore–herbivore",
      TRUE ~ "Other"
    )
  ) %>%
  drop_na(d13C_use, d15N_use)

# ---- Sample sizes per site × species ----
sizes <- df %>% count(Site_code, Fish_species, name = "n")

# For plotting ellipses: allow n ≥ 3
df_for_ellipse <- df %>%
  inner_join(filter(sizes, n >= 3), by = c("Site_code","Fish_species"))

# For SEAc: require n ≥ 5 to avoid eigen errors
df_for_seac <- df %>%
  inner_join(filter(sizes, n >= 5), by = c("Site_code","Fish_species"))

message("Cells available for ellipses (n ≥ 3):")
print(arrange(filter(sizes, n >= 3), Site_code, Fish_species))
message("Cells available for SEAc (n ≥ 5):")
print(arrange(filter(sizes, n >= 5), Site_code, Fish_species))

# ---- Build SIBER object on the safer subset (n ≥ 5) ----
if (nrow(df_for_seac) > 0) {
  species_key <- df_for_seac %>%
    distinct(Fish_species) %>% arrange(Fish_species) %>%
    mutate(group_id = row_number())
  
  site_key <- df_for_seac %>%
    distinct(Site_code) %>% arrange(Site_code) %>%
    mutate(comm_id = row_number())
  
  df_id <- df_for_seac %>%
    left_join(species_key, by = "Fish_species") %>%
    left_join(site_key,   by = "Site_code")
  
  siber_df <- df_id %>%
    transmute(
      iso1      = d13C_use,
      iso2      = d15N_use,
      group     = as.integer(group_id),
      community = as.integer(comm_id)
    ) %>% as.data.frame()
  
  siber_obj <- createSiberObject(siber_df)
  
  SEAc <- siberEllipses(siber_obj)
  cat("\nSEAc (rows = sites per site_key; cols = species per species_key):\n")
  print(SEAc)
  cat("\nSpecies key:\n"); print(species_key)
  cat("\nSite key:\n");    print(site_key)
} else {
  warning("No site × species cells have n ≥ 5; skipping SEAc to avoid numerical errors.")
}

# ---- Plot: points for all; ellipses where n ≥ 3; legend starts with L. altianalis ----
# Custom colors with desired legend order
sp_levels <- c("Labeobarbus altianalis", "Labeo victorianus")
sp_cols   <- c("Labeobarbus altianalis" = "#E41A1C",  # red first
               "Labeo victorianus"      = "#1F78B4")  # blue second

df <- df %>%
  mutate(
    Site_code    = factor(Site_code, levels = paste0("M", 4:9)),
    Fish_species = factor(Fish_species, levels = sp_levels)
  )

df_for_ellipse <- df_for_ellipse %>%
  mutate(
    Site_code    = factor(Site_code, levels = paste0("M", 4:9)),
    Fish_species = factor(Fish_species, levels = sp_levels)
  )

p <- ggplot() +
  geom_point(data = df,
             aes(x = d13C_use, y = d15N_use, color = Fish_species),
             size = 2.2, alpha = 0.9) +
  stat_ellipse(data = df_for_ellipse,
               aes(x = d13C_use, y = d15N_use, color = Fish_species),
               type = "norm", level = 0.40, linewidth = 0.9, linetype = "dashed") +
  facet_wrap(~ Site_code, nrow = 2) +
  scale_color_manual(values = sp_cols, name = "Species") +
  labs(
    x = expression(paste(delta^13, "C (‰)")),
    y = expression(paste(delta^15, "N (‰)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid      = element_blank(),
    legend.position = "bottom",
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(hjust = 0.5, face = "bold"),
    strip.background= element_rect(fill = "grey92"),
    strip.text      = element_text(face = "bold")
  )

print(p)

# ---- Descriptive summaries (all points) ----
means_by_trophic <- df %>%
  group_by(trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups   = "drop"
  )
cat("\nMeans by trophic group (all points):\n"); print(means_by_trophic)

means_by_trophic_site <- df %>%
  group_by(Site_code, trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups   = "drop"
  )
cat("\nMeans by trophic group and site (all points):\n"); print(means_by_trophic_site)
############################################################################
# ==============================
# δ13C box + dot plot by species
# ==============================

remove(list = ls())
# renv::restore()

# ---- Packages ----
library(tidyverse)
library(readr)

# ---- Load data ----
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code   %in% target_sites)

# ---- Helper to pick the δ13C column robustly ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (!length(hit)) stop("None of these columns were found: ",
                         paste(candidates, collapse = " | "))
  hit[[1]]
}

c13_name <- pick_first_col(df, c(
  "d13C (permil, vs VPDB)",
  "Normalized d13C",
  "d13C (‰, vs VPDB)",
  "d13C"
))

# ---- Prepare data for plotting ----
sp_levels <- c("Labeobarbus altianalis", "Labeo victorianus")

df <- df %>%
  mutate(
    d13C_use    = .data[[c13_name]],
    Fish_species = factor(Fish_species, levels = sp_levels)
  ) %>%
  drop_na(d13C_use)

# Optional: quick check of sample sizes and means
df %>%
  group_by(Fish_species) %>%
  summarise(
    n       = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ---- Plot: horizontal boxplot + jitter + mean point ----
p <- ggplot(df, aes(x = d13C_use, y = Fish_species)) +
  # boxplot for distribution
  geom_boxplot(
    width = 0.5,
    alpha = 0.4,
    outlier.shape = NA,
    fill = "grey80",
    color = "black"
  ) +
  # jittered points (individual fish)
  geom_jitter(
    height = 0.12,
    size   = 2,
    alpha  = 0.7
  ) +
  # mean as a distinct dot
  stat_summary(
    fun   = mean,
    geom  = "point",
    shape = 21,
    size  = 3.5,
    fill  = "black",
    color = "white"
  ) +
  labs(
    x = expression(paste(delta^13, "C (‰, vs VPDB)")),
    y = "Fish species"
  ) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.y       = element_text(face = "bold"),
    axis.title.x       = element_text(face = "bold")
  )

print(p)
############################################################################
# ==============================
# δ15N box + dot plot by species (vertical)
# ==============================

remove(list = ls())
# renv::restore()

# ---- Packages ----
library(tidyverse)
library(readr)

# ---- Load data ----
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code   %in% target_sites)

# ---- Helper to pick the δ15N column robustly ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (!length(hit)) stop("None of these columns were found: ",
                         paste(candidates, collapse = " | "))
  hit[[1]]
}

n15_name <- pick_first_col(df, c(
  "d15N (permil, vs AIR)",
  "d15N (‰, vs AIR)",
  "Normalized d15N",
  "d15N"
))

# ---- Prepare data for plotting ----
sp_levels <- c("Labeobarbus altianalis", "Labeo victorianus")

df <- df %>%
  mutate(
    d15N_use     = .data[[n15_name]],
    Fish_species = factor(Fish_species, levels = sp_levels)
  ) %>%
  drop_na(d15N_use)

# ---- Quick check of sample sizes and means ----
df %>%
  group_by(Fish_species) %>%
  summarise(
    n        = n(),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# ---- Plot: vertical boxplot + jitter + mean point ----
p <- ggplot(df, aes(x = Fish_species, y = d15N_use)) +
  
  # boxplot
  geom_boxplot(
    width = 0.55,
    alpha = 0.4,
    outlier.shape = NA,
    fill = "grey80",
    color = "black"
  ) +
  
  # jittered individual points
  geom_jitter(
    width  = 0.12,
    size   = 2,
    alpha  = 0.7
  ) +
  
  # mean point
  stat_summary(
    fun   = mean,
    geom  = "point",
    shape = 21,
    size  = 3.8,
    fill  = "black",
    color = "white"
  ) +
  
  # labels
  labs(
    x = "Fish species",
    y = expression(paste(delta^15, "N (‰, vs AIR)"))
  ) +
  
  # theme
  theme_bw(base_size = 13) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title.x       = element_text(face = "bold"),
    axis.title.y       = element_text(face = "bold")
  )

print(p)
############################################################################