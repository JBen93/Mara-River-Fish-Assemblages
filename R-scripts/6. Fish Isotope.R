# ===============================================
# Fish Isotope using SIBER (non-Bayesian)
# ===============================================
# clear everything in memory (of R)
remove(list=ls())
renv::restore()
# load the the required packages
library(tidyverse)
library(readr)
library(SIBER)

#data URL source if you need to inspect for the whole dataset
#browseURL("https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pubhtml")

# Load data from Google Sheets
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)
# ---- Filter species & sites ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code %in% target_sites)

# ---- Choose δ13C column (prefer normalized) ----
has_norm_c <- "Normalized d13C" %in% names(df)
df <- df %>%
  mutate(
    d13C_use = if (has_norm_c) `Normalized d13C` else `d13C (permil, vs VPDB)`,
    d15N_use = `d15N (permil, vs AIR)`
  ) %>%
  drop_na(d13C_use, d15N_use)

# ---- Trophic group tags (optional) ----
df <- df %>%
  mutate(trophic_group = case_when(
    Fish_species == "Labeobarbus altianalis" ~ "Omnivore–benthivore",
    Fish_species == "Labeo victorianus"      ~ "Detritivore–herbivore",
    TRUE ~ "Other"
  ))

# ---- Remove groups with n < 3 (fixes SIBER eigen error) ----
grp_sizes <- df %>%
  count(Site_code, Fish_species, name = "n")

df_ok <- df %>%
  inner_join(grp_sizes %>% filter(n >= 3),
             by = c("Site_code","Fish_species"))

# ---- Prepare SIBER object (exact column order) ----
species_key <- df_ok %>%
  distinct(Fish_species) %>% arrange(Fish_species) %>%
  mutate(group_id = row_number())

site_key <- df_ok %>%
  distinct(Site_code) %>% arrange(Site_code) %>%
  mutate(comm_id = row_number())

df_id <- df_ok %>%
  left_join(species_key, by = "Fish_species") %>%
  left_join(site_key,   by = "Site_code")

siber_df <- df_id %>%
  transmute(
    iso1      = d13C_use,        # δ13C
    iso2      = d15N_use,        # δ15N
    group     = as.integer(group_id),   # species id
    community = as.integer(comm_id)     # site id
  ) %>%
  as.data.frame()

siber_obj <- createSiberObject(siber_df)

# ---- SEAc (non-Bayesian) ----
SEAc <- siberEllipses(siber_obj)
cat("\nSEAc (rows = communities/sites in site_key order; cols = groups/species in species_key order):\n")
print(SEAc)
cat("\nSpecies key:\n"); print(species_key)
cat("\nSite key:\n");    print(site_key)

# ---- Publication figure (ggplot) with clean legend & axes ----
# Colours for the two species:
sp_cols <- c("Labeo victorianus" = "#1f78b4",    # blue
             "Labeobarbus altianalis" = "#e31a1c")  # red

# Order facets M4 -> M9
df_ok <- df_ok %>% mutate(Site_code = factor(Site_code, levels = paste0("M",4:9)))

p <- ggplot(df_ok, aes(x = d13C_use, y = d15N_use, color = Fish_species)) +
  geom_point(size = 2.2, alpha = 0.9) +
  # 40% normal ellipse ~ SIBER "standard ellipse" visual analogue
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
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text = element_text(face = "bold")
  )

print(p)

# ---- Summary means (your earlier tables) ----
means_by_trophic <- df_ok %>%
  group_by(trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups = "drop"
  )
cat("\nMeans by trophic group (filtered n>=3):\n")
print(means_by_trophic)

means_by_trophic_site <- df_ok %>%
  group_by(Site_code, trophic_group) %>%
  summarise(
    n         = n(),
    mean_d13C = mean(d13C_use, na.rm = TRUE),
    sd_d13C   = sd(d13C_use,   na.rm = TRUE),
    mean_d15N = mean(d15N_use, na.rm = TRUE),
    sd_d15N   = sd(d15N_use,   na.rm = TRUE),
    .groups = "drop"
  )
cat("\nMeans by trophic group and site (filtered n>=3):\n")
print(means_by_trophic_site)

