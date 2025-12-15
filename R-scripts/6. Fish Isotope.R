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

##############################################################################
# ===============================================
# Fish isotopes in the Mara River
# Isotopic niche shown as convex hull polygons
# (no SEAc, no ellipses, no Layman metrics)
# ===============================================

rm(list = ls())
graphics.off()

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
  filter(
    Fish_species %in% target_species,
    Site_code    %in% target_sites
  )

# ---- Robust column chooser (handles exact header + variants) ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) {
    stop("None of these columns were found: ",
         paste(candidates, collapse = " | "))
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

# ---- Ensure minimum n per species × site (n >= 3 is sensible) ----
df_ok <- df %>%
  group_by(Site_code, Fish_species) %>%
  filter(n() >= 3) %>%
  ungroup()

if (nrow(df_ok) == 0) stop("No species × site groups with n >= 3 after filtering.")

# ---- Build convex hull polygons for each species at each site ----
hulls <- df_ok %>%
  group_by(Site_code, Fish_species) %>%
  # chull returns row indices of convex hull vertices
  slice(chull(d13C_use, d15N_use)) %>%
  ungroup()

# ---- Plot δ13C–δ15N with convex hulls (isotopic niche) ----
sp_cols <- c(
  "Labeo victorianus"      = "#1F78B4",
  "Labeobarbus altianalis" = "#E41A1C"
)

df_ok <- df_ok %>%
  mutate(Site_code = factor(Site_code, levels = paste0("M", 4:9)))

p_niche <- ggplot() +
  # hull polygons show niche breadth per species × site
  geom_polygon(
    data = hulls,
    aes(x = d13C_use, y = d15N_use,
        group = interaction(Site_code, Fish_species),
        fill  = Fish_species),
    alpha = 0.20, colour = NA
  ) +
  # individual fish
  geom_point(
    data = df_ok,
    aes(x = d13C_use, y = d15N_use, color = Fish_species),
    size = 2.2, alpha = 0.9
  ) +
  facet_wrap(~ Site_code, nrow = 2) +
  scale_color_manual(values = sp_cols, name = "Species") +
  scale_fill_manual(values  = sp_cols, name = "Species") +
  labs(
    title = expression(paste(delta^13, "C vs ", delta^15, "N for two fish species (M4–M9)")),
    x = expression(paste(delta^13, "C (‰)")),
    y = expression(paste(delta^15, "N (‰)"))
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )

print(p_niche)
##############################################################################

# ===============================================================
# Fish isotopes in the Mara River
# Bayesian SEA (SEAb) density plots for isotopic niche width
# - NO SEAc
# - NO ellipses on the isotope plot
# - NO Layman metrics
# ===============================================================

# Optional: clean environment
rm(list = ls())
graphics.off()

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
  filter(
    Fish_species %in% target_species,
    Site_code    %in% target_sites
  )

# ---- Robust column chooser (handles exact header + variants) ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) {
    stop("None of these columns were found: ",
         paste(candidates, collapse = " | "))
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
    d15N_use = .data[[n15_name]]
  ) %>%
  drop_na(d13C_use, d15N_use)

# ---- Ensure minimum sample size per species × site (n >= 3) ----
df_ok <- df %>%
  group_by(Site_code, Fish_species) %>%
  filter(n() >= 3) %>%
  ungroup()

if (nrow(df_ok) == 0) stop("No species × site groups with n >= 3 after filtering.")

# ---- Encode for SIBER object ----
# group = species, community = site
df_siber <- df_ok %>%
  mutate(
    group     = as.numeric(as.factor(Fish_species)),
    community = as.numeric(as.factor(Site_code))
  ) %>%
  transmute(
    iso1      = d13C_use,
    iso2      = d15N_use,
    group,
    community
  )

siber_obj <- createSiberObject(as.data.frame(df_siber))

# ---- Optional: ML group metrics just to get clean group labels ----
group_ML <- groupMetricsML(siber_obj)
# Columns of group_ML correspond to groups (species × site) in the same order
group_labels <- colnames(group_ML)

# ---- Bayesian SIBER model settings ----
set.seed(1)

parms <- list(
  n.iter   = 2 * 10^4,
  n.burnin = 1 * 10^3,
  n.thin   = 10,
  n.chains = 2
)

priors <- list(
  R      = diag(2),
  k      = 2,
  tau.mu = 1.0E-3
)

# ---- Fit Bayesian ellipses (posterior) ----
ellipses.posterior <- siberMVN(siber_obj, parms, priors)

# ---- Bayesian SEAb: posterior Standard Ellipse Area per group ----
# This is the CORRECT function for SEA.B in SIBER
SEA_B <- siberEllipses(ellipses.posterior)
# SEA_B = matrix: rows = posterior draws, cols = groups
# give columns informative names
colnames(SEA_B) <- group_labels

# ---- Tidy SEA_B for density plots ----
SEA_df <- SEA_B %>%
  as.data.frame() %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols      = -draw,
    names_to  = "group_id",
    values_to = "SEA"
  )

# ---- Build a mapping from group_id to Site_code & Fish_species ----
# We use the original df_ok with numeric group/community codes
group_map <- df_ok %>%
  mutate(
    group     = as.numeric(as.factor(Fish_species)),
    community = as.numeric(as.factor(Site_code))
  ) %>%
  distinct(community, group, Site_code, Fish_species) %>%
  arrange(community, group)

# SIBER's internal order is "community.group"
# e.g. "1.1", "1.2", "2.1", ...
expected_group_ids <- group_map %>%
  mutate(group_id = paste(community, group, sep = ".")) %>%
  pull(group_id)

# Check: these should match the column names of SEA_B / group_labels
# print(group_labels)
# print(expected_group_ids)

# Join mapping into SEA_df
SEA_df <- SEA_df %>%
  left_join(
    group_map %>%
      mutate(group_id = paste(community, group, sep = ".")),
    by = "group_id"
  )

# ---- Summarise SEA posterior: median, mode, 95% credible interval ----
mode_est <- function(x) {
  d <- density(x)
  d$x[which.max(d$y)]
}

SEA_summary <- SEA_df %>%
  group_by(Site_code, Fish_species) %>%
  summarise(
    n_draws = n(),
    SEA_mode   = mode_est(SEA),
    SEA_median = median(SEA),
    SEA_lower95 = quantile(SEA, 0.025),
    SEA_upper95 = quantile(SEA, 0.975),
    .groups = "drop"
  )

cat("\nBayesian SEA (SEAb) summary per species × site:\n")
print(SEA_summary)

# ---- Plot: SEA posterior density curves (isotopic niche width) ----
SEA_df_plot <- SEA_df %>%
  filter(!is.na(Site_code), !is.na(Fish_species)) %>%
  mutate(
    Site_code    = factor(Site_code,    levels = paste0("M", 4:9)),
    Fish_species = factor(Fish_species)
  )

p_SEA <- ggplot(SEA_df_plot,
                aes(x = SEA, fill = Fish_species, colour = Fish_species)) +
  geom_density(alpha = 0.30) +
  facet_wrap(~ Site_code, scales = "free") +
  labs(
    title = "Bayesian posterior distribution of Standard Ellipse Area (SEAb)",
    x     = expression(paste("Standard Ellipse Area (", "\u2030"^2, ")")),
    y     = "Posterior density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )

print(p_SEA)

sp_cols <- c(
  "Labeo victorianus"      = "#1F78B4",
  "Labeobarbus altianalis" = "#E41A1C"
)

p_SEA_col <- ggplot(SEA_df_plot,
                    aes(x = SEA, fill = Fish_species, colour = Fish_species)) +
  geom_density(alpha = 0.30) +
  scale_fill_manual(values = sp_cols) +
  scale_colour_manual(values = sp_cols) +
  facet_wrap(~ Site_code, scales = "free") +
  labs(
    title = "Bayesian posterior distribution of Standard Ellipse Area (SEAb)",
    x     = expression(paste("Standard Ellipse Area (", "\u2030"^2, ")")),
    y     = "Posterior density"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    strip.background = element_rect(fill = "grey92"),
    strip.text       = element_text(face = "bold")
  )
print(p_SEA_col)
#arrange the Leged to start with Labeobarbus altianalis
p_SEA_col + guides(fill = guide_legend(reverse = TRUE), colour = guide_legend(reverse = TRUE))
