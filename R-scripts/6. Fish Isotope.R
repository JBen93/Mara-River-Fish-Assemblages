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


# ===============================================
# Isotopic niche overlap across sites using SIBER
# (1) maxLikOverlap (point estimate)
# (2) bayesianOverlap (posterior distribution; optional)
# ===============================================

remove(list = ls())

library(tidyverse)
library(readr)
library(SIBER)

# ---- Load data ----
raw <- read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code %in% target_sites)

# ---- Robust column chooser ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) stop("None of these columns were found: ", paste(candidates, collapse = " | "))
  hit[[1]]
}

c13_name <- pick_first_col(df, c("d13C (permil, vs VPDB)", "Normalized d13C", "d13C (‰, vs VPDB)", "d13C"))
n15_name <- pick_first_col(df, c("d15N (permil, vs AIR)",  "d15N (‰, vs AIR)",  "d15N"))

# ---- Clean + enforce numeric ----
df <- df %>%
  mutate(
    d13C_use = suppressWarnings(as.numeric(.data[[c13_name]])),
    d15N_use = suppressWarnings(as.numeric(.data[[n15_name]]))
  ) %>%
  filter(is.finite(d13C_use), is.finite(d15N_use))

# ---- Check sample sizes per site x species ----
grp_sizes <- df %>% count(Site_code, Fish_species, name = "n")
print(grp_sizes)

# Strong recommendation: n >= 5 for stable covariance / ellipses
MIN_N <- 4

df_ok <- df %>%
  inner_join(grp_sizes %>% filter(n >= MIN_N),
             by = c("Site_code", "Fish_species"))

if (nrow(df_ok) == 0) stop("No site × species groups with n >= ", MIN_N, ". Lower MIN_N to 3, but expect instability.")

# ---- Keys for SIBER ----
species_key <- df_ok %>%
  distinct(Fish_species) %>% arrange(Fish_species) %>%
  mutate(group_id = row_number())

site_key <- df_ok %>%
  distinct(Site_code) %>% arrange(Site_code) %>%
  mutate(comm_id = row_number())

df_id <- df_ok %>%
  left_join(species_key, by = "Fish_species") %>%
  left_join(site_key, by = "Site_code")

# ---- Build SIBER object ----
siber_df <- df_id %>%
  transmute(
    iso1      = d13C_use,
    iso2      = d15N_use,
    group     = as.integer(group_id),  # species
    community = as.integer(comm_id)    # site
  ) %>%
  as.data.frame()

siber_obj <- createSiberObject(siber_df)

cat("\nSpecies key:\n"); print(species_key)
cat("\nSite key:\n");    print(site_key)
cat("\nSample sizes used by SIBER:\n"); print(siber_obj$sample.sizes)


# ==========================================================
# 2) BAYESIAN OVERLAP PER SITE (uncertainty; requires JAGS)
# ==========================================================
# This section uses siberMVN() + bayesianOverlap() like your screenshot.
# If you don't have JAGS installed, skip this section.

RUN_BAYES <- FALSE  # set TRUE if JAGS is installed and you want posterior overlap

if (RUN_BAYES) {
  
  # These priors are the standard SIBER defaults/examples (reasonable starting point)
  # You can tune these later, but keep them as-is first.
  parms <- list(
    n.iter = 2 * 10^4,
    n.burnin = 5 * 10^3,
    n.thin = 10,
    n.chains = 2
  )
  
  priors <- list(
    R = diag(2),
    k = 2,
    tau.mu = 1.0E-3
  )
  
  # Fit posterior MVN ellipses across all communities/groups
  ellipses_posterior <- siberMVN(siber_obj, parms, priors)
  
  # Per-site Bayesian overlap draws
  DRAWS <- 1000
  
  bayes_overlap_by_site <- site_key %>%
    mutate(
      label_LV = paste0(comm_id, ".", g_LV),
      label_LB = paste0(comm_id, ".", g_LB)
    ) %>%
    rowwise() %>%
    mutate(
      bayes_obj = list(
        tryCatch(
          bayesianOverlap(label_LV, label_LB,
                          ellipses_posterior,
                          draws = DRAWS,
                          p.interval = P_ELLIPSE,
                          n = N_POLY),
          error = function(e) NA
        )
      )
    ) %>%
    ungroup()
  
  # bayesianOverlap returns a table; you can summarise overlap distribution:
  # (the output object structure can vary by SIBER version, so print first)
  print(bayes_overlap_by_site$bayes_obj[[1]])
  
  # If it returns a vector of overlap draws, summarise like:
  bayes_summary <- bayes_overlap_by_site %>%
   mutate(
  overlap_draws = map(bayes_obj, ~ .x[,"overlap"])  # adjust to match your output
  ) %>%
  mutate(
  overlap_med = map_dbl(overlap_draws, median, na.rm = TRUE),
  overlap_lo  = map_dbl(overlap_draws, ~ quantile(.x, 0.025, na.rm = TRUE)),
  overlap_hi  = map_dbl(overlap_draws, ~ quantile(.x, 0.975, na.rm = TRUE))
 )
 print(bayes_summary)
}


# ===============================================
# Posterior distribution of niche overlap (%) using SIBER
# Labeobarbus altianalis vs Labeo victorianus, per site
# ===============================================

remove(list = ls())

library(tidyverse)
library(readr)
library(SIBER)

# ---- Load data ----
raw <- read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---- Targets ----
target_species <- c("Labeobarbus altianalis", "Labeo victorianus")
target_sites   <- paste0("M", 4:9)

df <- raw %>%
  filter(Fish_species %in% target_species,
         Site_code %in% target_sites)

# ---- Robust column chooser ----
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) stop("None of these columns were found: ", paste(candidates, collapse = " | "))
  hit[[1]]
}

c13_name <- pick_first_col(df, c("d13C (permil, vs VPDB)", "Normalized d13C", "d13C (‰, vs VPDB)", "d13C"))
n15_name <- pick_first_col(df, c("d15N (permil, vs AIR)",  "d15N (‰, vs AIR)",  "d15N"))

# ---- Enforce numeric + drop NAs ----
df <- df %>%
  mutate(
    d13C_use = suppressWarnings(as.numeric(.data[[c13_name]])),
    d15N_use = suppressWarnings(as.numeric(.data[[n15_name]]))
  ) %>%
  filter(is.finite(d13C_use), is.finite(d15N_use))

# ---- Sample sizes and filter (Bayesian can run with n>=3, but n>=5 is more stable) ----
MIN_N <- 3
grp_sizes <- df %>% count(Site_code, Fish_species, name = "n")
print(grp_sizes)

df_ok <- df %>%
  inner_join(grp_sizes %>% filter(n >= MIN_N),
             by = c("Site_code","Fish_species"))

if (nrow(df_ok) == 0) stop("No site × species groups with n >= ", MIN_N)

# ---- Keys ----
species_key <- df_ok %>%
  distinct(Fish_species) %>% arrange(Fish_species) %>%
  mutate(group_id = row_number())

site_key <- df_ok %>%
  distinct(Site_code) %>% arrange(Site_code) %>%
  mutate(comm_id = row_number())

df_id <- df_ok %>%
  left_join(species_key, by = "Fish_species") %>%
  left_join(site_key, by = "Site_code")

# ---- Build SIBER object ----
siber_df <- df_id %>%
  transmute(
    iso1      = d13C_use,
    iso2      = d15N_use,
    group     = as.integer(group_id),
    community = as.integer(comm_id)
  ) %>%
  as.data.frame()

siber_obj <- createSiberObject(siber_df)

cat("\nSpecies key:\n"); print(species_key)
cat("\nSite key:\n");    print(site_key)
cat("\nSample sizes used by SIBER:\n"); print(siber_obj$sample.sizes)

# ---- Bayesian MVN fit (JAGS) ----
parms <- list(
  n.iter   = 20000,
  n.burnin = 5000,
  n.thin   = 10,
  n.chains = 2
)

priors <- list(
  R      = diag(2),
  k      = 2,
  tau.mu = 1.0E-3
)

ellipses_posterior <- siberMVN(siber_obj, parms, priors)

# ---- Overlap settings ----
P_ELLIPSE <- 0.40     # use 0.40 if you want to match your plotted 40% ellipses
N_POLY    <- 250      # polygon resolution for overlap
DRAWS     <- 90     # number of posterior draws to use (reduce if slow)

ls()

# group IDs for your two fishes
g_LV <- species_key %>% filter(Fish_species == "Labeo victorianus") %>% pull(group_id)
g_LB <- species_key %>% filter(Fish_species == "Labeobarbus altianalis") %>% pull(group_id)

# helper: run bayesianOverlap safely and return a tidy tibble of draws
run_bayes_overlap <- function(label_A, label_B, site_label,
                              ellipses_posterior,
                              draws = 250,
                              p.interval = 0.40,
                              n = 90) {
  
  # NOTE: call bayesianOverlap POSITIONALLY for compatibility
  bo <- bayesianOverlap(label_A, label_B,
                        ellipses_posterior,
                        draws = draws,
                        p.interval = p.interval,
                        n = n)
  
  bo_df <- as.data.frame(bo)
  nm <- names(bo_df)
  
  if (!("overlap" %in% nm)) stop("Expected column 'overlap' not found. Found: ", paste(nm, collapse=", "))
  
  if ("area1" %in% nm && "area2" %in% nm) {
    bo_df <- bo_df %>% dplyr::rename(area_A = area1, area_B = area2)
  } else if ("area.1" %in% nm && "area.2" %in% nm) {
    bo_df <- bo_df %>% dplyr::rename(area_A = `area.1`, area_B = `area.2`)
  } else {
    stop("bayesianOverlap output missing expected area columns. Found: ", paste(nm, collapse=", "))
  }
  
  bo_df %>%
    dplyr::mutate(
      Site = site_label,
      jaccard = overlap / (area_A + area_B - overlap),
      jaccard_pct = 100 * jaccard,
      p = p.interval
    ) %>%
    dplyr::select(Site, p, area_A, area_B, overlap, jaccard, jaccard_pct)
}

# ---- Compute posterior overlap draws per site ----
run_bayes_overlap(label_LV, label_LB, site_label = Site_code)
overlap_draws <- site_key %>%
  mutate(
    label_LV = paste0(comm_id, ".", g_LV),
    label_LB = paste0(comm_id, ".", g_LB)
  ) %>%
  purrr::pmap_dfr(function(Site_code, comm_id, label_LV, label_LB) {
    run_bayes_overlap(
      label_A = label_LV,
      label_B = label_LB,
      site_label = Site_code,
      ellipses_posterior = ellipses_posterior,
      draws = 250,
      p.interval = 0.40,
      n = 90
    )
  })


library(ggplot2)

ggplot(overlap_draws, aes(x = jaccard_pct)) +
  geom_density() +
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    x = "Posterior niche overlap (Jaccard, %)",
    y = "Density",
    title = "Posterior distribution of isotopic niche overlap by site"
  ) +
  theme_bw(base_size = 12)


ggplot(overlap_draws, aes(x = Site, y = jaccard_pct)) +
  geom_violin(color = "black", trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    title = paste0("Posterior niche overlap (%) by site (p = ", P_ELLIPSE, ")"),
    x = "Site",
    y = "Niche overlap (Jaccard, %)"
  ) +
  theme_minimal(base_size = 13)



overlap_summary <- overlap_draws %>%
  group_by(Site) %>%
  summarise(
    median = median(jaccard_pct, na.rm = TRUE),
    lo95   = quantile(jaccard_pct, 0.025, na.rm = TRUE),
    hi95   = quantile(jaccard_pct, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Site)

print(overlap_summary)


library(ggplot2)

ggplot(overlap_draws, aes(x = Site, y = jaccard_pct)) +
  geom_violin(
    fill = "grey80",
    color = "black",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.15,
    fill = "white",
    color = "black",
    outlier.shape = NA
  ) +
  labs(
    x = "Sampling site",
    y = "Isotopic niche overlap (Jaccard, %)",
    title = "Posterior isotopic niche overlap between fish species across sites"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid = element_blank()
  )
library(ggplot2)

ggplot(overlap_draws, aes(x = jaccard_pct)) +
  geom_density(color = "black", linewidth = 0.9) +
  
  # Median
  geom_vline(
    data = overlap_summary,
    aes(xintercept = median),
    linetype = "solid",
    linewidth = 0.8
  ) +
  
  # 95% credible interval
  geom_vline(
    data = overlap_summary,
    aes(xintercept = lo95),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  geom_vline(
    data = overlap_summary,
    aes(xintercept = hi95),
    linetype = "dashed",
    linewidth = 0.6
  ) +
  
  facet_wrap(~ Site, scales = "free_y") +
  labs(
    x = "Posterior niche overlap (Jaccard, %)",
    y = "Density",
    title = "Posterior distribution of isotopic niche overlap (core niche, p = 0.40)"
  ) +
  theme_bw(base_size = 12)


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

