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
# ============================================================
# SIBER vignette workflow ADAPTED to YOUR data
# Communities = Sites (M4–M9)
# Groups      = Species (2 fish species)
# ============================================================

rm(list = ls())
graphics.off()
set.seed(1)

library(tidyverse)
library(readr)
library(SIBER)
library(hdrcde)

# ---------------------------
# 1) LOAD YOUR DATA
# ---------------------------
raw <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=698972139&single=true&output=csv",
  show_col_types = FALSE
)

# ---------------------------
# 2) TARGETS (match your setup)
# ---------------------------
species_levels <- c("Labeobarbus altianalis", "Labeo victorianus")  # group order
site_levels    <- paste0("M", 4:9)                                  # community order

# ---------------------------
# 3) HELPERS
# ---------------------------
pick_first_col <- function(dat, candidates) {
  hit <- intersect(candidates, names(dat))
  if (length(hit) == 0) stop("None of these columns were found: ", paste(candidates, collapse = " | "))
  hit[[1]]
}

# Fix for SIBER sometimes returning SEA.B with NULL colnames:
# Create "community.group" labels in the SAME order as SIBER output.
label_seab_cols_if_missing <- function(SEA.B, siber_obj, min_n) {
  if (!is.null(colnames(SEA.B))) return(as.matrix(SEA.B))
  
  ss <- siber_obj$sample.sizes
  present <- which(!is.na(ss) & ss >= min_n, arr.ind = TRUE)
  
  if (nrow(present) == 0) stop("No valid community×group combos with n >= ", min_n)
  
  # SIBER order is community-major then group
  present <- present[order(present[, 1], present[, 2]), , drop = FALSE]
  labels  <- apply(present, 1, function(ix) paste(ix[1], ix[2], sep = "."))
  
  if (ncol(SEA.B) != length(labels)) {
    cat("\nDEBUG mismatch in SEA.B labeling\n")
    cat("ncol(SEA.B) =", ncol(SEA.B), "\n")
    cat("length(labels) =", length(labels), "\n")
    cat("labels:\n"); print(labels)
    stop("Cannot label SEA.B: mismatch between SEA.B columns and present community×group combos.\n",
         "Try raising min_n (e.g., 6–10) or check sample sizes.")
  }
  
  SEA.B <- as.matrix(SEA.B)
  colnames(SEA.B) <- labels
  SEA.B
}

# Lookup: "community.group" -> Site + Species labels
make_lookup <- function(col_ids, site_levels, species_levels) {
  tmp <- strsplit(col_ids, "\\.")
  community <- as.integer(vapply(tmp, `[`, "", 1))
  group     <- as.integer(vapply(tmp, `[`, "", 2))
  
  tibble(
    group_id     = col_ids,
    community    = community,
    group        = group,
    Site_code    = site_levels[community],
    Fish_species = species_levels[group]
  )
}

# ---------------------------
# 4) FILTER + CLEAN ISOTOPES
# ---------------------------
df0 <- raw %>%
  dplyr::filter(Fish_species %in% species_levels,
                Site_code %in% site_levels)

c13_name <- pick_first_col(df0, c(
  "d13C (permil, vs VPDB)",
  "Normalized d13C",
  "d13C (‰, vs VPDB)",
  "d13C", "d13C_corrected", "C13"
))
n15_name <- pick_first_col(df0, c(
  "d15N (permil, vs AIR)",
  "d15N (‰, vs AIR)",
  "d15N", "N15"
))

df <- df0 %>%
  transmute(
    Site_code    = as.character(Site_code),
    Fish_species = as.character(Fish_species),
    d13C_use     = suppressWarnings(as.numeric(.data[[c13_name]])),
    d15N_use     = suppressWarnings(as.numeric(.data[[n15_name]]))
  ) %>%
  drop_na(Site_code, Fish_species, d13C_use, d15N_use)

# ---------------------------
# 5) MINIMUM SAMPLE SIZE FILTER
# ---------------------------
min_n <- 5  # recommended; can drop to 3 if absolutely necessary

df_ok <- df %>%
  group_by(Site_code, Fish_species) %>%
  filter(n() >= min_n) %>%
  ungroup()

if (nrow(df_ok) == 0) stop("No Site×Species groups with n >= ", min_n, ". Try min_n <- 3 if needed.")

# stable coding (communities=sites; groups=species)
df_ok <- df_ok %>%
  mutate(
    Site_code    = factor(Site_code, levels = site_levels),
    Fish_species = factor(Fish_species, levels = species_levels),
    community    = as.integer(Site_code),
    group        = as.integer(Fish_species)
  )

# SIBER input frame
siber_df <- df_ok %>%
  transmute(
    iso1      = d13C_use,
    iso2      = d15N_use,
    group     = group,
    community = community
  ) %>%
  as.data.frame()

# ---------------------------
# 6) CREATE SIBER OBJECT (this replaces siber.example)
# ---------------------------
siber.example <- createSiberObject(siber_df)

cat("\nSample sizes (rows=community/site index, cols=group/species index):\n")
print(siber.example$sample.sizes)
cat("\nCommunity index -> Site:\n"); print(tibble(community = seq_along(site_levels), Site = site_levels))
cat("\nGroup index -> Species:\n");   print(tibble(group = seq_along(species_levels), Species = species_levels))

# ============================================================
# 7) PLOTS (same as vignette, using your siber.example object)
# ============================================================

community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.95, lty = 1, lwd = 2)
group.hulls.args     <- list(lty = 2, col = "grey20")

par(mfrow = c(1,1))
plotSiberObject(
  siber.example,
  ax.pad = 2,
  hulls = FALSE, community.hulls.args = community.hulls.args,
  ellipses = TRUE, group.ellipses.args = group.ellipses.args,
  group.hulls = TRUE, group.hulls.args = group.hulls.args,
  bty = "L",
  iso.order = c(1,2),
  xlab = expression({delta}^13*C~"‰"),
  ylab = expression({delta}^15*N~"‰")
)

# smaller points plot (like vignette)
group.hull.args <- list(lty = 2, col = "grey20")
par(mfrow = c(1,1))
plotSiberObject(
  siber.example,
  ax.pad = 2,
  hulls = FALSE, community.hulls.args,
  ellipses = FALSE, group.ellipses.args,
  group.hulls = FALSE, group.hull.args,
  bty = "L",
  iso.order = c(1,2),
  xlab = expression({delta}^13*C~"‰"),
  ylab = expression({delta}^15*N~"‰"),
  cex  = 0.5
)

# ============================================================
# 8) >>> THIS IS THE BLOCK YOU SAID WAS MISSING <<<
#    GROUP ML METRICS + ELLIPSES + COMMUNITY METRICS
# ============================================================

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
cat("\nGroup-level ML metrics (TA, SEA, SEAc). Columns = community.group:\n")
print(group.ML)

# Add a prediction ellipse
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  lty = 1, lwd = 2)

# Add CI around bivariate means
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95, ci.mean = TRUE,
                  lty = 1, lwd = 2)

# Plot convex hulls (community level), like vignette
par(mfrow = c(1,1))
plotSiberObject(
  siber.example,
  ax.pad = 2,
  hulls = TRUE, community.hulls.args,
  ellipses = FALSE, group.ellipses.args,
  group.hulls = FALSE, group.hull.args,
  bty = "L",
  iso.order = c(1,2),
  xlab = expression({delta}^13*C~"‰"),
  ylab = expression({delta}^15*N~"‰"),
  cex  = 0.5
)

# Optionally add CI ellipses on top of hull plot
plotGroupEllipses(siber.example, n = 100, p.interval = 0.95,
                  ci.mean = TRUE, lty = 1, lwd = 2)

# Community-level Layman metrics
community.ML <- communityMetricsML(siber.example)
cat("\nCommunity-level ML Layman metrics (per community/site index):\n")
print(community.ML)

# ============================================================
# 9) BAYESIAN SETTINGS + POSTERIOR (same as vignette)
#    (You can keep this if you want the calculations; plots are optional)
# ============================================================

parms <- list()
parms$n.iter   <- 2 * 10^4
parms$n.burnin <- 1 * 10^3
parms$n.thin   <- 10
parms$n.chains <- 2

priors <- list()
priors$R      <- 1 * diag(2)
priors$k      <- 2
priors$tau.mu <- 1.0E-3

ellipses.posterior <- siberMVN(siber.example, parms, priors)

# SEA.B (posterior draws)
SEA.B <- siberEllipses(ellipses.posterior)
SEA.B <- label_seab_cols_if_missing(SEA.B, siber.example, min_n)

# Create nice readable labels (Site | Species) in the same order as SEA.B
lookup <- make_lookup(colnames(SEA.B), site_levels, species_levels)
xticks <- paste0(lookup$Site_code, " | ", lookup$Fish_species)

# ---- OPTIONAL PLOT (you said you don't need the curves; so it’s optional) ----
# siberDensityPlot(SEA.B, xticklabels = xticks,
#                  xlab = "Site | Species",
#                  ylab = expression("Standard Ellipse Area " ("‰"^2)),
#                  bty = "L", las = 2,
#                  main = "SIBER ellipses on each group (SEA.B)")

# Add red x's for ML SEAc (matched safely)
# (This still computes even if you don’t plot)
seac_vec <- group.ML["SEAc", ]
seac_vec <- seac_vec[match(colnames(SEA.B), colnames(group.ML))]

# if you plot, then uncomment:
# points(1:ncol(SEA.B), seac_vec, col = "red", pch = "x", lwd = 2)

# Credible intervals and modes (CALCULATIONS)
cr.p <- c(0.95, 0.99)

SEA.B.credibles <- lapply(
  as.data.frame(SEA.B),
  function(x, ...) { hdrcde::hdr(x)$hdr },
  prob = cr.p
)

SEA.B.modes <- lapply(
  as.data.frame(SEA.B),
  function(x, ...) { hdrcde::hdr(x)$mode },
  prob = cr.p, all.modes = TRUE
)

# Posterior means (needed for bayesianLayman calcs)
mu.post <- extractPosteriorMeans(siber.example, ellipses.posterior)

# Bayesian Layman metric distributions
layman.B <- bayesianLayman(mu.post)

# ---- OPTIONAL: if you want to plot Layman.B (curves), uncomment ----
# for (i in seq_along(layman.B)) {
#   siberDensityPlot(layman.B[[i]], xticklabels = colnames(layman.B[[i]]),
#                    bty = "L", ylim = c(0, 20),
#                    main = paste0("Layman metrics: Site ", site_levels[i]))
# }

# ---- OPTIONAL: TA compare first two communities only (if they exist) ----
# if (length(layman.B) >= 2) {
#   par(mfrow=c(1,1))
#   siberDensityPlot(cbind(layman.B[[1]][,"TA"], layman.B[[2]][,"TA"]),
#                    xticklabels = c(paste0(site_levels[1]), paste0(site_levels[2])),
#                    bty="L", ylim=c(0,20), las=1,
#                    ylab="TA - Convex Hull Area", xlab="")
# }

# ============================================================
# 10) OUTPUT TABLES WITH HUMAN-READABLE LABELS (very useful!)
# ============================================================

# group.ML columns are community.group -> convert to Site + Species
groupML_tbl <- as.data.frame(t(group.ML)) %>%
  rownames_to_column("community_group") %>%
  separate(community_group, into = c("community", "group"), sep = "\\.", convert = TRUE) %>%
  mutate(
    Site_code    = site_levels[community],
    Fish_species = species_levels[group]
  ) %>%
  select(Site_code, Fish_species, TA, SEA, SEAc) %>%
  arrange(Site_code, Fish_species)

cat("\nGroup ML metrics table (Site × Species):\n")
print(groupML_tbl)

# community.ML columns are community index -> attach Site labels
communityML_tbl <- as.data.frame(t(community.ML)) %>%
  rownames_to_column("community") %>%
  mutate(
    community = as.integer(community),
    Site_code = site_levels[community]
  ) %>%
  select(Site_code, everything(), -community) %>%
  arrange(Site_code)

cat("\nCommunity ML Layman metrics table (Site):\n")
print(communityML_tbl)

