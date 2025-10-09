#Labeo victorianus distribution analysis in the Mara River
# clear everything in memory (of R)
remove(list=ls())

renv::restore()
# =========================================================
# Labeo victorianus distribution & abundance (M1–M9, 2021–2022)
# =========================================================

# (optional) start fresh
# rm(list = ls())

# ---- Packages ----
libs <- c("tidyverse","readr","dplyr","ggplot2","rstatix",
          "MASS","emmeans","performance","DHARMa")
to_install <- libs[!libs %in% installed.packages()[,"Package"]]
if(length(to_install)) install.packages(to_install, dependencies = TRUE)
lapply(libs, require, character.only = TRUE)

# ---- Data ----
currentfish <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) |>
  filter(Location_ID %in% paste0("M", 1:9),
         sampling_year %in% 2021:2022)

# Keep only Labeo victorianus
currentfish_lv <- currentfish |>
  filter(Species == "Labeo victorianus")

# ---- Choose the abundance column (edit if needed) ----
# Try to detect a likely abundance/count column
cand_cols <- c("Abundance","Count","Counts","N","Number")
abundance_col <- cand_cols[cand_cols %in% names(currentfish_lv)][1]
if (is.na(abundance_col)) {
  stop("Couldn't find an abundance column. Rename your abundance field to one of: ",
       paste(cand_cols, collapse=", "))
}

# Optional: account for sampling effort if available (e.g., time or area)
# Set `effort_col <- "EffortMinutes"` (or similar) if your data has it; else NA
effort_col <- NA

# Coerce factors for modeling
dat <- currentfish_lv |>
  mutate(
    Site = factor(Location_ID, levels = paste0("M", 1:9)),
    Year = factor(sampling_year)
  ) |>
  select(Site, Year, all_of(abundance_col), any_of(effort_col)) |>
  rename(Abund = all_of(abundance_col))

# Basic check
stopifnot(all(levels(dat$Site) %in% paste0("M", 1:9)))
stopifnot(all(levels(dat$Year) %in% c("2021","2022")))

# =========================================================
# 1) Descriptives & mean ± SE by Site × Year
# =========================================================
sum_by_site_year <- dat |>
  group_by(Site, Year) |>
  summarise(
    n_samples = n(),
    mean_abund = mean(Abund, na.rm = TRUE),
    sd_abund   = sd(Abund, na.rm = TRUE),
    se_abund   = sd_abund / sqrt(n_samples),
    .groups = "drop"
  )

# =========================================================
# 2) Model: Negative binomial GLM (robust for overdispersion)
#    Abundance ~ Site * Year  (+ optional offset for effort)
# =========================================================
if (!is.na(effort_col) && effort_col %in% names(dat)) {
  dat <- dat |> mutate(offset_term = log(.data[[effort_col]]))
  m_nb <- MASS::glm.nb(Abund ~ Site * Year + offset(offset_term), data = dat)
} else {
  m_nb <- MASS::glm.nb(Abund ~ Site * Year, data = dat)
}

# Model checks
check_overdispersion(m_nb)
check_zeroinflation(m_nb)
# Residual simulation
sim <- DHARMa::simulateResiduals(m_nb, n = 1000)
plot(sim)

# Type II/III tests via car::Anova if desired
# car::Anova(m_nb, type = 2)

# =========================================================
# 3) Post-hoc: 2021 vs 2022 within each site (EMMs)
# =========================================================
emm <- emmeans(m_nb, ~ Year | Site, type = "response")
contrast_site <- contrast(emm, "pairwise", adjust = "holm")  # 2021 vs 2022 per Site
print(contrast_site)

# Also get model-predicted means for plotting
pred <- emmeans(m_nb, ~ Site * Year, type = "response") |>
  as.data.frame() |>
  rename(pred_mean = response,
         lcl = asymp.LCL,
         ucl = asymp.UCL)

# =========================================================
# 4) Plots
#    A) Empirical means ± SE
#    B) Model-predicted means ± 95% CI
# =========================================================

# A) Sample means ± SE by Site × Year
p_empirical <- ggplot(sum_by_site_year,
                      aes(x = Site, y = mean_abund, fill = Year)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_errorbar(aes(ymin = mean_abund - se_abund,
                    ymax = mean_abund + se_abund),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "Labeo victorianus – Mean Abundance (±SE) by Site & Year",
       x = "Site (M1–M9)", y = "Mean abundance (samples)") +
  theme_minimal(base_size = 12) +
  scale_fill_brewer(palette = "Set2")

# B) Model-predicted means ± 95% CI
p_model <- ggplot(pred, aes(x = Site, y = pred_mean, color = Year, group = Year)) +
  geom_point(position = position_dodge(width = 0.4), size = 2.8) +
  geom_errorbar(aes(ymin = lcl, ymax = ucl),
                position = position_dodge(width = 0.4), width = 0.2) +
  labs(title = "Labeo victorianus – Model-Predicted Abundance (95% CI)",
       x = "Site (M1–M9)", y = "Predicted mean abundance (NB-GLM)") +
  theme_minimal(base_size = 12) +
  scale_color_brewer(palette = "Dark2")

print(p_empirical)
print(p_model)

# =========================================================
# 5) Optional: quick per-site test (2021 vs 2022) using NB GLM
# =========================================================
per_site_tests <- lapply(levels(dat$Site), function(s) {
  d <- dat |> filter(Site == s)
  fit <- try(MASS::glm.nb(Abund ~ Year, data = d), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  broom::tidy(fit, conf.int = TRUE)
})
# per_site_tests is a list of tidy summaries by site
