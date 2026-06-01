remove(list = ls())
# Setup renv (optional)
renv::restore()
# 1. Load packages
library(readxl)
library(tidyverse)
library(janitor)
library(vegan)
library(lme4)
library(DHARMa)
library(FSA)
library(indicspecies)

# 2. Import data
beetles <- read_excel("Ground beetles_Gloria.xlsx", sheet = "Raw data-1") %>%
  clean_names()

# 3. Clean data
beetles_clean <- beetles %>%
  filter(!is.na(species), !is.na(abundance)) %>%
  mutate(
    species = str_replace_all(species, "\\.", "_"),
    species = str_replace_all(species, " ", "_"),
    zone = as.factor(zone),
    season = as.factor(season),
    plot = as.factor(plot),
    plot2 = as.factor(plot2),
    method_used = as.factor(method_used)
  )

# Check data structure
str(beetles_clean)
summary(beetles_clean)