#Labeo victorianus distribution analysis in the Mara River
# clear everything in memory (of R)
remove(list=ls())

renv::restore()
# Load necessary libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)  # For Levene's test
library(rstatix)  # For normality and variance tests

#database source
#browseURL("https://docs.google.com/spreadsheets/d/1J78CV0XTYJPNbSXnmilWQUgDeqAm1A32LR4A2z9D4vw/edit?usp=sharing")

# Read data from Google Sheets
currentfish<- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  filter(Location_ID %in% paste0("M", 1:9), sampling_year %in% 2021:2022)
#filter on Labeo victorianus species for the current fish data for the sites M1-M9 and years 2021-2022.
currentfish_lv <- currentfish %>%
  filter(Species == "Labeo victorianus")


