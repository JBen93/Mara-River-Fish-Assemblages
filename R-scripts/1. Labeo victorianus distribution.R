#Labeo victorianus distribution analysis in the Mara River
# Load necessary libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)  # For Levene's test
library(rstatix)  # For normality and variance tests

#datasource

# Read data from Google Sheets
currentfish<- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=152464398&single=true&output=csv",
  show_col_types = FALSE
) %>%
  filter(Location_ID %in% paste0("M", 1:9), sampling_year %in% 2021:2023)
