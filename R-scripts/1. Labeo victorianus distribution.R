#Fish assemblages in the Mara River
# Load necessary libraries
library(tidyverse)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(car)  # For Levene's test
library(rstatix)  # For normality and variance tests

# Read data from Google Sheets
watermetals <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vSFJjla0jcvd2MqIFU0oPlIMGqLM5QQJnN0C9ognk6WRrqxtYnZgCRZ-14J2S1mtl4NTUqb7y16bWYm/pub?gid=1141386352&single=true&output=csv",
  show_col_types = FALSE
) %>%
  filter(Location_ID %in% paste0("M", 1:9), sampling_year %in% 2021:2023)
