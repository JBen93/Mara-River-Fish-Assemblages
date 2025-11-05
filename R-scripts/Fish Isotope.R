#Fish Isotope using SIBER R package
# clear everything in memory (of R)
remove(list=ls())
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