# Compare the historical and current macroinvertebrate community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.
# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)