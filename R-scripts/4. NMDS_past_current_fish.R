# Compare the historical and current fish community structure in the Mara River using NMDS, DCA, and PCoA to determine the temporal changes in the community structure.

# clear everything in memory (of R)
remove(list=ls())

renv::restore()

library(vegan) 
library(psych) 
library(tidyverse)
library(ggpubr)
library(janitor) # for clean_names
library(stringr) # for str_remove

# Past fish data (2013, 2014, 2016)
# --------------------------------------------------------------
# Past fish composition
# --------------------------------------------------------------
#NMDS analysis of past fish data
#database source 
#data URL source if you need to inspect for the whole dataset
#browseURL("https://docs.google.com/spreadsheets/d/1WsfU7zcpAl_Kwg9ZxoSvGHgYrnjZlQVs4zzcHqEHkaU/edit?usp=sharing")

#load data filter 2008,2009 and also group by Location_ID, month, year,River_reach and Family
# ---- Load past data ----
#set working directory
setwd("~/Desktop/Git hub/Mara-River-Fish-Assemblages")
file <- "/mnt/data/Historical fish.csv"   # <- replace if needed
pastfish <- readr::read_csv(
  "https://docs.google.com/spreadsheets/d/e/2PACX-1vRDo5laGSxF444O2xpHBPq4papf5IJd5VQ6BOFoUKGZIZZRqAp5gHsWrWfv-P3A2OBeJUH16Gn4N_ng/pub?gid=983226609&single=true&output=csv",
  show_col_types = FALSE
) %>%
  janitor::clean_names()   # -> location_id, sampling_year, fish_species, fish_weight, total_length, standard_length, ...
head(pastfish)

pastfish2<-pastfish %>% select(-c(Location_ID, month, year, Reach)) #remove the columns that are not needed for the analysis
head(pastfish2)

#PerMANOVA test to determine if the Location_ID, month, year, and River_reach are significant factors in explaining the variation in the macroinvertebrate community structure.
set.seed(123) #set seed for reproducibility
permmacros <- adonis2(historicalmacros2 ~ Location_ID + month + year + Reach, 
                      data = historicalmacros, 
                      method = "euclidean", permutations = 999,
                      by = "margin")
permmacros

#the PerMANOVA test for the historical data shows that the Location_ID, month, year, and River_reach are not significant factors in explaining the variation in the macroinvertebrate community structure.

#plot the NMDS ordination of the historical macroinvertebrate community structure
histomacros<-vegdist(historicalmacros2, "bray") #we use bray curtis distance matrix since it's count data
histomacros

#You are going to use the metaMDS function in the vegan package.
#K =2 because we are interested in only two dimension (which is common for NMDS)
#Trymax=100 because we have a small data set
nmdshisto<-metaMDS(histomacros,k=2, trace=T, trymax=100)
nmdshisto

stressplot(nmdshisto) #plot the stress plot
#What do the stress value and the fit (R2) of the monotonic regression tell you about the NMDS plot?
#The stress value is  0.092  which implies fairly good fit. The stress value is a measure of how well the data are represented in the NMDS plot.
#The closer the stress value is to 0, the better the representation.
#The R2 value is 0.992, which is also very good. The R2 value is a measure of how well the data are represented in the NMDS plot.
#The closer the R2 value is to 1, the better the representation.

#First, we need to create a grouping variable for the historical groups. #We will use the ca package to do this. Identify the river reach, month, year and site as groups: #River reach as a grouping variable
#First, we need to extract the site scores from the historical NMDS object.
#We will use the scores function in the vegan package to do this.
#We will also convert the site scores to a data frame.
histodata.scores <- as.data.frame(scores(nmdshisto)) #Using the scores function from vegan to extract the site scores and convert to a data.frame
histodata.scores$sites <- rownames(historicalmacros) # create a column of River reach from the original data frame macros
histodata.scores$Reach<-(historicalmacros$Reach) # create a column of 
head(histodata.scores)  #look at the data


#Plot NMDS for the past fish community structure

ggplot(histodata.scores, aes(x = NMDS1, y = NMDS2, col = Reach)) + 
  geom_point() +
  stat_ellipse(linetype = "dashed") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-2.5, 2.5) +
  ylim(-2, 3) +
  labs() +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", x = 2, y = 3, label = paste("Stress = ", round(nmdshisto$stress, 3)))

