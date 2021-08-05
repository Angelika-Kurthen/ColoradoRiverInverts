
install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
library(tidyr)
library(ggplot2)
drift <- readDB(gear = "Drift", type = "Sample", updater = T)

drift <- readDB(gear = "Drift", type = "Sample", updater = T)
# data from the upper thirds of the river
LFdriftUpper <- drift[drift$RiverMile <= 75, ]
LFdriftUpper <- sampspec(samp = LFdriftUpper)
# data from lower third of the river
LfdriftLower <- drift[drift$RiverMile >= 150, ]
LfdriftLower <- sampspec(samp = LfdriftLower)

# pull NZMS data
LFdriftUpper_NZMS <- sampspec(samp = LFdriftUpper, species = "NZMS", stats = T)
# get species data
LFdriftUpper_NZMS_SpecDel <- LFdriftUpper_NZMS$SpecDel
# get rid of NAs
LFdriftUpper_NZMS_SpecDel <- LFdriftUpper_NZMS_SpecDel[-which((is.na(LFdriftUpper_NZMS_SpecDel$CountTotal == "NA"))== T),] 

# make columns rows and rows columns
LFdriftUpper_NZMS_size <- gather(LFdriftUpper_NZMS_SpecDel[ , 3:17], key = "Size")

# remove counts of 0 (since they don't need to be in out distribution)
LFdriftUpper_NZMS_size <- LFdriftUpper_NZMS_size[which(LFdriftUpper_NZMS_size$value != 0), ]

# do this for all 13 size classes
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B0"), 0)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B1"), 1)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B2"), 2)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B3"), 3)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B4"), 4)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B5"), 5)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B6"), 6)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B7"), 7)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B8"), 8)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B9"), 9)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B10"), 10)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B11"), 11)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B12"), 12)
LFdriftUpper_NZMS_size$Size <- replace(LFdriftUpper_NZMS_size$Size, which(LFdriftUpper_NZMS_size$Size == "B13"), 13)

size_upper <- vector()
for (i in 1:length(LFdriftUpper_NZMS_size$Size)){
  a <- as.vector(rep(LFdriftUpper_NZMS_size$Size[i], times = LFdriftUpper_NZMS_size$value[i]))
  size_upper <- c(size_upper, a)
}

location <- rep("RM -15 to 75", length(size_upper))
df <- as.data.frame(cbind(size_upper, location))
df$size_upper <- as.numeric(as.character(df$size_upper))

ggplot(data = df, aes(x = location, y = size_upper))+
  geom_violin(adjust = 4, trim = F, draw_quantiles = c(.25, .50, .75))

#want size proportion - how to standardize by Volume? will not get a whole number?

# so in each row, we want to get the proportion of 

LfdriftLower_NZMS <- sampspec(samp = LfdriftLower, species = "NZMS", stats = T)
