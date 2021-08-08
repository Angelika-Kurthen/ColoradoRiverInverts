
install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
library(tidyr)
library(ggplot2)
drift <- readDB(gear = "Drift", type = "Sample", updater = T)

# data from the upper thirds of the river

# NZMS = New Zealand Mudsnail
# GAMM = Gammarus lacustris
# 
# SIML = Simuliidae
# 
species <- c("NZMS","GAMM", "SIML")


upper <- drift[drift$RiverMile <= 75, ]
upper <- sampspec(samp = upper)

# data from lower third of the river
lower <- drift[drift$RiverMile >= 150, ]
lower <- sampspec(samp = lower)

# pull species specific data
upper <- sampspec(samp = upper, species = "NZMS", stats = T)
# get species data
upper_sp <-  upper$SpecDel
# get rid of NAs
upper_sp <- upper_sp[-which((is.na(upper_sp$CountTotal == "NA"))== T),] 

# make columns rows and rows columns
LFdriftUpper_NZMS_size <- gather(upper_sp[ , 3:17], key = "Size")

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
