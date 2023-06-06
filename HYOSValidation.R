######################################
## Code to Validate 
######################################

library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#install.packages("devtools")
library(devtools)
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "GCgage")
library(GCgage)

source("1spFunctions.R")
# load Water Temperature data from above Diamond Creek Confluence (RM226)
temp <- read.delim("CRaboveDC_Temp.tsv", header=T)
colnames(temp) <- c("Date", "Temperature")
temp$Date <- as.Date(temp$Date, format = "%Y-%m-%d")
temp <- aggregate(temp$Temperature, by = list(temp$Date), FUN = mean)
colnames(temp) <- c("Date", "Temperature")
temps <- TimestepTemperature(temp) # intermittantly missing data until  2000-12-26 
temps <- temps[-c(1:82),]
# load discharge data from above Diamond Creek Confluence (RM226)
discharge <- readNWISdv("09404200", "00060", "2000-12-26", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)

source("HYOS_1sp.R")
out <- HYOSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 10000, Qmin = 0.25, extinct = 500, iteration = 1, peaklist = 0.27, peakeach = length(temps$Temperature))
