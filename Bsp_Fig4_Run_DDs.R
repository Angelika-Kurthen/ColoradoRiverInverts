##Baetid Figure 4: Importance of Disturbance Timing
# Importance of Disturbance Timing is just comparing delta(max) to delta(max)
# could also be a value (min(max)/max(max)) <- the larger the number, the bigger the difference between max values
#max(jday.df$jday_max)-min(jday.df$jday_max)
#vs
#min(jday.df$jday_max)/max(jday.df$jday_max)

# Code for HPC - tidyverse has some issues on our HPC because one of the packages is deprecated
# We have to manually load all tidyverse packages
library(purrr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(tibble, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(tidyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(readr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(stringr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(forcats, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(plyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(dplyr, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(ggplot2, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(dataRetrieval, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

source("B_1sp_Model.R")
source("1spFunctions.R")
Time <- c(1:36500)
Date <- rep(1:365, times = 100)
Day <- seq(as.Date("2022-01-01"), as.Date("2121-12-31"), by="days")
# find and remove leap days
find_leap = function(x){
  day(x) == 29 & month(x) == 2 
}
Day <- Day[which(find_leap(Day) == F)]


Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]

# we can now iterate (similar to a sensitivity analysis) and see how importance of disturbance timing impacts 
# what if we iterated fecundities 
dds <- seq(375,625, by = 12.5)
# # Initializes the progress bar
 pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(s]),  # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")   # Character used to create the bar
all.dates <- unique(format(temp$dts, "%m-%d"))[order(unique(format(temp$dts, "%m-%d")))]
IDT <- vector()
IDT.ratio <- vector()

for (f in 1:length(dds)) {
  # Sets the progress bar to the current state
   setTxtProgressBar(pb, f)
  # loop to select a date from a Week-Month combo from each unique year
  means <- list()
  for (d in 1:length(all.dates)){ # 30 reps takes 60 mins
    sample_dates <- temp$dts[which(format(temp$dts, "%m-%d")== all.dates[d])]
    samp <- which(temp$dts == sample(sample_dates[which(sample_dates > temp$dts[300] & sample_dates < temp$dts[2508])], size = 1))
    dates <- temp[(samp-300):(samp+100),]
    discharge <- rep(0.1, time = length(dates$dts)) # create a list of non-disturbance discharges
    discharge[match(samp, dates$dts)] <- 0.3 # from that list of dates from above, assign a disturbance discharge to that date
    print(length(dates$dts))
    source("B_1sp_Model.R")
    # run model
    out <- Bmodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature), dds = dds[f])
    # create summary dataframe 
    m <- mean.data.frame(out, burnin = 250, iteration = 9)
    means[d] <- list(m)
  }
  jday_max <- unlist(lapply(means, function(x)
    return(max(x$mean.abund)))) # would mean +/- se or just maximum (and maybe minimum values) be valuable
  jday.df <- as.data.frame(cbind(jday_max, all.dates))
  jday.df$jday_max <- as.numeric(jday.df$jday_max)
  jday.df$all.dates <- yday(as.Date(jday.df$all.dates, "%m-%d"))
  IDT[f] <- max(jday.df$jday_max)-min(jday.df$jday_max)
  IDT.ratio[f] <- min(jday.df$jday_max)/max(jday.df$jday_max)
  close(pb) # close progress bar
}

IDT.df <- as.data.frame(cbind(IDT, dds))
IDT.df1 <- ggplot(data = IDT.df, aes(x = dds, y = IDT))+
  geom_line(linewidth = 1)+
  ylab('Importance of Disturbance Timing') +
  xlab("Mean Degree Days")+
  theme_bw()

#IDT.ratio.df1 <- IDT.ratio.df <- as.data.frame(cbind(IDT.ratio, dds))
#ggplot(data = IDT.ratio.df, aes(x = dds, y = IDT.ratio))+
#  geom_line(linewidth = 1)+
#  ylab('Importance of Disturbance Timing')+
#  xlab("Mean Fecundity")+
#  theme_bw()
  
pdf("IDT_DegreeDay")
print(IDT.df1)

dev.off()
