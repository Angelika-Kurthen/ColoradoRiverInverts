###########################################
## Code to create a ramping temp increase
###########################################




library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)

flow.magnitude <- TimestepDischarge(flow, 85000) #discharge data from USGS is in cfs, bankfull dischage (pers comm TK) 85000 cfs

temp <- temp <- readNWISdv("09380000", "00010", "2007-10-01", "2021-09-30")

temp$Date <- as_datetime(temp$Date)
temp$Date <- yday(temp$Date)


# option 1: for each date, calculate mean and standard error so we can create a new dataset
temp <- temp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(X_00010_00003), 
                                                     sd = sd(X_00010_00003), 
                                                     count = n(), 
                                                     se = sd(X_00010_00003)/count)
temp







# Make an index to be used for aggregating
ID <- as.numeric(as.factor(temp$Date)) 
# want it to be every 14 days, hence the 14
ID <- ID %/% 14

ID[which(ID ==26)] <- 25
# aggregate over ID and TYPEall numeric data.
outs <- aggregate(temp[sapply(temp,is.numeric)],
                  by=list(ID),
                  FUN=mean)



# format output
names(outs)[2:3] <-c("dts","Temperature")
# add the correct dates as the beginning of every period
outs$dts <- strptime(round(outs$dts), "%j") ###Note need to subtract 365 if larger than 365

# order by date in chronological order#
#outs <- outs[order(outs$dts),]
outs$dts <- as_date(outs$dts)

n <- 13

# repeat this data frame for 13 years
temp_seq <- do.call("rbind", replicate(n, outs, simplify = FALSE))
temp_seq$dts <- as.Date(temp_seq$dts)
# now adjust the years so time can proceed in chronological order
year(temp_seq$dts[1:26]) <- 2000
year(temp_seq$dts[27:52]) <- 2001
year(temp_seq$dts[53:79]) <- 2002
year(temp_seq$dts[80:104]) <- 2003
year(temp_seq$dts[105:130]) <- 2004
year(temp_seq$dts[131:156]) <- 2005
year(temp_seq$dts[157:182]) <- 2006
year(temp_seq$dts[183:208]) <- 2007
year(temp_seq$dts[209:234]) <- 2008
year(temp_seq$dts[235:260]) <- 2009
year(temp_seq$dts[261:286]) <- 2010
year(temp_seq$dts[287:312]) <- 2011
year(temp_seq$dts[313:338]) <- 2012


# adjust the year temps

# first five years remain at mean temp for 2007 - 2022, which we will call baseline

# years 6 and 7  represent a 1 C increase from baseline
temp_seq$Temperature[131:182] <- temp_seq$Temperature[131:182] + 1

# years 8 and 9 represent a 2.5 C increase from baseline
temp_seq$Temperature[183:234] <- temp_seq$Temperature[183:234] + 2.5

# years 10 and 11 represent at 5 C increase from baseline
temp_seq$Temperature[235:286] <- temp_seq$Temperature[235:286] + 5

# year 11 and 12 represent a 7.5 C increase from baseline
temp_seq$Temperature[287:338] <- temp_seq$Temperature[287:338] + 7.5


arrows <- tibble(
  x1 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  x2 = c("2005-01-07", "2007-01-07", "2009-01-07", "2011-01-07"),
  y1 = c(20, 20, 20, 20), 
  y2 = c(7, 7, 7, 7)
)
arrows$x1 <- as.Date(arrows$x1)
arrows$x2 <- as.Date(arrows$x2)

ggplot(data = temp_seq, mapping = aes(x = dts, y = Temperature ))+
  geom_line()+
  ylab("Water Temperature (°C)")+
  xlab(" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13))+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")+
  annotate("segment", x = arrows$x1, y = arrows$y1, xend = arrows$x2, yend = arrows$y2, color = "red" )+
  annotate("text", x = arrows$x1[1]+7, y = 21, label = "+1°C", size = 5)+
  annotate("text", x = arrows$x1[2]+7, y = 21, label = "+2.5°C", size = 5)+
  annotate("text", x = arrows$x1[3]+7, y = 21, label = "+5°C", size = 5)+
  annotate("text", x = arrows$x1[4]+7, y = 21, label = "+7.5°C", size = 5 )





