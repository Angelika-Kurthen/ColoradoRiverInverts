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

source("1spFunctions.R")
source("HYOS_1sp.R")
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

out <- HYOSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 5000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))

drift.data.total <- readDB(gear = "LightTrap", type = "Sample", updater = TRUE)

drift.CR <- drift.data.total[which(drift.data.total$Region == "GrandCanyon" & drift.data.total$RiverMile >= 119 & drift.data.total$RiverMile <= 225),]

specieslist <- c("CHEP","HYOS", "HYSP")
HYOS.samp.CR <- sampspec(samp = drift.CR, stats = T)
HYOS.samp <- merge(HYOS.samp.CR$Statistics, HYOS.samp.CR$Samples, by = "BarcodeID", all = T)
HYOS.samp$Density <- HYOS.samp$CountTotal/HYOS.samp$TimeElapsed
#NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
#NZMS.samp$Density <- (NZMS.samp$Density/(9e-15 *(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))
HYOS.samp <- aggregate(HYOS.samp$Density, list(HYOS.samp$Date), FUN = mean)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Group.1 >= temps$dts[i] & HYOS.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na((d$x) == T))) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means[length(temps$dts)] <- HYOS.samp$x[length(HYOS.samp$x)]
#
# means.list.HYOS <- mean.data.frame(out, burnin = 286, iteration= 9)
# means.list.HYOS <- cbind(means.list.HYOS, temps$dts[286:last(means.list.HYOS$timesteps)])
# means.list.HYOS$`temps$dts` <- as.Date(means.list.HYOS$`temps$dts`)

means.list.HYOS <- rowSums(out[,3,])/9
means.list.HYOS <- as.data.frame(cbind(means.list.HYOS[150:571], temps$dts[150:571]))
colnames(means.list.HYOS) <- c("mean.abund", "Date")
means.list.HYOS$Date <- as.Date(as.POSIXct(means.list.HYOS$Date, origin = "1970-01-01"))


HYOS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.HYOS$Date), means[149:570])))
HYOS.samp.sum$V1 <- as.Date(HYOS.samp.sum$V1, origin = "1970-01-01")
HYOS.samp.sum <- HYOS.samp.sum[which(HYOS.samp.sum$V1 < "2022-01-01"),]

cor.df <- left_join(HYOS.samp.sum, means.list.HYOS, by=c('V1'="Date"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)

summary(cor.lm)
ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")+
  geom_text(x = 3, y = 300, label = "y = 7.77e-11x, R^2 = 0.22")+
  labs(y = "Hydropsychidae Model Output", x = "Hydropsychidae Empirical Data")

cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")

ggplot(data = means.list.HYOS, aes(x = Date,
                                                        y = mean.abund*10, group = 1, color = "Model")) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
              # alpha = .15,
              # show.legend = T) +
  geom_line(show.legend = T, linewidth = 0.7) +
  geom_line(data =HYOS.samp.sum, aes(x = as.Date(V1, origin = "1970-01-01"), y = V2, color = "Empirical"), show.legend = T)+
  #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  #coord_cartesian(ylim = c(0,2000000)) +
  ylab('Hydropsyichidae S3 Density inds/(m2)') +
  xlab("")+
  labs(colour=" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")

# plot(means.list.HYOS$Date, means.list.HYOS$mean.abund, type = "l")
# lines(as.Date(temps$dts), temps$Temperature * 5, col = "red")
# lines(as.Date(flow.magnitude$dts), flow.magnitude$Discharge * 1000, col = "blue")
# lines(as.Date(HYOS.samp.sum$V1), HYOS.samp.sum$V2*2, col = "darkgreen")

