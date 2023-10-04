##############################
## Code to validate BAET model
##############################
library(purrr)
library(tidyverse)
library(lubridate)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)


# load Baet model
source("BAET_1sp_Model.R")
# pull discharge and temps from below flaming gorge dam
discharge <- readNWISdv("09234500", "00060", "1986-10-01", "1999-10-06")
# Bankfull discharge for Green River from https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=7604&context=etd
flow.magnitude <- TimestepDischarge(discharge, 22424.813)
temp <- readNWISdv("09234500", "00010", "2004-02-05", "2023-05-01")
temps <- average.yearly.temp(temp, "X_00010_00003","Date")
temps <- rep.avg.year(temps, 15, change.in.temp = 0, years.at.temp = 15)
# align dates
temps <- temps[20:359,2:3]
temps$dts <- flow.magnitude$dts

out <- BAETmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 40000, baselineK = 7000, Qmin = 0.12, extinct = 50, iteration = 9, peaklist = 0.13, peakeach = length(temps$Temperature))

# upload larval baet data from Flaming Gorge Dam 
bugdata <- read_delim("~/ColoradoRiverInverts/APPL_11_13/bugdata.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
bugdata <- as.data.frame(bugdata[-c(1:6, 3732:3740),])
names(bugdata) <- c("Sample", "Location", "Date", "Citation", "Method", "Area", "Density", "Phylum", "Class", "Order", "Family", "Subfamily", "Genus", "Species")
bugdata <- bugdata[which(bugdata$Location == "0.8KDD" | bugdata$Location == "6KDD" | bugdata$Location == "12KDD"),]
bugdata$Date <- as.Date(bugdata$Date, "%m/%d/%Y")
bugdata$Density <- as.numeric(bugdata$Density)
BAETdata <- bugdata[which(bugdata$Date >= "1986-10-01" & bugdata$Family == "Baetidae"), ]

BAET.samp <- aggregate(BAETdata$Density, list(BAETdata$Date), FUN = sum)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- BAET.samp[which(BAET.samp$Group.1 >= temps$dts[i] & BAET.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na(d$x) == T)) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means[length(temps$dts)] <- BAET.samp$x[length(BAET.samp$x)]

means.list.BAET <- mean.data.frame(out, burnin = 188, iteration= 9)
means.list.BAET <- cbind(means.list.BAET, temps$dts[188:341])
means.list.BAET$`temps$dts` <- as.Date(means.list.BAET$`temps$dts`)

BAET.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.BAET$`temps$dts`), means[188:340])))
BAET.samp.sum$V1 <- as.Date(BAET.samp.sum$V1, origin = "1970-01-01")

#BAET.samp.sum <- BAET.samp.sum[which(BAET.samp.sum$V1 %in% sample(BAET.samp.sum$V1, size = 30)),]

cor.df <- left_join(BAET.samp.sum, means.list.BAET, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$V2)
cor.test((cor.df$V2), (cor.df$mean.abund), method = "spearman")


summary(cor.lm)
ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")+
  geom_text(x = 1000, y = 3250, label = "y = 0.00021x, R^2 = 0.37")+
  labs(y = "Baetidae Model Output", x = "Baetidae Emprical Data")

abund.trends.NZMS <- ggplot(data = means.list.BAET, aes(x = `temps$dts`,
                                                        y = mean.abund, group = 1, color = "Model")) +
  geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
                  ymax = mean.abund + 1.96 * se.abund),
              colour = 'transparent',
              alpha = .15,
              show.legend = T) +
  geom_line(show.legend = T, linewidth = 0.7) +
  geom_line(data =BAET.samp.sum, aes(x = as.Date(V1, origin = "1970-01-01"), y = V2, color = "Empirical"), show.legend = T)+
  #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  #coord_cartesian(ylim = c(0,6000)) +
  ylab('Baetidae S1 and S2 Density inds/(m2)') +
  xlab("")+
  labs(colour=" ")+
theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
      axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")


