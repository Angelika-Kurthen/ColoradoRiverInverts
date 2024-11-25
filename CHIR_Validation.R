##################
# CHIR N Mixture
#################
# data retrieval tool from USGS
#install.packages("dataRetrieval")
library(dataRetrieval)
#install.packages("devtools")
#install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(devtools)
library(foodbase)
library(lubridate)
library(AICcmodavg)
library(unmarked)
library(ubms)

source("1spFunctions.R")
source("CHIR_1sp_Model.R")
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
# pull drift data from DB
drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)
# get drift data from between Lees Ferry and RM -6
drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
#specify
CHIR.samp.LF <- sampspec(samp = drift.LF, species = "CHIL", stats = T)
# pull stats and merge with sample info
CHIR.samp <- merge(CHIR.samp.LF$Statistics, CHIR.samp.LF$Samples, by = "BarcodeID", all = T)
# make sure we are using the same gear
CHIR.samp <- CHIR.samp[which(CHIR.samp$GearID == 4),] 
CHIR.samp <- CHIR.samp[which(CHIR.samp$FlagStrange == 0), ]
# calculate density
CHIR.samp$Density <- CHIR.samp$CountTotal/CHIR.samp$Volume
CHIR.samp <- merge(CHIR.samp, discharge[, 3:4], by = "Date")
# for densities, apply equation to adjust for flow (but discharge is in cubic meter seconds not cubic feet seconds)
# C(#/m^3) ~ aB(m^2)Q(m^3/s)^g 
# B(#/m^2) ~ C(#/m^3)/aQ(m^3/s)^g
# right now we have C(#/m^3)aQ(cfs)^g
# must convert density from #/cfs to #/m^3/s and Q from cfs to m^s/s
# 1 cfs = 0.0283168 m^/s
CHIR.samp$Density <- (CHIR.samp$Density)/(1.3 *(CHIR.samp$X_00060_00003 * 0.0283168)^4.1)
CHIR.samp <- aggregate(CHIR.samp$Density, list(CHIR.samp$Date), FUN = mean)

source("CHIR_1sp_Model.R")
out <- CHIRmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 100000, baselineK = 10000, Qmin = 0.2, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))


means <- vector()
for (i in 1:length(temps$dts)){
  d <- CHIR.samp[which(CHIR.samp$Group.1 >= temps$dts[i] & CHIR.samp$Group.1 < temps$dts[i+1]),]
  if (any(is.nan(mean(d$x))) == T || any(is.na((d$x) == T))) {
    s = NA
  } else {
    s <- mean(d$x)}
  means <- append(means, s)
  # we know the last value doesn't fit into interval and needs to be manually added
}
means <- as.data.frame(cbind(means, as.Date(temps$dts)))
means$V2 <- as.Date(means$V2, origin = "1970-01-01")

means.list.CHIR <- rowSums(out[, 1:2, ])/9
means.list.CHIR <- as.data.frame(cbind(means.list.CHIR[1:404], temps$dts))
colnames(means.list.CHIR) <- c("mean.abund", "Date")
means.list.CHIR$Date <- as.Date(as.POSIXct(means.list.CHIR$Date, origin = "1970-01-01"))
means.list.CHIR <- as.data.frame(cbind(means.list.CHIR, means))
means.list.CHIR <- means.list.CHIR[-(1:199), ] # use first 150 as burn in 
means.list.CHIR <- means.list.CHIR[which(means.list.CHIR$Date < "2019-02-06"),] # we have some one off data points
cor.df <- na.omit(means.list.CHIR)
cor.lm <- lm(cor.df$mean.abund ~ cor.df$means)

# CHIR.tems <- left_join(CHIR.samp.sum, temps, by = c("V1" = "dts"), copy = T )
# cor.test(CHIR.tems$Temperature, CHIR.tems$means, method = "spearman")
# # 
# ggplot(data = CHIR.tems, aes(x = V1, y = means*10000000000000))+
#   geom_line()+
#   geom_line(data= CHIR.tems, aes(x = V1, y = Temperature, col = "red"))

# plot
# flow.magnitude$dts <- as.Date(flow.magnitude$dts)
# CHIR.flows <- left_join(CHIR.samp.sum, flow.magnitude, by = c("V1" = "dts"), copy = T)
# cor.test(CHIR.flows$means, CHIR.flows$Discharge)


cor.test((cor.df$means), (cor.df$mean.abund), method = "spearman")

ggplot(data = cor.df, aes(x = Date,
                                   y = mean.abund/10000, group = 1, color = "Model")) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
  # alpha = .15,
  # show.legend = T) +
  geom_line(show.legend = T, linewidth = 0.7) +
  geom_line(data =cor.df, aes(x = Date, y = means*10000000000, color = "Empirical"), show.legend = T)+
  #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
  #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
  #coord_cartesian(ylim = c(0,2000000)) +
  ylab('Hydropsyichidae S3 Density inds/(m2)') +
  xlab("")+
  labs(colour=" ")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%B", date_breaks  ="6 months")



