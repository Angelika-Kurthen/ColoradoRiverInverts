#########################
## Code to Validate NZMS model
###########################
library(purrr)
library(tidyverse)
library(lubridate)
library(readxl)
library(plyr)
library(dplyr)
library(ggplot2)
# data retrieval tool from USGS
library(dataRetrieval)
#install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)

source("1spFunctions.R")
source("NZMS_1sp_Model.R")

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)



#2000, baselineK = 5000, Qmin = 0.25
#disturbanceK = 9000, baselineK = 10000, Qmin = 0.3
out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 9000, baselineK = 5000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))


# adults<-as.data.frame(cbind(as.Date(temps$dts), out[1:length(temps$dts),2:3,1]))
# colnames(adults) <- c("Time","Adult")
# adults$Time <- as.Date(adults$Time, origin = "1970-01-01")

means.list.NZMS <- mean.data.frame(out, burnin = 50, iteration= 9)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[49:length(temps$dts)])
colnames(means.list.NZMS) <- c("timesteps", "mean.abund", "sd.abund", "se.abund", "date")
means.list.NZMS$date <- as.Date(means.list.NZMS$date)
# plot abundance over time

# drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)
# drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]
# 
# NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
# NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
# NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
# NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume
# 
# NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
# NZMS.samp <- aggregate(NZMS.samp$Density, list(NZMS.samp$Date), FUN = mean)

NZMS.samp <- NZMSsamp
NZMS.samp$Density <- (NZMS.samp$Density/(0.028*(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))

means <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date %within% interval(as.Date(temps$dts[i]), as.Date(temps$dts[i+1])) == T),]
  if (is.nan(mean(d$x)) == T) {
    s = NA
  } else {
    s<- mean(d$Density, na.rm = T)}
  means <- append(means, s)
}


NZMS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(temps$dts, format = "%Y-%m-%d"), means)))
NZMS.samp.sum <- NZMS.samp.sum[NZMS.samp.sum$means >= 0, ]
NZMS.samp.sum$V1 <- as.Date(NZMS.samp.sum$V1, origin = "1970-01-01")
 
# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  
ncor <- left_join(lam, means.list.NZMS, by = c("V1" = "date"), copy = T) 
cor.test(ncor$V2, ncor$mean.abund, method = "spearman")




cor.df <- left_join(NZMS.samp.sum, means.list.NZMS, by=c('V1'="date"), copy = T)
cor.lm <- lm((cor.df$mean.abund) ~ (cor.df$means))
cor.test((cor.df$means), (cor.df$mean.abund), method = "spearman")
cor.df <- cor.df[-125,]

NZMS.samp.sum1 <- NZMS.samp.sum %>% slice(which(row_number() %% 3 == 0))
NZMS.samp.sum2 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 1))
NZMS.samp.sum3 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 2))

cor.df1 <- left_join(NZMS.samp.sum1, means.list.NZMS, by=c('V1'="date"), copy = T)
cor.lm1 <- lm(cor.df1$mean.abund ~ cor.df1$V2)
summary(cor.lm1)
plot(cor.df1$means, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$means), (cor.df1$mean.abund), method = "spearman")

cor.df2 <- left_join(NZMS.samp.sum2, means.list.NZMS, by=c('V1'="date"), copy = T)
cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho2 <- cor.test((cor.df2$means), (cor.df2$mean.abund), method = "spearman")

cor.df3 <- left_join(NZMS.samp.sum3, means.list.NZMS, by=c('V1'="date"), copy = T)
cor.lm3 <- lm(cor.df3$mean.abund ~ cor.df3$means)
summary(cor.lm3)
rho3 <- cor.test((cor.df3$means), (cor.df3$mean.abund), method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))


summary(cor.df)
ggplot(data = cor.df, aes(x = (means) , y = (mean.abund)))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")+
  geom_text(x = 3e+05, y = 6100, label = "y = 0.0047x, R^2 = 0.05")+
  labs(y = "NZMS Model Output", x = "NZMS Emprical Data")


hist(cor.df$V2, xlab = "Empirical NZMS Density (#/m3)", col = "#CADBD7")
xlab = ("NZMS Abundance")

colors <- c("black", "#FF7F00")
ggplot(data = means.list.NZMS, aes(x = date, y = scale(mean.abund), group = 1, color = "Model")) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data = NZMS.samp.sum[-125,], aes(x =V1, y = scale(means), color = "Empirical"), linewidth = 1, alpha = 0.8, show.legend = T)+
  geom_point(data = NZMS.samp.sum[125,], aes(x = V1, y = scale(means), color = "Empirical"), show.legend = T)+
  ylab('Log New Zealand Mudsnail Abundance (#/m2)') +
  xlab("")+
  labs(colour=" ")+
  theme_bw()+
  scale_color_manual(values = colors)+
  scale_y_continuous(
    sec.axis = sec_axis(~., name="Log New Zealand Mudsnail Abundance (#/m3)")
  )+ 
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%Y")



# attempting to use zinb data
# est <- vector()
# for (i in 1:length(unique(vals))){
#   est[i] <- mean(std_index[[i]])
# }
# 
# est <- as.data.frame(cbind(temps$dts[unique(vals)], est))
# est$V1 <- as.Date(as.POSIXct(est$V1))
# 
# cor.df <- left_join(est, means.list.NZMS, by=c('V1'="date"), copy = T)
# 
# values <- means.list.NZMS[which(as.Date(temps$dts[unique(vals)]) %in% means.list.NZMS$date),]
# nmix <- est[!is.na(est$V2),]
# 
# # because data starts in 2007, we use first few years as burn in so we need to also 
# cor.test(values$mean.abund, est)
# fo