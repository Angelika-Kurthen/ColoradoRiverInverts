#########################
## Code to Validate NZMS model
###########################
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


source("NZMS_1sp_Model.R")

discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
flow.magnitude <- TimestepDischarge(discharge, 85000)
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
temps <- TimestepTemperature(temp)

set.seed(333)
out <- NZMSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 9000, baselineK = 5000, Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))


# adults<-as.data.frame(cbind(as.Date(temps$dts), out[1:length(temps$dts),2:3,1]))
# colnames(adults) <- c("Time","Adult")
# adults$Time <- as.Date(adults$Time, origin = "1970-01-01")

means.list.NZMS <- mean.data.frame(out,burnin = 50, iteration= 1000)
means.list.NZMS <- cbind(means.list.NZMS, temps$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temps$dts` <- as.Date(means.list.NZMS$`temps$dts`)
# plot abundance over time

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = F)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]

NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume

NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
NZMS.samp$Density <- (NZMS.samp$Density/(9e-15 *(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))
NZMS.samp <- aggregate(NZMS.samp$Density, list(NZMS.samp$Date), FUN = mean)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Group.1 %within% interval(temps$dts[i], temps$dts[i+1]) == T),]
  if (is.nan(mean(d$x)) == T) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
}

# test for temporal autocorrelation
#acf(na.omit(means)) #we have some

NZMS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.NZMS$`temps$dts`, format = "%Y-%m-%d"), means[51:406])))
NZMS.samp.sum$V1 <- as.Date(NZMS.samp.sum$V1, origin = "1970-01-01")

# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  

cor.df <- left_join(NZMS.samp.sum, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm <- lm((cor.df$mean.abund) ~ (cor.df$V2))
cor.test((cor.df$V2+1), (cor.df$mean.abund+1), method = "spearman")


# culling data by lags - Temporal autocorrelation: a neglected factor in the study of behavioral repeatability and plasticity  
NZMS.samp.sum1 <- NZMS.samp.sum %>% slice(which(row_number() %% 3 == 0))
NZMS.samp.sum2 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 1))
NZMS.samp.sum3 <- NZMS.samp.sum %>%  slice(which(row_number() %% 3 == 2))

cor.df1 <- left_join(NZMS.samp.sum1, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm1 <- lm(cor.df1$mean.abund ~ cor.df1$V2)
summary(cor.lm1)
plot(cor.df1$V2, cor.df1$mean.abund)
rho1 <- cor.test((cor.df1$V2+1), (cor.df1$mean.abund+1), method = "spearman")

cor.df2 <- left_join(NZMS.samp.sum2, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm2 <- lm(cor.df2$mean.abund ~ cor.df2$V2)
rho2 <- cor.test((cor.df2$V2+1), (cor.df2$mean.abund+1), method = "spearman")

cor.df3 <- left_join(NZMS.samp.sum3, means.list.NZMS, by=c('V1'="temps$dts"), copy = T)
cor.lm3 <- lm(cor.df3$mean.abund ~ cor.df3$V2)
summary(cor.lm3)
rho3 <- cor.test((cor.df3$V2+1), (cor.df3$mean.abund+1), method = "spearman")

rho <- mean(c(rho1$estimate, rho2$estimate, rho3$estimate))
sd <- sd(c(rho1$estimate, rho2$estimate, rho3$estimate))

summary(cor.df)
ggplot(data = cor.df, aes(x = (V2) , y = (mean.abund)))+
  geom_point()+
  stat_smooth(method = "lm",
              formula = y ~ x,
              geom = "smooth")+
  geom_text(x = 3e+05, y = 6100, label = "y = 0.0047x, R^2 = 0.05")+
  labs(y = "NZMS Model Output", x = "NZMS Emprical Data")

colors <- c("#4477AA", "black")
linetypes <- c("solid", "twodash")
NZMSts <- ggplot(data = cor.df, aes(x = V1, y = scale(mean.abund), group = 1, color = "Model", linetype = "Model")) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data = cor.df, aes(x =V1, y = scale(V2), color = "Empirical", linetype = "Empirical"), linewidth = 1, alpha = 0.8, show.legend = T)+
  #geom_point(data = NZMS.samp.sum[125,], aes(x = V1, y = scale(means), color = "Empirical"), show.legend = T)+
  labs(y=expression(paste(italic("P. antipodarum"), " Abund.")))+
  ylim(c(-4,7))+
  geom_text(mapping = aes(x = as.Date("2018-06-01"), y =5, label = paste('rho', "==", 0.59)), parse = T, color = "black", size = 4.5)+
  xlab("")+
  labs(colour=" ")+
  theme_bw()+
  scale_color_manual(values = colors)+
  scale_linetype_manual(values = linetypes)+
  guides(linetype=guide_legend(" "), color = "none")+
  theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%Y")

