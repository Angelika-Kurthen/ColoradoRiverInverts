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
temp <- subset(temp, Date >= "2004-05-23")
temp$Date <- as.Date(temp$Date, format = "%Y-%m-%d")
temp <- aggregate(temp$Temperature, by = list(temp$Date), FUN = mean)
colnames(temp) <- c("Date", "Temperature")
simtemp <- temp
# all dates
all_dates <- as.data.frame(seq.Date(from = as.Date("2004-05-23"), to = as.Date("2024-08-23"), by = "days"))
names(all_dates) <- "Date"
# join them
temp <- full_join(temp, all_dates)


simtemp$Date <- yday(simtemp$Date)
simtemp <- simtemp %>% group_by(Date) %>% dplyr::summarise(Temperature = mean(Temperature))


temp$Temperature[which(is.na(temp$Temperature)==T)] <- simtemp$Temperature[yday(temp$Date[which(is.na(temp$Temperature)==T)])]
temp <-temp[order(temp$Date),]
temps <- TimestepTemperature(temp) # intermittantly missing data until  2000-12-26 
# load discharge data from above Diamond Creek Confluence (RM226)
# now for the issues of NAs, need to add in averages 

# that is our filled in temperature data

discharge <- readNWISdv("09404200", "00060", "2004-05-22", "2024-08-23")
discharge <- full_join(discharge, all_dates)
flow.magnitude <- TimestepDischarge(discharge, 85000)

out <- HYOSmodel(flow.data = flow.magnitude$Discharge, temp.data = temps, disturbanceK = 20000, baselineK = 10000, Qmin = 0.2, extinct = 5, iteration = 9, peaklist = 0.17, peakeach = length(temps$Temperature))

drift.data.total <- readDB(gear = "LightTrap", type = "Sample", updater = F)

drift.CR <- drift.data.total[which(drift.data.total$Region == "GrandCanyon" & drift.data.total$RiverMile >= 219 & drift.data.total$RiverMile <= 225),]

specieslist <- c("CHEP","HYOS", "HYSP")
HYOS.samp.CR <- sampspec(samp = drift.CR, stats = T, species = specieslist)
HYOS.samp <- merge(HYOS.samp.CR$Statistics, HYOS.samp.CR$Samples, by = "BarcodeID", all = T)
HYOS.samp <- HYOS.samp[which(HYOS.samp$FlagStrange == F),]
HYOS.samp$Density <- HYOS.samp$CountTotal/HYOS.samp$TimeElapsed
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


hist(cor.df$V2, xlab = "Hydropsychidae Adults (#/hour)", col = "#CADBD7")
colors <- c("black", "#FF7F00")
ggplot(data = means.list.HYOS[100:571,], aes(x = Date, y = mean.abund, group = 1, color = "Model")) +
  # geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
  #                 ymax = mean.abund + 1.96 * se.abund),
  #             colour = 'transparent',
              # alpha = .15,
              # show.legend = T) +
  geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
  geom_line(data =HYOS.samp.sum, aes(x = as.Date(V1, origin = "1970-01-01"), y = V2/10, color = "USGS Light Trap Data"), linewidth = 1, alpha = 0.8,show.legend = T)+
  # geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = Discharge), color = "blue") +
  # geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*0.1), color = "green")+
  # #coord_cartesian(ylim = c(0,2000000)) +
  ylab('Hydropsychidae Adults (#/m2)') +
  scale_y_continuous(
    # Features of the first axis
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./10, name="Hydropsychidae Adults (#/hour)")
  ) + 
  xlab("Year")+
  labs(colour=" ")+
  theme_bw()+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), )+
  scale_x_date(date_labels="%Y")+
  scale_color_manual(values = colors)

# plot(means.list.HYOS$Date, means.list.HYOS$mean.abund, type = "l")
# lines(as.Date(temps$dts), temps$Temperature * 5, col = "red")
# lines(as.Date(flow.magnitude$dts), flow.magnitude$Discharge * 1000, col = "blue")
# lines(as.Date(HYOS.samp.sum$V1), HYOS.samp.sum$V2*2, col = "darkgreen")

drift.data.total <- readDB(gear = "LightTrap", type = "Sample", updater = F)

drift.CR <- drift.data.total[which(drift.data.total$Region == "GrandCanyon" & drift.data.total$RiverMile >= 219 & drift.data.total$RiverMile <= 225),]

specieslist <- c("CHEP","HYOS", "HYSP")
HYOS.samp.CR <- sampspec(samp = drift.CR, stats = T, species = specieslist)
HYOS.samp <- merge(HYOS.samp.CR$Statistics, HYOS.samp.CR$Samples, by = "BarcodeID", all = T)
HYOS.samp <- HYOS.samp[which(HYOS.samp$FlagStrange == F),]

vals <- vector()
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# add to data frame
HYOS.samp <- cbind(HYOS.samp, vals)

# now we need into include mean water temperature for each 
HYOS.samp <- cbind(HYOS.samp, temps$Temperature[HYOS.samp$vals])
HYOS.samp <- left_join(HYOS.samp, discharge, by = c("Date", "Date"))
max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- HYOS.samp[which(HYOS.samp$Date >= temps$dts[i] & HYOS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means - this is just for checking NAs
}
# make data frame with the timestep, the mean
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]
# phenology may also play into dynamics, so include month column as well
#month <- month(HYOS.samp$Date)
install.packages("aspace")
library(aspace)
df$circdate <- sin(as_radians((lubridate::yday(df$`temps$dts`)/365)*360))

# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)
obs_intercept <- matrix(data = 1, nrow = R, ncol = J)
# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
flows <- vector()
temperature <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
time <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- HYOS.samp[which(HYOS.samp$Date >= temps$dts[i] & HYOS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$X_00060_00003)
  #date <- df[]
  temperature[i] <- mean(d$`temps$Temperature[HYOS.samp$vals]`)
  time[i,] <- c(d$TimeElapsed, rep(NA, times = (J- length(d$CountTotal))))
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  dens_mat[i, ] <- c(as.integer(d$Density), rep (NA, times = (J - length(d$Density))))
  volumes[i, ] <- c((d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,]) # count data
#dens_mat <- as.matrix(dens_mat[-nodata, ]) # density data
obs_intercept <- as.matrix(obs_intercept[-nodata,]) # intercept for obs cov
time <- as.matrix(scale(time[-nodata,])) # duration in H20 obs cov
time[is.na(time)] <- mean(time, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors

flows <- as.data.frame(scale(flows[-nodata])) # site cov flow 
temperature <- as.data.frame(scale(temperature[-nodata])) # site cov temp
circdate <- as.data.frame(df$circdate)
# dimnames(time) <- list(temps$dts[-nodata], seq(1:48))
# time <- list(time)
# names(time) <- c("time")
site_intercept <- rep(1, times = length(flows$V1)) 
site_covs<- as.matrix(cbind(site_intercept, flows, circdate)) #flows,temperature, circdate)
#obs_covs <- array(data= NA, dim = c(length(flows$V1),J,2))
obs_covs <- array(data= NA, dim = c(length(flows$V1),J,1))
obs_covs[,,1] <- obs_intercept                                  
#obs_covs[,,2] <- time
#offset
offset <- as.matrix(scale(log(volumes[-nodata, ])))
offset[is.na(offset)] <- mean(offset, na.rm = TRUE) # replace NAs with mean duration time since NAs not allowed in predictors or offsets
# dimnames(volumes) <- list(temps$dts[-nodata], seq(1:48))
# volumes <- list(volumes)
# names(volumes) <- c("vol")


