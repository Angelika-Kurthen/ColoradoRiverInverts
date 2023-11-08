
#############################
# Code to run Species A Model
#############################



Time <- c(1:365)
library(lattice)
library(latticeExtra)
library(RColorBrewer)

## Plot setup
clrs <- colorRampPalette(brewer.pal(9, "YlOrRd"))
trellis.par.set("axis.line", list(col = NA, lty = 1, lwd = 1))
theme.novpadding <- list(
  layout.heights = list(top.padding = 0, bottom.padding = 0),
  layout.widths = list(left.padding = 0, right.padding = 0, zlab.axis.padding = 2)
)

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

# pull one random day in either January, April, July, or Octorber for each year
Year <- year(temp$dts) # make a list of years
uYear <- unique(Year) # unique years
Month <- month(temp$dts) # make a list of months
season <- as.numeric(unlist(c(1, 4, 7, 10))) # list of months we want to investigate 
Week <- week (temp$dts) # make a list of the week of the month
source("A_1sp_Model.R")



## Figure 3: Calculate Julian Date Effect - in each iteration, we want to create a disturbance at a specific date
# a few ways to do this- for every date, run out 300 timesteps prior (burn-in + some change) --> (Date - 300:Date+100)
# pull one random day in either January, April, July, or Octorber for each year
all.dates <- unique(format(temp$dts, "%m-%d"))[order(unique(format(temp$dts, "%m-%d")))]
# loop to select a date from a Week-Month combo from each unique year
means <- list()
jday_data <- data.frame(date = NULL, post_dist = NULL, abund = NULL)
jday_max <- vector()
for (d in 1:length(all.dates)){ # 30 reps takes 60 mins
  sample_dates <- temp$dts[which(format(temp$dts, "%m-%d")== all.dates[d])]
  samp <- which(temp$dts == sample(sample_dates[which(sample_dates > temp$dts[300] & sample_dates < temp$dts[2508])], size = 1))
  dates <- temp[(samp-300):(samp+100),]
  discharge <- rep(0.1, time = length(dates$dts)) # create a list of non-disturbance discharges
  discharge[match(temp$dts[samp], dates$dts)] <- 0.3 # from that list of dates from above, assign a disturbance discharge to that date
  
  source("A_1sp_Model.R")
  # run model
  out <- Amodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature))
  # create summary dataframe 
  m <- mean.data.frame(out, burnin = 250, iteration = 9)
  m <- cbind(m, discharge[250:402])
  means[d] <- list(m)
  abund <- means[[d]][(which(discharge == "0.3") - 249):(which(discharge == "0.3") - 240),2 ]
  post <- seq(0, 10, by = 1)
  dat <- rep(d, times = length(post))
  jday_data <- as.data.frame(rbind(jday_data, cbind(dat, post, abund)))
}
# jday_max <- unlist(lapply(means, function(x)
#         return(max(x$mean.abund[)])))) # would mean +/- se or just maximum (and maybe minimum values) be valuable
# }

#jday_max[d] <- max(means[[d]][(which(discharge == "0.3") - 250):(which(discharge == "0.3") - 239),2 ])


jday.df <- as.data.frame(cbind(jday_max, all.dates))
jday.df$jday_max <- as.numeric(jday.df$jday_max)
jday.df$all.dates <- yday(as.POSIXct(jday.df$all.dates, format = "%m-%d"))
jday.plot <- ggplot(data = jday.df, aes(x = all.dates, y = jday_max/10000, group = 1))+
  geom_line(linewidth = 1)+
  coord_cartesian(ylim = c(0,5)) +
  theme_bw()+
  ylab('Max Values of Mayfly Abundance Relative to Baseline K') +
  xlab("Julian Date")
theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
      axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
  scale_x_continuous(n.breaks = 11)

wireframe(abund ~ dat + post, data = jday_data, # 
          #distance =c(2, 5, 8),
          aspect = c(1, .4),
          drape = TRUE,
          shade = F,
          colorkey = FALSE,
          col = alpha('#ffeda0', 0.08),
          xlab = "Julian Date",
          ylab = list(label = "Time Post-Disturbance", rot = 33),
          zlab = list(label = "Abundance", rot = 90),
          zlim = c(6600, 40000),
          scales = list(arrows = FALSE),
          screen = list(z = -40, x = -70),
          #panel.3d.wireframe = panel.3d.contour, # can't use because Rtools is deprecated
          par.settings = theme.novpadding,
          col.regions = clrs(1000),
          main = 'A spp. Abund')



out <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0.13, peakeach = 
         length(temps$Temperature))
means.list.NZMS <- mean.data.frame(out,burnin = 50, iteration= 1)
means.list.NZMS <- cbind(means.list.NZMS, temp$dts[1:length(means.list.NZMS$timesteps)])
means.list.NZMS$`temp$dts` <- as.Date(means.list.NZMS$`temp$dts`)

plot(means.list.NZMS$`temp$dts`, means.list.NZMS$mean.abund, type = "l", xlab = "Time", ylab = "Type A species")
