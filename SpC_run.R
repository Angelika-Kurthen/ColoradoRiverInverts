############################
# Code to run Species C code
############################

# Time <- c(1:365)
# Date <- rep(c(1:365), times = 100)
# Day <- seq(as.Date("2022-01-01"), as.Date("2222-12-31"), by="days")
# Day <- Day[-which(Day == "2024-02-29" & Day = "2028-02-29" & Day = "20")]
# 
# Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243
# 
# temp <- as.data.frame(cbind(Time, Day, Temperature))
# temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
# colnames(temp) <- c("Time", "Date", "Temperature")
# temp <- TimestepTemperature(temp)
# temp <- temp[c(1,3)]
# discharge <- rep(0, times = 131)
# flow.magnitude <- as.data.frame(cbind(temp$dts, discharge))
# 
library(lubridate)


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


source("C_1sp_Model.R")
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

Year <- year(temp$dts)
uYear <- unique(Year)
Month <- month(temp$dts)


## Figure 1: single disturbance time - looking at how disturbance affects population dynamics
discharge <- rep(0.1, time = length(temp$dts)) # vector of discharge
discharge[367] <- 0.3 # create one discharge at 2036-01-15
flow.magnitude <- as.data.frame(cbind(temp$dts, discharge)) #match to dates

source("C_1sp_Model.R")

# run model
out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 99, peaklist = 0, peakeach = length(temp$Temperature))
# create summary dataframe
m <- mean.data.frame(out, burnin = 250, iteration = 99)
# add dates
m <- as.data.frame(cbind(m, temp$dts[249:2608]))
m$`temp$dts[249:2608]`<- as.Date(m$`temp$dts[249:2608]`)
# only want to look a little subset
m <- m[100:165,]
# plot
JanDist <- ggplot(data = m, aes(x = `temp$dts[249:2608]`, y = mean.abund/10000, group = 1))+
  geom_line(linewidth = 1, mapping = aes(color = "January Disturbance +/- S.E."))+
  # geom_ribbon(aes(ymin = (mean.abund - 1.96 * se.abund)/10000,
  #                 ymax = (mean.abund + 1.96 * se.abund)/10000),
  #             colour = 'transparent',
  #             fill = "#F5793A",
  #             alpha = .15, show.legend = T) +
  labs(colour = " ")+
  ylab('Sp C Relative to Baseline Carrying Capacity') +
  xlab("Time")+
  theme_bw()+
  theme(
    legend.title = element_text(margin=margin(b=10)),
    legend.spacing.y = unit(-3,'pt'),
    legend.position = c(0.76, 0.76)
  )+
  labs(colour = "One-Time Disturbance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5),
        axis.text.y = element_text(size = 13))+
  scale_color_manual(values=c("#F5793A"))+
  scale_x_date(date_labels="%B", date_breaks  ="2 month")

# season <- as.numeric(unlist(c(1, 4, 7, 10)))
# 
# for (s in season){
#   results <- lapply(uYear, function(x){
#     sample_dates <- temp$dts[Month == s & Year == x]
#     return(sample(sample_dates, size = 1))
#   })
#   
#   discharge <- rep(0.1, time = length(temp$dts))
#   results <- do.call("c", results)
#   discharge[match(results, temp$dts)] <- 0.3
#   flow.magnitude <- as.data.frame(cbind(temp$dts, discharge))
# 
# system.time( 
#   
# out <- Cmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 99, peaklist = 0, peakeach = length(temp$Temperature))
# 
# )
# m <- mean.data.frame(out, burnin = 250, iteration = 99)
# assign(paste0("means.list.", s), m)
# }
# 
# 
# 
# means.list <- as.data.frame(cbind(means.list.1$mean.abund, means.list.4$mean.abund, means.list.7$mean.abund, means.list.10$mean.abund, temp$dts[250:length(means.list.1$timesteps)]))
# means.list[,5] <- as.Date(as.POSIXct(means.list[,5], origin = "1970-01-01"))
# colnames(means.list) <- c("January", "April", "July", "October", "Time")
# x11()
# means.list.abridged <- means.list[106:158,]
# ggplot(data = means.list.abridged, aes(x = Time,
#                               y = January/10000, group = 1, color = "January")) +
#   #geom_ribbon(aes(ymin = mean.abund - 1.96 * se.abund,
#   #                ymax = mean.abund + 1.96 * se.abund),
#   #            colour = 'transparent',
#   #            alpha = .15,
#   #            show.legend = T) +
#   geom_line(show.legend = T) +
#   geom_line(data = means.list.abridged, aes(x = Time, y = April/10000, color = "April"), show.legend = T)+
#   geom_line(data = means.list.abridged, aes(x = Time, y = July/10000, color = "July"), show.legend = T)+
#   geom_line(data = means.list.abridged, aes(x = Time, y = October/10000, color = "October"), show.legend = T)+
#   coord_cartesian(ylim = c(0,2.5)) +
#   ylab('Caddisfly Abundance Relative to Baseline Recuritment Limit') +
#   xlab("")+
#   labs(colour=" ")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), )+
#   scale_x_date(date_labels="%B", date_breaks  ="4 months")
# 
# 
# means.list.C <- mean.data.frame(out,burnin = 50, iteration= 1)
# means.list.C <- cbind(means.list.C, temp$dts[1:length(means.list.C$timesteps)])
# means.list.C$`temp$dts` <- as.Date(means.list.C$`temp$dts`)
# 
# plot(means.list.C$`temp$dts`, means.list.C$mean.abund, type = "l", xlab = "Time", ylab = "Type C species")

## Figure 2: Calculate Julian Date Effect - in each iteration, we want to create a disturbance at a specific date
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
  
  source("C_1sp_Model.R")
  # run model
  out <- Cmodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature))
  # create summary dataframe 
  m <- mean.data.frame(out, burnin = 250, iteration = 9)
  m <- cbind(m, discharge[250:402])
  means[d] <- list(m)
  abund <- means[[d]][(which(discharge == "0.3") - 248):(which(discharge == "0.3") - 238),2 ]
  post <- seq(0, 10, by = 1)
  dat <- rep(d, times = length(post))
  jday_data <- as.data.frame(rbind(jday_data, cbind(dat, post, abund)))
}

jday.df <- as.data.frame(cbind(jday_max, all.dates))
jday.df$jday_max <- as.numeric(jday.df$jday_max)
jday.df$all.dates <- yday(as.POSIXct(jday.df$all.dates, format = "%m-%d"))
jday.plot <- ggplot(data = jday.df, aes(x = all.dates, y = jday_max/10000, group = 1))+
  geom_line(linewidth = 1)+
  coord_cartesian(ylim = c(0,5)) +
  theme_bw()+
  ylab('Max Values of Sp C Abundance Relative to Baseline K') +
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
          main = 'C spp. Abund')




# jday_max <- unlist(lapply(means, function(x)
#   return(max(x$mean.abund)))) # would mean +/- se or just maximum (and maybe minimum values) be valuable
# jday.df <- as.data.frame(cbind(jday_max, all.dates))
# jday.df$jday_max <- as.numeric(jday.df$jday_max)
# jday.df$all.dates <- yday(as.POSIXct(jday.df$all.dates, format = "%Y-%m-%d"))
# ggplot(data = jday.df, aes(x = all.dates, y = jday_max/10000, group = 1))+
#   geom_line(linewidth = 1)+
#   coord_cartesian(ylim = c(0,2.75)) +
#   theme_bw()+
#   ylab('Max Values of Caddisfly Abundance Relative to Baseline Recuritment Limit') +
#   xlab("Julian Date")+
#   theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))+
#   scale_x_continuous(n.breaks = 11)
# 
