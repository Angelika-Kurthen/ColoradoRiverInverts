#########################
# Code to run Species B Model
#######################
library(lubridate, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(ggplot2, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")
library(data.table, lib.loc = "/home/ib/kurthena/R_libs/4.2.1")

library(lubridate)
library(ggplot2)
library(data.table)


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

## Figure 1: single disturbance time - looking at how disturbance affects population dynamics
discharge <- rep(0.1, time = length(temp$dts)) # vector of discharge
discharge[367] <- 0.3 # create one discharge at 2036-01-15
flow.magnitude <- as.data.frame(cbind(temp$dts, discharge)) #match to dates

source("B_1sp_Model.R")

# run model
out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 99, peaklist = 0, peakeach = length(temp$Temperature))
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
  ylab('Mayfly Abundance Relative to Baseline Carrying Capacity') +
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

#pdf("JanDist")
#print("JanDist")
#dev.off()

## Figure 2: Yearly Averages --> this plot is ugly and not very informative
# pull one random day in either January, April, July, or Octorber for each year
Year <- year(temp$dts) # make a list of years
uYear <- unique(Year) # unique years
Month <- month(temp$dts) # make a list of months
season <- as.numeric(unlist(c(1, 4, 7, 10))) # list of months we want to investigate 
Week <- week (temp$dts) # make a list of the week of the month
# loop to select a date from a Week-Month combo from each unique year
#for (s in season){
#results <- lapply(uYear[2:101], function(x){
#  sample_dates <- temp$dts[Month == 1 & Year == x & Week == 2 | Week == 3]
#  return(sample(sample_dates, size = 1))
#})


#discharge <- rep(0.1, time = length(temp$dts)) # create a list of non-disturbance discharges
#results <- do.call("c", results)
#discharge[match(results, temp$dts)] <- 0.3 # from that list of dates from above, assign a disturbance discharge to that date
#flow.magnitude <- as.data.frame(cbind(temp$dts, discharge)) #match to dates

#source("B_1sp_Model.R")

# run model
#out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature))
# summarise iterations for each model run 
#m <- mean.data.frame(out, burnin = 250, iteration = 9)
# m <- as.data.frame(cbind(m, temp$dts[249:length(temp$dts)]))
#colnames(m) <- c("timesteps", "mean.abund", "sd.abund", "se.abund", "Time")
# reframe to get yearly averages
#m <- m %>% group_by(format(Time, "%Y")) %>% reframe(year.mean.abund = mean(mean.abund), 
#                                               count = n(),
#                                               sd = sd(mean.abund, na.rm = T),
#                                               se = sd/count)
#assign(paste0("means.list.", s), m)
#}

#means.list <- as.data.frame(cbind(means.list.1[, c(1:2,5)], means.list.4[, c(2,5)], means.list.7[,c(2,5)], means.list.10[,c(2,5)]))
#means.list[,1] <- as.Date(means.list[,1], "%Y")
#colnames(means.list) <- c("Year", "January", "Jase", "April", "Ase", "July", "Juse", "October", "Ose")
#x11()

#
#tem <- temp %>% group_by(format(dts, "%Y")) %>% reframe(mean.temp = mean(Temperature))
#tem <- as.data.frame(tem)
#tem$`format(dts, "%Y")` <- as.Date(tem$`format(dts, "%Y")`, "%Y")
#jplot <- ggplot(data = means.list, aes(x = Year,
#                                   y = January/10000, group = 1, col = "January")) +
#  geom_line(linewidth = 1, show.legend = T)+
#  geom_ribbon(aes(ymin = January/10000 - 1.96 * Jase/10000,
#                 ymax = January/10000 + 1.96 * Jase/10000),
#             colour = 'transparent',
#             fill = "#F5793A",
#             alpha = .15,
#             show.legend = T) +
#  geom_line(data = means.list, aes(x = Year, y = April/10000, col = "April"), show.legend = T, linewidth = 1)+
#  geom_ribbon(aes(ymin = April/10000 - 1.96 * Ase/10000,
#                 ymax = April/10000 + 1.96 * Ase/10000),
#             colour = 'transparent',
#             fill = "#A95AA1",
#             alpha = .15,
#             show.legend = T) +
#  geom_line(data = means.list, aes(x = Year, y = July/10000, col = "July"), show.legend = T, linewidth = 1)+
#    geom_ribbon(aes(ymin = July/10000 - 1.96 * Juse/10000,
#                   ymax = July/10000 + 1.96 * Juse/10000),
#               colour = 'transparent',
#               fill = "#85C0F9",
#               alpha = .15,
 #              show.legend = T) +
#  geom_line(data = means.list, aes(x = Year, y = October/10000, col = "October"), show.legend = T, linewidth = 1)+
#    geom_ribbon(aes(ymin = October/10000 - 1.96 * Ose/10000,
#                   ymax = October/10000 + 1.96 * Ose/10000),
#               colour = 'transparent',
#                fill = "#0F2080",
#                alpha = .15,
#                show.legend = T) +
#  coord_cartesian(ylim = c(0,1.25)) +
#  geom_line(data = tem, aes(x = `format(dts, "%Y")`-249, y = mean.temp/9, col = "Temperature"), linewidth = 1) +
#  ylab('Mayfly Abundance Relative to Baseline Recuritment Limit') +
#  xlab("")+
#  labs(colour=" ")+
#  theme_bw()+
#  theme(text = element_text(size = 14), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent") )+
#  scale_color_manual(values=c("#F5793A", "#A95AA1", "#85C0F9", "#0F2080", "#DDCC77"))+
#scale_x_date(date_labels="%Y", date_breaks  ="5 years") 

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

    source("B_1sp_Model.R")
    # run model
    out <- Bmodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature))
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

jday_event <- jday_data[c(TRUE,rep(FALSE,10)), ]

plot(data = jday_event, aes(x = dat, y = abund/1000))+
  geom_line(linewidth = 1)+
  theme_bw()+
  ylab("Spp B Immediate Disturbance Abundance Relative to Rb")+
  xlab("Julian Date")



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
          zlim = c(6600, 60000),
          scales = list(arrows = FALSE),
          screen = list(z = -40, x = -70),
          #panel.3d.wireframe = panel.3d.contour, # can't use because Rtools is deprecated
          par.settings = theme.novpadding,
          col.regions = clrs(1000),
          main = 'B spp. Abund')





pdf("B_Julian")
print(jday.plot)
dev.off()

plot(jday.df, jday.df$jday_max)
##Figure 4: Importance of Disturbance Timing
# Importance of Disturbance Timing is just comparing delta(max) to delta(max)
# could also be a value (min(max)/max(max)) <- the larger the number, the bigger the difference between max values
#max(jday.df$jday_max)-min(jday.df$jday_max)
#vs
#min(jday.df$jday_max)/max(jday.df$jday_max)
#source("B_1sp_Model.R")
#source("1spFunctions.R")
#Time <- c(1:36500)
#Date <- rep(1:365, times = 100)
#Day <- seq(as.Date("2022-01-01"), as.Date("2121-12-31"), by="days")
# find and remove leap days
#find_leap = function(x){
#  day(x) == 29 & month(x) == 2 
#}
#Day <- Day[which(find_leap(Day) == F)]


#Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 10.956243

#temp <- as.data.frame(cbind(Time, Day, Temperature))
#temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
#colnames(temp) <- c("Time", "Date", "Temperature")
#temp <- TimestepTemperature(temp)
#temp <- temp[c(1,3)]

# we can now iterate (similar to a sensitivity analysis) and see how importance of disturbance timing impacts 
# what if we iterated fecundities 
#fecs <- seq(675,1125, by = 9)
# # Initializes the progress bar
# pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
#                      max = length(fecs),  # Maximum value of the progress bar
#                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                      width = 50,   # Progress bar width. Defaults to getOption("width")
#                      char = "=")   # Character used to create the bar
#all.dates <- unique(format(temp$dts, "%m-%d"))[order(unique(format(temp$dts, "%m-%d")))]
#IDT <- vector()
#IDT.ratio <- vector()

#for (f in 1:length(fecs)) {
  # Sets the progress bar to the current state
  # setTxtProgressBar(pb, f)
  # loop to select a date from a Week-Month combo from each unique year
#  means <- list()
#    for (d in 1:length(all.dates)){ # 30 reps takes 60 mins
#      sample_dates <- temp$dts[which(format(temp$dts, "%m-%d")== all.dates[d])]
#      samp <- which(temp$dts == sample(sample_dates[which(sample_dates > temps$dts[300] & sample_dates < temps$dts[2508])], size = 1))
#      dates <- temp[(samp-300):(samp+100),]
#      discharge <- rep(0.1, time = length(dates$dts)) # create a list of non-disturbance discharges
#      discharge[match(samp, dates$dts)] <- 0.3 # from that list of dates from above, assign a disturbance discharge to that date
#      print(length(dates$dts))
#      source("B_1sp_Model.R")
      # run model
#      out <- Bmodel(discharge, dates, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 9, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[f])
      # create summary dataframe 
#      m <- mean.data.frame(out, burnin = 250, iteration = 9)
#      means[d] <- list(m)
#    }
#  jday_max <- unlist(lapply(means, function(x)
#    return(max(x$mean.abund)))) # would mean +/- se or just maximum (and maybe minimum values) be valuable
#  jday.df <- as.data.frame(cbind(jday_max, all.dates))
#  jday.df$jday_max <- as.numeric(jday.df$jday_max)
#  jday.df$all.dates <- yday(as.Date(jday.df$all.dates, "%m-%d"))
#  IDT[f] <- max(jday.df$jday_max)-min(jday.df$jday_max)
#  IDT.ratio[f] <- min(jday.df$jday_max)/max(jday.df$jday_max)
  # close(pb) # close progress bar
#}

#IDT.df <- as.data.frame(cbind(IDT, fecs))
#ggplot(data = IDT.df, aes(x = fecs, y = IDT))+
#  geom_line(linewidth = 1)+
#  ylab('Importance of Disturbance Timing') +
#  xlab("Mean Fecundity")+
#  theme_bw()

#IDT.ratio.df <- as.data.frame(cbind(IDT.ratio, fecs))
#ggplot(data = IDT.ratio.df, aes(x = fecs, y = IDT.ratio))+
#  geom_line(linewidth = 1)+
#  ylab('Importance of Disturbance Timing')+
#  xlab("Mean Fecundity")+
#  theme_bw()
