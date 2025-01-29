############################
# Hydropeaking Sensitivity BAET
#############################

source("1spFunctions.R")
source("BAET_1sp_Model.R")


# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")

# calculate average yearly flows
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date" )
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)
# create a timeseries of average flows 100 years long
flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))
# match dates
flows$dts <- as.Date(temps$dts)
# get discharge magnitude by dividing by bankfull discharge 
flows$Discharge <- flows$Discharge/85000

# create sequence of hydropeaking intensities
temp_seq <- seq(0, 10, by = 1)

# makes some vectors for data to go into
means <- vector()
sd <- vector()
sizemeans <- vector()
sizesd <- vector()
S3Yrprod <- vector()
S3Yrprod_sd <- vector()
# cycle though hydropeaking scenarios
for (te in 1:length(temp_seq)){
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature + temp_seq[te]
  
  sizes <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 2, peaklist = 0.17, peakeach = length(temps$Temperature), stage_output = "size")
  # model abundances
  out <- BAETmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 40000 , Qmin = 0.1, extinct = 50, iteration = 2, peaklist = 0.17, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature - temp_seq[te]
  
  # for each stage, calculate mean biomass
  s1s <- mean(out[-c(1:260), 1, ]) * (0.0053 * (mean(sizes[-c(1:260)])/2)^2.875)
  s2s <- mean(out[-c(1:260), 2,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  s3s <- mean(out[-c(1:260), 3,]) * (0.0053 * (mean(sizes[-c(1:260)]))^2.875)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizemeans[te] <- sum(c(s1s, s2s, s3s))
  # produce standard devation
  sizesd[te] <- sd(c(s1s, s2s, s3s))
  
  # now instead of getting the mean, calculate biomass at every timestep
  s1ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 1, ]) * (0.0053 * ((sizes[-c(1:260)])/2)^2.875), year(temps$dts[-c(1:259)])))
  s2ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 2, ]) * (0.0053 * ((sizes[-c(1:260)]))^2.875), year(temps$dts[-c(1:259)])))
  s3ss <- as.data.frame(cbind(rowMeans(out[-c(1:260), 3, ]) * (0.0053 * ((sizes[-c(1:260)]))^2.875), year(temps$dts[-c(1:259)])))
  
  # find annual stage specific biomass based on year
  s1sYr <- aggregate(V1 ~ V2, data = s1ss, FUN = sum, na.rm = TRUE)
  s2sYr <- aggregate(V1 ~ V2, data = s2ss, FUN = sum, na.rm = TRUE)
  s3sYr <- aggregate(V1 ~ V2, data = s3ss, FUN = sum, na.rm = TRUE)
  
  # add all stages together to get average annual biomass (aka secondary production)
  #Yrprod[hyd] <- sum(mean(c(s1sYr$V1, s2sYr$V1, s3sYr$V1), na.rm = T))
  
  # average annual biomass of only stage 3 (emergent adults)
  S3Yrprod[te] <- mean( s3sYr$V1, na.rm = T)
  S3Yrprod_sd[te] <- sd(s3sYr$V1, na.rm = T)
  # calculate mean abundances at each timestep
  means.list.BAET <- mean.data.frame(out, burnin = 260, iteration = 2)
  # calculate the average of mean abundances at each hydropeaking intensity
  means[te] <- mean(means.list.BAET$mean.abund)  
  # calculate the standard deviation of mean abundances at each hydropeaking intensity
  sd[te] <- sd(means.list.BAET$mean.abund, na.rm = T)
}

# compile abundance data
BAET_te_means <- as.data.frame(cbind(temp_seq, means, sd, rep("BAET", length(means))))
# compile timeseries biomass data
BAET_te_size <- as.data.frame(cbind(temp_seq, sizemeans, sizesd, rep("BAET", length(sizemeans))))
# compile annual biomass data
BAET_te_yrprod <- as.data.frame(cbind(temp_seq, S3Yrprod, S3Yrprod_sd, rep("BAET", length(sizemeans))))


# ggplot(data = BAET_hyd_means, aes(x = hydropeak,  y = means, group = 1, color = "BAET")) +
#   geom_ribbon(aes(ymin = means - sd,
#                   ymax = means + sd),
#               colour = 'transparent',
#               alpha = .15,
#               show.legend = T) +
#   geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
#   #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
#   #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
#   #coord_cartesian(ylim = c(0,6000)) +S1 and S2 (inds/m2)
#   ylab('BAET spp. Abund.') +
#   xlab("")+
#   labs(colour=" ")+
#   theme_bw()+
#   scale_color_manual(values = colors)+
#   #scale_y_continuous(
#   # sec.axis = sec_axis(~., name="Baetidae Larvae (inds/m2)"
#   # ))+
#   theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), )
# 
# 
# ggplot(data = BAET_hyd_size, aes(x = hydropeak,  y = sizemeans, group = 1, color = "BAET")) +
#   geom_ribbon(aes(ymin = sizemeans - sizesd,
#                   ymax = sizemeans + sizesd),
#               colour = 'transparent',
#               alpha = .15,
#               show.legend = T) +
#   geom_line(show.legend = T, linewidth = 1, alpha = 0.8) +
#   #geom_line(data = flow.magnitude, aes(x = as.Date(dts), y = X_00060_00003), color = "blue") +
#   #geom_line(data = temps, aes(x = as.Date(dts), y = Temperature*1000), color = "green")+
#   #coord_cartesian(ylim = c(0,6000)) +S1 and S2 (inds/m2)
#   ylab('BAET spp. Abund.') +
#   xlab("")+
#   labs(colour=" ")+
#   theme_bw()+
#   scale_color_manual(values = colors)+
#   #scale_y_continuous(
#   # sec.axis = sec_axis(~., name="Baetidae Larvae (inds/m2)"
#   # ))+
#   theme(text = element_text(size = 13), axis.text.x = element_text(angle=45, hjust = 1, size = 12.5), 
#         axis.text.y = element_text(size = 13), )
# 