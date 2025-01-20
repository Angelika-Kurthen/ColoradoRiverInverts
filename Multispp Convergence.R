#############
# Convergence 
#############

library(ggpubr)
source("Multispp.R")
source("MultisppFunctions.R")
# read in Lees Ferry temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")

# calculate average yearly flows
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date" )
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 50, change.in.temp = 0, years.at.temp = 0)
# create a timeseries of average flows 100 years long
flows <- do.call("rbind", replicate(50, flow, simplify = FALSE))
# match dates
flows$dts <- as.Date(temps$dts)
# get discharge magnitude by dividing by bankfull discharge 
flows$Discharge <- flows$Discharge/85000

# create sequence of hydropeaking intensities
abunds <- list()
biomass <- list()
annual.biomass <- list()
its <- seq(2, 199, by = 15)
for(it in 1:length(its)){
  #set.seed(123) # make reproducible
  # model rel abundances 
  print(its[it])
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = its[it], peaklist = 0.2, peakeach = length(temps$Temperature), stage_output = "all")
  means.abund <- multispp.data.frame(out, burnin = 260, iteration = its[it], value = "abund")
  
  # model rel biomass 
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = its[it], peaklist = 0.2, peakeach = length(temps$Temperature), stage_output = "biomass")
  # for each stage, calculate mean biomass
  means.biomass <- multispp.data.frame(out, burnin = 260, iteration = its[it], value = "biomass")
  
  # model annual avg biomass
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = its[it], peaklist =0.2,  peakeach = length(temps$Temperature), stage_output = "biomass")
  means.s3.biomass <- multispp.data.frame(out, burnin = 260, iteration = its[it], value = "S3.biomass")
  
  # calculate the average of mean abundances at each hydropeaking intensity
  abunds[[it]] <- means.abund %>% 
    dplyr::group_by(taxa) %>%
    dplyr::summarise(avg.abund = mean(mean.rel.abund, na.rm = T)) %>%
    ungroup()
  abunds[[it]] <- cbind(abunds[[it]], rep(its[it], times = 5))
  
  # calculate the average mean biomass per timestep
  biomass[[it]] <- means.biomass %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(avg.biomass = mean(mean.rel.biomass, na.rm = T)) %>%
    ungroup()
  biomass[[it]] <- cbind(biomass[[it]], rep(its[it], times = 5))
  
  # annual S3 biomasses
  # make sure we can align date info with timestep
  datelist <- as.data.frame(cbind(temps$dts, seq(2, length(temps$dts+1), by = 1)))
  datelist$V1 <- as.factor(year(as.Date(datelist$V1))) # only want year
  datelist$V2 <- as.factor(datelist$V2)
  # merge data frames by timestep and year
  means.s3.biomass <- merge(means.s3.biomass, datelist, by.x = "timesteps", by.y =  "V2")
  # calculate average annual biomass for each year
  annual.means <- aggregate(means.s3.biomass$mean.S3.biomass ~ taxa + V1, data = means.s3.biomass, FUN = sum, na.rm = TRUE)
  #average the averages
  annual.avg <- as.data.frame(aggregate(`means.s3.biomass$mean.S3.biomass` ~ taxa , data =annual.means, FUN = mean, na.rm = TRUE))
  annual.biomass[[it]] <- annual.avg
  
  annual.biomass[[it]] <- cbind(annual.biomass[[it]], rep(its[it], times = 5))
  
  
}
iteration.biomass2 <- do.call("rbind", biomass)
iteration.abunds2<- do.call("rbind", abunds)
iteration.avg.annual2 <- do.call("rbind", annual.biomass)

iteration.biomass <- rbind(iteration.biomass, iteration.biomass2)
iteration.abunds <- rbind(iteration.abunds, iteration.abunds2)
iteration.avg.annual <- rbind(iteration.avg.annual, iteration.avg.annual2)

iter.abund<- ggplot(data = iteration.abunds, aes(`rep(its[it], times = 5)`, avg.abund, group = taxa, color = taxa))+ 
  geom_point(alpha = 0.8)+
  # geom_ribbon(data = Hydro_Abund, aes(ymin = scale(means - sd),
  #                 ymax = scale(means + sd), fill = V4),
  #                        alpha = .15,
  #                        show.legend = T) +
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linetype = "dotted", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.17 ), linetype="dashed", 
             linewidth=1)+
  geom_vline(aes(xintercept = 0.55 ), linetype = "dotdash", 
             linewidth=1)+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme_bw()+
  xlab("Iterations")+
  ylab("Relative Abundance")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

iter_ts_multi <- ggplot(data = iteration.biomass, aes(`rep(its[it], times = 5)`,avg.biomass , color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  xlab("Iterations")+
  ylab("Relative Biomass")+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

iter_yr_multi <- ggplot(data = iteration.avg.annual, aes(`rep(its[it], times = 5)`, log(`means.s3.biomass$mean.S3.biomass`), group = taxa, color = taxa)) + 
  # geom_ribbon(aes(ymin = sizemeans - sizesd,
  #                   ymax = sizemeans + sizesd),
  #               colour = V3,
  #               alpha = .15,
  #               show.legend = T) +
  geom_line(linewidth = 1, alpha = 0.8)+
  scale_color_manual(name = "Taxon", labels=c("Hydropsyche spp.", "Baetidae spp.", "P. anitpodarum", "Chironomidae spp.", "G. lacustris"), values=c("#AA3377","#66CCEE","#4477AA","#228833", "#CCBB44"))+
  geom_vline(aes(xintercept = 0.01), linewidth = 1, linetype = "dotted")+
  geom_vline(aes(xintercept = 0.17 ), linewidth = 1, linetype="dashed")+ 
  geom_vline(aes(xintercept = 0.55 ), linewidth = 1, linetype = "dotdash")+
  theme_bw()+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 4))+
  xlab("Hydropeaking Intensity")+
  ylab(" Annual S3 Biomass (mg)")+
  theme(text = element_text(size = 14), axis.text.x = element_text(hjust = 1, size = 12.5), 
        axis.text.y = element_text(size = 13), legend.key = element_rect(fill = "transparent"))

ggarrange(iter.abund, iter_ts_multi , iter_yr_multi,
          labels = c("a", "b", "c"),
          ncol = 3, nrow = 1, common.legend = T)

# going to run it 99x to converge on a mean - abundance converges around 60, Biomass, emegence about 50

