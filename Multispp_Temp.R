##########
# Multispecies Temp increase Simulation
#############

# Load required scripts with custom functions and models
source("Multispp.R")
source("MultisppFunctions.R")

# Load necessary libraries
library(doParallel)
library(foreach)

# Read in Lees Ferry temperature and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")  # Water temperature data
discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")  # River discharge data

# Calculate average yearly flows from discharge data
flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date")

# Calculate average yearly temperatures from temperature data
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")

# Create a time series of average temperatures extending for 100 years with no temperature change
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Create a time series of average flows for 100 years by repeating the original yearly data
flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))

# Match the dates between temperature and flow datasets
flows$dts <- as.Date(temps$dts)

# Normalize discharge values by dividing by the bankfull discharge (85,000 cfs)
flows$Discharge <- flows$Discharge / 85000

# Create a sequence of hydropeaking intensity levels from 0.00 to 0.70 in increments of 0.05
temp_seq <- c(1,1.1,1.2,1.5)

# Initialize lists to store results for abundance, biomass, and annual biomass calculations
Multispp_temp_abund <- data.frame(timesteps= numeric(), taxa = factor(), abundance=numeric(), sd.abund = numeric(), se.abund = numeric(), temperature=factor())
Multispp_temp_biomass <- data.frame(timesteps= numeric(), taxa = factor(), biomass=numeric(), sd.biomass = numeric(), se.biomass = numeric(), temperature=factor())
MultisppS3_temp_biomass <- data.frame(timesteps= numeric(), taxa = factor(), S3biomass=numeric(), sd.biomass = numeric(), se.biomass = numeric(), temperature=factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  
  # Run multispecies model for abundance under given hydropeaking intensity
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, 
                  disturbanceK = 100000, Qmin = 0.25, extinct = 50, iteration = 10 , 
                  peaklist = 0, peakeach = length(temps$Temperature), stage_output = c("all", "biomass"))
  
  # Calculate mean abundances at each timestep
  means.abund <- multispp.data.frame(out$N, burnin = 260, iteration = 10, value = "abund")

    # Calculate mean biomass at each timestep
  means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 10, value = "biomass")
  
  # Extract mean S3 biomass
  means.s3.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 10, value = "S3.biomass")
  
  temps$Temperature <- temps$Temperature / te
  
  # Append hydropeaking intensity value to abundance data
  average_means <- cbind(means.abund, rep(te, times = length(means.abund$mean.rel.abund)))
  colnames(average_means) <- colnames(Multispp_temp_abund)

  # Calculate the average mean biomass per timestep
  average_biomass <- cbind(means.biomass, rep(te, times = length(means.biomass$mean.rel.biomass)))
  colnames(average_biomass) <- colnames(Multispp_temp_biomass)
  
  # Append only stage 3s
  s3biomass <- cbind(means.s3.biomass, rep(te, times = length(means.s3.biomass$mean.S3.biomass)))
  colnames(s3biomass) <- colnames(MultisppS3_temp_biomass)
  
#   # Create a reference data frame to align time steps with actual years
#   datelist <- as.data.frame(cbind(temps$dts, seq(2, length(temps$dts) + 1, by = 1)))
#   datelist$V1 <- as.factor(year(as.Date(datelist$V1)))  # Extract year
#   datelist$V2 <- as.factor(datelist$V2)
#   
#   # Merge biomass data with reference dates by time step
#   means.s3.biomass <- merge(means.s3.biomass, datelist, by.x = "timesteps", by.y = "V2")
#   
#   # Aggregate S3 biomass by year
#   annual.means <- aggregate(means.s3.biomass$mean.S3.biomass ~ taxa + V1, data = means.s3.biomass, FUN = sum, na.rm = TRUE)
#   
#   # Compute the mean annual biomass for each species across years
#   annual.avg <- as.data.frame(aggregate(`means.s3.biomass$mean.S3.biomass` ~ taxa , data = annual.means, FUN = mean, na.rm = TRUE))
#   
#   # Store annual biomass results
#   annual.biomass[[hyd]] <- annual.avg
#   annual.biomass[[hyd]] <- cbind(annual.biomass[[hyd]], rep(hydropeak[hyd], times = 5))
# }

  # Append results to the main data storage
  Multispp_temp_abund <- rbind(Multispp_temp_abund, average_means)
  Multispp_temp_biomass <- rbind(Multispp_temp_biomass, average_biomass)
  MultisppS3_temp_biomass <- rbind(MultisppS3_temp_biomass, s3biomass)
  # Return results as a list
  return(list(Multispp_temp_hyd_abund = average_means, Multispp_temp_hyd_biomass = average_biomass, MultisppS3_temp_biomass = s3biomass))
}, mc.cores =   1)
#  }, mc.cores =  detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
Multispp_temp_abund <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_abund"))
Multispp_temp_biomass <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_biomass"))
MultisppS3_temp_biomass <- do.call(rbind, lapply(results, `[[`, "MultisppS3_temp_biomass"))
# Write results to CSV files
write.csv(Multispp_temp_abund, "Multispp_temp_abund.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass, "Multispp_temp_biomass.csv", row.names = FALSE)
write.csv(MultisppS3_temp_biomass, "MultisppS3_temp_biomass.csv", row.names = FALSE)



## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:22] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Initialize lists to store results for abundance, biomass, and annual biomass calculations
Multispp_temp_abund_spike <- data.frame(timesteps= numeric(), taxa = factor(), abundance=numeric(), sd.abund = numeric(), se.abund = numeric(), temperature=factor())
Multispp_temp_biomass_spike <- data.frame(timesteps= numeric(), taxa = factor(), biomass=numeric(), sd.biomass = numeric(), se.biomass = numeric(), temperature=factor())
MultisppS3_temp_biomass_spike <- data.frame(timesteps= numeric(), taxa = factor(), S3biomass=numeric(), sd.biomass = numeric(), se.biomass = numeric(), temperature=factor())

results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  # Run multispecies model for abundance under given hydropeaking intensity
  out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, 
                  disturbanceK = 100000, Qmin = 0.25, extinct = 50, iteration = 1000 , 
                  peaklist = 0, peakeach = length(temps$Temperature), stage_output = c("all", "biomass"))
  
  # Calculate mean abundances at each timestep
  means.abund <- multispp.data.frame(out$N, burnin = 260, iteration = 1000, value = "abund")
  
  # Calculate mean biomass at each timestep
  means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "biomass")
  
  # Extract mean S3 biomass
  means.s3.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "S3.biomass")
  
  temps$Temperature <- temps$Temperature / te
  
  # Append hydropeaking intensity value to abundance data
  average_means <- cbind(means.abund, rep(te, times = length(means.abund$mean.rel.abund)))
  colnames(average_means) <- colnames(Multispp_temp_abund_spike)
  
  # Calculate the average mean biomass per timestep
  average_biomass <- cbind(means.biomass, rep(te, times = length(means.biomass$mean.rel.biomass)))
  colnames(average_biomass) <- colnames(Multispp_temp_biomass_spike)
  
  # Append hydropeaking intensity value to biomass data
  s3biomass <- cbind(means.s3.biomass, rep(te, times = length(mean.S3.biomass$S3biomass)))
  colnames(s3biomass) <- colnames(MultisppS3_temp_biomass_spike)
  
  #   # Create a reference data frame to align time steps with actual years
  #   datelist <- as.data.frame(cbind(temps$dts, seq(2, length(temps$dts) + 1, by = 1)))
  #   datelist$V1 <- as.factor(year(as.Date(datelist$V1)))  # Extract year
  #   datelist$V2 <- as.factor(datelist$V2)
  #   
  #   # Merge biomass data with reference dates by time step
  #   means.s3.biomass <- merge(means.s3.biomass, datelist, by.x = "timesteps", by.y = "V2")
  #   
  #   # Aggregate S3 biomass by year
  #   annual.means <- aggregate(means.s3.biomass$mean.S3.biomass ~ taxa + V1, data = means.s3.biomass, FUN = sum, na.rm = TRUE)
  #   
  #   # Compute the mean annual biomass for each species across years
  #   annual.avg <- as.data.frame(aggregate(`means.s3.biomass$mean.S3.biomass` ~ taxa , data = annual.means, FUN = mean, na.rm = TRUE))
  #   
  #   # Store annual biomass results
  #   annual.biomass[[hyd]] <- annual.avg
  #   annual.biomass[[hyd]] <- cbind(annual.biomass[[hyd]], rep(hydropeak[hyd], times = 5))
  # }
  
  # Append results to the main data storage
  Multispp_temp_abund_spike <- rbind(Multispp_temp_abund_spike, average_means)
  Multispp_temp_biomass_spike <- rbind(Multispp_temp_biomass_spike, average_size)
  MultisppS3_temp_biomass_spike <- rbind(MultisppS3_temp_biomass_spike, s3biomass)
  # Return results as a list
  return(list(NZMS_temp_hyd_abund_spike = average_means, NZMS_temp_hyd_biomass_spike = average_size, MultisppS3_temp_biomass_spike = s3biomass))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
Multispp_temp_abund_spike <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_abund_spike"))
Multispp_temp_biomass_spike <- do.call(rbind, lapply(results, `[[`, "Multispp_temp_biomass_spike"))
MultisppS3_temp_biomass_spike <- do.call(rbind, lapply(results, `[[`, "MultisppS3_temp_biomass_spike"))
# Write results to CSV files
write.csv(Multispp_temp_abund_spike, "Multispp_temp_abund_spike.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass_spike, "Multispp_temp_biomass_spike.csv", row.names = FALSE)
write.csv(MultisppS3_temp_biomass_spike, "MultisppS3_temp_biomass_spike.csv", row.names = FALSE)

 # ##########
# # Multispecies Hydropeaking
# #############
# 
# source("Multispp.R")
# source("MultisppFunctions.R")
# # read in Lees Ferry temp and discharge data from 2007 to 2023
# temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# discharge <- readNWISdv("09380000", "00060", "2007-10-01", "2023-05-01")
# 
# # calculate average yearly flows
# flow <- average.yearly.flows(flowdata = discharge, flow.column_name = "X_00060_00003", date.column_name = "Date" )
# # calculate average yearly temperatures
# temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# # create a timeseries of average temperatures 100 years long
# temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)
# # create a timeseries of average flows 100 years long
# flows <- do.call("rbind", replicate(100, flow, simplify = FALSE))
# # match dates
# flows$dts <- as.Date(temps$dts)
# # get discharge magnitude by dividing by bankfull discharge 
# flows$Discharge <- flows$Discharge/85000
# 
# # create sequence of hydropeaking intensities
# abunds <- list()
# biomass <- list()
# annual.biomass <- list()
# hydropeak <- seq(0.00, 0.7, by = 0.05)
# for(hyd in 1:length(hydropeak)){
#   #set.seed(123) # make reproducible
#   # model rel abundances 
#   out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = 2, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature), stage_output = "all")
#   # calculate mean abundances at each timestep
#   means.abund <- multispp.data.frame(out, burnin = 260, iteration = 2, value = "abund")
#   
#   # model rel biomass 
#   out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = 2, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature), stage_output = "biomass")
#   # for each stage, calculate mean biomass
#   means.biomass <- multispp.data.frame(out, burnin = 260, iteration = 2, value = "biomass")
#   
#   # model annual avg biomass
#   out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, disturbanceK = 100000 , Qmin = 0.25, extinct = 50, iteration = 2, peaklist = hydropeak[hyd], peakeach = length(temps$Temperature), stage_output = "biomass")
#   means.s3.biomass <- multispp.data.frame(out, burnin = 260, iteration = 2, value = "S3.biomass")
#   
#   # calculate the average of mean abundances at each hydropeaking intensity
#   abunds[[hyd]] <- means.abund %>% 
#     dplyr::group_by(taxa) %>%
#     dplyr::summarise(avg.abund = mean(mean.rel.abund, na.rm = T)) %>%
#                        ungroup()
#   abunds[[hyd]] <- cbind(abunds[[hyd]], rep(hydropeak[hyd], times = 5))
#   
#   # calculate the average mean biomass per timestep
#   biomass[[hyd]] <- means.biomass %>%
#     dplyr::group_by(taxa) %>%
#     dplyr::summarise(avg.biomass = mean(mean.rel.biomass, na.rm = T)) %>%
#     ungroup()
#   biomass[[hyd]] <- cbind(biomass[[hyd]], rep(hydropeak[hyd], times = 5))
#   
#   # annual S3 biomasses
#   # make sure we can align date info with timestep
#   datelist <- as.data.frame(cbind(temps$dts, seq(2, length(temps$dts+1), by = 1)))
#   datelist$V1 <- as.factor(year(as.Date(datelist$V1))) # only want year
#   datelist$V2 <- as.factor(datelist$V2)
#   # merge data frames by timestep and year
#   means.s3.biomass <- merge(means.s3.biomass, datelist, by.x = "timesteps", by.y =  "V2")
#   # calculate average annual biomass for each year
#   annual.means <- aggregate(means.s3.biomass$mean.S3.biomass ~ taxa + V1, data = means.s3.biomass, FUN = sum, na.rm = TRUE)
#   #average the averages
#   annual.avg <- as.data.frame(aggregate(`means.s3.biomass$mean.S3.biomass` ~ taxa , data =annual.means, FUN = mean, na.rm = TRUE))
#   annual.biomass[[hyd]] <- annual.avg
#   
#   annual.biomass[[hyd]] <- cbind(annual.biomass[[hyd]], rep(hydropeak[hyd], times = 5))
# 
# 
# }
# hydropeak.biomass <- do.call("rbind", biomass)
# hydropeak.abunds<- do.call("rbind", abunds)
# hydropeak.avg.annual <- do.call("rbind", annual.biomass)

