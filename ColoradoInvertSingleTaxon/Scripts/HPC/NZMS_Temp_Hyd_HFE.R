############################
# Temperature % increase + Hyd + HFE NZMS
#############################

# Load required libraries
library(parallel)  # For parallel processing
library(dataRetrieval)
# Load required functions and model
source("1spFunctions.R")
source("NZMS_1sp_Model.R")


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
# Get row indices for high flow events
HFE_rows <- which(month(flows$dts) == 3 & day(flows$dts) == 31)
HFE_rows <- HFE_rows[seq(1, length(HFE_rows), by = 3)]  # Keep every 3rd occurrence


# Update Discharge to 0.45 for those dates
flows$Discharge[HFE_rows] <- 0.45

# Expand to include the next 6 rows for each event
#expanded_HFE_rows <- unlist(lapply(HFE_rows, function(row) seq(row, row + 6)))

#if we want to look at everything
expanded_HFE_rows <- seq(261, length(temps$dts))

# Ensure we don't exceed the number of rows in the dataset
expanded_HFE_rows <- expanded_HFE_rows[expanded_HFE_rows <= nrow(flows) & expanded_HFE_rows > 260]

# Ensure we don't exceed the number of rows in the dataset
expanded_HFE_rows <- expanded_HFE_rows[expanded_HFE_rows <= nrow(flows)]

# create sequence of temperature modifiers
temp_seq <- c(1, 1.1, 1.2, 1.5)

# makes some vectors for data to go into
NZMS_temp_hyd_abund_HFE <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
NZMS_temp_hyd_biomass_HFE <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # model sizes
  temps$Temperature <- temps$Temperature * te
  out <- NZMSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 2000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  # calculate mean abundances at each timestep
  means.list.NZMS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.NZMS$mean.abund[which(means.list.NZMS$timesteps %in% expanded_HFE_rows)]
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("NZMS",times = length(means)))
  colnames(average_means) <- colnames(NZMS_temp_hyd_abund_HFE)
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.02 * mean(c(0.5, 3.2))^2.4315)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.02 * mean(c(3.2, 4))^2.4315)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.02 * mean(c(4, 5.5))^2.4315)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("NZMS",times = length(sizes_list)))
  colnames(average_size) <- colnames(NZMS_temp_hyd_biomass_HFE)

  # Append results to the main data storage
  NZMS_temp_hyd_abund_HFE <- rbind(NZMS_temp_hyd_abund_HFE, average_means)
  NZMS_temp_hyd_biomass_HFE <- rbind(NZMS_temp_hyd_biomass_HFE, average_size)
  # Return results as a list
  return(list(NZMS_temp_hyd_abund_HFE = average_means, NZMS_temp_hyd_biomass_HFE = average_size))
}, mc.cores = detectCores() - 1)

# Combine results from all temperature scenarios into final dataframes
NZMS_temp_hyd_abund_HFE <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_abund_HFE"))
NZMS_temp_hyd_biomass_HFE <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_biomass_HFE"))

# Write results to CSV files
write.csv(NZMS_temp_hyd_abund_HFE, "NZMS_temp_hyd_abund_HFE.csv", row.names = FALSE)
write.csv(NZMS_temp_hyd_biomass_HFE, "NZMS_temp_hyd_biomass_HFE.csv", row.names = FALSE)

# tabula rasa
rm(temp)
rm(temps)
#rm(NZMS_temp_hyd_abund_HFE)
#rm(NZMS_temp_hyd_biomass_HFE)

## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)


# makes some vectors for data to go into
NZMS_temp_hyd_abund_HFE_spike <- data.frame(abundance=numeric(), temperature=factor(), taxa=factor())
NZMS_temp_hyd_biomass_HFE_spike <- data.frame(biomass=numeric(), temperature = factor(), taxa = factor())
results <- mclapply(temp_seq, function(te) {
  set.seed(123) # make reproducible
  # modify temperature regim
  temps$Temperature <- temps$Temperature * te
  # add temp spike in September
  out <- NZMSmodel(flow.data = flows$Discharge, temp.data = temps, baselineK = 5000, disturbanceK = 2000 , Qmin = 0.3, extinct = 50, iteration = 1000, peaklist = 0.17, peakeach = length(temps$Temperature))
  temps$Temperature <- temps$Temperature / te
  # calculate mean abundances at each timestep
  means.list.NZMS <- mean.data.frame(out, burnin = 260, iteration = 1000)
  means <- means.list.NZMS$mean.abund[which(means.list.NZMS$timesteps %in% expanded_HFE_rows)]
  # Store abundance data in a dataframe
  average_means <- cbind(means, rep(te, times = length(means)), rep("NZMS",times = length(means)))
  colnames(average_means) <- colnames(NZMS_temp_hyd_abund_HFE_spike)
  # for each stage, calculate mean biomass
  s1s <- colMeans(out[-c(1:260), 1, ]) * (0.02 * mean(c(0.5, 3.2))^2.4315)
  s2s <- colMeans(out[-c(1:260), 2, ]) * (0.02 * mean(c(3.2, 4))^2.4315)
  s3s <- colMeans(out[-c(1:260), 3, ]) * (0.02 * mean(c(4, 5.5))^2.4315)
  # sum the mean biomass of each stage to get mean timestep biomass
  sizes_list <- as.vector(s1s + s2s + s3s)
  # Store biomass data in a dataframe
  average_size <- cbind(sizes_list, rep(te, times = length(sizes_list)), rep("NZMS",times = length(sizes_list)))
  colnames(average_size) <- colnames(NZMS_temp_hyd_biomass_HFE_spike)
  
  # Append results to the main data storage
  NZMS_temp_hyd_abund_HFE_spike <- rbind(NZMS_temp_hyd_abund_HFE_spike, average_means)
  NZMS_temp_hyd_biomass_HFE_spike <- rbind(NZMS_temp_hyd_biomass_HFE_spike, average_size)
  # Return results as a list
  return(list(NZMS_temp_hyd_abund_HFE_spike = average_means, NZMS_temp_hyd_biomass_HFE_spike = average_size))
}, detectCores() -1 )

# Combine results from all temperature scenarios into final dataframes
NZMS_temp_hyd_abund_HFE_spike <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_abund_HFE_spike"))
NZMS_temp_hyd_biomass_HFE_spike <- do.call(rbind, lapply(results, `[[`, "NZMS_temp_hyd_biomass_HFE_spike"))

# Write results to CSV files
write.csv(NZMS_temp_hyd_abund_HFE_spike, "NZMS_temp_hyd_abund_HFE_spike.csv", row.names = FALSE)
write.csv(NZMS_temp_hyd_biomass_HFE_spike, "NZMS_temp_hyd_biomass_HFE_spike.csv", row.names = FALSE)

