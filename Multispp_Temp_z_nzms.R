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

# Create a sequences
temp_seq <- c(1,1.1,1.2,1.5)
q_seq <- c(1, 2, 4, 8)

# Initialize final storage data frames
Multispp_temp_abund <- data.frame(timesteps= numeric(), taxa = factor(), abundance=numeric(), temperature=factor(), q = factor())
Multispp_temp_biomass <- data.frame(timesteps= numeric(), taxa = factor(), biomass=numeric(), temperature=factor(), q = factor())
MultisppS3_temp_biomass <- data.frame(timesteps= numeric(), taxa = factor(), S3biomass=numeric(), temperature=factor(), q = factor())

# Parallel processing across temp_seq
results <- mclapply(temp_seq, function(te) {
  set.seed(123)

  # Loop over q_seq inside temp_seq
  temp_results <- lapply(q_seq, function(q_val) {

    # Modify inputs
    qs <- list(HYOS = 1, BAET = 1, NZMS = q_val, CHIR = 1, GAMM = 1)
    temps$Temperature <- temps$Temperature * te

    # Run model
    out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000,
                    disturbanceK = 100000, Qmin = 0.25, extinct = 50, iteration = 1000,
                    peaklist = 0, peakeach = length(temps$Temperature),
                    stage_output = c("all", "biomass"), q = qs)

    # Process outputs
    means.abund <- multispp.data.frame(out$N, burnin = 260, iteration = 1000, value = "abund")
    means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "biomass")
    means.s3.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "S3.biomass")

    # Reset input scaling
    temps$Temperature <- temps$Temperature / te

    # Convert to data frame & add metadata
    average_means <- data.frame(means.abund, temperature = rep(te, nrow(means.abund)), q = rep(q_val, nrow(means.abund)))
    colnames(average_means) <- colnames(Multispp_temp_abund)

    average_biomass <- data.frame(means.biomass, temperature = rep(te, nrow(means.biomass)), q = rep(q_val, nrow(means.biomass)))
    colnames(average_biomass) <- colnames(Multispp_temp_biomass)

    s3biomass <- data.frame(means.s3.biomass, temperature = rep(te, nrow(means.s3.biomass)), q = rep(q_val, nrow(means.s3.biomass)))
    colnames(s3biomass) <- colnames(MultisppS3_temp_biomass)

    # Ensure factor consistency
    average_means$taxa <- as.factor(average_means$taxa)
    average_means$temperature <- as.factor(average_means$temperature)
    average_means$q <- as.factor(average_means$q)

    average_biomass$taxa <- as.factor(average_biomass$taxa)
    average_biomass$temperature <- as.factor(average_biomass$temperature)
    average_biomass$q <- as.factor(average_biomass$q)

    s3biomass$taxa <- as.factor(s3biomass$taxa)
    s3biomass$temperature <- as.factor(s3biomass$temperature)
    s3biomass$q <- as.factor(s3biomass$q)

    return(list(Multispp_temp_abund = average_means,
                Multispp_temp_biomass = average_biomass,
                MultisppS3_temp_biomass = s3biomass))
  })

  return(temp_results)
}, mc.cores = 1)

}, mc.cores = detectCores() - 1)

# Flatten the nested results (convert list of lists to a single list)
flat_results <- do.call(c, results)

# Combine results from all temperature and q scenarios into final data frames
Multispp_temp_abund <- do.call(rbind, lapply(flat_results, `[[`, "Multispp_temp_abund"))
Multispp_temp_biomass <- do.call(rbind, lapply(flat_results, `[[`, "Multispp_temp_biomass"))
MultisppS3_temp_biomass <- do.call(rbind, lapply(flat_results, `[[`, "MultisppS3_temp_biomass"))
# Write results to CSV files
write.csv(Multispp_temp_abund, "Multispp_temp_abund_z_nzms_full.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass, "Multispp_temp_biomass_z_nzms_full.csv", row.names = FALSE)
write.csv(MultisppS3_temp_biomass, "MultisppS3_temp_biomass_z_nzms.csv", row.names = FALSE)

Multispp_temp_abund <- Multispp_temp_abund %>%
  group_by(temperature, taxa, q) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE),
            sd_abund = sd(abundance, na.rm = TRUE), .groups = "drop")

Multispp_temp_biomass <- Multispp_temp_biomass %>%
  group_by(temperature, taxa, q) %>%
  summarise(mean_biomass = mean(biomass, na.rm = TRUE),
            sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")

# Write results to CSV files
write.csv(Multispp_temp_abund, "Multispp_temp_abund_z_nzms.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass, "Multispp_temp_biomass_z_nzms.csv", row.names = FALSE)

## Now we add the summer spike to temperatures 
# read in LF temp and discharge data from 2007 to 2023
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")
# calculate average yearly temperatures
temps <- average.yearly.temp(tempdata = temp, temp.column_name = "X_00010_00003", date.column_name = "Date")
# create summertime spike (up to 21 C, then scale from there)
temps$Temperature[16:23] <- c(14, 16, 18, 21, 21, 18, 16, 14)

# create a timeseries of average temperatures 100 years long
temps <- rep.avg.year(temps, n = 100, change.in.temp = 0, years.at.temp = 0)

# Initialize lists to store results for abundance, biomass, and annual biomass calculations
Multispp_temp_abund_spike <- data.frame(timesteps= numeric(), taxa = factor(), abundance=numeric(), temperature=factor(), q = factor())
Multispp_temp_biomass_spike <- data.frame(timesteps= numeric(), taxa = factor(), biomass=numeric(),  temperature=factor(), q = factor())
MultisppS3_temp_biomass_spike <- data.frame(timesteps= numeric(), taxa = factor(), S3biomass=numeric(), temperature=factor(), q = factor())
# Parallel processing across temp_seq
results <- mclapply(temp_seq, function(te) {
  set.seed(123)
  
  # Loop over q_seq inside temp_seq
  temp_results <- lapply(q_seq, function(q_val) {

    # Modify inputs
    qs <- list(HYOS = 1, BAET = 1, NZMS = q_val, CHIR = 1, GAMM = 1)
    temps$Temperature <- temps$Temperature * te
    
    # Run model
    out <- Multispp(flow.data = flows$Discharge, temp.data = temps, baselineK = 10000, 
                    disturbanceK = 100000, Qmin = 0.25, extinct = 50, iteration = 1000, 
                    peaklist = 0, peakeach = length(temps$Temperature), 
                    stage_output = c("all", "biomass"), q = qs)
    
    # Process outputs
    means.abund <- multispp.data.frame(out$N, burnin = 260, iteration = 1000, value = "abund")
    means.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration =1000, value = "biomass")
    means.s3.biomass <- multispp.data.frame(out$biomass, burnin = 260, iteration = 1000, value = "S3.biomass")
    
    # Reset input scaling_
    temps$Temperature <- temps$Temperature / te
    
    # Convert to data frame & add metadata
    average_means <- data.frame(means.abund, temperature = rep(te, nrow(means.abund)), q = rep(q_val, nrow(means.abund)))
    colnames(average_means) <- colnames(Multispp_temp_abund_spike)
    
    average_biomass <- data.frame(means.biomass, temperature = rep(te, nrow(means.biomass)), q = rep(q_val, nrow(means.biomass)))
    colnames(average_biomass) <- colnames(Multispp_temp_biomass_spike)
    
    s3biomass <- data.frame(means.s3.biomass, temperature = rep(te, nrow(means.s3.biomass)), q = rep(q_val, nrow(means.s3.biomass)))
    colnames(s3biomass) <- colnames(MultisppS3_temp_biomass_spike)
    
    # Ensure factor consistency
    average_means$taxa <- as.factor(average_means$taxa)
    average_means$temperature <- as.factor(average_means$temperature)
    average_means$q <- as.factor(average_means$q)
    
    average_biomass$taxa <- as.factor(average_biomass$taxa)
    average_biomass$temperature <- as.factor(average_biomass$temperature)
    average_biomass$q <- as.factor(average_biomass$q)
    
    s3biomass$taxa <- as.factor(s3biomass$taxa)
    s3biomass$temperature <- as.factor(s3biomass$temperature)
    s3biomass$q <- as.factor(s3biomass$q)
    
    return(list(Multispp_temp_abund_spike = average_means, 
                Multispp_temp_biomass_spike = average_biomass, 
                MultisppS3_temp_biomass_spike = s3biomass))
  })
  
  return(temp_results)
}, mc.cores = 1)

}, mc.cores = detectCores() - 1)

# Flatten the nested results (convert list of lists to a single list)
flat_results <- do.call(c, results)

# Final binding
Multispp_temp_abund_spike <- do.call(rbind, lapply(flat_results, `[[`, "Multispp_temp_abund_spike"))
Multispp_temp_biomass_spike <- do.call(rbind, lapply(flat_results, `[[`, "Multispp_temp_biomass_spike"))
MultisppS3_temp_biomass_spike <- do.call(rbind, lapply(flat_results, `[[`, "MultisppS3_temp_biomass_spike"))

write.csv(Multispp_temp_abund_spike, "Multispp_temp_abund_spike_z_nzms_full.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass_spike, "Multispp_temp_biomass_spike_z_nzms_full.csv", row.names = FALSE)
write.csv(MultisppS3_temp_biomass_spike, "MultisppS3_temp_biomass_spike_z_nzms.csv", row.names = FALSE)

Multispp_temp_abund_spike <- Multispp_temp_abund_spike %>%
  group_by(temperature, taxa, q) %>%
  summarise(mean_abund = mean(abundance, na.rm = TRUE),
            sd_abund = sd(abundance, na.rm = TRUE), .groups = "drop")


Multispp_temp_biomass_spike <- Multispp_temp_biomass_spike %>%
  group_by(temperature, taxa, q) %>%
  summarise(mean_biomass = mean(biomass, na.rm = TRUE),
            sd_biomass = sd(biomass, na.rm = TRUE), .groups = "drop")

write.csv(Multispp_temp_abund_spike, "Multispp_temp_abund_spike_z_nzms.csv", row.names = FALSE)
write.csv(Multispp_temp_biomass_spike, "Multispp_temp_biomass_spike_z_nzms.csv", row.names = FALSE)

