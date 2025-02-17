###################
## Multispp Sensitivity
##################
library(ggplot2)
library(dplyr)
library(dataRetrieval)
library(parallel)

source("1spFunctions.R")
source("HYOSSurvivorship.R")
source("BAETSurvivorship.R")
source("NZMS_shell_length_fecundity.R")
source("NZMSSurvivorship.R")
source("CHIRSurvivorship.R")
source("GAMMSurvivorship.R")
source("MultisppFunctions.R")
source("Multispp.R")

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

# Define your species-stage parameters
parameters <- c("G1_HYOS", "G1_BAET", "G1_NZMS", "G1_GAMM", "G1_CHIR", 
                "G2_HYOS", "G2_BAET", "G2_NZMS", "G2_GAMM", "G2_CHIR")

# Sensitivity increments from -0.01 to +0.01 in 0.001 steps
increments <- seq(-0.01, 0.01, by = 0.001)



# Create a grid of all parameter/increment combinations
param_grid <- expand.grid(Parameter = parameters, Increment = increments, 
                          stringsAsFactors = FALSE)

# Function to modify the parameter and run the model
run_model_with_modification <- function(param, inc) {
  # Modify the model based on param and increment (Replace this with your actual model function)
  modified_model <- Multispp(flow.data = flows$Discharge, temp.data = temps, 
                            baselineK =10000 , disturbanceK = 100000, Qmin = 0.25,
                            extinct = 50, iteration = 1000, peaklist = 0, 
                            peakeach = length(temps$Temperature), stage_output = "all", 
                            modify_parameter = param, increment =inc)

    average_abundance <- apply(modified_model$N[-c(1:199),,,], c(2, 4), function(x) mean(x, na.rm = TRUE))
    # Create a data frame from the average abundances
    average_abundance_df <- as.data.frame(as.table(average_abundance))
    # Combine stage and taxon into a single column
    average_abundance_df$StageGroup <- paste(average_abundance_df$Var1,"_",average_abundance_df$Var2, sep = "")
    # Rename columns for clarity
    average_abundance_df <- as.data.frame(cbind(average_abundance_df[ , c("StageGroup", "Freq")], 
                                                rep(param, times = length(average_abundance_df$Var1)),
                                                rep(inc,times = length(average_abundance_df$Var1))))
    colnames(average_abundance_df) <- c("StageGroup", "Abundance", "Parameter", "SensitivityIncrement")
    # Extract relevant values (ensure your model function returns these)
    data.frame(
      average_abundance_df
  )
}

# Determine the number of cores to use (leave one core free)
numCores <- detectCores() - 1

# Use mclapply to loop over each row in the grid in parallel.
# Note: mclapply returns a list that we then combine into a single data frame.
sensitivity_results_list <- mclapply(1:nrow(param_grid), function(i) {
  param <- as.character(param_grid$Parameter[i])
  inc <- param_grid$Increment[i]
  run_model_with_modification(param, inc)
}, numCores)

# Combine the list into a data frame
sensitivity_results <- do.call(rbind, sensitivity_results_list)

# Save the data frame as a CSV file
write.csv(sensitivity_results, "Multispp_sensitivity_vitalrates.csv", row.names = FALSE)


