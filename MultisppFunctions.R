###############
#Multispp Functions
###############

library(data.table)

multispp.data.frame <- function(data, burnin = NULL, iteration, value) {
  # Check if the input data is a 4D array
  if (length(dim(data)) != 4) {
    stop("Input data must be a 4-dimensional array.")
  }
  
  # Extract dimensions
  dims <- dim(data)
  n_timesteps <- dims[1]
  n_stages <- dims[2]
  n_reps <- dims[3]
  n_taxa <- dims[4]
  
  # Reshape the 4D array into a 2D matrix
  data_matrix <- matrix(data, nrow = n_timesteps * n_stages * n_reps, ncol = n_taxa)
  
  # Create a data.table from the matrix
  dt <- as.data.table(data_matrix)
  
  # Add identifier columns
  dt[, timesteps := rep(1:n_timesteps, each = n_stages * n_reps)]
  dt[, stage := rep(rep(1:n_stages, each = n_reps), times = n_timesteps)]
  dt[, rep := rep(1:n_reps, times = n_timesteps * n_stages)]

  # Melt the data.table to long format
  dt_long <- melt(dt, id.vars = c("timesteps", "stage", "rep"), variable.name = "taxa", value.name = "value")
  
  
  # Assign the desired taxa names
  dt_long[, taxa := factor(taxa, levels = paste0("V", 1:n_taxa), labels = c("HYOS", "BAET", "NZMS", "CHIR", "GAMM"))]
  # Calculate total biomass for each combination of timesteps, stage, and rep
  dt_long[, total := sum(value), by = .(timesteps, rep)]
  
  # Calculate relative biomass
  dt_long[, rel.value := value / total]
  
  # Apply burnin if specified
  if (!is.null(burnin)) {
    dt_long <- dt_long[timesteps > burnin]
  }
  
  # Calculate mean and standard deviation of relative biomass by timesteps and taxa
  if (value == "biomass") {
  means.list <- dt_long[, .(
    mean.rel.biomass = mean(rel.value, na.rm = TRUE), #mean
    sd.rel.biomass = sd(rel.value, na.rm = TRUE) #sd
  ), by = .(timesteps, taxa)]
  
  } else if (value == "abund") {
    means.list <- dt_long[, .(
      mean.rel.abund = mean(rel.value, na.rm = TRUE), #mean
      sd.rel.abund = sd(rel.value, na.rm = TRUE) #sd
    ), by = .(timesteps, taxa)]
    
  } else if (value == "S3.biomass") {
    # Filter for stage 3
    dt3 <- dt_long[stage == "3"]
    
      means.list <- dt3[, .(
      mean.S3.biomass = mean(rel.value, na.rm = TRUE), #mean
      sd.abund = sd(rel.value, na.rm = TRUE) #sd
      ), by = .(timesteps, taxa)]
  }

  return(means.list)
}


#dat <- multispp.data.frame(out$biomass, burnin = 260, 1000, value = "S3.biomass")
