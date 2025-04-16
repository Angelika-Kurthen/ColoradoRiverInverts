###############
#Multispp Functions
###############
multispp.data.frame <- function(data, burnin, iteration, value = "biomass"){
  
  pop_matrix_perm <- aperm(data, c(1, 3, 4, 2))  # Now: [2601, 1000, 5, 3]
  # Sum over stages directly with rowSums (preserving dimensions)
  taxon_stage_sums <- rowSums(pop_matrix_perm, dims = 3)
  
  # Total sum across all taxa
  total_stage_sums <- rowSums(taxon_stage_sums, dims = 2)
  
  # Compute relative abundance
  relative_abundance <- sweep(taxon_stage_sums, c(1, 2), total_stage_sums, "/")
  
  
  dimnames(relative_abundance) <- list(
    timesteps = 1:dim(relative_abundance)[1],
    rep = 1:dim(relative_abundance)[2],
    taxa = 1:dim(relative_abundance)[3]
  )

  df <- as.data.frame(as.table(relative_abundance))  
  
  if (value == "abund"){
  
  # Rename the columns
  colnames(df) <- c('timesteps', 'rep', 'taxa', 'rel_abund')
  }
  
  else if (value == "biomass"){
    colnames(df) <- c('timesteps', 'rep', 'taxa', 'rel_biomass')
  }
  
  else if (value == "S3.biomass"){
      stage_3_data <- data[, 3, , ] 

      total_abundance <- apply(stage_3_data, c(1, 2), sum)  #
      relative_abundance <- sweep(stage_3_data, c(1, 2), total_abundance, "/")
      
      dimnames(relative_abundance) <- list(
        timesteps = 1:dim(relative_abundance)[1],
        rep = 1:dim(relative_abundance)[2],
        taxa = 1:dim(relative_abundance)[3]
      )
      
      # Convert the array to a data frame
      df <- as.data.frame(as.table(relative_abundance))
      
      # Rename the columns
      colnames(df) <- c('timesteps', 'rep', 'taxa', 'S3.biomass')
      
  }
  
  # Convert 'timesteps', 'rep', and 'taxa' to numeric (they are factors by default)
  # Convert 'timesteps' and 'rep' to numeric
  df$timesteps <- as.numeric(as.character(df$timesteps))
  df$rep <- as.numeric(as.character(df$rep))
  
  # Define taxa names
  taxa_names <- c("HYOS", "BAET", "NZMS", "CHIR", "GAMM")
  
  # Assign these names to the 'taxa' column
  df$taxa <- factor(df$taxa, levels = 1:5, labels = taxa_names)

  if (is.null(burnin)== F){
    df <- df <- df[df$timesteps > burnin, ]
  }
  # Remove the 'rep' column
  df$rep <- NULL
  
  return(df)
}

partial.eta.data.frame <- function(data, burnin, iteration){
  
  pop_matrix_perm <- aperm(data, c(1, 3, 4, 2))  # Now: [2601, 1000, 5, 3]
  # Sum over stages directly with rowSums (preserving dimensions)
  taxon_stage_sums <- rowSums(pop_matrix_perm, dims = 3)
  
  # Total sum across all taxa
  total_stage_sums <- rowSums(taxon_stage_sums, dims = 2)
  
  # Compute relative abundance
  relative_abundance <- sweep(taxon_stage_sums, c(1, 2), total_stage_sums, "/")
  
  dimnames(relative_abundance) <- list(
    timesteps = 1:dim(relative_abundance)[1],
    rep = 1:dim(relative_abundance)[2],
    taxa = 1:dim(relative_abundance)[3]
  )
  
  df <- as.data.frame(as.table(relative_abundance))  

  colnames(df) <- c('timesteps', 'rep', 'taxa', 'rel_biomass')
  df$timesteps <- as.numeric(as.character(df$timesteps))
  
  if (is.null(burnin)== F){
    df <- df[df$timesteps > burnin, ]
  }
  
  # Convert 'timesteps', 'rep', and 'taxa' to numeric (they are factors by default)
  # Convert 'timesteps' and 'rep' to numeric
  df$timesteps <- as.numeric(as.character(df$timesteps))
  df <- df %>% 
    group_by(rep, taxa) %>% 
    summarise(across(where(is.numeric), ~mean(.x, na.rm = T))) %>% 
    select(timesteps, taxa, rel_biomass)
  # Define taxa names
  taxa_names <- c("HYOS", "BAET", "NZMS", "CHIR", "GAMM")
  
  # Assign these names to the 'taxa' column
  df$taxa <- factor(df$taxa, levels = 1:5, labels = taxa_names)
  

  # Remove the 'rep' column
  df$rep <- NULL
  
  return(df)
}


