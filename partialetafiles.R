# code to produce files for compulation for partial eta square

library(readr)
library(dplyr)
library(stringr)
library(fs)

# Define file list and metadata
# files <- list(
#   "Multispp_temp_biomass_z_hyos_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
#   "Multispp_temp_biomass_hyd_z_hyos_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_hyd_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
#   "Multispp_temp_biomass_hfe_z_hyos_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_z_hyos_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
#   "Multispp_temp_hyd_HFE_biomass_spike_z_hyos_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
#   )
# 
# # Loop over each file
# for (fname in names(files)) {
#   cat("Processing file:", fname, "\n")
#   meta <- files[[fname]]
#   base_filename <- tools::file_path_sans_ext(basename(fname))
#   
#   read_csv_chunked(
#     file = fname,
#     callback = SideEffectChunkCallback$new(function(data, pos) {
#       # Filter only q values of interest
#       data <- data %>%
#         filter(q %in% c(1, 2, 4, 8))
#       
#       # Add metadata
#       data <- data %>%
#         mutate(
#           source = meta$source,
#           temp_spike = meta$temp_spike,
#           hydropeaking = meta$hydropeaking,
#           HFE = meta$HFE
#         )
#       
#       # Write a file per taxon
#       unique_taxa <- unique(data$taxa)
#       for (taxon in unique_taxa) {
#         taxon_data <- data %>% filter(taxa == taxon)
#         taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
#         out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
#         
#         write_csv(taxon_data, out_fname, append = file.exists(out_fname))
#       }
#     }),
#     chunk_size = 100000,
#     col_types = cols(
#       timesteps = col_double(),
#       taxa = col_character(),
#       biomass = col_double(),
#       temperature = col_double(),
#       q = col_double()
#     )
#   )
#   
#   # Remove the original large file after processing
#   file_delete(fname)
#   cat("Deleted original file:", fname, "\n\n")
# }

# Define file list and metadata
files <- list(
    "Multispp_temp_biomass_z_baet_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_hyd_z_hyos_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    "Multispp_temp_hyd_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    "Multispp_temp_biomass_hfe_z_baet_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
    "Multispp_temp_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_z_baet_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_spike_z_baet_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
)

# Loop over each file
for (fname in names(files)) {
  cat("Processing file:", fname, "\n")
  meta <- files[[fname]]
  base_filename <- tools::file_path_sans_ext(basename(fname))
  
  read_csv_chunked(
    file = fname,
    callback = SideEffectChunkCallback$new(function(data, pos) {
      # Filter only q values of interest
      data <- data %>%
        filter(q %in% c(1, 2, 4, 8))
      
      # Add metadata
      data <- data %>%
        mutate(
          source = meta$source,
          temp_spike = meta$temp_spike,
          hydropeaking = meta$hydropeaking,
          HFE = meta$HFE
        )
      
      # Write a file per taxon
      unique_taxa <- unique(data$taxa)
      for (taxon in unique_taxa) {
        taxon_data <- data %>% filter(taxa == taxon)
        taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
        out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
        
        write_csv(taxon_data, out_fname, append = file.exists(out_fname))
      }
    }),
    chunk_size = 100000,
    col_types = cols(
      timesteps = col_double(),
      taxa = col_character(),
      biomass = col_double(),
      temperature = col_double(),
      q = col_double()
    )
  )
  
  # Remove the original large file after processing
  file_delete(fname)
  cat("Deleted original file:", fname, "\n\n")
}

# Define file list and metadata
files <- list(
    "Multispp_temp_biomass_z_chir_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_hyd_z_chir_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    "Multispp_temp_hyd_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    "Multispp_temp_biomass_hfe_z_chir_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
    "Multispp_temp_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_z_chir_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_spike_z_chir_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
  
)

# Loop over each file
for (fname in names(files)) {
  cat("Processing file:", fname, "\n")
  meta <- files[[fname]]
  base_filename <- tools::file_path_sans_ext(basename(fname))
  
  read_csv_chunked(
    file = fname,
    callback = SideEffectChunkCallback$new(function(data, pos) {
      # Filter only q values of interest
      data <- data %>%
        filter(q %in% c(1, 2, 4, 8))
      
      # Add metadata
      data <- data %>%
        mutate(
          source = meta$source,
          temp_spike = meta$temp_spike,
          hydropeaking = meta$hydropeaking,
          HFE = meta$HFE
        )
      
      # Write a file per taxon
      unique_taxa <- unique(data$taxa)
      for (taxon in unique_taxa) {
        taxon_data <- data %>% filter(taxa == taxon)
        taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
        out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
        
        write_csv(taxon_data, out_fname, append = file.exists(out_fname))
      }
    }),
    chunk_size = 100000,
    col_types = cols(
      timesteps = col_double(),
      taxa = col_character(),
      biomass = col_double(),
      temperature = col_double(),
      q = col_double()
    )
  )
  
  # Remove the original large file after processing
  file_delete(fname)
  cat("Deleted original file:", fname, "\n\n")
}

# Define file list and metadata
files <- list(
    "Multispp_temp_biomass_z_gamm_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_hyd_z_gamm_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    "Multispp_temp_hyd_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    "Multispp_temp_biomass_hfe_z_gamm_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
    "Multispp_temp_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_z_gamm_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_spike_z_gamm_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
  
)

# Loop over each file
for (fname in names(files)) {
  cat("Processing file:", fname, "\n")
  meta <- files[[fname]]
  base_filename <- tools::file_path_sans_ext(basename(fname))
  
  read_csv_chunked(
    file = fname,
    callback = SideEffectChunkCallback$new(function(data, pos) {
      # Filter only q values of interest
      data <- data %>%
        filter(q %in% c(1, 2, 4, 8))
      
      # Add metadata
      data <- data %>%
        mutate(
          source = meta$source,
          temp_spike = meta$temp_spike,
          hydropeaking = meta$hydropeaking,
          HFE = meta$HFE
        )
      
      # Write a file per taxon
      unique_taxa <- unique(data$taxa)
      for (taxon in unique_taxa) {
        taxon_data <- data %>% filter(taxa == taxon)
        taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
        out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
        
        write_csv(taxon_data, out_fname, append = file.exists(out_fname))
      }
    }),
    chunk_size = 100000,
    col_types = cols(
      timesteps = col_double(),
      taxa = col_character(),
      biomass = col_double(),
      temperature = col_double(),
      q = col_double()
    )
  )
  
  # Remove the original large file after processing
  file_delete(fname)
  cat("Deleted original file:", fname, "\n\n")
}

# Define file list and metadata
files <- list(
    "Multispp_temp_biomass_z_nzms_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_hyd_z_nzms_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    "Multispp_temp_hyd_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    "Multispp_temp_biomass_hfe_z_nzms_full.csv" = list(source = "Temperature & HFE", temp_spike = 0, hydropeaking = 0, HFE = 1),
    "Multispp_temp_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_z_nzms_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_spike_z_nzms_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
  
)

# Loop over each file
for (fname in names(files)) {
  cat("Processing file:", fname, "\n")
  meta <- files[[fname]]
  base_filename <- tools::file_path_sans_ext(basename(fname))
  
  read_csv_chunked(
    file = fname,
    callback = SideEffectChunkCallback$new(function(data, pos) {
      # Filter only q values of interest
      data <- data %>%
        filter(q %in% c(1, 2, 4, 8))
      
      # Add metadata
      data <- data %>%
        mutate(
          source = meta$source,
          temp_spike = meta$temp_spike,
          hydropeaking = meta$hydropeaking,
          HFE = meta$HFE
        )
      
      # Write a file per taxon
      unique_taxa <- unique(data$taxa)
      for (taxon in unique_taxa) {
        taxon_data <- data %>% filter(taxa == taxon)
        taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
        out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
        
        write_csv(taxon_data, out_fname, append = file.exists(out_fname))
      }
    }),
    chunk_size = 100000,
    col_types = cols(
      timesteps = col_double(),
      taxa = col_character(),
      biomass = col_double(),
      temperature = col_double(),
      q = col_double()
    )
  )
  
  # Remove the original large file after processing
  file_delete(fname)
  cat("Deleted original file:", fname, "\n\n")
}

# Define file list and metadata
files <- list(
    "Multispp_temp_biomass_z_full.csv" = list(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
     "Multispp_temp_biomass_spike_z_full.csv" = list(source = "Temperature & spike", temp_spike = 1, hydropeaking = 0, HFE = 0),
    "Multispp_temp_biomass_hyd_z_full.csv" = list(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    "Multispp_temp_hyd_biomass_spike_z_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    "	Multispp_temp_biomass_hfe_z_full.csv" = list(source = "Temperature & spike & hydropeaking", temp_spike = 0, hydropeaking = 0, HFE = 1),
    "Multispp_temp_HFE_biomass_spike_z_full.csv" = list(source = "Temperature & spike & HFE", temp_spike = 1, hydropeaking = 0, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_z_full.csv" = list(source = "Temperature & hydropeaking & HFE", temp_spike = 0, hydropeaking = 1, HFE = 1),
    "Multispp_temp_hyd_HFE_biomass_spike_z_full.csv" = list(source = "Temperature & spike & hydropeaking & HFE", temp_spike = 1, hydropeaking = 1, HFE = 1)
  
)

# Loop over each file
for (fname in names(files)) {
  cat("Processing file:", fname, "\n")
  meta <- files[[fname]]
  base_filename <- tools::file_path_sans_ext(basename(fname))
  
  read_csv_chunked(
    file = fname,
    callback = SideEffectChunkCallback$new(function(data, pos) {
      # Filter only q values of interest
      data <- data %>%
        filter(q %in% c(1, 2, 4, 8))
      
      # Add metadata
      data <- data %>%
        mutate(
          source = meta$source,
          temp_spike = meta$temp_spike,
          hydropeaking = meta$hydropeaking,
          HFE = meta$HFE
        )
      
      # Write a file per taxon
      unique_taxa <- unique(data$taxa)
      for (taxon in unique_taxa) {
        taxon_data <- data %>% filter(taxa == taxon)
        taxon_clean <- str_replace_all(taxon, "[^A-Za-z0-9_]", "_")
        out_fname <- paste0(base_filename, "_", taxon_clean, ".csv")
        
        write_csv(taxon_data, out_fname, append = file.exists(out_fname))
      }
    }),
    chunk_size = 100000,
    col_types = cols(
      timesteps = col_double(),
      taxa = col_character(),
      biomass = col_double(),
      temperature = col_double(),
      q = col_double()
    )
  )
  
  # Remove the original large file after processing
  file_delete(fname)
  cat("Deleted original file:", fname, "\n\n")
}





