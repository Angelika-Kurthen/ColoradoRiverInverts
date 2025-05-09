# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# zi = 1 
temp8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z8.rds")
temp_spike8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Spike_Z8.rds")
hyd8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Z8.rds")
hyd_spike8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Spike_Z8.rds")
hfe8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Z8.rds")
hfe_spike8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Spike_Z8.rds")
hyd_hfe8 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_HFE_Z8.rds")
hyd_hfe_spike8 <-readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z8.rds")


partialeta <- bind_rows(
    temp8  %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    temp_spike8 %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
    hyd8 %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    hyd_spike8 %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    hfe8 %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
    hfe_spike8 %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
    hyd_hfe8 %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
    hyd_hfe_spike8 %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
  )
  
rm(temp8)
rm(temp_spike8)
rm(hyd8)
rm(hyd_spike8)
rm(hfe8)
rm(hfe_spike8)
rm(hyd_hfe8)
rm(hyd_hfe_spike8)

label.names <- c(
  "Temperature", 
  "Temperature +\nSpike", 
  "Temperature +\nHydropeaking",
  "Temperature +\nSpike + Hydropeaking",
  "Temperature +\nHFE",
  "Temperature +\nSpike + HFE",
  "Temperature +\nHydropeaking + HFE",
  "Temperature +\nSpike + Hydropeaking + HFE"
)

partialeta$source <- ordered(partialeta$source, levels = c(
  "Temperature", 
  "Temperature & spike",
  "Temperature & hydropeaking",
  "Temperature & spike & hydropeaking",
  "Temperature & HFE",
  "Temperature & spike & HFE",
  "Temperature & hydropeaking & HFE",
  "Temperature & spike & hydropeaking & HFE"
))  

# Function to fit ANOVA and extract partial eta squared
compute_eta_sq <- function(df, response_var) {
  model <- aov(as.formula(paste(response_var, "~ (temperature + temp_spike + hydropeaking + HFE)^4")), data = df)
  eta_squared(model, partial = TRUE) %>%
    as.data.frame()
}

# Run analysis per taxon for biomass
eta_sq_taxon_biomass <- partialeta %>%
  group_by(taxa) %>%
  group_split() %>%
  map_df(~ compute_eta_sq(.x, "rel_biomass"), .id = "taxon")

write.csv(eta_sq_taxon_biomass, "Multispp_partialetasquare_biomass_z8.csv")
