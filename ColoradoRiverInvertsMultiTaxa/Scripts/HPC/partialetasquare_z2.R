# ## Multispecies Partial Eta Square
library(data.table)
library(purrr)
library(effectsize)
library(dplyr)
library(parallel)

# zi = 1 
temp2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z2.rds")
temp_spike2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Spike_Z2.rds")
hyd2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Z2.rds")
hyd_spike2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_Spike_Z2.rds")
hfe2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Z2.rds")
hfe_spike2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_HFE_Spike_Z2.rds")
hyd_hfe2 <- readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Hyd_HFE_Z2.rds")
hyd_hfe_spike2 <-readRDS("/home/kurthena/ColoradoRiverInverts/Multispp_Temp_Z2.rds")


partialeta <- bind_rows(
    temp2  %>%  mutate(source = "Temperature", temp_spike = 0, hydropeaking = 0, HFE = 0),
    temp_spike2 %>% mutate(source = "Temperature", temp_spike = 1, hydropeaking = 0, HFE = 0),
    hyd2 %>%  mutate(source = "Temperature & hydropeaking", temp_spike = 0, hydropeaking = 1, HFE = 0),
    hyd_spike2 %>% mutate(source = "Temperature & spike & hydropeaking", temp_spike = 1, hydropeaking = 1, HFE = 0),
    hfe2 %>% mutate(source = "Temperature & HFE",  temp_spike = 0, hydropeaking = 0, HFE = 1),
    hfe_spike2 %>% mutate(source = "Temperature & spike & HFE",  temp_spike = 1, hydropeaking = 0, HFE = 1),
    hyd_hfe2 %>% mutate(source = "Temperature & hydropeaking & HFE",  temp_spike = 0, hydropeaking = 1, HFE = 1),
    hyd_hfe_spike2 %>% mutate(source = "Temperature & spike & hydropeaking & HFE",  temp_spike = 1, hydropeaking = 1, HFE = 1)
  )
  
rm(temp2)
rm(temp_spike2)
rm(hyd2)
rm(hyd_spike2)
rm(hfe2)
rm(hfe_spike2)
rm(hyd_hfe2)
rm(hyd_hfe_spike2)

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

write.csv(eta_sq_taxon_biomass, "Multispp_partialetasquare_biomass_z2.csv")
