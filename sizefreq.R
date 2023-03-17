####################################
## Size Freqency Plot Function - adapted from Mark Novak
## Angelika Kurthen
#####################################

# need to use devtools to load github repo
install.packages("devtools")
library(devtools)
# jeff's github repo with USGS data
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
library(tidyverse)
#library(ggplot2)

# use for mark's code
library(gam)
library(lubridate)


# pull drift data from Jeff's repo
drift <- readDB(gear = "Drift", type = "Sample", updater = T)

species <- c("NZMS","GAMM", "QUAG", "CHIA", "CHIP", "CHIL", "SIMA", "SIMP", "SIML", "HYDE", "BAEL")
locations <- c("upper", "middle", "lower")


for (sp in 1:length(species)) {
  for (loc in 1:length(locations)){
    a <- TaxaLocationData(species[sp], locations[loc])
    nam <- paste0(species[sp], locations[loc])
    assign(nam, a)
    write.csv(nam, 
              file = 
                paste0("ColoradoRiverInverts/CleanData/", paste0(species[sp], locations[loc])))
    
  }
}
NZMSupper <- TaxaLocationData(spp = "NZMS", loc = "upper")
NZMSmiddle <- TaxaLocationData(spp = "NZMS", loc = "middle")
NZMSlower <- TaxaLocationData(spp = "NZMS", loc = "lower")





















