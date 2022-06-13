###################################################################
## Stage Structured Model for Colorado River below Glen Canyon Dam
###################################################################
# By Angelika Kurthen
# Partially adapted from Rogosch et al., 2019)

#---------------------------------------------------------------------
source("EggMortalityIntegration.R") # this code includes 'h' (hydropeaking factor for invertes)

# Packages needed
install.packages("tidyverse")
install.packages("purrr")
library(plyr)
library(purrr)
library(tidyverse)

# in one year, each fortnight is one timestep
# we will assign the first timestep later, but need a sequence with 27 timesteps

time = seq(2,27, 1) 

species <- c("BAET","CHIRO", "SIMU")

iterations = 1  # can be changed
# set k (will be added later)

#k = 10000


# N ---------------------------------------------------------------------------
# # Output of biomass and no. ind. for each age class for each year projected  
## Array w/ 3 cols (stage classes) and however many rows there are yrs projected
## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated (adapted from Rogosch et al. 2019)

output.N.array <- array(0, dim = c(length(time) + 1, length(species)))

## Repeating the array 7 times 
output.N.list <- rep(list(output.N.array), length(species))
## Assigning names to each array from sppnames vector
names(output.N.list) <- species

# Total. N ---------------------------------------------------------------
## Takes all stages 2 and 3 and sums them for each year and  (adpated from Rogosch et al. 2019
Total.N <- array(0,
                 dim  <-c((length(time) +1 ), iterations),
                 dimnames <- list(1:(length(time) + 1), 1:iterations)
)

## Creating a list of 7 arrays to fill in. One for each spp. 
## Create an array to be repeated
reparray <- array(0,## replist. List of arrays w/ abundance data for each spp ----------------------

                  dim = c(length(time) + 1, 3, iterations),
                  dimnames = list(1:27, c("S1", "S2", "S3"), 1:iterations)
)

## Repeating the array 7 times 
replist <- rep(list(reparray), 3)
names(replist) <- species

# need to assign starting value
# in the future, we can pull these #s from a randomly selected date in the Colorado River data
# for now, will start with 100 S1 individuals for each species
for (sp in species){
  output.N.list[[sp]][1,1] <- 100
}

## Iteration Loop ------------------------------------------------------------------------------
# need to loop through all the iterations

for (iter in iteration) {

## Time Loop -----------------------------------------------------------------------------------
# loop through each timestep (2 week period)

for (t in time){
  
  # pull fecundity from normal distribution of egg amounts, multiply by hydropeaking modifier
  F_BAET = rnorm(1, mean = 1104.5, sd = 42.75) * H_BAET #Baetidae egg minima and maxima from Degrange, 1960
  F_CHIRO = rnorm(1, mean = 1675, sd = 522.5) * H_CHIRO # Chironomidae egg minima, and maxima from Thienemann 1954, Jones 1968, and Risbec 1951
  F_SIMU = rnorm(1, mean = 350, sd = 75) * H_SIMU # Simulliidae egg minima, and maxima from Hinton 1981
  
  # assign probablities for transitioning into next stage
  # at present, these are all random
    # Baetidae 
    P1_BAET <- 0.3
    P2_BAET <- 0.3
    P3_BAET <- 0.3 # P3 is probability adult will remain in that stage and not die
  
    # Chironomidae
    P1_CHIRO <- 0.476
    P2_CHIRO <- 0.476
    P3_CHIRO <- 0.476
  
    #Simulliidae
    P1_SIMU <- 0.2
    P2_SIMU <- 0.2
    P3_SIMU <- 0.2
    
  # assign probabilities for staying in same stage
    # at present, these are all random
    # Baetidae
    G1_BAET <- 0.3
    G2_BAET <- 0.3
    
    # Chironomidae
    G1_CHIRO <- 0.3
    G2_CHIRO <- 0.3
    
    # Simulliidae
    G1_SIMU <- 0.2
    G2_SIMU <- 0.2
    
  # create Leslie/Lefkovitch matrix for each species
    
    #Baetidae
    BAET1 <- c(P1_BAET, 0, F_BAET)
    BAET2 <- c(G1_BAET, P2_BAET, 0)
    BAET3 <- c(0, G2_BAET, P3_BAET)
    
    #Chironomidae
    CHIRO1 <- c(P1_CHIRO, 0, F_CHIRO)
    CHIRO2 <- c(G1_CHIRO, P2_CHIRO, 0)
    CHIRO3 <- c(0, G2_CHIRO, P3_CHIRO)
    
    #Simulliidae
    SIMU1 <- c(P1_SIMU, 0, F_SIMU)
    SIMU2 <- c(G1_SIMU, P2_SIMU, 0)
    SIMU3 <- c(0, G2_SIMU, P3_SIMU)
    
    for (sp in species) {
        assign(paste0('A', sp), rbind(
          get(paste0(sp, '1')),
          get(paste0(sp, '2')),
          get(paste0(sp, '3'))
        ))
    }
    # multiply matrix by N
    # sort N values into each df/array in the list       
    for(sp in species) {
      output.N.list[[sp]][t,1:3] <- output.N.list[[sp]][t-1, 1:3] %*% get(paste0('A', sp))
       }
   
    #include in iteration wide list so we can converge on average at the end
    for(sp in species) {
      replist[[sp]][,,iter] <- output.N.list[[sp]]    
    }
}
    
    ## Total.N -----------------------------------------------------------------
    ## Caculating Total.N for each year, and adding it to total.N data frame
    ## includes all individals from all stages and species

    Total.N[,iter] <- output.N.list %>%
      do.call('cbind', .) %>%
      apply(1, sum)
    
    ## turning replist into a df
    repdf <- ldply(replist, function(x) {
      adply(x, c(1,2,3))
    })
    
    names(repdf) <- c('Species', 'Timestep', 'Stage', 'Rep', 'Abund')
    repdf$Timestep <- as.numeric(as.character(repdf$Timestep))
    
    totn <- adply(Total.N, c(1,2))
    names(totn) <- c('Timestep', 'Rep', 'Tot.abund')
    totn$Timestep <- as.numeric(as.character(totn$Timestep))
    
    ## joining totn and repdf together
    repdf <- left_join(totn, repdf)
    
    ## calculating relative abundance
    repdf <- mutate(repdf, rel.abund = Abund/Tot.abund)

}

## Plotting  
ggplot(repdf, aes(Timestep, rel.abund, colour = Species, linetype= Stage))+
  geom_line()+
  facet_wrap(~Species, ncol = 2) 
theme_classic_facet() +
  coord_cartesian(ylim = c(0,1)) +
  ylab('Relative abundance') +
  xlab('Timestep')


  


#Next steps
# how does K interact with species? no reproduction? Mass mortality at beginning? 
#running multiple iterations, getting means, getting an average
#temperature and flow fluctuations
     
