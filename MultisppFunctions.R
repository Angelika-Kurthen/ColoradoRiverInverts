###############
#Multispp Functions
###############
multispp.data.frame <- function(data, burnin, iteration, value = "biomass"){
  
  repdf <- plyr::adply(data, c(1,2,3,4))
  #repdf <- plyr::adply(data, c(1,2))
  names(repdf) <- c('timesteps', 'stage', 'rep', "taxa", value)
  repdf$timesteps <- as.numeric(as.character(repdf$timesteps))
  
  if (value == "biomass"){
  # joining totn and repdf together
  # repdf <- dplyr::left_join(totn, repdf)
  tot.biomass <- rowSums(data)/iteration
  ## calculating relative abundance
  repdf <- mutate(repdf, rel.biomass = biomass/tot.biomass)
  repdf$timesteps <- as.factor(repdf$timesteps)
  ## Taking mean results to cf w/ observed data'
  
  means.list<- repdf %>%
    # select(-tot.abund) %>%
    dplyr::group_by(timesteps, rep, taxa) %>% # combining stages and taxa
    dplyr::summarise(rel.biomass = sum(rel.biomass)) %>%
    ungroup() %>%
    dplyr::group_by(timesteps, taxa) %>%
    dplyr::summarise(mean.rel.biomass = mean(rel.biomass),
                     sd.abund = sd(rel.biomass),
                     se.abund = sd(rel.biomass)/sqrt(iteration)) %>%
    ungroup()
  }
  
  if (value == "abund"){
    # repdf <- dplyr::left_join(totn, repdf)
    tot.n <- rowSums(data)/iteration
    ## calculating relative abundance
    repdf <- mutate(repdf, rel.n = abund/tot.n)
    repdf$timesteps <- as.factor(repdf$timesteps)
    ## Taking mean results to cf w/ observed data'
    
    means.list<- repdf %>%
      # select(-tot.abund) %>%
      dplyr::group_by(timesteps, rep, taxa) %>% # combining stages and taxa
      dplyr::summarise(rel.n = sum(rel.n)) %>%
      ungroup() %>%
      dplyr::group_by(timesteps, taxa) %>%
      dplyr::summarise(mean.rel.abund = mean(rel.n),
                       sd.abund = sd(rel.n),
                       se.abund = sd(rel.n)/sqrt(iteration)) %>%
      ungroup()
  }
    if (value == "S3.biomass"){

      ## subset stage 3s
      repdf3 <- subset(repdf, stage == "S3")
      repdf3$timesteps <- as.factor(repdf3$timesteps)
      ## Taking mean results to cf w/ observed data'
      
      means.list<- repdf3 %>%
        dplyr::group_by(timesteps, taxa) %>% # combining taxa, reps
        dplyr::summarise(mean.S3.biomass = mean(S3.biomass),
                         sd.abund = sd(S3.biomass),
                         se.abund = sd(S3.biomass)/sqrt(iteration)) %>%
        ungroup()
      
    }
  
  if (is.null(burnin)== F){
    means.list <- means.list[burnin:length(means.list$timesteps), ]
  }
  return(means.list)
}