##################
# Code dump
####################

# place to dump old code that I haven't used in a minute but might need

# we can also create a random flow scenario by sampling flows
#out_sample <- sample(out$Discharge,length(out$Discharge), replace=TRUE)
#Q <- out_sample

# another option is to keep flow the same each timestep (to test temp effects)
#out_same <- rep(10000, length(out$Discharge))
#Q <- out_same


# need to assign starting value
# we can start with 10 S1 individuals for each species
#for (sp in species){
#  output.N.list[[sp]][1,1] <- 10
#}