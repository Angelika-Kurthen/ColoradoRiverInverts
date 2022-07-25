#########################################
# Temperature Theshold Equation Graphing
#########################################

# make a sequence of temperatures
temp <- seq(7, 16, by = 0.05)

# empty lists for the probabilities to go into
devlist <- vector()
growthlist <- vector()

# cycle through all the temperatures and calculate the different probability values
for(t in temp){
  # development measures (basically, if below 10C, no development, if between 10 and 12, follows a function, if above 12, prob of transition to next stage is 0.6395)
  if (10 > t) d <- 0
  if (t > 13) d <- 0.6  
  if (10 <= t & t <= 13) d <- (0.2 * t) -2
  
  devlist <- append(devlist, d)

  # growth (if below 10C, no growth can occur - everything basically freezes, if between 10 and 11, prob of remaining in same stage = 0.6395, if above 13, prob of transition to next stage is 0 )
  if (10 >= t) g <- 0.6
  if (t > 13) g <- 0
  if (10 < t & t <=  13) g <- (-0.2 * t) + 2.6
  
  growthlist <- append(growthlist, g)
}

# create dataframe to plot with
tempdf <- as.data.frame(cbind(temp, devlist, growthlist))

# reset graphing parameters (in case something was messed up)
par(mfrow = c(1,1))

#plot
plot(temp, growthlist, type = "l", col = "red", ylab = "Matrix Probability Value", xlab = "Temperature (C)")
lines(temp, devlist, type = "l", col = "blue")  
legend(13, 0.3, legend=c("Growth (remain)", "Development (transition)"),
       col=c("Red", "Blue"), lty=1, cex=0.8)


     