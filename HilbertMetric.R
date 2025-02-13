#######################
## Testing the weak ergodic theorem
##########################

# we can use Hilberts metric (Caswell p370)
# Hilberts Distance d(x, y) = ln(maxi(ni/mi)/mini(ni/mi))

source("A_1sp_Model.R")
source("B_1sp_Model.R")
#different starting values
source("1spFunctions.R")
Time <- c(1:36500)
Date <- rep(1:365, times = 100)
Day <- seq(as.Date("2022-01-01"), as.Date("2121-12-31"), by="days")
# find and remove leap days
find_leap = function(x){
  day(x) == 29 & month(x) == 2 
}
Day <- Day[which(find_leap(Day) == F)]


Temperature <-13.95

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]
discharge <- rep(0.1, times = length(temp$dts))

# output1 <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))
# out <- as.vector(rowSums(output1)/2)
# out <- mean.data.frame(output1, burnin = 1, iteration = 2)
# identical(output1[,,1], output1[,,2])
# identical(tail(output1[,,1]), tail(output2[,,2]))
# # do not coverge on set values
# # that is, they are sensitive to initial conditions
# 
# 
# # testing H distance
# 
# n1 <- output1[,1,1]
# n3 <- output1[,3,1]
# m1 <- output1[,1,2]
# m3 <- output1[,3,2]
# 
# 
# log((1/.01)/(.01/1))
# log((n1/m1)/(n3/m3))
# 
# # does not converge 
# plot(log((n1/n3)/(m1/m3)))


# 
install.packages("fNonlinear")
library(fNonlinear)
install.packages("nonlinearTseries")
install.packages("tseriesChaos")
library(nonlinearTseries)
library(tseriesChaos)

le <- vector()
fecs <- seq(1200, 100000, by = 2000)
for (fec in 1:length(fecs)){
out1 <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[fec])
out2 <- mean.data.frame(out1, 250, 2)
out3 <- as.vector(out2$mean.abund)
lyapun <- lyap_k(out3, m=1, d=2, s=200, t=40, ref=1700, k=2, eps=4)
le[fec] <- max(lyapun)
}

# same results as chaos 0-1 test 

# 
# 
# 
# lyap(lyapun, -0.3131367,9.180819 )
# maxLyapunov(time.series=out$mean.abund,
#               min.embedding.dim=1,
#               max.embedding.dim=3,
#               time.lag=1,
#               radius=0.001,theiler.window=4,
#               min.neighs=2,min.ref.points=500,
#                max.time.steps=40,do.plot=T)
# 
# 
# 
# class(LE) <- "lyapunov"
# return(LE)
# install.packages("DChaos")
# library(DChaos)
# 
# le <- vector()
# fecs <- seq(1200, 100000, by = 2000)
# 
# for (fec in 1:length(fecs)){
#   out1 <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[fec])
#   out2 <- mean.data.frame(out1, 250, 2)
#   out3 <- as.vector(out2$mean.abund)
#   
#   # Check if out3 has enough data points
#   if (length(out3) < 100) {
#     warning(paste("Not enough data points for fecundity", fecs[fec], "Skipping..."))
#     next
#   }
#   
#   # Check for NA or infinite values
#   if (any(is.na(out3)) || any(is.infinite(out3))) {
#     warning(paste("NA or infinite values in out3 for fecundity", fecs[fec], "Skipping..."))
#     next
#   }
#   
#   # Compute Lyapunov exponent
#   exponent <- tryCatch({
#     DChaos::lyapunov(out3, timelapse="FIXED", B=100, doplot=T)
#   }, error = function(e) {
#     warning(paste("Error computing Lyapunov exponent for fecundity", fecs[fec], ":", e$message))
#     return(NULL)
#   })
#   
#   if (!is.null(exponent)) {
#     le[fec] <- exponent$exponent.median[1]
#   } else {
#     le[fec] <- NA  # Assign NA if there was an error
#   }
# }
# 
# plot(fecs, le)
# install.packages("Chaos01")
# library(Chaos01)

chaos1 <- vector()
fecs <- seq(1200, 100000, by = 2000)
for (fec in 1:length(fecs)){
  out <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = fecs[fec])
  out <- mean.data.frame(out, 250, 2)
  chaos <- testChaos01(out$mean.abund)
  chaos1[fec] <- chaos
}

chaostestdf <- as.data.frame(cbind(fecs, chaos1))

outchaos <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = 5200 )
outchaos2 <- mean.data.frame(outchaos, 250, 2)
length(unique(outchaos2$mean.abund))
#chaosn <- plot(outchaos2$mean.abund[1:2359], outchaos2$mean.abund[2:2360], type = "b")
chaosdf <- as.data.frame(cbind(outchaos2$mean.abund[1:2359], outchaos2$mean.abund[2:2360]))


outstable <- Bmodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature), fecundity = 1200)
outstable2 <- mean.data.frame(outstable, 250, 2)
length(unique(outstable2$mean.abund))
#stablen <- plot(outstable2$mean.abund[1:2359], outstable2$mean.abund[2:2360], type = "b")
stabledf <- as.data.frame(cbind(outstable2$mean.abund[1:2359], outstable2$mean.abund[2:2360]))


