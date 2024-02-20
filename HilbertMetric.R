#######################
## Testing the weak ergodic theorem
##########################

# we can use Hilberts metric (Caswell p370)
# Hilberts Distance d(x, y) = ln(maxi(ni/mi)/mini(ni/mi))

source("A_1sp_Model.R")

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


Temperature <-  -7.374528  * (cos(((2*pi)/365)*Date))  +  (-1.649263* sin(2*pi/(365)*Date)) + 13.956243

temp <- as.data.frame(cbind(Time, Day, Temperature))
temp$Day <- as.Date(temp$Day, origin= "1970-01-01")
colnames(temp) <- c("Time", "Date", "Temperature")
temp <- TimestepTemperature(temp)
temp <- temp[c(1,3)]
discharge <- rep(0.1, times = length(temp$dts))

output1 <- Amodel(discharge, temp, baselineK = 10000, disturbanceK = 40000, Qmin = 0.25, extinct = 50, iteration = 2, peaklist = 0, peakeach = length(temp$Temperature))

identical(output1[,,1], output1[,,2])
identical(tail(output1[,,1]), tail(output2[,,2]))
# do not coverge on set values
# that is, they are sensitive to initial conditions


# testing H distance

n1 <- output1[,1,1]
n3 <- output1[,3,1]
m1 <- output1[,1,2]
m3 <- output1[,3,2]


log((1/.01)/(.01/1))
log((n1/m1)/(n3/m3))

# does not converge 
plot(log((n1/n3)/(m1/m3)))
