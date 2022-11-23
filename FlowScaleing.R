##############################
## Flow Re-Scaleing
##############################
library(scales)
library(geoR)
library(forecast)
# We are working with a variety of different flow scenarios in wildly different environments
# We need to create a relativized flow magnitude scale
# Flows will be scaled from 0 to 1000 (like McMullen model) but removed from m^3/s units
# Qmin will be =  5


# allows us to artificially input a max and min if we want to 
# 
flow.scale <- function(flow, min = NA, max = NA) {
  if (is.na(min)==F){
    newflow <- append(flow, min)
    resc <- rescale(newflow, to = c(1,1000))
  }
  if (is.na(max)==F){
    newflow <- append(flow, max)
    resc <- rescale(newflow, to = c(1,1000))
  }
  else { resc <- rescale(flow, to = c(1,1000))
  return(resc)
  }}

flow <- readNWISdv("09380000", "00060", "1985-10-01", "2021-09-30")
flows <- as.vector(flow$X_00060_00003)
?boxcox(
)

# unhinged option - rescale data to normal distribution, then artificially skew data to the right so majority of points are below 5 (Qmin) 

lambda <- BoxCox.lambda(flow$X_00060_00003)
fit <- BoxCox(flow$X_00060_00003,lambda)
hist(fit)
resfit <- rescale(fit, to = c(0,1))
hist(resfit)
skew <- 1/resfit
resckew <- rescale(skew, to = c(0,1000))
hist(resckew)


# doesn't like 0 values, can also use 2 parameter lambda from geoR package
#https://www.researchgate.net/post/How-to-perform-a-two-parameter-Box-Cox-transformation-to-a-data-series-in-R


oflow <- readNWISdv("09510200", "00060", "1982-10-01", "1987-09-30")

l1 <- boxcoxfit(oflow$X_00060_00003, lambda2= T)
fit1 <- BoxCox(oflow$X_00060_00003, l1)
hist(fit1)
resfit1 <- rescale(T.norm, to = c(0,1))
hist(resfit1)
skew1 <- 1/resfit
resckew1 <- rescale(skew, to = c(0,1000))
hist(resckew1)
