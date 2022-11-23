##############################
## Flow Re-Scaleing
##############################
library(scales)
library(geoR)
library(forecast)
library(dataRetrieval)
library(MASS)
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

# unhinged option - rescale data to normal distribution, then artificially skew data to the right so majority of points are below 5 (Qmin) 
# OR just have different Q min standardss? 

flow.scale.distribution <- function(flow){ #requires geoR package
  fit <- boxcoxfit(flow, lambda2 = T)
  lambda = fit$lambda[1]
  lambda2 = fit$lambda[2]
  if(lambda==0){norm = log(flow + lambda2)}
  if(lambda!=0){norm = ((flow + lambda2) ^ lambda - 1) / lambda} # create normal distribution
  resc <- rescale(norm, to = c(0, 1000))
  return(resc)
}

# doesn't like 0 values, can also use 2 parameter lambda from geoR package
#https://www.researchgate.net/post/How-to-perform-a-two-parameter-Box-Cox-transformation-to-a-data-series-in-R
