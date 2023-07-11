##############################
# Code to calculate negative exponential fit of survivorships
##############################


library(readxl)
library(minpack.lm)
library(tidyverse)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  #df <- rbind(df, c(Qmin, 1.11))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function
  return(nls.fit)
}

mag <- c(0.5, 1, 3)
mort <- c(0.85, 0.999, 0.9999999)
high <- flow.surv.fit(mag, mort, 0.25)

mag <- c(0.75, 1, 3)
mort <- c(0.6, 0.65, 75)
med <- flow.surv.fit(mag, mort, 0.25)

mag <- c(0.75, 1, 2, 3)
mort <- c(0.4, 0.45, 0.55, 0.6)
low <- flow.surv.fit(mag, mort, 0.25)