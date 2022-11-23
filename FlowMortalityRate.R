#####################################
## Flow Mortality Rate graph
###########################
library(tidyverse)
library(minpack.lm)
# follows exp(-h*Q[t-1])
h =  0.0075
Q <- seq(from = 1, to = 400, by = 2)
surv <- exp(-h*Q)

flow.surv.rate <- function(h, k, max, min, interval) {
  Q <- seq(min, max, by = interval)
  surv <- k*(exp(-h*Q))
  surv.df <- as.data.frame(cbind(Q, surv))
  return(surv.df)
}

surv.df.McMullen <- flow.surv.rate(0.02, 1, 100, 1, 1)

ggplot(surv.df.McMullen, aes(x = Q, y = surv))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival (%)') +
  xlab('McMullen Model Discharge')

# so we want to scale the max and min values to the McMullen model h values
# using a neg exponential eq y = k*exp(-h*x)
# can scale by creating data frame with min and max values we want and min and max survival for those values from McMullen

mcmullen.scale <- function(min, max, h){
  scale <- data.frame(seq(min, max, by = ((max-min)/99)), surv.df.McMullen[,2])
  colnames(scale) <- c("x", "y")
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = scale, start = c(k = 0, h = 0))
  return(nls.fit)
}

syc.fit <- mcmullen.scale(1, 400, 0.02)
surv.df.Sycamore <- flow.surv.rate(syc.fit$m$getPars()[2] , syc.fit$m$getPars()[1], 400, 1, (399/99))

ggplot(surv.df.Sycamore, aes(x = Q, y = surv))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Mortality (%)') +
  xlab('Sycamore Creek Discharge (cfs)')

col.fit <- mcmullen.scale(10000, 100000, 0.02)
surv.df.Colorado <- flow.surv.rate(1, 0.02, 100000, 10000, (90000/99))

ggplot(surv.df.Colorado, aes(x = Q, y = surv))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Mortality (%)') +
  xlab('Sycamore Creek Discharge (cfs)')

