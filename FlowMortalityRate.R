######################################
## Flow Mortality Rate scale and graph
#######################################
library(tidyverse)
library(minpack.lm)
# follows exp(-h*Q[t-1])


flow.surv.rate <- function(h, k, max, min, interval) {
  Q <- seq(min, max, by = interval)
  surv <- k*(exp(-h*Q))
  surv.df <- as.data.frame(cbind(Q, surv))
  return(surv.df)
}


# original McMullen model runs with discharges from 1 - 1000 m/s with 100% mortality after around 250 m/s
# but that seems very steep, considering that only about 96% of individuals are wiped out post max disturbance (in Sycamore Creek)
# going to run it between 1 and 100 (and scale to that because at 100 m/s about 86% die)

surv.df.McMullen <- flow.surv.rate(0.02, 1, 1000, 1, 1) 

ggplot(surv.df.McMullen, aes(x = Q, y = surv))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  theme_bw()+
  ylab('Immediate Post-Disturbance Survival (%)') +
  xlab('McMullen Model Discharge (m^3/s)')

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
  ylab('Immediate Post-Disturbance Survival (%)') +
  theme_bw()+
  xlab('Sycamore Creek Discharge (cfs)')

col.fit <- mcmullen.scale(10000, 100000, 0.02)
surv.df.Colorado <- flow.surv.rate(1, 0.02, 100000, 10000, (90000/99))

ggplot(surv.df.Colorado, aes(x = Q, y = surv))+
  geom_line()+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival (%)') +
  theme_bw()+
  xlab('Colorado River Discharge (cfs)')

