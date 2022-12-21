#############################################################
## Code to fit curve to Baetidae spp mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
BAETVitalRates <- read_excel("VitalRates.xlsx", sheet = "Baetid Mortality Rates")
BAETVitalRates <- as.data.frame(BAETVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1))  # make sure we specify that we have 100% survival at Qmin (0.25)
  nls.fit <- nlsLM(formula = (y ~ k* exp(-h*x)), data = df, start = c(k = 1, h = 1)) # fit to negative exponential function
  return(nls.fit)
}

flow.surv.rate <- function(h, k, max, min, interval, Qmin) {
  Q <- seq(min, max, by = interval)
  surv <- k*(exp(-h*Q))
  surv.df <- as.data.frame(cbind(Q, surv))
  surv.df$surv[which(surv.df$Q <= Qmin)] <- 1
  return(surv.df)
}

surv.fit.BAET <- flow.surv.fit(BAETVitalRates$`Max Event Discharge/Bankfull Discharge`, BAETVitalRates$Mortality, 0.25)
surv.df.BAET <- flow.surv.rate(surv.fit.BAET$m$getPars()[2] , surv.fit.BAET$m$getPars()[1], 4, 0.001, 0.001, 0.25)

ggplot(surv.df.BAET, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = BAETVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

