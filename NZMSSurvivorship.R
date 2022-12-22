############################################################
## Code to fit curve to P. anitpodarum mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
NZMSVitalRates <- read_excel("VitalRates.xlsx", sheet = "NZMS Mortality Rates")
NZMSVitalRates <- as.data.frame(NZMSVitalRates)

flow.surv.fit <- function(magnitude, mortality, Qmin){
  x <- magnitude # rename x and y
  y <- 1 - mortality # convert to survival
  df <- as.data.frame(cbind(x,y))
  df <- rbind(df, c(Qmin, 1.375))  # make sure we specify that we have 100% survival at Qmin (0.25)
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

surv.fit.NZMS <- flow.surv.fit(NZMSVitalRates$`Max Event Discharge/Bankfull Discharge`, NZMSVitalRates$Mortality, 0.25)
surv.df.NZMS <- flow.surv.rate(surv.fit.NZMS$m$getPars()[2] , surv.fit.NZMS$m$getPars()[1], 2, 0.001, 0.001, 0.25)

ggplot(surv.df.NZMS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = NZMSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

