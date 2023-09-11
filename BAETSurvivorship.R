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
  df <- rbind(df, c(Qmin, 1.11))  # make sure we specify that we have 100% survival at Qmin (0.25)
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


#Sweeney et al 2017, Baetidae development rate follows 4th degree polynomial eq
# data from Supplementary Material
dev_days <- c(110.3, 72.6, 46.9, 29, 25.5, 20.6, 16, 13.5, 16.9)
temperature <- c(12.1, 14.3, 16.2, 20.2, 21.2, 23.9, 27.8, 30, 31.7)
df <- as.data.frame(cbind(dev_days, temperature))
polyfit <- nlsLM(dev_days ~ a*temperature^4 + b*temperature^3 + c*temperature^2 + d*temperature + e, data = df, start = list(a = 1, b = 1, c = 1, d = 1, e = 1))

# equation=== maturation rate/day = -2.246e-06*temperature^4 + 1.778e-04*temperature^3 -5.062e-03*temperature^2 + 6.483e-02*temperature -3.020e-01 
# equation ==== development time = 2.525e-03*temperature^4 -2.508e-01*temperature^3+  9.379e+00*temperature^2 -1.580e+02*temperature+  1.040e+03 

temps <- c(5, 6, 7, 8, 6, 4)
m_t <- function(temps){
  -2.246e-06*temps^4 + 1.778e-04*temps^3 -5.062e-03*temps^2 + 6.483e-02*temps -3.020e-01 
}
d_t <- function(temps){
  2.525e-03*temps^4 -2.508e-01*temps^3+  9.379e+00*temps^2 -1.580e+02*temps +  1.040e+03 
}



originaltaus <- 32.47879 27.05311 22.40065 18.44206 27.05311 38.76137

d(t)  = d(t - tau(t))*(m(T)/m(T(t-tau(t)))

d(t)  = d(t - tau(t))*(m(T)/m(T(t-tau(t))) 
round(1/d(t))/14.

ggplot(surv.df.BAET, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = BAETVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')


