#####################################
## Code to fit curve to Hydropsyche spp mortality rates to data
####################################
library(readxl)
library(minpack.lm)
library(tidyverse)
library(car)
library(boot)
library(data.table)
HYOSVitalRates <- read_excel("VitalRates.xlsx", sheet = "Hydropsyche Mortality Rates")
HYOSVitalRates <- as.data.frame(HYOSVitalRates)

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

surv.fit.HYOS <- flow.surv.fit(HYOSVitalRates$`Max Event Discharge/Bankfull Discharge`, HYOSVitalRates$Mortality, 0.25)
surv.df.HYOS <- flow.surv.rate(surv.fit.HYOS$m$getPars()[2] , surv.fit.HYOS$m$getPars()[1], 2, 0.001, 0.001, 0.25)

ggplot(surv.df.HYOS, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = HYOSVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# original rates and temps from McCarty et al 2022, critical thermal max from Kremer and Caldwell 2022
rate <- c(0.0037, 0.0164,0.0293,0.0293,0.00956,0.00321,0.00857,0.00539, 0.0037)
temps <- c(6.9,7.66,8.67,12,2.86,3.28,1.35,1.09, 34.4)
df <- as.data.frame(cbind(temps, rate))
plot(temps, rate)

polyfit <- nlsLM(1/rate ~ a * temps^2 + b * temps + c, start = c(a = 1, b = 1, c = 1))

devtime <- function(x){
  y = 0.598*x^2-19.224*x+219.578
  return(y)
}

# from Gaufin et al 1972, table 2
s <- c(0.0001, 1,1,1,0.9, 0.6, 0.79, 0.45, 0.00001, 0.0001)
temp <- c(0, 18.8, 20.4, 22.1, 24, 28, 29.2, 31.1, 32.7, 36.5)
fit <- nlsLM(logit(s) ~ a*temp^2 + b*temp + c, start = c(a = 1, b = 1, c = 1))

TempSurv <- function(temps){
  a <- -0.0230*temps^2 + 0.8086*temps -3.5611
  return(inv.logit(a))
}

