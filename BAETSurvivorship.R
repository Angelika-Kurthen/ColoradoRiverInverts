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
polyfit <- nlsLM(dev_days ~ ((l^temperature)/factorial(temperature)) * exp(-l), start = c(l = 1))

# equation=== maturation rate/day = -2.246e-06*temperature^4 + 1.778e-04*temperature^3 -5.062e-03*temperature^2 + 6.483e-02*temperature -3.020e-01 
# equation ==== development time = 2.525e-03*temperature^4 -2.508e-01*temperature^3+  9.379e+00*temperature^2 -1.580e+02*temperature+  1.040e+03 

MaturationRate <- function(x){
  a <- 1/x
  return(a)
}
devtime <- function(temps){
  a <- 2.525e-03*temps^4 -2.508e-01*temps^3+  9.379e+00*temps^2 -1.580e+02*temps +  1.040e+03 
}

TempSurv <- function(temps){
  a <-0.914*exp(-((temps-21)^2)/(2*8.27^2))
return(a)
}

TempSurv <- function(temps){
  a <-  -0.05795*temp^2 + 2.40764*temp -21.94990 
  return(inv.logit(a))
}

devtime(temp)
#survivorship over time
temp<- c(12.1, 14.3, 16.2, 20.2, 21.2, 23.9, 27.8, 30, 31.7, 33.5)
temp <- temp/100
s <- c(43.1, 62.9, 85.4, 87.1, 87, 79, 70.5, 74.3, 57.3,0.001)
s <- s/100
s <- logit(s)
# fit <- nlsLM(s ~ a*temp^4+b*temp^3+ c*temp^2 + d*temp + e, start = c(a =1,b=1,c=1, d = 1, e = 1))
# fit1 <- nlsLM(s ~ a*temp^2+b*temp + c, start = c(a =1, b = 1, c = 1))

ggplot(surv.df.BAET, aes(x = Q, y = surv))+
  geom_line()+
  geom_point(data = BAETVitalRates, aes(x = `Max Event Discharge/Bankfull Discharge` , y = 1-(Mortality), color = Citation))+
  #coord_cartesian(ylim = c(0,1)) +
  ylab('Immediate Post-Disturbance Survival') +
  theme_bw()+
  xlab('`Max Event Discharge/Bankfull Discharge`')

# fit.betalogit <- betareg(s~ temp)
# fit.betalogit2 <- betareg(s ~ a*temp^2 + b*temp + c, start = c(a=1, c=1, c=1))
# predict(fit.betalogit)
# exp(-0.117*(temp^2)+2.736)/(1+(exp((-0.117*temp^2) + 2.736)))

# fit4 <- nlsLM(logit(s) ~ a*temp^4 + b*temp^3 + c*temp^2 + d*temp + e, start = c(a = 1, b = 1, c = 1, d = 1, e = 1) )
fit2 <- nlsLM(logit(s) ~ a*temp^2 + b*temp + c, start = c(a = 1, b = 1, c = 1))
# inv.logit(predict(fit4))
inv.logit(predict(fit2))



plot(temp, s, xlab = "Temperature C", ylab = "Survival", col = "red", pch = 16, cex = 1.5, xlim = c(0,40), ylim = c(0,1))
points(temp, predict(fit.betalogit), col = "blue", pch = 1)
points(temp, inv.logit(predict(fit4)), col = "hotpink", pch = 2)
points(temp, inv.logit(predict(fit2)), col = "green", pch = 5)     
legend(-0.5, 1, legend=c("Data", "Beta Regression 2D Poly", "Logit 4D Poly", "Logit 2D Poly"),
       col=c("red", "blue", "hotpink", "green"), pch=c(16, 1, 2, 5), cex=0.8)     
