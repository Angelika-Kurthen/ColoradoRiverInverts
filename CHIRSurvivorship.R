#############################################################
## Code to fit curve to Chirnomus spp mortality rates to data
#############################################################

library(readxl)
library(minpack.lm)
library(tidyverse)
library(car)
library(boot)
CHIRVitalRates <- read_excel("VitalRates.xlsx", sheet = "Chiro Mortality Rates")
CHIRVitalRates <- as.data.frame(CHIRVitalRates)

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

surv.fit.CHIR <- flow.surv.fit(CHIRVitalRates$`Max/Bankful`, CHIRVitalRates$Mortality, 0.25)
surv.df.CHIR <- flow.surv.rate(surv.fit.CHIR$m$getPars()[2] , surv.fit.CHIR$m$getPars()[1], 4, 0.001, 0.001, 0.25)

# Calculate Development Time Based on Literature
# BAETDevRates <- read_excel("VitalRates.xlsx", sheet = "Baetid Development")
# BAETDevRates <- as.data.frame(BAETDevRates)
# polyfit <- nlsLM(DevelopmentTime ~ a*Temperature^4 + b* Temperature ^3 + c*Temperature ^2 +d*Temperature  + e, data = BAETDevRates, start = c(a = 1, b =1, c = 1, d = 1, e = 1))
# 
# devtime <- function(x){
#   a <- 2.525e-03*x^4 -2.508e-01*x^3+  9.379e+00*x^2 -1.580e+02*x +  1.040e+03 
#   return(a)
# }
# 
# MaturationRate <- function(x){
#   a <- 1/x
#   return(a)
# }

# Calculate Temperature Dependent Survival
CHIRSurvRate <- read_excel("VitalRates.xlsx", sheet = "Chiro Survival")
CHIRSurvRate <- as.data.frame(CHIRSurvRate)
#fit <- nlsLM(logit(Survival) ~ a*Temp^4 + b*Temp^3 + c*Temp^2 + d*Temp + e, data = CHIRSurvRate, start = c(a = 1, b = 1, c = 1, d = 1, e = 1))
fit <- nlsLM(logit(Survival) ~ a*Temp^2 + b*Temp + c, data = CHIRSurvRate, start = c(a = 1, b = 1, c = 1))
inv.logit(predict(fit))
#a = -1.016e-04  b = 9.412e-03 c = -3.121e-01 d = 4.317e+00 e =-2.032e+01 
# TempSurv <- function(n){
#   a <-  -1.016e-04*n^4 +  9.412e-03*n^3 -3.121e-01*n^2 + 4.317*n -2.032e+01
#   return(inv.logit(a))
# }
# 
TempSurv <- function(n){
  a <- -0.01462*n^2+  0.67162*n -6.72790
  return(inv.logit(a))
}
# # 
# min.RSS <- function(par){
#   mod <- dnbinom(as.integer(-CHIRSurvRate$Temp + 37.5), size = par[2], prob = par[1])
#   a <- sum(CHIRSurvRate$Survival - (mod*(max(CHIRSurvRate$Survival)/max(mod))))^2
# }
# params <- optim(par = c(0.2, 4), fn = min.RSS, method = "BFGS")
# 
# TempSurv <- function(n){
#   if (n <= 0){
#     a <- 0
#   }else {
#     a <-  dnbinom(as.integer(-n + 37.5), size = params$par[2] , prob = params$par[1])*(max(CHIRSurvRate$Survival)/max(dnbinom(as.integer(-CHIRSurvRate$Temp + 37.5), size =params$par[2], prob = params$par[1])))
#   }
#   return((a))
# }


# ggplot(surv.df.CHIR, aes(x = Q, y = surv))+
#   geom_line()+
#   geom_point(data = CHIRVitalRates, aes(x = `Max/Bankful` , y = 1-(Mortality), color = Citation))+
#   #coord_cartesian(ylim = c(0,1)) +
#   ylab('Immediate Post-Disturbance Survival') +
#   theme_bw()+
#   xlab('`Max Event Discharge/Bankfull Discharge`')

# 
# tem <- seq(0, 40, by = 1)
# plot(CHIRSurvRate$Temp, CHIRSurvRate$Survival, col = "red", pch = 16, xlab = "Temperature", ylab = "Survival", xlim = c(0,40), ylim = c(0, 1))
# lines(tem, TempSurv(tem))
# lines(tem,  dnbinom(as.integer(-tem + 37.5), size = params$par[2] , prob = params$par[1])*(max(CHIRSurvRate$Survival)/max(dnbinom(as.integer(-CHIRSurvRate$Temp + 37.5), size =params$par[2], prob = params$par[1]))))
# plot(temp, s, xlab = "Temperature C", ylab = "Survival", col = "red", pch = 16, cex = 1.5, xlim = c(0,40), ylim = c(0,1))
# points(temp, predict(fit.betalogit), col = "blue", pch = 1)
# points(temp, inv.logit(predict(fit4)), col = "hotpink", pch = 2)points(temp, inv.logit(predict(fit2)), col = "green", pch = 5)     
# legend(-0.5, 1, legend=c("Data", "Beta Regression 2D Poly", "Logit 4D Poly", "Logit 2D Poly"),
#        col=c("red", "blue", "hotpink", "green"), pch=c(16, 1, 2, 5), cex=0.8)     
