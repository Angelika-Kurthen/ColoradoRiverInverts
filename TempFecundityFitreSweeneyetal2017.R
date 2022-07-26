##################################
# Fitting the Sweeney et al. Temperature v Fecundity Data
##################################

# Values taken from Sweeney et al., 2017 (DOI: 10.1086/696611)  Supplemental Table 1
# https://www-journals-uchicago-edu.ezproxy.proxy.library.oregonstate.edu/doi/suppl/10.1086/696611/suppl_file/TableS1.pdf

temp <- c(12.1, 14.3, 16.2, 20.2, 21.2, 23.9, 27.8, 30, 31.7)
fecundity <- c(1473, 1446, 1231, 935, 1025, 787, 616, 426, 57) 

tempdf <- as.data.frame(cbind(temp, fecundity))

# in the paper, the author say they fit this data to a 3rd degree polynomial (so we will too)

#linear_model1 <- lm(y~x, data=sample_data)
#linear_model2 <- lm(y~poly(x,2,raw=TRUE), data=sample_data)
linear_model3 <- lm(fecundity~poly(temp,3,raw=TRUE), data=tempdf)
#linear_model4 <- lm(y~poly(x,4,raw=TRUE), data=sample_data)
#linear_model5 <- lm(y~poly(x,5,raw=TRUE), data=sample_data)

summary(linear_model3)

# construct model using the coefficient values 
# y = -379.8021x + 16.4664x^2 -0.2684x^3 + 4196.8608

# no we can make a series of temperature and plot both the collected data and the fitted model

tempseq <- seq(7, 32, by = 0.1)
fecfit <- vector()

for (t in tempseq){
  f <- -379.8021*(t) + 16.4664*(t^2) -0.2684*(t^3) + 4196.8608
  fecfit <- append(fecfit, f)
}

plot(tempseq, fecfit, type = "l", col = "gray", lwd = 2, xlab = "Temperature (C)", ylab = "Predicted Fecundity per female")
lines(temp, fecundity, type = "l", col = "red", lwd = 2)
legend(7, 500, legend=c("Fitted Line (y = -379.8021x + 16.4664x^2 - 0.2684x^3 + 4196.8608)", "Sweeney et al. 2017 Data"),
      col=c("gray", "red"), lty=1, cex=0.8)
