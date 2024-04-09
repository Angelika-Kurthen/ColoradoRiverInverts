##################
# NZMS N Mixture
#################

source()
library(unmarked)
library(ubms)

# plot abundance over time

drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -6 & drift.data.total$RiverMile <= 0),]

NZMS.samp.LF <- sampspec(samp = drift.LF, species = "NZMS", stats = T)
NZMS.samp <- merge(NZMS.samp.LF$Statistics, NZMS.samp.LF$Samples, by = "BarcodeID", all = T)
NZMS.samp <- NZMS.samp[which(NZMS.samp$GearID == 4),] 
NZMS.samp$Density <- NZMS.samp$CountTotal/NZMS.samp$Volume

NZMS.samp <- merge(NZMS.samp, discharge[, 3:4], by = "Date")
NZMS.samp$Density <- (NZMS.samp$Density/(9e-15 *(NZMS.samp$X_00060_00003 * 0.02831683)^4.1))
write.csv(NZMS.samp, "NZMSsamp.csv")


NZMS.samp <- aggregate(NZMS.samp$Density, list(NZMS.samp$Date), FUN = mean)

means <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Group.1 %within% interval(temps$dts[i], temps$dts[i+1]) == T),]
  if (is.nan(mean(d$x)) == T) {
    s = NA
  } else {
    s<- mean(d$x)}
  means <- append(means, s)
}


NZMS.samp.sum <- na.omit(as.data.frame(cbind(as.Date(means.list.NZMS$`temps$dts`, format = "%Y-%m-%d"), means[51:406])))
NZMS.samp.sum$V1 <- as.Date(NZMS.samp.sum$V1, origin = "1970-01-01")


temps <- read.csv(file = "temps.csv", header = T)

NZMS.samp <- read.csv("NZMSsamp.csv", header = T)

# need to fit our data to an unmarked frame

# observations need to be in RxJ matrix, where R is # of sites and J is max number of obs per site
# we measured over temps$dts (those are our timesteps)

max_visits <- vector() # vector to put max obs per site
means <- vector()
for (i in 1:length(temps$dts)){  # cycle through all the 14 day timesteps that we have model output for
  # pull abundances between each each timestep
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  max_visits[i] <- length(d$Date) # number of observations
  means[i] <- mean(d$CountTotal) # means
}

# make data frame with the timestep, the mean, and standard error
df <- cbind.data.frame(temps$dts, means)
# remove anythig where there are NA values of Density (means that volume or count is missing)
df <- df[!is.na(df$means), ]

# define our RxJ matrix
R <- length(temps$dts)
J <- max(max_visits)
site_mat <- matrix(data = NA, nrow = R, ncol = J)
dens_mat <- matrix(data = NA, nrow = R, ncol = J)

# make vector for flows at each timestep
# make RxJ matrix full of densities
# make RxJ matrix full of raw counts
# make RxJ matrix full of volumes sampled for each abundance
flows <- vector()
volumes <- matrix(data = NA, nrow = R, ncol = J)
for (i in 1:length(temps$dts)){
  d <- NZMS.samp[which(NZMS.samp$Date >= temps$dts[i] & NZMS.samp$Date < temps$dts[i+1]), ]
  flows[i] <- mean(d$X_00060_00003)
  site_mat[i, ] <- c(d$CountTotal, rep(NA, times = (J- length(d$CountTotal))))
  dens_mat[i, ] <- c(as.integer(d$Density), rep (NA, times = (J - length(d$Density))))
  volumes[i, ] <- c(scale(d$Volume),rep(NA, times = (J- length(d$CountTotal))))
}

# we need to remove all timesteps that are just NAs
nodata <- which(is.na(site_mat[,1]))
# first identify all the timsteps that don't have data (so we can match them up later)
site_mat <- as.matrix(site_mat[-nodata,])
dens_mat <- as.matrix(dens_mat[-nodata, ])
flows <- as.data.frame(flows[-nodata])
volumes <- as.matrix(volumes[-nodata, ])
dimnames(volumes) <- list(temps$dts[-nodata], seq(1:48))
volumes <- list(volumes)
names(volumes) <- c("vol")

# convert organized data into unmarked dataframe format, with flow as site covariate and volume of sample as observation covariate
dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)

# pcount is N-mixture model to estimate counts at each site
# we are doing space-time approximation

# we want to offset our count data by volume (since volumes vary slighlty by each sample)
# these are basically our three different hypotheses
#Poisson mixture
mod_um_vol_P <- pcount(~ 1 + offset(vol) ~1 , dat, K = 11000, mixture = "P")
# Neg Binom mixture
mod_um_vol_NB <- pcount(~ 1 + offset(vol) ~1 , dat, K = 11000, mixture = "NB")
# Zero Inflated Poisson mixture
mod_um_vol_ZIP <- pcount(~ 1 + offset(vol) ~1 , dat, K =11000, mixture = "ZIP")

# now we want to compare the models
fm <- fitList(mod_um_vol_P, mod_um_vol_NB, mod_um_vol_ZIP)
modSel(fm)
# Model selectionfit
# suggests that NB is best model


# use parboot function to imulate datasets from a fitted model, refit the model, 
#and generate a sampling distribution for a user-specified fit-statistic

# Function returning three fit-statistics
fitstats <- function(fm, na.rm=TRUE) {
  observed <- getY(fm@data)
  expected <- fitted(fm)
  resids <- residuals(fm)
  sse <- sum(resids^2, na.rm=na.rm)
  chisq <- sum((observed - expected)^2 / expected, na.rm=na.rm)
  freeTuke <- sum((sqrt(observed) - sqrt(expected))^2, na.rm=na.rm)
  out <- c(SSE=sse, Chisq=chisq, freemanTukey=freeTuke)
  return(out)
}

pb <- parboot(mod_um_vol_NB, fitstats, nsim=25, report=1)
pb

pb_plot_mod_vol <- plot(pb, main="")
plot(pb_plot_mod_vol)

pb_plot_mod_dens <- plot(pb_dens, main = "")

# Compare to empirical Bayes confidence intervals
colSums(confint(ranef(mod_um_vol_NB, K=11000)))
# Estimates of conditional abundance distribution at each site
re <- ranef(mod_um_vol_NB)
# Best Unbiased Predictors
est <- bup(re) # estimated population size at each timestep
conf <- as.data.frame(confint(re, level=0.95)) # 90% CI
# get matching dates that match with non-NA estimations
date <- df$`temps$dts`[!is.nan(est)]
conf <- (conf[!is.nan(est),])


# remove the NAs from out estimation 
est <- est[!is.nan(est)]
# check for temporal autocorrelation, there is some at 8 and 9, and 17?
acf(est)


estimate <- cbind.data.frame(date, est, conf)
estimate$date <- as.Date(estimate$date)

ggplot(data = estimate, aes(x = (date), y = est))+
    geom_line(size = 1, alpha = 0.8)+
    geom_point(size = 2, alpha = 0.5)
    #geom_ribbon(aes(ymin = `2.5%`, ymax = `97.5%` ), alpha = 0.2)




#


#dens_mod_mull <- stan_pcount(~1 ~1, dat_dens, K = 1176244, chains = 3, iter = 2000, seed = 123, log_lik = T)

# occupancy to get global p 
# site_mat[site_mat >0] <- 1
# occu_dat <- unmarkedFramePCount(y = site_mat, siteCovs = flows, obsCovs = volumes)
# 
# occu_mod_null <- stan_occu(~1 ~1, occu_dat, chains = 3, iter = 2000)
# occu_mod_vol <- stan_occu(~scale(vol)~1, occu_dat, chains = 3, iter = 2000)
# 
# #posterior
# names(mod_um)
# ## [1] "beta_state[(Intercept)]" "beta_det[(Intercept)]"
# occ_intercept <- extract(mod_um, "det")[[1]]
# hist(occ_intercept, freq=FALSE)
# lines(density(occ_intercept), col='red', lwd=2)
# 
# #Compare the models
# #First we combine the models into a fitList:
# mods <- fitList(mod_um, mod_um_vol)
# #Then we generate a model selection table:
# round(modSel(mods), 3)
# #the model with the largest elpd performed best
# plot_residuals(mod_um_vol, submodel="det")
# fit_top_gof <- gof(mod_um_vol, raws=100, quiet=TRUE)
# 
# prob_det <- predict(mod_um_vol, submodel="det")
# # from predicted detection probabilities, get range
# range(prob_det$Predicted, na.rm = T)
# 
# 
