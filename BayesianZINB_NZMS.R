###########################################################
# Fit ZI NB model to invasive NZ mudsnail data from the Colorado River, AZ
# code originally from https://github.com/songsqian/BeesRStan/blob/main/SFS/2024/SFS%202024%20Short%20Course-2-3.Rmd
############################################################

# load libraries needed
library(dataRetrieval)
library(lubridate)
require(rstan)
library(rv)
library(car)

# function to index and summarize temperature data over timesteps length
TimestepTemperature <- function(temp){
  # Make an index to be used for aggregating
  ID <- as.numeric(as.factor(temp$Date))-1
  # want it to be every 14 days, hence the 14
  ID <- ID %/% 14
  # aggregate over ID and TYPE for all numeric data.
  outs <- aggregate(temp[sapply(temp,is.numeric)],
                    by=list(ID),
                    FUN=mean)
  # format output
  names(outs)[1:2] <-c("dts","Temperature")
  # add the correct dates as the beginning of every period
  outs$dts <- as.POSIXct(temp$Date[((outs$dts*14)+1)], origin = "1970-01-01")
  # order by date in chronological order
  temps <- outs[order(outs$dts),]
  
  return(temps)
}


# make sure rstan can run in parallel
rstan_options(auto_write = TRUE)
options(mc.cores = min(c(parallel::detectCores(), 8)))

# specify run/iterations
nchains <-  min(c(parallel::detectCores(), 8))
niters <- 5000
nkeep <- 2500
nthin <- ceiling((niters/2)*nchains/nkeep)

# load my data
NZMSsamp <- read.csv("~/ColoradoRiverInverts/NZMSsamp.csv")

# make sure dates are dates
NZMSsamp$Date <- as.Date(NZMSsamp$Date)

# remove any samples that are flagged as "strange"
NZMSsamp <- NZMSsamp[which(NZMSsamp$FlagStrange == F),]

# read in temperature data - optional
temp <- readNWISdv("09380000", "00010", "2007-10-01", "2023-05-01")

# put temp data into 2 week timestep bins
temps <- TimestepTemperature(temp)

# we have uneven sample distribution (some days, 11 samples taken, others only 1)
# therefore, we use 2 week timesteps
# we now will index each sample based on which timestep it falls into
vals <- vector()
for (i in 1:length(temps$dts)){
  d <- NZMSsamp[which(NZMSsamp$Date %within% interval(temps$dts[i], temps$dts[i+1]-1) == T),]
  if (length(d$CountTotal) > 0) {
    s<- rep(i, times = length(d$CountTotal))
    vals <- append(vals, s)}
}

# phenology may also play into dynamics, so include month column as well
month <- month(NZMSsamp$Date)

# add to data frame
NZMSsamp <- cbind(NZMSsamp, vals, month)

## Zero-inflated negative binomial stan model with covariates and year hierarchy
## from Mark DuFour's example of juvenile sturgen with yearly hierachy
zinb_sturg2 <-"
data {
  int<lower=1> N0;      // number of zero counts
  int<lower=1> Np;      // number of positive counts
  int<lower=0> Yp[Np];  // response variable (+ counts only)

  real Zp[Np];          // zero model design matrix
  real Z0[N0];          // zero model design matrix
  
  int<lower=1> Ka;      // number of count model coefficients
  matrix[Np, Ka] Xp;    // count model design matrix
  matrix[N0, Ka] X0;    // count model design matrix

  vector[Np] offsetP;   // exposure variable (Volume)
  vector[N0] offset0;   // exposure variable (Volume)
  
  int<lower=1> Nts;                   // number of timesteps
  int<lower=1,upper=Nts> tsP[Np];   // timestep hierarchy
  int<lower=1,upper=Nts> ts0[N0];   // timestep hierarchy
  }

parameters {
  real alpha0;          // zero model intercept
  real alpha1;          // zero model coefficient
  vector[Ka] beta;      // count model coefficients
  real<lower=0> r;      // inverse disperson parameter
  
  vector[Nts] b0;       // count model timestep specific intercept
  real mu_b;            // hyper mean for b0
  real<lower=0> tau_b;  // hyper sd for b0
  }
  
transformed parameters{
  vector[Np] zi_p;          // initialize postive count predictor term for zi
  vector[N0] zi_0;          // initialize zero predictor term for zi
  vector[Np] eta_p;         // initialize positive count predictor term for eta
  vector[N0] eta_0;         // initialize zero count predictor term for eta

  for (i in 1:Np){          // linear predictor for postivie counts
    eta_p[i] = b0[tsP[i]]+ Xp[i] * beta + offsetP[i];
    zi_p[i]  = alpha0 + Zp[i] * alpha1;
  }
  
  for (i in 1:N0){          // linear predictor for true zeros
    eta_0[i] = b0[ts0[i]]+ X0[i] * beta + offset0[i];
    zi_0[i]  = alpha0 + Z0[i] * alpha1;
  }
 
  }
  
model {
  // priors
  alpha0 ~ normal(0,5);         // zero model intercept
  alpha1 ~ normal(0,5);         // zero model coefficients
  beta ~ normal(0,5);           // count model coefficients
  r ~ gamma(0.001, 0.001);        // inverse dispersion parameter
  
  b0 ~ normal(mu_b, tau_b);     // timestep specific count model intercepts
  mu_b ~ normal(0,5);           // hyper mean for b0
  tau_b ~ normal(0,5);          // hyper sd for b0

  // likelihood
  for(i in 1:N0){    
    target += log_sum_exp(bernoulli_logit_lpmf(1 | zi_0[i]), 
                       bernoulli_logit_lpmf(0 | zi_0[i]) + 
                       neg_binomial_2_log_lpmf(0 | eta_0[i], r)); 
    } 
    
  for(i in 1:Np){ 
    target += bernoulli_logit_lpmf(0 | zi_p[i]) +  
             neg_binomial_2_log_lpmf(Yp[i] | eta_p[i], r); 
    } 
}
"

##############################################
# put the model into stan
fit_zinb_sturg2 <- stan_model(model_code = zinb_sturg2)

## bundle data, initial values, n.chains, and parameters to monitor
ZINB_in <- function(data=NZMSsamp, n.chains=nchains){
  y <- data$CountTotal # count data
  tmp <- !is.na(y) # isn't NAs
  data <- data[tmp,] #remove NAs
  N <- dim(data)[1] #get length (N)
  Y <- y[tmp] # count data that isn't NA
  # log(mu) = X(beta) --> predictor variables for mean
  X <- cbind(data$TimeElapsed - mean(data$TimeElapsed), 
             data$X_00060_00003 - mean(data$X_00060_00003), #flow (cfs)
             data$month - mean(data$month)) # month
  Kbeta <- dim(X)[2] # number of columns of predictor variable
  # logit(theta) = Z(alpha) --> predictor variables for true 0
  Z <- data$X_00060_00003-mean(data$X_00060_00003) # just flow
  offsets <- log(data$Volume) # offset by log Volume filtered
  timestp <- as.numeric(ordered(data$vals)) #ordered list of timesteps
  nts <- length(unique(vals)) #number of timesteps
  tmp <- y==0 # which catch numbers are 0
  np <- sum(!tmp) # number of non 0 values
  n0 <- sum(tmp) # number of 0 values
  data <- list(Np=np, N0=n0, Nts=nts, Ka=Kbeta,
               Xp=X[!tmp,], X0=X[tmp,], Yp=y[!tmp],
               Zp=Z[!tmp], Z0=Z[tmp],
               offset0=offsets[tmp], offsetP=offsets[!tmp],
               tsP=timestp[!tmp], ts0=timestp[tmp])
  inits <- list()
  for (i in 1:n.chains)
    inits[[i]] <- list(alpha0=rnorm(1), alpha1=rnorm(1), b0=rnorm(nts), beta=rnorm(Kbeta), r=runif(1),   tau_b=runif(1)) # mu_b=rnorm(1),
  paras <- c("alpha0","alpha1","b0","beta","r","tau_b", "mu_b")
  return(list(data=data, init=inits, nchains=n.chains, para=paras ))
}

# put into stan
input.to.stan <- ZINB_in()

## fun model --> run 
keep_zinb_sturg2 <- sampling(fit_zinb_sturg2, data=input.to.stan$data,
                             init=input.to.stan$init,
                             pars=input.to.stan$para,
                             iter=niters,thin=nthin,
                             chains=input.to.stan$nchains)

##save(keep_zinb_sturg2,file="zinb_sturg2.RData")

## view pairs plots 
pairs(keep_zinb_sturg2, pars = c("mu_b","tau_b","r"))
#pairs(keep_zinb_sturg2, pars = c("b0","beta","r"))
pairs(keep_zinb_sturg2, pars = c("alpha0","alpha1","r"))

#Usually runs with no problem, sometimes has 1 diverging chain
# pairs for mu_b and tau_b and r look good
# pairs for alpha0 and alpha1 and r do not look good
#The model runs quickly without producing any warnings. Looking over the pair plots, there are no concerning patterns in the upper level parameters. Correlations between `r` and the `b0` and `beta[1:12]` look good. However, the non-identifiability issue in the zero model parameters (`alpha0` and `alpha1`) remains.


#```{R zinb_sturg2 standardized index, tidy=TRUE}

## extract MCMC draws and summarize parameters
fitcoef <- rvsims(as.matrix(as.data.frame(
  rstan::extract(keep_zinb_sturg2, permute = T))))
options(digits = 3)
summary(fitcoef)


## un-center intercepts
real_zi_intercept <- fitcoef[1] - (fitcoef[2]*mean(NZMSsamp$X_00060_00003)) 
real_Intercept_12 <- fitcoef[3:139] - (fitcoef[140]*mean(NZMSsamp$TimeElapsed) + fitcoef[141]*mean(NZMSsamp$X_00060_00003) + fitcoef[142]*mean(NZMSsamp$month))


## standardized index - sample at survey mean volume
std_index <- exp(real_Intercept_12 + mean(log(NZMSsamp$Volume)))


#tikz(file=paste(plotDIRch5, "zinbPred.tex", sep="/"),
#     height=2.75, width=3.25
#, standAlone=F)

# plot
par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.5, 0.125,0), las=1, tck=0.01)
x.dtsf <- seq(0,1,by=0.01)
plot(jitter(as.numeric(as.factor(NZMSsamp$vals))),
     NZMSsamp$CountTotal/NZMSsamp$Volume,pch=1,col="gray", cex=0.5, xaxt="n",
     xlab="Sampling Time", ylab="#/m3",ylim = c(0, 1000))
nom_index <- aggregate(CountTotal/Volume~vals, data=NZMSsamp,FUN="mean")
points(nom_index$'CountTotal/Volume', col="black", typ="b",pch=4)
points(std_index,cex = 0.5)

# this plot is hard to read - will take a different approach 
# just the b0s 

est <- vector()
for (i in 1:length(unique(vals))){
  est[i] <- mean(std_index[[i]])
}


par(mfrow=c(1,1), mar=c(3,3,3,1), mgp=c(1.5, 0.125,0), las=1, tck=0.01)
x.dtsf <- seq(0,1,by=0.01)
plot(jitter(as.numeric(as.factor(NZMSsamp$vals))),
     (NZMSsamp$CountTotal/NZMSsamp$Volume),pch=1,col="gray", cex=0.5, xaxt="n",
     xlab="Sampling Time", ylab="#/m3")

nom_index <- aggregate(CountTotal/Volume~vals, data=NZMSsamp,FUN="mean")
points(nom_index$'CountTotal/Volume', col="black", typ="b",pch=4)
points(est, col = "red", cex = 0.5)


# real points probably fall within 95% CI of all predictions but how can we tell for sure?

# what are other ways to test how well this model fits? And how can we extract values
# of mu_b (model mean for each timestep)