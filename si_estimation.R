library("fitdistrplus")
library(MASS)
library(survival)
library(cubature)
library(schoolmath)
library(tidyverse)
library(rstan)
library(loo)
library(sn)
library(readr)
library(lme4)

#time-varying estimation of serial interval

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
data <- read.csv("delta.csv", colClasses = "character")

data$serial_infector <- as.numeric(data$serial_infector)
data1 <- data[data$serial == 1,] # select the 93 transmission pairs
data1 <- data1[data1$serial_infector <= 1,] # select cohort by symptom onset of infectors


data1$expol1 <- as.numeric(data1$expol1)
data1$expou1 <- as.numeric(data1$expou1)
data1$shedl1 <- as.numeric(data1$shedl1)
data1$shedu1 <- as.numeric(data1$shedu1)
data1$onset_infectee1 <- as.numeric(data1$onset_infectee1)
data1$onset_infector1 <- as.numeric(data1$onset_infector1)

input.data1 <- list(
  N = nrow(data1),
  tStartExposure = data1$expol1,
  tEndExposure = data1$expou1,
  tStartshed = data1$shedl1,
  tSymptomOnset = data1$onset_infectee1,
  tSymptomOnsetor = data1$onset_infector1,
  tEndshed = data1$shedu1)
input.data1

## Simulated serial interval using weibull distribution
model.serial <- stan(data = input.data1, 
                        chains = 0, 
                        iter = 0,
                        model_code = "
data{
  int <lower = 1> N;
  vector[N] tSymptomOnset;
  vector[N] tSymptomOnsetor;
}

parameters{
  real shape; 	// Shape parameter of skew normal distribution
  real scale;
}

transformed parameters{
  vector[N] sE; 	// infectious moment
  real minvalue;
  sE = tSymptomOnset - tSymptomOnsetor;
  minvalue = min(sE) - 1;
}

model{
  // Contribution to likelihood of incubation period
  target += weibull_lpdf(sE - minvalue | shape, scale);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(sE[i] - minvalue | shape, scale);
  }
}
"
)

## weibull distribution fitting for serial interval
stanfitG_serial <- stan(fit = model.serial, data = input.data1, 
                           init = "random",
                           warmup = 4000,
                           iter = 14000, 
                           chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG_serial)

shape_serial <- rstan::extract(stanfitG_serial)$shape
scale_serial <- rstan::extract(stanfitG_serial)$scale
xiI_serial <- rstan::extract(stanfitG_serial)$minvalue

quantile(scale_serial*gamma(1+1/shape_serial) + xiI_serial, probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(sqrt(scale_serial^2*(gamma(1+2/shape_serial)-(gamma(1+1/shape_serial))^2)), probs = c(0.025,0.5,0.975)) #sd



