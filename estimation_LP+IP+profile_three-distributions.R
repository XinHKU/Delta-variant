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

### simulating latent, incubation, infectiousness profile for lognormal, gamma and weibull distribution

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
data <- read.csv("delta.csv", colClasses = "character")



data <- data[data$serial == 1,]

data$expol1 <- as.numeric(data$expol1)
data$expou1 <- as.numeric(data$expou1)
data$shedl1 <- as.numeric(data$shedl1)
data$shedu1 <- as.numeric(data$shedu1)
data$onset_infectee1 <- as.numeric(data$onset_infectee1)
data$onset_infector1 <- as.numeric(data$onset_infector1)

input.data <- list(
  N = nrow(data),
  tStartExposure = data$expol1,
  tEndExposure = data$expou1,
  tStartshed = data$shedl1,
  tSymptomOnset = data$onset_infectee1,
  tSymptomOnsetor = data$onset_infector1,
  tEndshed = data$shedu1)


## Simulated using lognormal distribution
model.lognormal <- stan(data = input.data, 
                        chains = 0, 
                        iter = 0,
                        model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tStartshed;
  vector[N] tEndshed;
  vector[N] tSymptomOnset;
  vector[N] tSymptomOnsetor;
}

parameters{
  real<lower = 0> alphaP;
  real<lower = 0> betaP;
  real<lower = 0> alphaI;
  real<lower = 0> betaI;
  real<lower = 0> alphaL;
  real<lower = 0> betaL; 
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
  vector<lower = 0, upper = 1>[N] uS;
}

transformed parameters{
  vector[N] tE; 	// infection moment
  vector[N] sE; 	// infectious moment
  vector[N] ma;
  real minfec;
  vector[N] Infprofile;
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
  for (h in 1:N) {
    ma[h] = fmax(tE[h], tStartshed[h]);
  }  
  sE = ma + uS .* (tEndshed - ma);
  Infprofile = tE - tSymptomOnsetor;
  minfec = min(Infprofile) - 1;
}

model{
  // Contribution to likelihood of incubation period
  target += lognormal_lpdf(tSymptomOnset - tE  | alphaI, betaI);
  target += lognormal_lpdf(sE - tE  | alphaL, betaL);
  target += lognormal_lpdf(Infprofile - minfec | alphaP, betaP);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = lognormal_lpdf(Infprofile[i] - minfec | alphaP, betaP) + lognormal_lpdf(tSymptomOnset[i] - tE[i]  | alphaI, betaI) + lognormal_lpdf(sE[i] - tE[i]  | alphaL, betaL);
  }
}
"
)

#lognormal distribution fitting
stanfitG_lognormal <- stan(fit = model.lognormal, data = input.data, 
                           init = "random",
                           warmup = 4000,
                           iter = 14000, 
                           chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG_lognormal) # check results and convergence


LLG_lognormal = extract_log_lik(stanfitG_lognormal, parameter_name = 'log_lik') # Pareto Loo PSIS
loo2_lognormal <- loo(LLG_lognormal)
loo2_lognormal

#extract parameters
xiI <- rstan::extract(stanfitG_lognormal)$minfec 
alpha_infec <- rstan::extract(stanfitG_lognormal)$alphaP
beta_infec <- rstan::extract(stanfitG_lognormal)$betaP
alpha_incu <- rstan::extract(stanfitG_lognormal)$alphaI 
beta_incu <- rstan::extract(stanfitG_lognormal)$betaI
alpha_latent <- rstan::extract(stanfitG_lognormal)$alphaL 
beta_latent <- rstan::extract(stanfitG_lognormal)$betaL

quantile(exp(alpha_infec + (beta_infec^2)/2) + xiI, probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(exp(alpha_incu + (beta_incu^2)/2), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(exp(alpha_latent + (beta_latent^2)/2), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for latent


quantile(sqrt((exp(beta_incu^2)-1) * exp(2*alpha_incu + beta_incu^2)), probs = c(0.025,0.5,0.975))
quantile(sqrt((exp(beta_latent^2)-1) * exp(2*alpha_latent + beta_latent^2)), probs = c(0.025,0.5,0.975))
quantile(sqrt((exp(beta_infec^2)-1) * exp(2*alpha_infec + beta_infec^2)), probs = c(0.025,0.5,0.975))

percentilesL_latent <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                              function(p) quantile(qlnorm(p = p, meanlog = alpha_latent, sdlog = beta_latent), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_latent) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesL_latent

percentilesL_incu <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                              function(p) quantile(qlnorm(p = p, meanlog = alpha_incu, sdlog = beta_incu), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_incu) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesL_incu

percentilesL_infec <- sapply(c(-7, -4, 0, 4), 
                       function(p) quantile(plnorm(p - xiI, meanlog = alpha_infec, sdlog = beta_infec), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_infec) <- c(-7, -4, 0, 4)
percentilesL_infec







##Gamma distribution

model.gamma <- stan(data = input.data, 
                    chains = 0, 
                    iter = 0,
                    model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tStartshed;
  vector[N] tEndshed;
  vector[N] tSymptomOnset;
  vector[N] tSymptomOnsetor;
}

parameters{
  real<lower = 0> alphaP;
  real<lower = 0> betaP;
  real<lower = 0> alphaI;
  real<lower = 0> betaI;
  real<lower = 0> alphaL;
  real<lower = 0> betaL; 
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
  vector<lower = 0, upper = 1>[N] uS;
}

transformed parameters{
  vector[N] tE; 	// infection moment
  vector[N] sE; 	// infectious moment
  vector[N] ma;
  real minfec;
  vector[N] Infprofile;
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
  for (h in 1:N) {
    ma[h] = fmax(tE[h], tStartshed[h]);
  }  
  sE = ma + uS .* (tEndshed - ma);
  Infprofile = tE - tSymptomOnsetor;
  minfec = min(Infprofile) - 1;
}

model{
  // Contribution to likelihood of incubation period
  target += gamma_lpdf(tSymptomOnset - tE  | alphaI, betaI);
  target += gamma_lpdf(sE - tE  | alphaL, betaL);
  target += gamma_lpdf(Infprofile - minfec | alphaP, betaP);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = gamma_lpdf(Infprofile[i] - minfec | alphaP, betaP) + gamma_lpdf(tSymptomOnset[i] - tE[i]  | alphaI, betaI) + gamma_lpdf(sE[i] - tE[i]  | alphaL, betaL);
  }
}
"
)

#gamma distribution fitting
stanfitG_gamma <- stan(fit = model.gamma, data = input.data, 
                       init = "random",
                       warmup = 4000,
                       iter = 14000, 
                       chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG_gamma) # check results and convergence


LLG_gamma = extract_log_lik(stanfitG_gamma, parameter_name = 'log_lik') # Pareto Loo PSIS
loo2_gamma <- loo(LLG_gamma)
loo2_gamma

xiIG <- rstan::extract(stanfitG_gamma)$minfec #extract parameters
alpha_infecG <- rstan::extract(stanfitG_gamma)$alphaP
beta_infecG <- rstan::extract(stanfitG_gamma)$betaP
alpha_incuG <- rstan::extract(stanfitG_gamma)$alphaI #extract parameters
beta_incuG <- rstan::extract(stanfitG_gamma)$betaI
alpha_latentG <- rstan::extract(stanfitG_gamma)$alphaL #extract parameters
beta_latentG <- rstan::extract(stanfitG_gamma)$betaL

quantile((alpha_infecG / beta_infecG) + xiIG, probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile((alpha_incuG / beta_incuG), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile((alpha_latentG / beta_latentG), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for latent


quantile(sqrt(alpha_infecG / beta_infecG^2), probs = c(0.025,0.5,0.975))
quantile(sqrt(alpha_incuG / beta_incuG^2), probs = c(0.025,0.5,0.975))
quantile(sqrt(alpha_latentG / beta_latentG^2), probs = c(0.025,0.5,0.975))






#weibull distribution

model.weibull <- stan(data = input.data, 
                      chains = 0, 
                      iter = 0,
                      model_code = "
data{
  int <lower = 1> N;
  vector[N] tStartExposure;
  vector[N] tEndExposure;
  vector[N] tStartshed;
  vector[N] tEndshed;
  vector[N] tSymptomOnset;
  vector[N] tSymptomOnsetor;
}

parameters{
  real<lower = 0> alphaP;
  real<lower = 0> betaP;
  real<lower = 0> alphaI;
  real<lower = 0> betaI;
  real<lower = 0> alphaL;
  real<lower = 0> betaL; 
  vector<lower = 0, upper = 1>[N] uE;	// Uniform value for sampling between start and end exposure
  vector<lower = 0, upper = 1>[N] uS;
}

transformed parameters{
  vector[N] tE; 	// infection moment
  vector[N] sE; 	// infectious moment
  vector[N] ma;
  real minfec;
  vector[N] Infprofile;
  tE = tStartExposure + uE .* (tEndExposure - tStartExposure);
  for (h in 1:N) {
    ma[h] = fmax(tE[h], tStartshed[h]);
  }  
  sE = ma + uS .* (tEndshed - ma);
  Infprofile = tE - tSymptomOnsetor;
  minfec = min(Infprofile) - 1;
}

model{
  // Contribution to likelihood of incubation period
  target += weibull_lpdf(tSymptomOnset - tE  | alphaI, betaI);
  target += weibull_lpdf(sE - tE  | alphaL, betaL);
  target += weibull_lpdf(Infprofile - minfec | alphaP, betaP);
}

generated quantities {
  // likelihood for calculation of looIC
  vector[N] log_lik;
  for (i in 1:N) {
    log_lik[i] = weibull_lpdf(Infprofile[i] - minfec | alphaP, betaP) + weibull_lpdf(tSymptomOnset[i] - tE[i]  | alphaI, betaI) + weibull_lpdf(sE[i] - tE[i]  | alphaL, betaL);
  }
}
"
)

#weibull distribution fitting
stanfitG_weibull <- stan(fit = model.weibull, data = input.data, 
                         init = "random",
                         warmup = 4000,
                         iter = 14000, 
                         chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG_weibull) # check results and convergence


LLG_weibull = extract_log_lik(stanfitG_weibull, parameter_name = 'log_lik') # Pareto Loo PSIS
loo2_weibull <- loo(LLG_weibull)
loo2_weibull

xiIW <- rstan::extract(stanfitG_weibull)$minfec #extract parameters
alpha_infecW <- rstan::extract(stanfitG_weibull)$alphaP
beta_infecW <- rstan::extract(stanfitG_weibull)$betaP
alpha_incuW <- rstan::extract(stanfitG_weibull)$alphaI #extract parameters
beta_incuW <- rstan::extract(stanfitG_weibull)$betaI
alpha_latentW <- rstan::extract(stanfitG_weibull)$alphaL #extract parameters
beta_latentW <- rstan::extract(stanfitG_weibull)$betaL

quantile(beta_infecW*gamma(1+1/alpha_infecW) + xiIW, probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(beta_incuW*gamma(1+1/alpha_incuW), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(beta_latentW*gamma(1+1/alpha_latentW), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for latent


quantile(sqrt(beta_infecW^2*(gamma(1+2/alpha_infecW)-(gamma(1+1/alpha_infecW))^2)), probs = c(0.025,0.5,0.975))
quantile(sqrt(beta_incuW^2*(gamma(1+2/alpha_incuW)-(gamma(1+1/alpha_incuW))^2)), probs = c(0.025,0.5,0.975))
quantile(sqrt(beta_latentW^2*(gamma(1+2/alpha_latentW)-(gamma(1+1/alpha_latentW))^2)), probs = c(0.025,0.5,0.975))

quantile(beta_infecW*(((alpha_infecW-1)/alpha_infecW)^(1/alpha_infecW)) + xiIW, probs = c(0.025,0.5,0.975))

percentilesL_latentW <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                              function(p) quantile(qweibull(p = p, shape = alpha_latentW, scale = beta_latentW), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_latentW) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesL_latentW

percentilesL_incuW <- sapply(c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99), 
                            function(p) quantile(qweibull(p = p, shape = alpha_incuW, scale = beta_incuW), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_incuW) <- c(0.025, 0.05, 0.5, 0.95, 0.975, 0.99)
percentilesL_incuW

percentilesL_infecW <- sapply(c(-7, -4, 0, 4), 
                             function(p) quantile(pweibull(p - xiI, shape = alpha_infecW, scale = beta_infecW), probs = c(0.025, 0.5, 0.975)))
colnames(percentilesL_infecW) <- c(-7, -4, 0, 4)
percentilesL_infecW


