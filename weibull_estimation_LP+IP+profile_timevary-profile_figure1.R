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




#weibull distribution jointly estimation of latent,incubation and infectiousness profile

##model construstion
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

##weibull distribution fitting
stanfitG_weibull <- stan(fit = model.weibull, data = input.data, 
                         init = "random",
                         warmup = 4000,
                         iter = 14000, 
                         chains = 4, control = list(adapt_delta = 0.99))
print(stanfitG_weibull) # check results and convergence


LLG_weibull = extract_log_lik(stanfitG_weibull, parameter_name = 'log_lik') # Pareto Loo PSIS
loo2_weibull <- loo(LLG_weibull)
loo2_weibull

xiIW <- rstan::extract(stanfitG_weibull)$minfec #extract parameters: minximum value in estimation of infectiousness profile 
alpha_infecW <- rstan::extract(stanfitG_weibull)$alphaP #extract parameters
beta_infecW <- rstan::extract(stanfitG_weibull)$betaP #extract parameters
alpha_incuW <- rstan::extract(stanfitG_weibull)$alphaI #extract parameters 
beta_incuW <- rstan::extract(stanfitG_weibull)$betaI #extract parameters
alpha_latentW <- rstan::extract(stanfitG_weibull)$alphaL #extract parameters
beta_latentW <- rstan::extract(stanfitG_weibull)$betaL #extract parameters

quantile(beta_infecW*gamma(1+1/alpha_infecW) + xiIW, probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for infectiousness profile
quantile(beta_incuW*gamma(1+1/alpha_incuW), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for incubation
quantile(beta_latentW*gamma(1+1/alpha_latentW), probs = c(0.025,0.5,0.975)) # posterior mean and 95%CI of mean for latent


quantile(sqrt(beta_infecW^2*(gamma(1+2/alpha_infecW)-(gamma(1+1/alpha_infecW))^2)), probs = c(0.025,0.5,0.975)) #SD
quantile(sqrt(beta_incuW^2*(gamma(1+2/alpha_incuW)-(gamma(1+1/alpha_incuW))^2)), probs = c(0.025,0.5,0.975)) #SD
quantile(sqrt(beta_latentW^2*(gamma(1+2/alpha_latentW)-(gamma(1+1/alpha_latentW))^2)), probs = c(0.025,0.5,0.975)) #SD

quantile(beta_infecW*(((alpha_infecW-1)/alpha_infecW)^(1/alpha_infecW))+xiIW, probs = c(0.025,0.5,0.975)) #mode of infectiousness profile


##forward estimate of time-varying infectiousness profile mode
Infprofile1 <- rstan::extract(stanfitG_weibull)$Infprofile #extract infectiousness profile for each case in the model
forward <- Infprofile1[,1:93] #select the cases
forward_vector <- c(forward[,1:93]) #vectorized
hist(forward_vector) # pictured it
min_value <- min(forward_vector)-1 # define the minimum value
fw <- fitdist(forward_vector-min_value, "weibull") # fit weibull distribution for the values
result <- summary(fw) # got shape and scale for the weibull distribution
result$estimate[2]*(((result$estimate[1]-1)/result$estimate[1])^(1/result$estimate[1])) + min_value # estimating mode


#bootstrap of forward infectiousness profile to get 95% CI
library(boot)
library("fitdistrplus")
library(scales)

x.loverall.delta <- rweibull(92, 8.280231, 15.668116)

y.loverall.delta <- data.frame(list(sample=x.loverall.delta))

air.fun.loverall.delta <- function(data) {
  fln.delta <- fitdist(data$sample, "weibull")
  par1.delta <- summary(fln.delta)$estimate[1]
  par2.delta <- summary(fln.delta)$estimate[2]
  dist.stat.fx = function(par.1, par.2) {
    alpha = par.1
    beta  = par.2
    mode1 = beta*(((alpha-1)/alpha)^(1/alpha)) - 16.40243
    stat = c("mode1" = mode1)
    stat = round(stat, 2)
    return(stat)
  }
  
  c(par1.delta,par2.delta,
    dist.stat.fx(par1.delta,par2.delta)[1])
}

air.rg.loverall.delta <- function(data, mle) {
  
  out.delta <- data
  out.delta$sample <- rweibull(nrow(out.delta),mle[1],mle[2])
  out.delta
}

air.boot.loverall.delta <- boot(y.loverall.delta, air.fun.loverall.delta, R = 1000, sim = "parametric",
                                ran.gen = air.rg.loverall.delta,mle = c(8.280231, 15.668116))
boot.ci(boot.out=air.boot.loverall.delta,index = 3,type = "perc")





#draw figure 1


par(mfrow=c(3,1))

## plot latent period
dataset_latent <- data.frame(shape_latent=alpha_latentW,scale_latent=beta_latentW)
dataset_latent <- dataset_latent[sample(nrow(dataset_latent),100),]


mean_latent_delta <- dweibull(3.9, 1.56, 4.33)

percentile_latent_delta <- dweibull(8.8 , 1.56, 4.33)

latent_delta <- dweibull(seq(0,21,0.01), 1.56, 4.33)


xax <- seq(0,21,0.01)
plot(c(0,21), c(0,0.25),cex.axis=1.5,cex.lab=1.5,xaxt="n",yaxs="i",xaxs="i",
     type = "n",xlab = "Latent period (days)",ylab = "Density")
axis(1,cex.axis=1.5)


for (i in 1:100) {
  latent_delta_1 <- dweibull(seq(0,21,0.01), dataset_latent$shape_latent[i], dataset_latent$scale_latent[i])
  
  lines(xax,latent_delta_1,col = adjustcolor("rosybrown1", alpha.f =0.2),lwd=1,lty=1)
  
}

lines(xax,latent_delta,col = "#F8766D",lwd=2,lty=1)

lines(3.9,mean_latent_delta,col = "#F8766D",lwd=2,lty=3,type="h",)
points(8.8,percentile_latent_delta,type="p", lwd=2,pch=2, col = "#F8766D", bg=NA, cex=1.5)


title("A: Latent period distribution",line = -2,adj=0.9)


#### plot Incubation period

dataset_incu <- data.frame(shape_incu=alpha_incuW,scale_incu=beta_incuW)
dataset_incu <- dataset_incu[sample(nrow(dataset_incu),100),]


mean_incuba_delta <- dweibull(5.8, 1.98, 6.52)

percentile_incuba_delta <- dweibull(11.3, 1.98, 6.52)

incuba_delta <- dweibull(seq(0,21,0.01), 1.98, 6.52)

xax <- seq(0,21,0.01)
plot(c(0,21), c(0,0.2),cex.axis=1.5,cex.lab=1.5,xaxt="n",yaxs="i",xaxs="i",
     type = "n",xlab = "Incubation period (days)",ylab = "Density")
axis(1,cex.axis=1.5)



for (i in 1:100) {
  incuba_delta_1 <- dweibull(seq(0,21,0.01), dataset_incu$shape_incu[i], dataset_incu$scale_incu[i])
  
  lines(xax,incuba_delta_1,col = adjustcolor("rosybrown1", alpha.f =0.2),lwd=1,lty=1)
  
}

lines(xax,incuba_delta,col = "#F8766D",lwd=2,lty=1)

lines(5.8,mean_incuba_delta,col = "#F8766D",lwd=2,lty=3,type="h",)
points(11.3,percentile_incuba_delta,type="p", lwd=2,pch=2, col = "#F8766D", bg=NA, cex=1.5)


title("B: Incubation period distribution",line = -2,adj=0.93)




## ## plot infectiousness profile
mean(xiIW)
dataset_infec <- data.frame(shape_infec=alpha_infecW,scale_infec=beta_infecW, minvalue=xiIW)
dataset_infec <- dataset_infec[sample(nrow(dataset_infec),100),]


mode_infect_delta <- dweibull(9.90789, 5.13,  10.38)

infect_delta <- dweibull(seq(0,40,0.01), 5.13,  10.38)

xax <- seq(0,40,0.01)


plot(c(-12,12), c(0,0.25),cex.axis=1.5,cex.lab=1.5,xaxt="n",yaxs="i",xaxs="i",
     type = "n",xlab = "Days after symptom onset (days)",ylab = "Density")
axis(1,at=c(-12,-8,-4,0,4,8,12),cex.axis=1.5)


for (i in 1:100) {
  infect_delta_1 <- dweibull(seq(0,40,0.01), dataset_infec$shape_infec[i], dataset_infec$scale_infec[i])
  
  lines(xax+dataset_infec$minvalue[i],infect_delta_1,col = adjustcolor("rosybrown1", alpha.f =0.2),lwd=1,lty=1)
  
}


lines(xax-11.20789,infect_delta,col = "#F8766D",lwd=2,lty=1)
lines(-1.3,mode_infect_delta,col = "#F8766D",lwd=2,lty=3,type="h",)


title("C: Infectiousness profile",line = -2,adj=0.9)
