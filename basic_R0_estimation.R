library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(tikzDevice)

##simulating time-varying R0 and got the initial R0

nboot <- 200

r1 <- as.numeric(0.3291696)

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
rr <- read.csv("delta.csv", colClasses = "character")
rr$serial1 <- as.numeric(rr$serial1)
rr$serial <- as.numeric(rr$serial)
rr$serial_infector <- as.numeric(rr$serial_infector)
rr <- rr[rr$serial == 1,]

tcut <- seq(1, 16, by=1)
serdata1 <- lapply(tcut, function(x) {
  set.seed(123)
  tmp <- rr %>%
    filter(serial_infector <= x)
  
  R0samp <- replicate(nboot, {
    ss <- tmp$serial1
    ss <- sample(ss, replace=TRUE)
    1/mean(exp(-r1*ss), na.rm=TRUE)
  })
  
  data.frame(
    t=x,
    R0=1/mean(exp(-r1*tmp$serial1), na.rm=TRUE),
    lwr=quantile(R0samp, 0.025),
    upr=quantile(R0samp, 0.975)
  )
}) %>%
  bind_rows
serdata1
