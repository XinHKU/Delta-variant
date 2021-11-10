library(EpiEstim)
library("fitdistrplus")

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
data1 <- read.csv("inciden.csv", colClasses = "character")

data1$I <- as.numeric(data1$I)
data1$dates <- as.Date(data1$dates)

R_Parametric <- estimate_R(data1,method =  "parametric_si",
                           config =  make_config(list(t_start = 28,t_end = 33,mean_si = 4.029449,std_si= 3.992959)))
R_Parametric$R
a <- R_Parametric$R[3]
a[,1]

