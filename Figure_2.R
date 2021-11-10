
library(scales)

## drawing Figure 2

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
incuba <- read.csv("serial_r_trend1026.csv")
r <- as.numeric(incuba$date)
lowerbound_si <- as.numeric(incuba$lower_mean_si)
upperbound_si <- as.numeric(incuba$upper_mean_si)
mean_si <- as.numeric(incuba$mean_si)

lowerbound_rt <- as.numeric(incuba$lower_mean_r)
upperbound_rt <- as.numeric(incuba$upper_mean_r)
mean_rt <- as.numeric(incuba$mean_r)

lowerbound_peak <- as.numeric(incuba$lower_peak)
upperbound_peak <- as.numeric(incuba$upper_peak)
mean_peak <- as.numeric(incuba$mean_peak)

lowerbound_ro <- as.numeric(incuba$lower_mean_ro)
upperbound_ro <- as.numeric(incuba$upper_mean_ro)
mean_ro <- as.numeric(incuba$mean_ro)

layout(mat = matrix(c(1,1,2,3,4,5),nrow = 3, ncol = 2, byrow = TRUE))

# A: barplot for epidemic curve
bar <- barplot(incuba$case,ylim = c(0,16),xlim = c(0,35.5),col = alpha("grey"),
               space = 0,xlab = "Date (year of 2021)",ylab = "Number of cases",cex.axis=1.5,cex.lab=1.5)

axis(side = 1,at = seq(0.5,35.5,5),
     labels = c("18 May","23 May","28 May","02 Jun","07 Jun","12 Jun","17 Jun","22 Jun"),cex.axis=1.5)
title("A",line = -2,adj=0.1,cex.main = 2)

# B: time-varying serial interval
plot(c(0,35.5), c(3,7.5),cex.axis=1.5,cex.lab=1.5,frame.plot = 0,xaxt = "n",
     type = "n",xlab = "Date (year of 2021)",ylab = "Daily estimate of forward serial interval (days)")
segments(r,lowerbound_si,r,upperbound_si)
axis(side = 1,at=c(0,5,10,15,20,25,30,35),
     labels = c("18 May","23 May","28 May","02 Jun","07 Jun","12 Jun","17 Jun","22 Jun"),cex.axis=1.5)
points(r,mean_si,col="black",lwd=2,type = "p",lty = 2,pch =19,cex=1.5)
segments(x0=-1.5,x1=35,y0=6.3,y1=6.3,lty=2,col="red",lwd=2)

title("B",line = -2,adj=0.1,cex.main = 2)

# C: time-varying Rt
plot(c(0,35), c(0,16),cex.axis=1.5,cex.lab=1.5,frame.plot = 0,xaxt = "n",
     type = "n",xlab = "Date (year of 2021)",ylab = "Daily estimate of Rt")
polygon(c(r,rev(r)),c(lowerbound_rt,rev(upperbound_rt)),col="thistle",border=NA)
axis(side = 1,at=c(0,5,10,15,20,25,30,35),
     labels = c("18 May","23 May","28 May","02 Jun","07 Jun","12 Jun","17 Jun","22 Jun"),cex.axis=1.5)
lines(r,mean_rt,col="black",lwd=2)
segments(x0=-1.5,x1=35,y0=1,y1=1,lty=3,col="red",lwd=2)

title("C",line = -2,adj=0.1,cex.main = 2)

# D: time-varying infectiousness profile
plot(c(0,35), c(-2.0,1.0),cex.axis=1.5,cex.lab=1.5,frame.plot = 0,xaxt = "n",
     type = "n",xlab = "Date (year of 2021)",ylab = "Shift in infectiousness peak (days)")
segments(r,lowerbound_peak,r,upperbound_peak)
axis(side = 1,at=c(0,5,10,15,20,25,30,35),
     labels = c("18 May","23 May","28 May","02 Jun","07 Jun","12 Jun","17 Jun","22 Jun"),cex.axis=1.5)
points(r,mean_peak,col="black",lwd=2,type = "p",lty = 2,pch =19,cex=1.5)
title("D",line = -2,adj=0.1,cex.main = 2)


# E: time-varying R0
plot(c(0,35), c(0,16),cex.axis=1.5,cex.lab=1.5,frame.plot = 0,xaxt = "n",
     type = "n",xlab = "Date (year of 2021)",ylab = "Daily estimate of R0")
polygon(c(r,rev(r)),c(lowerbound_ro,rev(upperbound_ro)),col="thistle",border=NA)
axis(side = 1,at=c(0,5,10,15,20,25,30,35),
     labels = c("18 May","23 May","28 May","02 Jun","07 Jun","12 Jun","17 Jun","22 Jun"),cex.axis=1.5)
lines(r,mean_ro,col="black",lwd=2)
segments(x0=-1.5,x1=35,y0=4.9,y1=4.9,lty=3,col="red",lwd=2)

title("E",line = -2,adj=0.1,cex.main = 2)


