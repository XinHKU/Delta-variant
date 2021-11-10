library(outbreaks)
library(ggplot2)
library(incidence)

## simulating the exponential growth rate for Guangdong outbreak

setwd("C:/Users/xinhualei/Dropbox (Personal)/Xin Hualei_research work/delta")
data <- read.csv("guangdong_growth.csv", colClasses = "character")
data <- data$date

i <- incidence(data,interval = "3 days")
plot(i)
best.fit <- fit(i[1:5])
best.fit$fit

best.fit1 <- fit_optim_split(i)
best.fit1$fit


p <- plot(i, fit=best.fit,border = "grey",stack = is.null(fit),color = "darkorange")
p <- p+geom_text(label="exponential growth phase",x=as.Date("2021-05-21"),y=20,size=3)+
       geom_text(label="fitted r =0.33",x=as.Date("2021-05-22"),y=15,size=3)+ylim(0,50)
p <- p+ theme_bw()+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p <- p+ theme(axis.text.x= element_text(size=10, colour="black"))+
  theme(axis.text.y= element_text(size=12,colour="black"))+scale_x_date(date_labels = "%B %d")+
  xlab("Date") +ylab("Daily incidence")

p



