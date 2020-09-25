setwd("~/Desktop/Research/Ordinal BUCP Simulation/Ordinal-BUCP-master/Real data analysis")
library(runjags)
source('plot_SSD.R')
source("fit_OBUCP.R")



data <- read.csv("Lambert.csv")
obs <- c(table(data[, 1]))
maxobs <- max(obs)
cums <- c(0, cumsum(obs))
redata <- matrix(NA, 4, maxobs)
for (i in 1:4){
  redata[i, (1: obs[i])] <- data[(cums[i] + 1):cums[i + 1], 2]
}
plot_SSD(redata)

res.mat <- matrix(NA, 3, 21)
colnames(res.mat) <- c("phase", "log.mean1.data", "log.mean2.data",
                       "rho.l95", "rho.u95","rho.mean", "rho.sd", 
                       "brr.l95", "brr.u95", "brr.mean", "brr.sd",
                       "cp.l95", "cp.median", "cp.u95", "cp.mean", "cp.sd", "cp.mode",
                       "mu1.mean.est", "mu2.mean.est", "sd1.mean.est", "sd2.mean.est")

for(i in 1:3) {
  tempdata <- c(na.omit(c(redata[i, ], redata[(i + 1),])))
  results <- fit_OBUCP(tempdata)
  x <- as.data.frame(summary(results))
  res.mat[i,] <- as.numeric(c(i,
                              log(mean(redata[i,], na.rm = T)),
                              log(mean(redata[(i + 1), ], na.rm = T)),
                              x["rho",c(1, 3:5)],
                              x["brr",c(1, 3:5)], 
                              x["cp", 1:6], x["mu[1]", 4],
                              x["mu[2]", 4], x["sigma[1]", 4], 
                              x["sigma[2]", 4]))
  
}
write.csv(res.mat, "Lambert-results.csv")
