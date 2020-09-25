fit_OBUCP <- function(data){
  library(runjags)
  n <- length(data)/2

model <- "model {
    x[1] <- 0
    lambda[1] <- exp(mu[1])
    y[1] ~ dpois(lambda[1])
    
    for (i in 2:T) {
    
      dummy[i] <- step(cp - i) # dummy == 1 if i in baseline phase, and dummy == 0 if i in intervention
      x[i] <- dummy[i] * mu[1] + (1 - dummy[i]) * mu[2] 
      lambda[i] <- ifelse(cp == i + 1, exp(x[i]), exp(x[i]) + rho * (y[i - 1] - lambda[i - 1]))   
      y[i] ~ dpois(lambda[i])
    }
    mu[1] ~ dnorm(1, prec[1])
    mu[2] ~ dnorm(1, prec[2])
    cp ~ dcat(pi)
    for (i in 1:P){

      prec[i] ~ dgamma(1, 1)
      sigma[i] <- 1/sqrt(prec[i])
      #prec is the epsilon
    }
  
    
  rho ~ dunif(-1, 1)
  brr <- exp(mu[2]- mu[1])
#change to 2 - 1 for cases 7-10

}"

results <- autorun.jags(
  model = model,
  data = list(y = data, T = 2*n, P = 2, 
              pi = c(rep(0, 3), rep(1/(2*n - 6),(2*n - 6)), rep(0, 3))),
  monitor = c("rho",  "mu", "sigma", "brr","cp"),
  n.chains = 4,
  startsample = 30000,
  raftery.options = FALSE,
  method = "rjparallel"
)
results$draws <- combine.mcmc(results$mcmc)
results
}