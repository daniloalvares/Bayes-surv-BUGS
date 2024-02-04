# install.packages("KMsurv", dep = TRUE)
library("KMsurv")
data("larynx")
str(larynx)

# Covariates
larynx$age <- as.numeric(scale(larynx$age))
larynx$diagyr <- as.numeric(scale(larynx$diagyr))
larynx$stage <- as.factor(larynx$stage)
X <- model.matrix(~ stage + age + diagyr, data = larynx)[,-1] # Remove intercept

# Survival information
time <- larynx$time
delta <- larynx$delta

# Time axis partition 
K <- 3 # number of intervals  
a <- seq(0, max(larynx$time) + 0.001, length.out = K + 1)

# intobs: vector that tells us at which interval each observation is
intobs <- matrix(data = NA, nrow = nrow(larynx), ncol = length(a) - 1)
d <- matrix(data = NA, nrow = nrow(larynx), ncol = length(a) - 1)

for(i in 1:nrow(larynx)){
  for(k in 1:(length(a) - 1)){
    d[i, k] <- ifelse(time[i] - a[k] > 0, 1, 0) * ifelse(a[k + 1] - time[i] > 0, 1, 0)
    intobs[i, k] <- d[i, k] * k
  }
}
intobs <- rowSums(intobs)


library("rstan")
m2 <- stan(file   = "PH.stan", 
           data   = list(n = length(time), m = length(a) - 1, Nbetas = ncol(X), 
                         time = time, delta = delta, X = X, a = a, intobs = intobs),
           warmup = 200,                 
           iter   = 500,
           chains = 3,
           seed   = 2,
           cores  = getOption("mc.cores",3)) 

print(m2)


# Diagnostics
library("shinystan")
launch_shinystan(m2)

pars <- c("beta", "lambda")
plot(m2, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m2, plotfun = "hist", pars = pars)

post2 <- extract(m2, pars)

hist(post2$beta[,1])
hist(post2$lambda)
