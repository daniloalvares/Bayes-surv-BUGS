# install.packages("p3state.msm", dep = TRUE)
library("p3state.msm")
data("heart2")
str(heart2)

# Covariates
X <- model.matrix(~ age + year + surgery, data = heart2)[,-1]

# Survival information
time3 <- heart2$times2
time3[time3 == 0] <- 0.0001
time <- matrix(c(heart2$times1, heart2$time, time3), ncol = 3)
delta <- matrix(c(heart2$delta, 
                  heart2$status * (1 - heart2$delta), 
                  heart2$delta * heart2$status), ncol = 3)


library("rstan")
m5 <- stan(file   = "IllDeath.stan", 
           data   = list(n = nrow(time), Ntrans = ncol(delta), Nbetas = ncol(X), 
                         time = time, delta = delta, X = X),
           warmup = 200,                 
           iter   = 500,
           chains = 3,
           seed   = 5,
           cores  = getOption("mc.cores",3)) 

print(m5)


# Diagnostics
library("shinystan")
launch_shinystan(m5)

pars <- c("beta", "lambda", "alpha")
plot(m5, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m5, plotfun = "hist", pars = pars)

post5 <- extract(m5, pars)

hist(post5$beta[,1,1])
hist(post5$lambda[,1])
