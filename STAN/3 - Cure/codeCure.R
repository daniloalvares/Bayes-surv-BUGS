# install.packages("smcure", dep = TRUE)
library("smcure")
data("bmt")
str(bmt)

# Covariates
XC <- model.matrix(~ TRT, data = bmt) # Reference = allogeneic
XU <- matrix(XC[,-1], ncol = 1) # Remove intercept

# Survival information
time <- bmt$Time
delta <- bmt$Status


library("rstan")
m3 <- stan(file   = "Cure.stan", 
           data   = list(n = length(time), Nbetas = ncol(XC), time = time, 
                         delta = delta, XC = XC, XU = XU),
           warmup = 200,                 
           iter   = 500,
           chains = 3,
           seed   = 3,
           cores  = getOption("mc.cores",3)) 

print(m3)


# Diagnostics
library("shinystan")
launch_shinystan(m3)

pars <- c("betaC", "betaU", "lambda", "alpha")
plot(m3, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m3, plotfun = "hist", pars = pars)

post3 <- extract(m3, pars)

hist(post3$betaC[,1])
hist(post3$lambda)
