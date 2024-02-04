# install.packages("KMsurv", dep = TRUE)
library("KMsurv")
data("larynx")
str(larynx)

# Covariates
larynx$age <- as.numeric(scale(larynx$age))
larynx$diagyr <- as.numeric(scale(larynx$diagyr))
larynx$stage <- as.factor(larynx$stage)
X <- model.matrix(~ stage + age + diagyr, data = larynx)

# Survival information
time <- larynx$time
delta <- larynx$delta


library("rstan")
m1 <- stan(file   = "AFT.stan", 
           data   = list(n = length(time), Nbetas = ncol(X), 
                         time = time, delta = delta, X = X),
           warmup = 200,                 
           iter   = 500,
           chains = 3,
           seed   = 1,
           cores  = getOption("mc.cores",3)) 

print(m1)


# Diagnostics
library("shinystan")
launch_shinystan(m1)

pars <- c("beta", "sigma")
plot(m1, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m1, plotfun = "hist", pars = pars)

post1 <- extract(m1, pars)

hist(post1$beta[,1])
hist(post1$sigma)
