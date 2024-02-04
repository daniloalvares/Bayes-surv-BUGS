# install.packages("compeir", dep = TRUE)
library("compeir")
data("okiss")
str(okiss)

# Covariates
X <- model.matrix(~ allo + sex, data = okiss)[,-1] # Reference = female

# Survival information
time <- as.vector(okiss$time)
delta <- matrix(c(as.integer(okiss$status == 1), 
                  as.integer(okiss$status == 2), 
                  as.integer(okiss$status == 7)), ncol = 3)


library("rstan")
m4 <- stan(file   = "CR.stan", 
           data   = list(n = length(time), Nrisks = ncol(delta), Nbetas = ncol(X), 
                         time = time, delta = delta, X = X),
           warmup = 200,                 
           iter   = 500,
           chains = 3,
           seed   = 4,
           cores  = getOption("mc.cores",3)) 

print(m4)


# Diagnostics
library("shinystan")
launch_shinystan(m4)

pars <- c("beta", "lambda", "alpha")
plot(m4, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m4, plotfun = "hist", pars = pars)

post4 <- extract(m4, pars)

hist(post4$beta[,1,1])
hist(post4$lambda[,1])
