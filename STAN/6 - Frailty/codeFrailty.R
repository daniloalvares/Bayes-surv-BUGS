# install.packages("frailtyHL", dep = TRUE)
library("frailtyHL")
data("kidney")
str(kidney)

# Number of patients and catheters
n <- length(unique(kidney$id))
J <- 2

# Covariates
sex <- kidney$sex[seq(1, 2 * n, 2)] - 1 # Reference = male
X <- model.matrix(~ sex)

# Survival information
time <- kidney$time
delta <- kidney$status
# Matrix format
time <- matrix(time, n, J, byrow = TRUE)
delta <- matrix(delta, n, J, byrow = TRUE)


library("rstan")
m6 <- stan(file   = "Frailty.stan", 
           data   = list(n = n, J = J, Nbetas = ncol(X), 
                         time = time, delta = delta, X = X),
           warmup = 1000,                 
           iter   = 2000,
           chains = 3,
           seed   = 6,
           cores  = getOption("mc.cores",3)) 

print(m6)


# Diagnostics
library("shinystan")
launch_shinystan(m6)

pars <- c("beta[2]", "alpha", "lambda", "psi")
plot(m6, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m6, plotfun = "hist", pars = pars)

post6 <- extract(m6, pars)

hist(post6$beta)
hist(post6$alpha)
