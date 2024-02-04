# install.packages("JMbayes", dep = TRUE)
library("JMbayes")
data("prothro")
data("prothros")
str(prothro); str(prothros)

# Covariates
treat <- as.numeric(prothros$treat) - 1 # Reference = placebo
X <- model.matrix(~ treat) # Fixed effects

# Longitudinal information
y <- log(prothro$pro)
tt <- prothro$time
N <- length(y)
M <- as.numeric(table(prothro$id))

# Survival information
time <- prothros$Time
delta <- prothros$death
n <- length(time)

# Gauss-Legendre quadrature (15 points)
# install.packages("statmod")
library("statmod")
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
K <- length(xk)   # K-points


library("rstan")
m7 <- stan(file   = "JM.stan", 
           data   = list(N = N, n = n, NbetasL = ncol(X)+1, NbetasS = ncol(X), Nb = 2,
                         y = y, tt = tt, ID = rep(1:n,M), time = time, delta = delta, 
                         X = X, K = K, xk = xk, wk = wk),
           warmup = 500,                 
           iter   = 1000,
           chains = 3,
           seed   = 7,
           cores  = getOption("mc.cores",3)) 

print(m7)


# Diagnostics
library("shinystan")
launch_shinystan(m7)

pars <- c("betaL", "betaS", "gamma", "sigma", "Sigma", "lambda", "alpha")
plot(m7, plotfun = "trace", pars = pars, inc_warmup = TRUE)
plot(m7, plotfun = "hist", pars = pars)

post7 <- extract(m7, pars)

hist(post7$betaL[,1])
hist(post7$gamma)
