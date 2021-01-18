#---------------------------------------------------------------#
#                       USEFUL FUNCTION                         #
#---------------------------------------------------------------#
summary_parameters <- function(table) {
  table <- cbind(apply(table, 1, mean), apply(table, 1, sd),
    t(apply(table, 1, quantile, probs = c(0.025, 0.5, 0.975))))
    colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%")
  return(table)
}
#---------------------------------------------------------------#

# install.packages("JMbayes", dep = TRUE)
library("JMbayes")
data("prothro")
data("prothros")
str(prothro); str(prothros)

# Number of patients and number of longitudinal observations per patient
n <- nrow(prothros)
M <- table(prothro$id)
# Survival and censoring times
Time <- prothros$Time
death <- prothros$death

# Longitudinal information in matrix format
time <- matrix(NA, n, max(M))
log.proth <- matrix(NA, n, max(M))
count <- 1
for(i in 1:n) {
  log.proth[i, 1:M[i]] <- log(prothro$pro[count:(M[i] + count - 1)])
  time[i, 1:M[i]] <- prothro$time[count:(M[i] + count - 1)]
  count <- count + M[i]
}

treat <- as.numeric(prothros$treat) - 1 # Reference = placebo
XS <- model.matrix(~ treat) # Fixed effects

XL <- array(1, dim = c(n, max(M), 3)) # Fixed effects
XL[, , 2] <- time
XL[, , 3] <- treat
ZL <- array(1, dim = c(n, max(M), 2)) # Random effects
ZL[, , 2] <- time

# Gauss-Legendre quadrature (15 points)
# install.packages("statmod")
library("statmod")
glq <- gauss.quad(15, kind = "legendre")
xk <- glq$nodes   # nodes
wk <- glq$weights # weights
K <- length(xk)   # K-points


d.jags <- list(n = n, M = M, Time = Time, XS = XS, log.proth = log.proth, XL = XL, ZL = ZL, 
  death = death, mub = rep(0, 2), V = diag(1, 2), Nb = 2, zeros = rep(0, n),
  NbetasL = dim(XL)[3], NbetasS = ncol(XS), K = length(xk), xk = xk, wk = wk)

i.jags <- function() {
  list(betaS = rnorm(ncol(XS)), gamma = rnorm(1), alpha = runif(1),
    betaL = rnorm(dim(XL)[3]), tau = runif(1), Omega = diag(runif(2)))
}

p.jags <- c("betaS", "gamma", "alpha", "lambda", "betaL", "sigma", "Sigma", "b")

# install.packages("rjags", dep = TRUE)
library(rjags)
m7 <- jags.model(data = d.jags, file = "JM.txt", inits = i.jags, n.chains = 3)
update(m7, 1000)
res <- coda.samples(m7, variable.names = p.jags, n.iter = 10000, thin = 10)

result <- as.mcmc(do.call(rbind, res))

Sigma2.11 <- result[, 1]
Sigma2.12 <- result[, 2]
Sigma2.22 <- result[, 4]
alpha <- result[, 5]
b1 <- result[, 6:(n + 5)]
b2 <- result[, (n + 6):(2 * n + 5)]
betaL1 <- result[, (2 * n + 6)]
betaL2 <- result[, (2 * n + 7)]
betaL3 <- result[, (2 * n + 8)]
betaS2 <- result[, (2 * n + 10)]
gamma <- result[, (2 * n + 11)]
lambda <- result[, (2 * n + 12)]
sigma <- result[, (2 * n + 13)]

table <- rbind(betaS2, gamma, lambda, alpha, betaL1, betaL2, betaL3, sigma, Sigma2.11, Sigma2.12, Sigma2.22)
table <- round(cbind(summary_parameters(table), c(mean(betaS2 > 0), mean(gamma > 0), 1, 1,
  mean(betaL1 > 0), mean(betaL2 > 0), mean(betaL3 > 0), 1, 1, mean(Sigma2.12 > 0), 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("betaS2", "gamma", "lambda", "alpha", "betaL1", "betaL2", "betaL3", "sigma", "Sigma2.11", "Sigma2.12", "Sigma2.22")
table
