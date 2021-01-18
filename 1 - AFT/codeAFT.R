# install.packages("KMsurv", dep = TRUE)
library("KMsurv")
data("larynx")
str(larynx)

# Covariates
larynx$age <- as.numeric(scale(larynx$age))
larynx$diagyr <- as.numeric(scale(larynx$diagyr))
larynx$stage <- as.factor(larynx$stage)
X <- model.matrix(~ stage + age + diagyr, data = larynx)
un <- which(larynx$delta == 1)

# Uncensored observations
time_un <- larynx$time[un]
cen_un <- larynx$time[un]
is.censored_un <- rep(0, length(cen_un))
X_un <- X[un,]

# Right-censored observations
cen_cen <- larynx$time[-un]
time_cen <- rep(NA, length(cen_cen))
is.censored_cen <- rep(1, length(cen_cen))
X_cen <- X[-un,]

# Uncensored and right-censored observations
time <- c(time_un, time_cen)
cen <- c(cen_un, cen_cen)
is.censored <- c(is.censored_un, is.censored_cen)
Xnew <- rbind(X_un, X_cen)


d.jags <- list(n = nrow(X), time = time , cen = cen, X = Xnew, is.censored = is.censored, Nbetas = ncol(X))

i.jags <- function(){
  list(beta = rnorm(ncol(X)), alpha = runif(1))
}

p.jags <- c("beta", "alpha")

# install.packages("rjags", dep = TRUE)
library("rjags")
m1 <- jags.model(data = d.jags, file = "AFT.txt", inits = i.jags, n.chains = 3)
update(m1, 1000)
set.seed(1)
res <- coda.samples(m1, variable.names = p.jags, n.iter = 50000, thin = 10)
summary(res)

par(mfrow = c(2,4))
traceplot(res)

gelman.diag(res)

par(mfrow = c(2,4))
densplot(res, xlab = "")

result <- as.mcmc(do.call(rbind, res))

alpha <- result[, 1]
beta1 <- result[, 2]; beta2 <- result[, 3]
beta3 <- result[, 4]; beta4 <- result[, 5]
beta5 <- result[, 6]; beta6 <- result[, 7]

RM.s3_s4 <- exp(beta3 - beta4)
summary(RM.s3_s4)
