#---------------------------------------------------------------#
#                       USEFUL FUNCTIONS                        #
#---------------------------------------------------------------#
summary_parameters <- function(table) {
  table <- cbind(apply(table, 1, mean), apply(table, 1, sd),
    t(apply(table, 1, quantile, probs = c(0.025, 0.5, 0.975))))
    colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%")
  return(table)
}

fk <- function(u.vect, lambda, alpha, beta, x, k) {
  res <- sapply(u.vect, function(u) {
    # Cause-specific hazard
    hk <- lambda[k] * alpha[k] * (u^(alpha[k] - 1)) * exp(sum(unlist(beta[, k]) * x))
    # Cumulative cause-specific hazard
    Hk <- lambda * (rep(u,length(lambda))^alpha) * exp((t(beta) %*% matrix(x, ncol = 1))[, 1])
    # Cause-specific hazard x Overall survival
    aux <- hk * exp(-sum(Hk))
    return(aux)
  })
  return(res)
}

# install.packages("pracma")
library("pracma")
cif <- function(tt, lambda, alpha, beta, x, k) {
  return(integral(fk, xmin = 0, xmax = tt, method = "Simpson", lambda = lambda,
    alpha = alpha, beta = beta, x = x, k = k))
}

library("parallel")
options(mc.cores = detectCores())
mcmc_cif <- function(obj, t.pred, x) {
  # Put all chains together
  var.names <- names(obj)
  # Indices of alpha's, beta's, and lambda's
  a.idx <- which(substr(var.names, 1, 5) == "alpha")
  b.idx <- which(substr(var.names, 1, 4) == "beta")
  l.idx <- which(substr(var.names, 1, 6) == "lambda")
  # Number of causes
  K <- length(a.idx)
  # Number of covariates
  n.b <- length(b.idx) / K
  # Sub-sample to speed up computations
  samples.idx <- sample(1:nrow(obj), 200)
  
  res <- lapply(1:K, function(k) {
    sapply(t.pred, function(tt) {
      aux <- mclapply(samples.idx, function(i) {
        cif(tt, alpha = unlist(obj[i, a.idx]), lambda = unlist(obj[i, l.idx]),
            beta = matrix(unlist(c(obj[i, b.idx])), nrow = n.b), x = x, k = k)
      })
      return(mean(unlist(aux)))
    })
  })
  return(res)
}
#---------------------------------------------------------------#

# install.packages("compeir", dep = TRUE)
library("compeir")
data("okiss")
str(okiss)

delta <- matrix(c(as.integer(okiss$status == 1), 
   as.integer(okiss$status == 2), as.integer(okiss$status == 7)), ncol = 3)
head(delta)

X <- model.matrix(~ allo + sex, data = okiss) # Reference = female
X <- X[, -1] # Remove intercept


d.jags <- list(n = nrow(okiss), t = as.vector(okiss$time),
  X = X, delta = delta, zeros = rep(0, nrow(okiss)),
  Nbetas = ncol(X), Nrisks = ncol(delta))

i.jags <- function() {
  list(beta = matrix(rnorm(ncol(X) * ncol(delta)), ncol = ncol(delta)),
    lambda = runif(ncol(delta)), alpha = runif(ncol(delta)))
}

p.jags <- c("beta", "alpha", "lambda")

# install.packages("rjags", dep = TRUE)
library("rjags")
m4 <- jags.model(data = d.jags, file = "CR.txt", inits = i.jags, n.chains = 3)
update(m4, 1000)
res <- coda.samples(m4, variable.names = p.jags, n.iter = 10000, n.thin = 10)

result <- as.mcmc(do.call(rbind, res))

alpha1 <- result[, 1]; alpha2 <- result[, 2]; alpha3 <- result[, 3]
beta11 <- result[, 4]; beta21 <- result[, 5]
beta12 <- result[, 6]; beta22 <- result[, 7]
beta13 <- result[, 8]; beta23 <- result[, 9]
lambda1 <- result[, 10]; lambda2 <- result[, 11]; lambda3 <- result[, 12]

table <- rbind(beta11, beta21, lambda1, alpha1, beta12, beta22, lambda2, alpha2, beta13, beta23, lambda3, alpha3)
table <- round(cbind(summary_parameters(table), c(mean(beta11 > 0), mean(beta21 > 0), 1, 1,
  mean(beta12 > 0), mean(beta22 > 0), 1, 1,
  mean(beta13 > 0), mean(beta23 > 0), 1, 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("beta11", "beta21", "lambda1", "alpha1", "beta12", "beta22",
  "lambda2", "alpha2", "beta13", "beta23", "lambda3", "alpha3")
table

res2 <- as.data.frame(result)
t.pred <- seq(0, 100, by = 2.5)
cum_inc <- mcmc_cif(res2, t.pred, c(1, 1))

# install.packages("ggplot2", dep = TRUE)
library("ggplot2")
df <- data.frame(cif = unlist(cum_inc), time = t.pred, 
  cause = rep(c("infection", "end of neutropenia", "death"), each = length(cum_inc[[1]])))

ggplot(data = df, aes(x = time, y = cif, group = cause)) + geom_line(aes(color = cause)) + 
  ylab("cumulative incidence") + ylim(c(0,1)) + theme_bw() + theme(legend.position = "top")
