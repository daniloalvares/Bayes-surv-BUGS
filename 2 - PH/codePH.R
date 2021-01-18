#---------------------------------------------------------------#
#                       USEFUL FUNCTION                         #
#---------------------------------------------------------------#
summary_parameters <- function(table){
  table <- cbind(apply(table, 1, mean), apply(table, 1, sd),
    t(apply(table, 1, quantile, probs = c(0.025, 0.5, 0.975))))
    colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%")
  return(table)
}
#---------------------------------------------------------------#

# install.packages("KMsurv", dep = TRUE)
library("KMsurv")
data("larynx")
str(larynx)

time <- larynx$time
larynx$age <- as.numeric(scale(larynx$age))
larynx$diagyr <- as.numeric(scale(larynx$diagyr))
larynx$stage <- as.factor(larynx$stage)
X <- model.matrix(~ -1 + stage + age + diagyr, data = larynx)

# Time axis partitition 
K <- 3 # number of intervals  
a <- seq(0, max(larynx$time) + 0.001, length.out = K + 1)

# int.obs: vector that tells us at which interval each observation is
int.obs <- matrix(data = NA, nrow = nrow(larynx), ncol = length(a) - 1)
d <- matrix(data = NA, nrow = nrow(larynx), ncol = length(a) - 1)

for(i in 1:nrow(larynx)){
  for (k in 1:(length(a) - 1)){
    d[i, k] <- ifelse(time[i] - a[k] > 0, 1, 0) * ifelse(a[k + 1] - time[i] > 0, 1, 0)
    int.obs[i, k] <- d[i, k] * k
  }
}
int.obs <- rowSums(int.obs)


d.jags <- list(n = nrow(larynx), m = length(a) - 1, delta = larynx$delta,
  time = larynx$time, X = X[, -1], a = a, int.obs = int.obs, Nbetas = ncol(X) - 1,
  zeros = rep(0, nrow(larynx)))

i.jags <- function(){
  list(beta = rnorm(ncol(X) - 1), lambda = runif(3, 0.1))
}

p.jags <- c("beta", "lambda")

# install.packages("rjags")
library("rjags")
m2 <- jags.model(data = d.jags, file = "PH.txt", inits = i.jags, n.chains = 3)
update(m2, 1000)
res <- coda.samples(m2, variable.names = p.jags, n.iter = 50000, thin = 10)
summary(res)

# Save results to vectors
result <- as.mcmc(do.call(rbind, res))
beta2 <- result[, 1]; beta3 <- result[, 2]; beta4 <- result[, 3]; beta5 <- result[, 4]
beta6 <- result[, 5]; lambda1 <- result[, 6]; lambda2 <- result[, 7]; lambda3 <- result[, 8]

table <- rbind(beta2, beta3, beta4, beta5, beta6, lambda1, lambda2, lambda3)
table <- round(cbind(summary_parameters(table), c(mean(beta2 > 0), mean(beta3 > 0),
         mean(beta4 > 0), mean(beta5 > 0), mean(beta6 > 0), 1, 1, 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("beta2", "beta3", "beta4", "beta5", "beta6", "lambda1", "lambda2", "lambda3")
table

HR.s3_s4 <- exp(beta3 - beta4)
summary(HR.s3_s4)
