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

# install.packages("frailtyHL", dep = TRUE)
library("frailtyHL")
data("kidney")
str(kidney)

# Number of patients and catheters
n <- length(unique(kidney$id))
J <- 2
# Survival and censoring times
time <- kidney$time
cens <- time
time[kidney$status == 0] <- NA # Censored
is.censored <- as.numeric(is.na(time))

# Matrix format
time <- matrix(time, n, J, byrow = TRUE)
cens <- matrix(cens, n, J, byrow = TRUE)
is.censored <- matrix(is.censored, n, J, byrow = TRUE)

sex <- kidney$sex[seq(1, 2 * n, 2)] - 1 # Reference = male
X <- model.matrix(~ sex)


d.jags <- list(n = n, J = J, time = time, cens = cens, X = X,
  is.censored = is.censored, Nbetas = ncol(X))

i.jags <- function(){
  list(beta = rnorm(ncol(X)), alpha = runif(1), psi = runif(1))
}

p.jags <- c("beta", "alpha", "lambda", "psi", "w")

# install.packages("rjags")
library("rjags")
m6 <- jags.model(data = d.jags, file = "Frailty.txt", inits = i.jags, n.chains = 3)
update(m6, 10000)
res <- coda.samples(m6, variable.names = p.jags, n.iter = 100000, thin = 100)

result <- as.mcmc(do.call(rbind, res))

alpha <- result[, 1]
beta2 <- result[, 3]
lambda <- result[, 4]
psi <- result[, 5]
w <- result[, 6:ncol(result)]

table <- rbind(beta2, alpha, lambda, psi)
table <- round(cbind(summary_parameters(table), c(mean(beta2 > 0), 1, 1, 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("beta2", "alpha", "lambda", "psi")
table

grid <- 1000
time <- seq(0, max(kidney$time), len = grid)
surv <- matrix(NA, n, grid) 
for(i in 1:n){
  for(k in 1:grid){
    surv[i, k] <- mean( exp(-w[i] * lambda * exp(beta2 * sex[i]) * time[k]^alpha) )
  }
}  

# install.packages("ggplot2")
library("ggplot2")
sex.colour <- sex
sex.colour[sex == 0] <- "male"
sex.colour[sex == 1] <- "female"
df <- data.frame(time = rep(time,n), survival = c(t(surv)), 
  patient = rep(1:n, each = grid), sex = rep(sex.colour, each = grid))

ggplot(data = df, aes(x = time, y = survival, group = patient, colour = sex)) + 
  geom_line() + theme_bw() + theme(legend.position = "top")
