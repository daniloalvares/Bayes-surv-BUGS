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

# install.packages("smcure", dep = TRUE)
library("smcure")
data("bmt")
str(bmt)

XC <- model.matrix(~ TRT, data = bmt) # Reference = allogeneic
XU <- model.matrix(~ TRT, data = bmt)
XU <- matrix(XU[,-1], ncol = 1) # Remove intercept


d.jags <- list(n = nrow(bmt), t = bmt$Time,
  XC = XC, XU = XU, delta = bmt$Status, zeros = rep(0,nrow(bmt)),
  NbetasC = ncol(XC), NbetasU = ncol(XU))

i.jags <- function() {
  list(betaC = rnorm(ncol(XC)), betaU = rnorm(ncol(XU)),
    lambda = runif(1), alpha = runif(1))
}

p.jags <- c("betaC", "betaU", "alpha", "lambda")

# install.packages("rjags", dep = TRUE)
library(rjags)
m3 <- jags.model(data = d.jags, file = "Cure.txt", inits = i.jags, n.chains = 3)
update(m3, 10000)
res <- coda.samples(m3, variable.names = p.jags, n.iter = 100000, n.thin = 100)
summary(res)

result <- as.mcmc(do.call(rbind, res))

alpha <- result[, 1]
betaC1 <- result[, 2]
betaC2 <- result[, 3]
betaU <- result[, 4]
lambda <- result[, 5]

table <- rbind(betaC1, betaC2, betaU, alpha, lambda)
table <- round(cbind(summary_parameters(table), c(mean(betaC1 > 0), mean(betaC2 > 0),
         mean(betaU > 0), 1, 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("betaC1", "betaC2", "betaU", "alpha", "lambda")
table

CP.allo <- exp(betaC1) / (1 + exp(betaC1))
CP.auto <- exp(betaC1 + betaC2) / (1 + exp(betaC1 + betaC2))
summary(cbind(CP.allo, CP.auto))

grid <- 100
time <- seq(0, max(bmt$Time), len = grid)
surv.allo <- vector(); surv.auto <- vector()
for(l in 1:grid) {
  surv.allo[l] <- mean(exp(-lambda * time[l]^alpha))
  surv.auto[l] <- mean(exp(-lambda * exp(betaU) * time[l]^alpha))
}

# install.packages("ggplot2", dep = TRUE)
library("ggplot2")
treat.colour <- rep(0:1, each = grid)
treat.colour[treat.colour == 0] <- "allogeneic"
treat.colour[treat.colour == 1] <- "autologous"
df <- data.frame(time = rep(time, 2), survival = c(surv.allo, surv.auto), treatment = treat.colour)

ggplot(data = df, aes(x = time, y = survival, group = treatment, colour = treatment)) +
    geom_line() + theme_bw() + theme(legend.position = "top")
