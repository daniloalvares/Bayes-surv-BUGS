#---------------------------------------------------------------#
#                       USEFUL FUNCTIONS                        #
#---------------------------------------------------------------#
summary_parameters <- function(table) {
  table <- cbind(apply(table, 1, mean), apply(table, 1, sd),
    t(apply(table, 1, quantile, probs = c(0.025, 0.5, 0.975))))
    colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%")
  return(table)
}

p11_s_t <- function(tt, ss, l1, a1, b1, l2, a2, b2, x) {
  # Cumulative hazard 1
  H1 <- l1 * exp(sum(b1 * x)) * (tt^a1 - ss^a1)
  # Cumulative hazard 2
  H2 <- l2 * exp(sum(b2 * x)) * (tt^a2 - ss^a2)
  return(exp(-H1 - H2))
}

int_p22 <- function(u, tt, ss, l1, a1, b1, l3, a3, b3, x) {
  # Cumulative distribution 1
  F1 <- 1 - exp(-l1 * exp(sum(b1 * x)) * ss^a1)
  # Hazard 1
  u1 <- ifelse(u > 0, u^(a1 - 1),0)
  h1 <- l1 * a1 * exp(sum(b1 * x)) * u1
  # Cumulative hazard 1
  H1 <- l1 * exp(sum(b1 * x)) * u^a1
  # Cumulative hazard 3
  H3 <- l3 * exp(sum(b3 * x)) * ((tt - u)^a3 - (ss - u)^a3)
  return(h1 * exp(-H1 - H3) / F1)
}

int_p12 <- function(u, tt, ss, l1, a1, b1, l2, a2, b2, l3, a3, b3, x) {
  # Hazard 1
  u1 <- ifelse(u > 0, u^(a1 - 1), 0)
  h1 <- l1 * a1 * exp(sum(b1 * x)) * u1
  # Cumulative hazard 1
  H1 <- l1 * exp(sum(b1 * x)) * (u^a1 - ss^a1)
  # Cumulative hazard 2
  H2 <- l2 * exp(sum(b2 * x)) * (u^a2 - ss^a2)
  # Cumulative hazard 3
  H3 <- l3 * exp(sum(b3 * x)) * (tt - u)^a3
  return(h1 * exp(-H1 - H2 - H3))
}

# install.packages("pracma", dep = TRUE)
library("pracma")
p22_s_t <- function(tt, ss, l1, a1, b1, l3, a3, b3, x) {
  return(integral(int_p22, xmin = 0, xmax = ss, method = "Kronrod", tt = tt, 
    ss = ss, l1 = l1, a1 = a1, b1 = b1, l3 = l3, a3 = a3, b3 = b3, x = x))
}

# install.packages("pracma")
library(pracma)
p12_s_t <- function(tt, ss, l1, a1, b1, l2, a2, b2, l3, a3, b3, x) {
  return(integral(int_p12, xmin = ss, xmax = tt, method = "Kronrod", tt = tt, 
    ss = ss, l1 = l1, a1 = a1, b1 = b1, l2 = l2, a2 = a2, b2 = b2,
    l3 = l3, a3 = a3, b3 = b3, x = x))
}

library("parallel")
options(mc.cores = detectCores())
mcmc_p11_s_t <- function(t.pred, s.pred, l1, a1, b1, l2, a2, b2, x) {
  # Sub-sample to speed up computations
  samples.idx <- sample(1:length(l1), 200)
  
  res <- sapply(t.pred, function(tt) {
    aux <- mclapply(samples.idx, function(i) {
      p11_s_t(tt, ss = s.pred, l1 = l1[i], a1 = a1[i], b1 = b1[i, ],
        l2 = l2[i], a2 = a2[i], b2 = b2[i, ], x = x)
    })
    return(mean(unlist(aux)))
  })
  return(res)
}

mcmc_p22_s_t <- function(t.pred, s.pred, l1, a1, b1, l3, a3, b3, x) {
  # Sub-sample to speed up computations
  samples.idx <- sample(1:length(l1), 200)
  
  res <- sapply(t.pred, function(tt) {
    aux <- mclapply(samples.idx, function(i) {
      p22_s_t(tt, ss = s.pred, l1 = l1[i], a1 = a1[i], b1 = b1[i, ],
        l3 = l3[i], a3 = a3[i], b3 = b3[i, ], x = x)
    })
    return(mean(unlist(aux)))
  })
  return(res)
}

mcmc_p12_s_t <- function(t.pred, s.pred, l1, a1, b1, l2, a2, b2, l3, a3, b3, x) {
  # Sub-sample to speed up computations
  samples.idx <- sample(1:length(l1), 200)
  
  res <- sapply(t.pred, function(tt) {
    aux <- mclapply(samples.idx, function(i) {
      p12_s_t(tt, ss = s.pred, l1 = l1[i], a1 = a1[i], b1 = b1[i, ],
        l2 = l2[i], a2 = a2[i], b2 = b2[i, ], l3 = l3[i], a3 = a3[i], b3 = b3[i, ], x = x)
    })
    return(mean(unlist(aux)))
  })
  return(res)
}
#---------------------------------------------------------------#

# install.packages("p3state.msm", dep = TRUE)
library("p3state.msm")
data("heart2")
str(heart2)

event <- matrix(c(heart2$delta, heart2$status * (1 - heart2$delta), heart2$delta * heart2$status), ncol = 3)
head(event)

X <- model.matrix(~ age + year + surgery, data = heart2)
X <- X[,-1] # Remove intercept
time3 <- heart2$times2
time3[time3 == 0] <- 0.0001


d.jags <- list(n = nrow(heart2), t1 = heart2$times1, t2 = heart2$time, t3 = time3,
  X = X, event = event, zeros = rep(0, nrow(heart2)), Nbetas = ncol(X))

i.jags <- function() { 
  list(beta = matrix(rnorm(3 * ncol(X)), ncol = 3), lambda = runif(3), alpha = runif(3)) 
}

p.jags <- c("beta", "alpha", "lambda")

# install.packages("rjags", dep = TRUE)
library("rjags")
m5 <- jags.model(data = d.jags, file = "IllDeath.txt", inits = i.jags, n.chains = 3)
update(m5, 1000)
res <- coda.samples(m5, variable.names = p.jags, n.iter = 10000, n.thin = 10)

result <- as.mcmc(do.call(rbind, res))

alpha1 <- result[, 1]; alpha2 <- result[, 2]; alpha3 <- result[, 3]
beta11 <- result[, 4]; beta21 <- result[, 5]; beta31 <- result[, 6]
beta12 <- result[, 7]; beta22 <- result[, 8]; beta32 <- result[, 9]
beta13 <- result[, 10]; beta23 <- result[, 11]; beta33 <- result[, 12]
lambda1 <- result[, 13]; lambda2 <- result[, 14]; lambda3 <- result[, 15]

table <- rbind(beta11, beta21, beta31, lambda1, alpha1, beta12, beta22, beta32, lambda2, alpha2, beta13, beta23, beta33, lambda3, alpha3)
table <- round(cbind(summary_parameters(table), c(mean(beta11 > 0), mean(beta21 > 0), mean(beta31 > 0), 1, 1,
  mean(beta12 > 0), mean(beta22 > 0), mean(beta32 > 0), 1, 1, 
  mean(beta13 > 0), mean(beta23 > 0), mean(beta33 > 0), 1, 1)), dig = 3)
colnames(table) <- c("Mean", "Sd", "2.5%", "Median", "97.5%", "P(.>0)")
rownames(table) <- c("beta11", "beta21", "beta31", "lambda1", "alpha1", "beta12", "beta22", "beta32", "lambda2", "alpha2", "beta13", "beta23", "beta33", "lambda3", "alpha3")
table

# p11
grid <- seq(0, 1000, length.out = 51)
x0 <- c(median(heart2$age), median(heart2$year), 0)
x1 <- c(median(heart2$age), median(heart2$year), 1)
p11_0_t_surgery0 <- mcmc_p11_s_t(t.pred = grid, s.pred = 0, l1 = lambda1, a1 = alpha1,
  b1 = cbind(beta11, beta21, beta31), l2 = lambda2, a2 = alpha2,
  b2 = cbind(beta12, beta22, beta32), x = x0) 
p11_0_t_surgery1 <- mcmc_p11_s_t(t.pred = grid, s.pred = 0, l1 = lambda1, a1 = alpha1,
  b1 = cbind(beta11, beta21, beta31), l2 = lambda2, a2 = alpha2,
  b2 = cbind(beta12, beta22, beta32), x = x1)

# p22
m <- median(heart2$times1[heart2$delta == 1])
grid2 <- grid + m
p22_m_t_surgery0 <- mcmc_p22_s_t(t.pred = grid2, s.pred = m, l1 = lambda1, a1 = alpha1,
  b1 = cbind(beta11, beta21, beta31), l3 = lambda3, a3 = alpha3,
  b3 = cbind(beta13, beta23, beta33), x = x0) 
p22_m_t_surgery1 <- mcmc_p22_s_t(t.pred = grid2, s.pred = m, l1 = lambda1, a1 = alpha1, b1 = cbind(beta11, beta21, beta31), l3 = lambda3, a3 = alpha3, b3 = cbind(beta13, beta23, beta33), x = x1)

# p12
p12_0_t_surgery0 <- mcmc_p12_s_t(t.pred = grid, s.pred = 0, l1 = lambda1, a1 = alpha1, b1 = cbind(beta11, beta21, beta31), l2 = lambda2, a2 = alpha2, b2 = cbind(beta12, beta22, beta32), l3 = lambda3, a3 = alpha3, b3 = cbind(beta13, beta23, beta33), x = x0) 
p12_0_t_surgery1 <- mcmc_p12_s_t(t.pred = grid, s.pred = 0, l1 = lambda1, a1 = alpha1, b1 = cbind(beta11, beta21, beta31), l2 = lambda2, a2 = alpha2, b2 = cbind(beta12, beta22, beta32), l3 = lambda3, a3 = alpha3, b3 = cbind(beta13, beta23, beta33), x = x1)

# p13
p13_0_t_surgery0 <- 1 - p11_0_t_surgery0 - p12_0_t_surgery0
p13_0_t_surgery1 <- 1 - p11_0_t_surgery1 - p12_0_t_surgery1

# p23
p23_m_t_surgery0 <- 1 - p22_m_t_surgery0
p23_m_t_surgery1 <- 1 - p22_m_t_surgery1


# install.packages("ggplot2", dep = TRUE)
library("ggplot2")
surg.colour <- rep(c("no", "yes"), each = length(grid))

# p11
df11 <- data.frame(time = rep(grid, 2), transition = c(p11_0_t_surgery0, p11_0_t_surgery1), surgery = factor(surg.colour))
p11 <- ggplot(data = df11, aes(x = time, y = transition, group = surgery, colour = surgery)) + ylim(0, 1) + theme_bw() +
  geom_line() + annotate("text", x = 880, y = 0.95, label=expression(1 %->% 1), size=5) + ylab("probability") + theme(legend.position = "top")

# p12
df12 <- data.frame(time = rep(grid, 2), transition = c(p12_0_t_surgery0, p12_0_t_surgery1), surgery = factor(surg.colour))
p12 <- ggplot(data = df12, aes(x = time, y = transition, group = surgery, colour = surgery)) + ylim(0, 1) + theme_bw() +
  geom_line() + annotate("text", x = 880, y = 0.95, label=expression(1 %->% 2), size=5) + ylab("probability") + theme(legend.position = "top")

# p13
df13 <- data.frame(time = rep(grid, 2), transition = c(p13_0_t_surgery0, p13_0_t_surgery1), surgery = factor(surg.colour))
p13 <- ggplot(data = df13, aes(x = time, y = transition, group = surgery, colour = surgery)) + ylim(0, 1) + theme_bw() +
  geom_line() + annotate("text", x = 880, y = 0.95, label=expression(1 %->% 3), size=5) + ylab("probability") + theme(legend.position = "top")

# p22
df22 <- data.frame(time = rep(grid, 2), transition = c(p22_m_t_surgery0, p22_m_t_surgery1), surgery = factor(surg.colour))
p22 <- ggplot(data = df22, aes(x = time, y = transition, group = surgery, colour = surgery)) + ylim(0, 1) + theme_bw() +
  geom_line() + annotate("text", x = 880, y = 0.95, label=expression(2 %->% 2), size=5) + ylab("probability") + theme(legend.position = "top")

# p23
df23 <- data.frame(time = rep(grid, 2), transition = c(p23_m_t_surgery0, p23_m_t_surgery1), surgery = factor(surg.colour))
p23 <- ggplot(data = df23, aes(x = time, y = transition, group = surgery, colour = surgery)) + ylim(0, 1) + theme_bw() +
  geom_line() + annotate("text", x = 880, y = 0.95, label=expression(2 %->% 3), size=5) + ylab("probability") + theme(legend.position="top")

# install.packages("gridExtra", dep = TRUE)
library("gridExtra")
grid.arrange(p11, p12, p13, p22, p23, nrow = 2, ncol = 3, layout_matrix = rbind(c(1, 2, 3), c(NA, 5, 6)))
