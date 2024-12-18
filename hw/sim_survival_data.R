# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-11-25
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Simulate interval survival data


# Parameters
# ----------
# n <- 100
# beta <- c(0.5, 1.5)
# shape <- 1
# scale <- 1
# cens_prob <- 0.3

sim_interval_censored_data <- function(n, beta = c(0.5, 1.5), shape = 1, scale = 1, cens_prob = 0.3) {

Z <- rnorm(n)
A <- rbinom(n,1,0.5)
X <- cbind(Z, A)
lambda_0 <- function(t) shape * scale * t^(shape-1) # baseline hazard pdf (Weibull) 1
Lambda_0 <- function(t) scale * t^shape # cumulative hazard function t


# Simulate failure times
# ----------------------
sim_T <- function(n, X, beta, Lambda_0) {
  
  U <- runif(n)
  T <- sapply(1:n, function(i) {
    inv_cum_hazard <- function(t) {
      Lambda_0(t) * exp(X[i,] %*% beta) - log(1 / U[i])
    }
    uniroot(function(t) inv_cum_hazard(t), lower = 0, upper = 100)$root
  })
  return(T)
}

T <- 10*sim_T(n, X, beta, Lambda_0)
T


# Generate interval-censored data (L, R)
# --------------------------------------
L <- R <- numeric(n)
for(i in 1:n) {
  L[i] <- max(0, T[i] - runif(1, 1, 5))
  R[i] <- T[i] + runif(1, 1, 5)
}


# Censored observations
# ---------------------
censored_indexes <- sample(1:n, size = n*cens_prob, replace = FALSE)
R[censored_indexes] <- Inf


# Pack into dataframe
# -------------------
data <- data.frame(T=T, L=L, R=R, Z=Z, A=A, 
                   C=as.numeric(is.infinite(R)), 
                   E=as.numeric(!is.infinite(R)))
data
}

