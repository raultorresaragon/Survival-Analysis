# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-12-12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
set.seed(1234)
library(survival)
library(dplyr)
library(tibble)
source("Fang_algo.R")
beta_true <- 1.5

# TESTING
### N <- 200
### beta_true <- 1.5
### cens_prob <- 0.3
### n <- 100
### data <- sim_interval_censored_data(n = n, beta_true = beta_true, cr = 1-cens_prob) 
### res <- iterative_algo(data = data)
### res$beta


run_sim <- function(N, n, beta_true, event_prob) {
    
    beta_hats <- numeric(N)
    for(i in 1:N) {
      data <- sim_interval_censored_data(n = n, beta_true = beta_true, cr = 1-event_prob)
      res <- iterative_algo(data = data)
      beta_hats[i] <- res$beta
    }
    res_tab <- data.frame(true_beta = beta_true, beta_hat = beta_hats)

    return <- tibble(mean_beta_hat = mean(res_tab$beta_hat, na.rm = TRUE), 
                     sd_beta_hat = sd(res_tab$beta_hat, na.rm = TRUE))
    return$coverage <- as.numeric(
      ((return$mean_beta_hat - 1.96*return$sd_beta_hat) <= beta_true) & 
      (beta_true <= (return$mean_beta_hat + 1.96*return$sd_beta_hat))
    )
    return
}
    

res_tab <- tibble(n=numeric(), 
                  cens_prob=numeric(), 
                  mean_beta_hat=numeric(), 
                  sd_beta_hat=numeric(), 
                  coverage=numeric())

ns <- c(30, 100, 300)
evs <- c(0.3, 0.5, 0.7)
for(ev in evs) {
  for(n in ns) {
      print(paste0("cr=", 1-ev, ";  n=", n))
      row <- run_sim(N=100, n=n, beta_true=beta_true, event_prob=ev) |> 
             as_tibble() |> 
             mutate(n=n, cens_rate=1-ev)
      res_tab <- rbind(res_tab, row)
  }
}
res_tab <- res_tab[, c("cens_rate","n","mean_beta_hat","sd_beta_hat")]


# Lambda_0
# --------
plot_Lambdas <- function(n, beta_true, event_rate = 0.3, N=50) {
    data <- sim_interval_censored_data(n = n, beta_true = beta_true, cr = 1-event_rate)
    res <- iterative_algo(data = data)
    
    Ts <- res$time_points
    Lambda_0_true <- 1 * res$time_points
    Lambda_0_estimated <- res$Lambda_0[1:length(Ts)]
    
    plot(res$time_points, Lambda_0_true, type = "l", col = "darkred", lwd = 3, 
         xlab = "Time", ylab = "Cumulative Baseline Hazard",
         main = paste0("True vs Estimated Lambda_0(t) \n","n = ",n, " cr = ", 1-event_rate))
    lines(res$time_points, Lambda_0_estimated, col = "lightblue", lwd = 1, lty = "solid")
    
    for(i in 1:N) {
      data <- sim_interval_censored_data(n = n, beta_true = beta_true, cr = 1-event_rate)
      res <- iterative_algo(data = data, verbose = FALSE)
      Ts <- res$time_points
      Lambda_0_estimated <- res$Lambda_0[1:length(Ts)]
      lines(Ts, Lambda_0_estimated, col = "lightblue", lwd = 1, lty = "solid")
    }
    legend("bottomright", legend = c("True Lambda_0(t)", "Estimated Lambda_0(t)"),
           col = c("darkred", "lightblue"), lty = c(1, 1), lwd = 2)
}

jpeg(file='Lambda_0_plot.jpeg')
par(mfrow = c(2,2))
plot_Lambdas(n= 50, event_rate=0.7, beta_true=1.5, N=25)
plot_Lambdas(n=100, event_rate=0.7, beta_true=1.5, N=25)
plot_Lambdas(n=500, event_rate=0.7, beta_true=1.5, N=25)
plot_Lambdas(n=1000,event_rate=0.7, beta_true=1.5, N=25)
par(mfrow = c(1,1))
dev.off()
save.image(file='survival_env.RData')




### DEPRECATED
### compute_Lambda_0 <- function(data, T) {
###   data$T <- T
###   data <- data[order(data$T), ]
###   unique_times <- unique(data$T[data$status == 1])  # Unique event times
###   cumulative_hazard <- numeric(length(unique_times))  # To store Lambda_0(t)
###   
###   # Calculate cumulative hazard using Nelson-Aalen estimator
###   risk_set <- n  # Start with all individuals at risk
###   cumsum_hazard <- 0
###   for (i in seq_along(unique_times)) {
###     event_time <- unique_times[i]
###     d_i <- sum(data$T == event_time & data$status == 1)  # Number of events at t_i
###     R_i <- sum(data$T >= event_time)                    # Risk set at t_i
###     delta_hazard <- d_i / R_i                           # Incremental hazard
###     cumsum_hazard <- cumsum_hazard + delta_hazard       # Update cumulative hazard
###     cumulative_hazard[i] <- cumsum_hazard
###   }
###   cumulative_hazard
### }


