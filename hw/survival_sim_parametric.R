# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-12-12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rm(list = ls())
set.seed(1234)
library(icenReg)
library(survival)
library(dplyr)
library(tibble)

source("sim_survival_data.R")
#source("survival_algo.R")
source("survival_algo_nonparam.R")


N <- 5e2
beta <- c(0.5, 1.5)
shape <- 1
scale <- 1

cens_prob <- 0.3
n<- 100

run_sim <- function(N, n, beta, shape, scale, cens_prob) {
    
    beta2v <- scalev <- shapev <- beta2sev <- beta2sigv <- numeric(N)
    for(i in 1:N) {
      data <- sim_interval_censored_data(n, beta, shape, scale, cens_prob)
      res <- get_estimates_algo(data = data, verbose = FALSE)
      beta2v[i] <- res[1,1]
      beta2sev[i] <- res[1,2]
      beta2sigv[i] <- res[1,3]
      scalev[i] <- res[1,4]
      shapev[i] <- res[1,5]
    }
    res_tab <- data.frame(true_beta2 = beta[2], beta2_hat = beta2v, beta2_se = beta2sev, beta2_sig = beta2sigv,
                          true_shape = shape, shape_hat = shapev,
                          true_scale = scale, scale_hat = scalev)
    
    coverage <- ((res_tab$beta2_hat - 1.96*res_tab$beta2_se) <= beta[2]) & 
                (beta[2] <= (res_tab$beta2_hat + 1.96*res_tab$beta2_se))

    return <- list(mean_beta_hat = mean(res_tab$beta2_hat, na.rm = TRUE), 
                   se_beta_hat = sd(res_tab$beta2_hat, na.rm = TRUE),
                   mean_shape_hat = mean(res_tab$shape_hat, na.rm = TRUE), 
                   se_shape_hat = sd(res_tab$shape_hat, na.rm = TRUE),
                   mean_scale_hat = mean(res_tab$scale_hat, na.rm = TRUE), 
                   se_scale_hat = sd(res_tab$scale_hat, na.rm = TRUE),
                   mean_coverage = mean(coverage, na.rm = TRUE),
                   mean_pvallt5 = mean(res_tab$beta2_sig, na.rm = TRUE))
}

res_tab <- tibble(n=numeric(), 
                  cens_prob=numeric(), 
                  mean_beta_hat=numeric(), 
                  se_beta_hat=numeric(), 
                  mean_shape_hat=numeric(), 
                  se_shape_hat=numeric(),
                  mean_scale_hat=numeric(), 
                  se_scale_hat=numeric(), 
                  mean_coverage=numeric(),
                  mean_pvallt5=numeric())

ns <- c(30, 100, 500)
cps <- c(0.3, 0.5, 0.7)
for(cens_prob in cps) {
  for(n in ns) {
      print(paste0("cp = ", cens_prob, "; n = ", n))
      row <- run_sim(N, n=n, beta, shape, scale, cens_prob=cens_prob) |> 
             as_tibble() |> 
             mutate(n=n, cens_prob=cens_prob)
      res_tab <- rbind(res_tab, row)
  }
}
res_tab[,c("cens_prob","n","mean_beta_hat","se_beta_hat","mean_coverage","mean_pvallt5")]

save.image(file='survival_env.RData')
