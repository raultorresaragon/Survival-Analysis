# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-12-12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~


library(icenReg)
library(survival)

# true_beta <- beta
# true_shape <- shape
# true_scale <- scale
# true_lambda_0 <- lambda_0
# true_T <- T
# keep <- c("data","true_beta","true_lambda_0", "true_T","true_shape","true_scale","sim_T")
# rm(list = ls()[!(ls() %in% keep)])

get_estimates_algo <- function(data, tol = 0.001) {


# Simulated interval-censored data from earlier
data_initial <- data
head(data_initial)
n <- nrow(data_initial)
X <- as.matrix(data[,c("Z","A")])
K <- 0


beta_2_hat_vec <- shape_new_vec <- scale_new_vec <- beta_2_se_vec <- beta_2_sig_vec <- numeric(0)
#for(j in 1:K) {
cont <- 1
j <- 0
while(cont == 1) {
  
  j <- j + 1
  
  # Step1
  # ~~~~~
    fit <- ic_par(Surv(L, R, type = "interval2") ~ Z + A, data = data)
    fit_sum <- summary(fit)
    beta_hat_new <- c(fit$coefficients[[3]], fit$coefficients[[4]])
    shape_new <- exp(fit$baseline[["log_shape"]])
    scale_new <- exp(fit$baseline[["log_scale"]])
    
    beta_2_hat_vec[j] <- beta_hat_new[2]
    beta_2_se_vec[j] <- fit_sum$summaryParameters["A","Std.Error"]
    beta_2_sig_vec[j] <- as.numeric(fit_sum$summaryParameters["A","p"] < 0.05)
    shape_new_vec[j] <- shape_new
    scale_new_vec[j] <- scale_new
    
    if (j > 2) {
      if (abs(beta_2_hat_vec[j-1] - beta_2_hat_vec[j]) > tol) {
        cont <- 1
      } else {
        cont <- 0
      }
    }
    
    
  # Step2: get new density function from estimates
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    lambda_0_new <- function(t) {
      ifelse(t > 0, shape_new * scale_new * t^(shape_new-1), 0) # hazard function from Weibull
    }
    Lambda_0_new <- function(t) {
      ifelse(t > 0, scale_new * t^shape_new, 0) # cumulative hazard function from Weibull
    }
    
    f_ <- function(t, X) {
      xbeta <- as.vector(X %*% beta_hat_new)
      lambda_0_t <- lambda_0_new(t)
      Lambda_0_t <- Lambda_0_new(t)
      lambda_0_t * exp(xbeta) * exp(Lambda_0_t * exp(xbeta))
    }
    
    
  # Step3: get T_new and data_new using f_
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    maximize_loglike <- function(L_i, R_i, X_i, f_){
      objective <- function(t) -log(f_(t,X_i))
      start_val <- runif(1,0,1)
      result <- optim(start_val, fn = objective, lower = L_i, upper = R_i, method = "Brent")
      result$par
    }
    
    T_new <- numeric(n)
    for(i in 1:n) {
      if(is.infinite(data$R[i])) {
        T_new[i] <- Inf #data$L[i] + 0.1
      } else {
        T_new[i] <- maximize_loglike(L_i = data$L[i], R_i = data$R[i], X_i = X[i,], f_ = f_)
      }
    }
    
    data$R <- T_new

}

data.frame(beta2 = beta_2_hat_vec, 
           beta2_se = beta_2_se_vec,
           beta2_sig = beta_2_sig_vec,
           shape = shape_new_vec, 
           scale = scale_new_vec)[j,]
}


