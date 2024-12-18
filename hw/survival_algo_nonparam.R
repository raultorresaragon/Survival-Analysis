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

get_estimates_algo <- function(data, tol = 0.001, verbose = FALSE) {


# Simulated interval-censored data from earlier
data_initial <- data
head(data_initial)
n <- nrow(data_initial)
X <- as.matrix(data[,c("Z","A")])
K <- 0


beta_2_hat_vec <- beta_2_se_vec <- beta_2_sig_vec <- numeric(0)
shape_new_vec <- scale_new_vec <- numeric(0)
#for(j in 1:K) {
cont <- 1
j <- 0
while(cont == 1) {
  
  j <- j + 1
  
  # Step1
  # ~~~~~
  
    if(j==1){
      fit <- ic_par(Surv(L, R, type = "interval2") ~ Z + A, data = data)
      fit_sum <- summary(fit)
      beta_hat_new <- c(fit$coefficients[[3]], fit$coefficients[[4]])
      beta_hat_se <- fit_sum$summaryParameters["A","Std.Error"]
      beta_sig <- as.numeric(fit_sum$summaryParameters["A","p"] < 0.05)
      shape_new <- exp(fit$baseline[["log_shape"]])
      scale_new <- exp(fit$baseline[["log_scale"]])
      
    } else {
      
      # non-parametric fit
      # fit <- coxph(Surv(T_new, E) ~ Z + A, data = data)
      # fit_sum <- summary(fit)
      # beta_hat_new <- coefficients(fit_sum)[2,1]
      # beta_hat_se <- coefficients(fit_sum)[2,3]
      # 
      # beta_2_hat_vec[j] <- beta_hat_new
      # beta_2_se_vec[j] <- beta_hat_se
      # beta_2_sig_vec[j] <- as.numeric(coefficients(fit_sum)[2,5] < 0.05)
      
      # parametric fit
      fit <- survreg(Surv(T_new, E) ~ Z + A, data = data, dist = "weibull")
      fit_sum <- summary(fit)
      beta_hat_new <- c(fit_sum$table[2,1] |> abs(), fit_sum$table[3,1] |> abs())
      beta_hat_se <- fit_sum$table[3,2] |> abs()
      beta_sig <- as.numeric(fit_sum$table[3,4] < 0.05)
      scale_new <- fit$scale
      shape_new <- 1 / scale_new
    
    }      
  
    beta_2_hat_vec[j] <- beta_hat_new[2]
    beta_2_se_vec[j] <- beta_hat_se
    beta_2_sig_vec[j] <- beta_sig
    shape_new_vec[j] <- shape_new
    scale_new_vec[j] <- scale_new
    
    
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
        T_new[i] <- data$L[i] + 0.1
      } else {
        T_new[i] <- maximize_loglike(L_i = data$L[i], R_i = data$R[i], X_i = X[i,], f_ = f_)
      }
    }
    
    data$T_new <- T_new
    
    
  # Assess whether to continue iterating based on set tolerance
    if (j > 2) {
        
        beta_new_minus_beta_old <-abs(beta_2_hat_vec[j-1] - beta_2_hat_vec[j])
        proceed1 <- ifelse(!is.na(beta_new_minus_beta_old), 
                           beta_new_minus_beta_old > tol,
                           FALSE)
        proceed2 <- ifelse(j<20, TRUE, FALSE)
        cont <- ifelse(proceed1 == TRUE & proceed2 == TRUE, 1, 0) 
    }
}


if(verbose==TRUE) {
  print(paste0("convergence reached on iteration: ", j))
}
data.frame(beta2 = beta_2_hat_vec, 
           beta2_se = beta_2_se_vec,
           beta2_sig = beta_2_sig_vec,
           shape = shape_new_vec, 
           scale = scale_new_vec)[j,]

}


