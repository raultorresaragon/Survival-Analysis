# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-12-12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#library(icenReg)
library(survival)


# Simulate a dataset
# ------------------
sim_interval_censored_data <- function(n, beta_true, cr) {
    Z <- rnorm(n)
    L <- rexp(n, rate = exp(Z * beta_true))
    R <- L + rexp(n, rate = 1)
    status <- ifelse(runif(n) < cr, 1, 0) 
    R[status == 0] <- Inf 
    data <- data.frame(L=L, R=R, Z=Z, status=status)
    data
}


# Iterative algorithm
# -------------------
#fit <- ic_par(Surv(L, R, type = "interval2") ~ Z, data = data) |> summary()
#beta <- fit$summaryParameters[3,1]

get_estimates_algo <- function(data, tol = 0.00001, max_iter = 50, verbose = FALSE) {
  
    L <-data$L
    R <- data$R
    status <- data$status
    Z <- data$Z
  
    # Set initial values with vanilla cox model
    # -----------------------------------------
    fit <- coxph(Surv(L, status) ~ Z, data = data) 
    fit_sum <- fit |> summary() |> coef()
    beta <- fit_sum[[1]]
    time_points <- sort(unique(c(L,R[R < Inf])))
    Lambda_0 <- rep(0.1,length(time_points))
    
    basehaz <- basehaz(initial_cox, centered = FALSE)
    time_points <- basehaz$time
    Lambda_0 <- basehaz$hazard
    
    
    
    # Log-likelihood function
    # -----------------------
    f_i <- function(beta, Lambda_0, L, R, Z, time_points) {
      sapply(1:n, function(i) {
        exp_Z_beta <- exp(Z[i] * beta)
        if (R[i] == Inf) {
          (-Lambda_0[time_points >= L[i]] * exp_Z_beta)[1]
        } else {
          Lambda_L <- Lambda_0[time_points >= L[i]][1]
          Lambda_R <- Lambda_0[time_points > R[i]][1]
          max(exp(-Lambda_L * exp_Z_beta) - exp(-Lambda_R * exp_Z_beta),.001)
        }
      }) |> sum(na.rm = TRUE)
    }
    
    
    
    # iterative algo
    # --------------
    for(iter in 1:max_iter){
      
      ### Step 2: for each subject we have density _f
      f_ <- numeric(n)
      f_ <- sapply(1:n, function(i) {
        exp_Z_beta <- exp(Z[i] * beta)
        if (R[i] == Inf) {
          f_[i] <- exp(-Lambda_0[time_points >= L[i]] * exp_Z_beta)[1]
        } else {
          Lambda_L <- Lambda_0[time_points >= L[i]][1]
          Lambda_R <- Lambda_0[time_points > R[i]][1]
          f_[i] <- max((exp(-Lambda_L * exp_Z_beta) - exp(-Lambda_R * exp_Z_beta)), 1e-10)
        }
      })
      
      ### Step 3: define new failure time
      Lambda_0_new <- Lambda_0
      
      for (j in seq_along(time_points)) {
        at_risk <- (L <= time_points[j]) & (time_points[j] < R)
        weights <- sum(f_[at_risk], na.rm = TRUE)
        if (weights > 0) {
          Lambda_0_new[j] <- sum(at_risk) / weights
        } else {
          Lambda_0_new[j] <- Lambda_0[j]  # No update if no at-risk weights
        }
      }
      
      beta_new <- optimize(function(b) -f_i(b, Lambda_0_new, L, R, Z, time_points), 
                           interval = c(-2,2))$minimum
      
      ### Assess whether to iterate or stop
      params_max <- max(abs(beta_new - beta), max(abs(log(Lambda_0_new) - log(Lambda_0))),na.rm = TRUE)
      if(params_max < tol) break
      beta <- beta_new
      Lambda_0 <- Lambda_0_new
      
    }
    
    if (verbose == TRUE) {
        print(paste0("number of iterations:", iter))
        print(paste0("beta hat = ", beta))
        print("Lambda_hat[1:5] = ")
        print(Lambda_0[1:5])
    }
    r <- list(beta = beta, Lambda_0 = Lambda_0)

}
