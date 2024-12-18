# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Author: Raul Torres Aragon
# Course: Advanced Survival 
# Instructor:
# Date: 2024-12-12
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(survival)

# Simulate interval-censored data
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

# Custom iterative algorithm
iterative_algo <- function(data, max_iter = 30, tol = 1e-4, verbose = FALSE) {
  L <- data$L
  R <- data$R
  Z <- data$Z
  status <- data$status
  n <- nrow(data)
  
  # Step 1: Initial Cox PH model to estimate beta and Lambda_0(t)
  initial_cox <- coxph(Surv(L, status) ~ Z, data = data)
  beta <- coef(initial_cox)
  basehaz <- basehaz(initial_cox, centered = FALSE)
  time_points <- basehaz$time
  Lambda_0 <- basehaz$hazard
  
  # Calculate likelihood f_i(t) for each subject
  compute_likelihood <- function(t, Lambda_0, beta, Z, L, R, i) {
    idx <- which.min(abs(t - time_points))  # Find closest Lambda_0(t)
    lambda_0_t <- ifelse(idx <= length(Lambda_0), Lambda_0[idx], 0)
    Lambda_0_t <- ifelse(idx <= length(Lambda_0), sum(Lambda_0[1:idx]), 0)
    exp_Z_beta <- exp(Z[i] * beta)
    lambda_0_t * exp_Z_beta * exp(-Lambda_0_t * exp_Z_beta)
  }
  
  # Step 2: Iterative updates of T, beta, and Lambda_0
  for (iter in 1:max_iter) {
    T <- numeric(n)  # New failure times
    
    for (i in 1:n) {
      # Define search range for T_i
      t_range <- if (R[i] < Inf) {
        time_points[time_points > L[i] & time_points <= R[i]]
      } else {
        time_points[time_points > L[i]]
      }
      
      # Calculate likelihood f_i(t) for all t in range
      # T_i = argmax(f_i(t)) if R < Inf
      # T_i > L_i for R = Inf, set 
      if (length(t_range) > 0) {
        f_t <- sapply(t_range, compute_likelihood, Lambda_0, beta, Z, L, R, i)
        T[i] <- t_range[which.max(f_t)]  
      } else {
        T[i] <- L[i] + .001  
      }
    }
    
    # Update data with new failure times
    new_data <- data
    new_data$status <- 1  # All events are considered uncensored since we predicted all Ts
    new_data$L <- T
    new_data$R <- T
    
    # Step 3: Estimate Cox PH model with updated failure times
    cox_model <- coxph(Surv(L, status) ~ Z, data = new_data)
    beta_new <- coef(cox_model)
    basehaz <- basehaz(cox_model, centered = FALSE)
    Lambda_0_new <- basehaz$hazard
    
    # Step 4: Convergence check
    beta_diff <- abs(beta_new - beta)
    if (beta_diff < tol) {
      if (verbose == TRUE) print(paste0("Converged after ", iter, " iterations"))
      break
    }
    
    # Update parameters for next iteration
    beta <- beta_new
    Lambda_0 <- Lambda_0_new
  }
  
  # Output results
  list(beta = beta, Lambda_0 = Lambda_0, time_points = basehaz$time, iterations = iter)
}

# Run the algorithm
# result <- iterative_algo(data)