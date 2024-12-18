#~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Student: Raul
# Date: 2024-08-29
# Course: Survival Analysis
# Instructor: Dr. Fang
#~~~~~~~~~~~~~~~~~~~~~~~~~~

library(survival)
library(survtools) 
library(ggsurvfit)
library(dplyr)

head(lung)
# time = survival time  
# status (=1) censored indicator.



# (2.1). 
# Use sex as group indicator, 
# Test if the survival curves for the two groups are the same.

survfit2(Surv(time, status) ~ sex, data = lung) |> 
  ggsurvfit() +
  labs(
    x = "Time (in days)",
    y = "Survival probability"
  ) + 
  add_confidence_interval()

myfit <- survdiff(Surv(time, status) ~ sex, data = lung)
(myfit)

# (2.2). 
# Give the KM estimator of the survival function and its 95% confidence intervals

summary(survfit(Surv(time, status) ~ sex, data = lung))


