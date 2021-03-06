library(data.table)
library(sl3)
library(tlverse)
library(dplyr)  # %>% 

# dependency
library(R6)  # R6Class
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
library(methods)  # is


home <- getwd()  # specify it if needed
source(file.path(home, "code", "basic_functions-202008.R"))
source(file.path(home, "code", "generate_Zheng_data-202008.R"))
source("./code/try_lmed3_to_import_202008.R")
code_list <- list.files("./R", full.names = T)
lapply(code_list, source)

sample_size <- 100
data_sim <- generate_Zheng_data(B = sample_size, tau = 2)
data_wide <- data_sim %>% data.frame

# ARZLY model
node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1", 
                  A_2 = "A_2", 
                  R_2 = "R_2", 
                  Z_2 = "Z_2", 
                  L_2 = "L1_2", 
                  Y_2 = "Y_2" 
                  )

middle_npsem(node_list)

middle_spec <- lmed_middle(
  treatment_level = 1,
  control_level = 0
)

tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)

tmle_task$npsem[[2]] %>% str
tmle_task$print()

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(lrnr_mean, lrnr_glm)
))
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates


initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)
print(initial_likelihood)

initial_likelihood$get_likelihoods(tmle_task)

# targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood)

# targeted_likelihood_no_cv <-
#   Targeted_Likelihood$new(initial_likelihood,
#                           updater = list(cvtmle = FALSE)
#   )


test <- Param_middle$new(initial_likelihood, treatment, control, outcome_node = last(temp_names))
test$estimates(tmle_task)$psi


# # try initial_likelihood substitution for now
# tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood)
# print(tmle_params)
# 
# 
# 
# tmle_fit_manual <- fit_tmle3(
#   tmle_task, targeted_likelihood, tmle_params,
#   targeted_likelihood$updater
# )
# print(tmle_fit_manual)
# 
