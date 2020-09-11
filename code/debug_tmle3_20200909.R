node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1"
                  # ,
                  # A_2 = "A_2",
                  # R_2 = "R_2",
                  # Z_2 = "Z_2",
                  # L_2 = "L1_2",
                  # Y_2 = "Y_2"
)

if_misspec <- F
data_sim <- generate_Zheng_data(B = sample_size, tau = 2, if_LY_misspec = if_misspec)
data_wide <- data.frame(data_sim)
middle_spec <- lmed_middle(
  treatment_level = 1,
  control_level = 0
)
tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)


# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
  learners = list(
    lrnr_mean,
    # lrnr_glm, 
    lrnr_glm_fast)
))
# learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates


initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)
# updater <- lmed3_Update$new(maxit = 5, convergence_type = "scaled_var",
#                             fluctuation_type = "standard", constrain_step = 1)
updater <- lmed3_Update$new(maxit = 1, convergence_type = "scaled_var",
                            fluctuation_type = "standard")
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
test <- fit_lmed3(tmle_task, targeted_likelihood, tmle_params, updater)
test$estimates[[1]]$psi
test$estimates[[1]]$IC %>% mean
test$initial_psi
targeted_likelihood$cache$get_update_step(likelihood_factor = targeted_likelihood$factor_list[[3]], tmle_task = tmle_task, fold_number = "full")
updater$check_convergence(tmle_task = tmle_task, fold_number = "full")

initial_psi <- sapply(
  tmle_params,
  function(tmle_param) {
    tmle_param$estimates(tmle_task)$psi
  }
)  
initial_psi
nontargeting$psi

# check updater$update;  if there has been 1 targeting step
updater$update(targeted_likelihood, tmle_task)
update_fold <- updater$update_fold
maxit <- updater$maxit
for (steps in seq_len(maxit)) {
  updater$update_step(targeted_likelihood, tmle_task, update_fold)
  if (updater$check_convergence(tmle_task, update_fold)) {
    break
  }
}
updater$steps
estimates <- lapply(
  tmle_params,
  function(tmle_param) {
    tmle_param$estimates(tmle_task, updater$update_fold)
  }
)
temp_density_sub
temp_density_tmle
estimates
ED_from_estimates(estimates)
targeted_likelihood$factor_list[[8]]
targeted_likelihood$cache$get_update_step(likelihood_factor = targeted_likelihood$factor_list[[8]], tmle_task = tmle_task, fold_number = "full")

# check updater$update_step and updater$check_convergence
# every time updater$update_step and updater$generate_submodel_data are called, update clever covariates obs in params
updater$update_step(targeted_likelihood, tmle_task, update_fold)
###
# consider cache clever covariate in param, everytime full lkd of targeted_likelihood is updated
###
updater$check_convergence(tmle_task, update_fold)
fold_number <- "full"
all_submodels <- updater$generate_submodel_data(
  targeted_likelihood, tmle_task,
  fold_number
)
all_submodels$R_1$observed
all_submodels$R_1$H
all_submodels$R_1$initial
new_epsilons <- updater$fit_submodels(all_submodels)
targeted_likelihood$update(new_epsilons, updater$step_number, fold_number)  # now full lkd list is updated too
# updater$step_number <- updater$step_number + 1

# check updater$generate_submodel_data
update_nodes <- updater$update_nodes
clever_covariates <- lapply(tmle_params, function(tmle_param) {
  tmle_param$clever_covariates(tmle_task, fold_number, update = T)  # this returns the list of H
})

# check updater$fit_submodels

# check likelihood$update
targeted_likelihood$update()
update_nodes <- updater$update_nodes
full_updates <- updater$apply_update_full(self$training_task, self, fold_number, new_epsilons)  # this returns updated (full) lkd