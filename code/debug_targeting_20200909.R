code_list <- list.files("./R", full.names = T)
lapply(code_list, source)


middle_spec <- lmed_middle(
  treatment_level = 1,
  control_level = 0
)
tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)
updater <- lmed3_Update$new(maxit = 1, convergence_type = "scaled_var",
                            fluctuation_type = "standard")
targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
updater$tmle_params <- tmle_params
test <- fit_lmed3(tmle_task, targeted_likelihood, tmle_params, updater)
test$estimates


updater$update(targeted_likelihood, tmle_task)
update_fold <- updater$update_fold
updater$maxit
updater$update_step(targeted_likelihood, tmle_task, update_fold)
if (self$check_convergence(tmle_task, update_fold)) {
  break
}

fold_number = "full"
all_submodels <- updater$generate_submodel_data(
  targeted_likelihood, tmle_task,
  fold_number
)
new_epsilons <- updater$fit_submodels(all_submodels)
updater$step_number

###
# problem
###
tmle_params[[1]]$intervention_list_treatment
targeted_likelihood$update(new_epsilons, updater$step_number, fold_number)  # now full lkd list is updated too

full_updates <- targeted_likelihood$updater$apply_update_full(targeted_likelihood$training_task, targeted_likelihood, fold_number, new_epsilons)

# apply_update_full
{
  update_nodes <- updater$update_nodes
  list_all_predicted_lkd <- targeted_likelihood$list_all_predicted_lkd
  tmle_task <- targeted_likelihood$training_task
  temp_node_names <- names(tmle_task$npsem)
  obs_data <- tmle_task$data
  obs_variable_names <- names(obs_data)
  
  tmle_params <- updater$tmle_params
  intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params$intervention_list_control))
  intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
  intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
  # apply update to all nodes
  updated_likelihoods <- lapply(update_nodes, function(update_node) {
    loc_node <- which(temp_node_names == update_node)
    
    current_newH <- get_current_newH(loc_node, 
                                     tmle_task, obs_data,
                                     intervention_variables, intervention_levels_treat, intervention_levels_control, 
                                     list_all_predicted_lkd
    )
    observed <- list_all_predicted_lkd[[update_node]][[update_node]] %>% as.vector()
    observed <- tmle_task$scale(observed, update_node)
    initial <- list_all_predicted_lkd[[update_node]]$output
    initial <- tmle_task$scale(initial, update_node)
    initial <- bound(initial, 0.005)
    submodel_data <- list(
      observed = observed,
      H = current_newH,
      initial = initial
    )
    # ZW todo: why H has to be a matrix here
    submodel_data_1 <- list(
      observed = observed[observed == 1],
      H = current_newH %>% as.matrix,
      initial = initial[observed == 1]
    )
    
    all_epsilon <- new_epsilons
    epsilon <- all_epsilon[[update_node]]
    updated_likelihood_1 <- updater$apply_submodel(submodel_data_1, epsilon)
    updated_likelihood <- submodel_data$initial
    updated_likelihood[submodel_data$observed == 1] <- updated_likelihood_1
    updated_likelihood[submodel_data$observed == 0] <- 1 - updated_likelihood_1  # list of probs are symmetric for now
    
    # un-scale to handle bounded continuous
    updated_likelihood <- tmle_task$unscale(
      updated_likelihood,
      update_node
    )
  })
  names(updated_likelihoods) <- update_nodes
  
}



tmle_params
intervention_nodes <- union(names(tmle_params$intervention_list_treatment), names(tmle_params$intervention_list_control))
intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)

tasks_at_step <- targeted_likelihood$cache$tasks
tasks_at_step %>% length
task_updates <- lapply(tasks_at_step, targeted_likelihood$updater$apply_update, targeted_likelihood, fold_number, new_epsilons)  # this returns updated (obs) lkd
###
# problem
###
targeted_likelihood$updater$apply_update(tmle_task = tasks_at_step[[1]], targeted_likelihood, fold_number, new_epsilons)
targeted_likelihood$updater$apply_update(tmle_task = tasks_at_step[[2]], targeted_likelihood, fold_number, new_epsilons)
targeted_likelihood$updater$apply_update(tmle_task = tasks_at_step[[3]], targeted_likelihood, fold_number, new_epsilons)

###
# debug
###
update_nodes <- updater$update_nodes
tmle_task <- tasks_at_step[[1]]
# tmle_task <- initial_likelihood$training_task

# get submodel data for all nodes
all_submodel_data <- updater$generate_submodel_data(
  likelihood, tmle_task,
  fold_number
)
# apply update to all nodes
updated_likelihoods <- lapply(update_nodes, function(update_node) {
  submodel_data <- all_submodel_data[[update_node]]
  epsilon <- all_epsilon[[update_node]]
  updated_likelihood <- updater$apply_submodel(submodel_data, epsilon)
  
  # un-scale to handle bounded continuous
  updated_likelihood <- tmle_task$unscale(
    updated_likelihood,
    update_node
  )
})
names(updated_likelihoods) <- update_nodes





update_nodes <- updater$update_nodes
all_submodel_data <- updater$generate_submodel_data(
  targeted_likelihood, tmle_task,
  fold_number
)



clever_covariates <- lapply(updater$tmle_params, function(tmle_param) {
  tmle_param$clever_covariates(tmle_task, fold_number)  # this returns the list of H
})
clever_covariates
observed_values <- lapply(update_nodes, tmle_task$get_tmle_node)
all_submodels <- lapply(update_nodes, function(update_node) {
  node_covariates <- lapply(clever_covariates, `[[`, update_node)
  covariates_dt <- do.call(cbind, node_covariates)  # there might be multiple targets

  observed <- tmle_task$get_tmle_node(update_node)  # raw data on this node
  ###
  # problem
  ###
  initial <- targeted_likelihood$get_likelihood(
    tmle_task, update_node,
    fold_number
  )  # observed (updated or initial) likelihood at this node
  
  # scale observed and predicted values for bounded continuous
  observed <- tmle_task$scale(observed, update_node)
  initial <- tmle_task$scale(initial, update_node)
  
  # protect against qlogis(1)=Inf
  initial <- bound(initial, 0.005)
  
  submodel_data <- list(
    observed = observed,
    H = covariates_dt,
    initial = initial
  )
})
names(all_submodels) <- update_nodes




# cache step without updating value; maybe happened in params
targeted_likelihood$cache$get_values(targeted_likelihood$factor_list[[5]], tmle_task, fold_number)
initial_likelihood$cache$get_values(initial_likelihood$factor_list[[5]], tmle_task, fold_number)


clever_covariates <- lapply(updater$tmle_params, function(tmle_param) {
  tmle_param$clever_covariates(tmle_task, fold_number)  # this returns the list of H
})
targeted_likelihood$cache$get_update_step(likelihood_factor = targeted_likelihood$factor_list[[6]], tmle_task = tmle_task, fold_number = "full")


names(all_submodels) <- update_nodes


self$generate_submodel_data(
  initial_likelihood, tmle_task,
  fold_number
)