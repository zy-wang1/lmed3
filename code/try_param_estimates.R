# self$observed_likelihood -> initial_likelihood
# self$ -> tmle_params$

if (is.null(tmle_task)) {
  tmle_task <- initial_likelihood$training_task
}

intervention_nodes <- union(names(tmle_params$intervention_list_treatment), names(tmle_params$intervention_list_control))

# # clever_covariates happen here (for this param) only, but this is repeated computation
# HA <- tmle_params$clever_covariates(tmle_task, fold_number)[[tmle_params$outcome_node]]
# 
# 
# # todo: make sure we support updating these params
# pA <- initial_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_treatment <- tmle_params$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_control <- tmle_params$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

# todo: extend for stochastic
cf_task_treatment <- tmle_params$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- tmle_params$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

Y <- tmle_task$get_tmle_node(tmle_params$outcome_node)

# all not A, not t=0 nodes
temp_node_names <- names(tmle_task$npsem)
loc_A <- grep("A", temp_node_names)
loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)

# get list of all possible predicted lkds
obs_data <- tmle_task$data
obs_variable_names <- colnames(obs_data)
list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
  if (loc_node > 1) {
    # currently only support univariate node for t>0
    current_variable <- tmle_task$npsem[[loc_node]]$variables
    temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
    temp_task <- lmed3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
    temp_output <- initial_likelihood$factor_list[[loc_node]]$get_likelihood(temp_task, fold_number = "full")  # corresponding outputs
    data.frame(temp_input, output = temp_output) %>% return
  }
})

intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
intervention_levels_treat <- map_dbl(tmle_params$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
intervention_levels_control <- map_dbl(tmle_params$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
# nodes to integrate out in the target identification
# only support univaraite node for now; assume treatment level is one
all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                     rule_variables = c(intervention_variables, 
                                                        last(obs_variable_names)), 
                                     rule_values = c(intervention_levels_treat, 1))
all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                     rule_variables = c(intervention_variables, 
                                                        last(obs_variable_names)), 
                                     rule_values = c(intervention_levels_control, 1))

# for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
unique_L0 <- obs_data[, tmle_task$npsem[[1]]$variables, with = F] %>% unique
library_L0 <- data.frame(unique_L0, output = 
                           map_dbl(1:nrow(unique_L0), function(which_row) {
                             temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_RZLY_0)
                             temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_RZLY_1)
                             # for all non-A, non-0 variables, calculate the variable by rule
                             # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                             # note that list_all_predicted_lkd is ordered by node
                             temp_list_0 <- lapply(loc_Z, 
                                                   function(each_t) {
                                                     left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                   })
                             temp_list_1 <- lapply(loc_RLY, 
                                                   function(each_t) {
                                                     left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                   })
                             temp_list <- c(temp_list_0, temp_list_1)
                             pmap_dbl(temp_list, prod) %>% sum %>% return
                           })
)
# substitution estimator
vec_est <- left_join(data_sim[[1]], library_L0)$output
psi <- mean(vec_est)

# tmle_params$.list_all_predicted_lkd <- list_all_predicted_lkd

# get true IC

# getting H's  
all_observed_1 <- all_observed_0 <- obs_data
for (temp_A in intervention_variables) {
  all_observed_1 <- all_observed_1 %>% mutate(!!temp_A := 1)
  all_observed_0 <- all_observed_0 %>% mutate(!!temp_A := 0)
}

list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = initial_likelihood, 
                    cf_task_treatment, cf_task_control, 
                    intervention_variables, intervention_levels_treat, intervention_levels_control)
# get a list of needed deltaQ
list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H, 
                      intervention_variables, intervention_levels_treat, intervention_levels_control, 
                      list_all_predicted_lkd, 
                      lt = 1)
list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H, 
                      intervention_variables, intervention_levels_treat, intervention_levels_control, 
                      list_all_predicted_lkd, 
                      lt = 0)

list_D <- list()
for (ind_var in 1:length(list_H)) {
  if(!is.null(list_H[[ind_var]])) {
    # ZW todo: for discretized variables
    current_ind <- (obs_data[[tmle_task$npsem[[ind_var]]$variables]] == 1)*1
    if (ind_var %in% loc_Z) temp_p <- initial_likelihood$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else 
      temp_p <- initial_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
    list_D[[ind_var]] <- (current_ind - temp_p) *
      list_H[[ind_var]] * (list_Q_1[[ind_var]] - list_Q_0[[ind_var]])
  }
}
list_D[[1]] <- vec_est - psi

vec_D <- list_D %>% compact %>% pmap_dbl(sum)
IC <- vec_D

result <- list(psi = psi, IC = IC)

return(result)