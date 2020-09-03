# self$observed_likelihood -> self$observed_likelihood
# self$ -> self$

if (is.null(tmle_task)) {
  tmle_task <- self$observed_likelihood$training_task
}

intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))

# # clever_covariates happen here (for this param) only, but this is repeated computation
# HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
# 
# 
# # todo: make sure we support updating these params
# pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

# todo: extend for stochastic
cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

Y <- tmle_task$get_tmle_node(self$outcome_node)

# EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
# EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
# EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)

# ZW todo: store all possible inputs and outputs as df, and use left_join to call for outputs
# R and L(t!=0) nodes: treated
# Z nodes: control
# the function bellow from Likelihood class might be used to get all wanted interventions
# one group of interventions for treated, one group of interventions for controls; each group contains all combos

# self$observed_likelihood is just a likelihood; it has get_possible_counterfactuals method






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
    temp_output <- self$observed_likelihood$factor_list[[loc_node]]$get_likelihood(temp_task, fold_number = "full")  # corresponding outputs
    data.frame(temp_input, output = temp_output) %>% return
  }
})


obs_data[, tmle_task$npsem[[1]]$variables, with = F]

intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
intervention_levels <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
# nodes to integrate out in the target identification
# only support univaraite node for now; assume treatment level is one
all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                    rule_variables = c(intervention_variables, 
                                                       last(obs_variable_names)), 
                                    rule_values = c(intervention_levels, 1))
all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                     rule_variables = c(intervention_variables, 
                                                        last(obs_variable_names)), 
                                     rule_values = c(intervention_levels, 1))


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

# # ZW todo: get true IC
# # vars_A <- lapply(loc_A, function(s) tmle_task$npsem[[s]]$variables) %>% unlist  # if raw variable name is inquired from obs_data
# 
# # get a list of corresponding H covariates; ordered by nodes, not variables
# list_H <- list()
# # calculate RLY nodes
# for (temp_ind in loc_RLY) {
#   loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
#   loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
#   # this is the At indicators for H_RLY; now 
#   A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == 
#     self$intervention_list_treatment[[ temp_node_names[last(loc_A_needed)] ]]$value
#   # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
#   # these A probs will be taken as product
#   part_A <- lapply(loc_A_needed, function(k) self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])) %>% pmap_dbl(prod)
#   part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)), 
#                        lapply(loc_Z_needed, function(k) {
#                          temp_p_0 <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
#                          temp_p_1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
#                          ifelse(cf_task_control$data[[tmle_task$npsem[[k]]$variables]] == 1, temp_p_0, 1 - temp_p_0)
#                          / ifelse(cf_task_treatment$data[[tmle_task$npsem[[k]]$variables]] == 1, temp_p_1, 1 - temp_p_1)
#                        }) %>% pmap_dbl(prod))
#   list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
# }
# # calculate Z nodes
# for (temp_ind in loc_Z) {
#   loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
#   loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
#   A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == 
#     self$intervention_list_control[[ temp_node_names[last(loc_A_needed)] ]]$value
#   part_A <- lapply(loc_A_needed, function(k) {
#     1 - self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])  # this it the prob of being control level (0)
#   } ) %>% pmap_dbl(prod)
#   part_RLY <- lapply(loc_RLY_needed, function(k) {
#     temp_p_1 <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
#     temp_p_0 <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
#     ifelse(cf_task_treatment$data[[tmle_task$npsem[[k]]$variables]] == 1, temp_p_1, 1 - temp_p_1)
#     / ifelse(cf_task_control$data[[tmle_task$npsem[[k]]$variables]] == 1, temp_p_0, 1 - temp_p_0)
#   }) %>% pmap_dbl(prod)
#   list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
# }
# 
# 
# 
# 
# 
# # get a list of needed deltaQ
# list_deltaQ <- list()
# for (ind_var in 1:length(list_H)) {
#   if(!is.null(list_H[[ind_var]])) {
#     
#     loc_Z_needed <- loc_Z[loc_Z > ind_var]  # only product children variables
#     loc_RLY_needed <- loc_RLY[loc_RLY > ind_var]
#     
#     # deltaQ part in EIC
#     # for the order of ARZLY, at each RZLY node, integrate out all children
#     # note that no A probs will be involved in the product
#     loc_in_nodes_to_combo <- which(nodes_to_combo == temp_node_names[ind_var])
#     # Y_tau = 1, current Lt = 1 or 0
#     if ( (loc_in_nodes_to_combo + 1) > (length(nodes_to_combo) - 1) ) {
#       temp_combos_current_1 <- data.frame(1, 
#                                           1)
#       colnames(temp_combos_current_1)[c(1, ncol(temp_combos_current_1))] <- c(nodes_to_combo[loc_in_nodes_to_combo], 
#                                                                               last(nodes_to_combo)
#       )
#       temp_combos_current_0 <- temp_combos_current_1
#       temp_combos_current_0[, 1] <- 0
#     } else {
#       temp_combos_current_1 <- self$observed_likelihood$get_possible_counterfactuals(nodes = nodes_to_combo[(loc_in_nodes_to_combo + 1) : (length(nodes_to_combo) - 1)])
#       temp_combos_current_1 <- data.frame(1, 
#                                           temp_combos_current_1, 
#                                           1)
#       colnames(temp_combos_current_1)[c(1, ncol(temp_combos_current_1))] <- c(nodes_to_combo[loc_in_nodes_to_combo], 
#                                                                               last(nodes_to_combo)
#       )
#       temp_combos_current_0 <- temp_combos_current_1
#       temp_combos_current_0[, 1] <- 0
#     }
#     
#     if (ind_var == length(list_H)) {
#       list_deltaQ[[ind_var]] <- rep(1, length(list_H[[ind_var]]))
#     } else {
#       # For Q1, use current Lt = 1 (_1) for sum
#       ancesters_needed <- lapply(nodes_to_combo[ifelse_vec(loc_in_nodes_to_combo > 1, 1:(loc_in_nodes_to_combo - 1), 0)], function(temp_name) tmle_task$npsem[[ temp_name ]]$variables) %>% unlist
#       Q1 <- lapply(1:nrow(temp_combos_current_1), function(s) {
#         each_row <- temp_combos_current_1[s, ]
#         if (length(ancesters_needed) == 0) {
#           loc_tasks <- left_join(each_row, 
#                                  temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
#           
#           # for all non-A, non-0 variables, calculate the variable by rule
#           # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
#           ifelse_vec(length(loc_Z_needed) == 0, 
#                      data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_RLY_needed])), 
#                      data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_Z_needed]), 
#                                 self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_RLY_needed])
#                      )  
#           ) %>% pmap_dbl(prod)
#         } else {
#           combos_to_search <- data.frame(obs_data %>% select(ancesters_needed), 
#                                          each_row)
#           # ZW todo: for multivariate nodes
#           names(combos_to_search) <- names(temp_combos)
#           loc_tasks <- left_join(combos_to_search, 
#                                  temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
#           unique_locs <- loc_tasks %>% unique
#           temp_p <- rep(0, length(loc_tasks))
#           for (l in unique_locs) {
#             temp_p[which(loc_tasks == l)] <- 
#               ifelse_vec(length(loc_Z_needed) == 0, 
#                          data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                              temp_node_names[loc_RLY_needed])[which(loc_tasks == l)]
#                          ), 
#                          data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
#                                                                              temp_node_names[loc_Z_needed])[which(loc_tasks == l)], 
#                                     self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                              temp_node_names[loc_RLY_needed])[which(loc_tasks == l)]
#                          )
#               )
#           }
#           
#           tmle_task$npsem[[k]]$variables
#           
#           
#           sapply(1:length(loc_tasks), function(i) {
#             # for all non-A, non-0 variables, calculate the variable by rule
#             # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
#             
#             
#             
#             
#             
#             
#           })
#         }
#       }) %>% pmap_dbl(sum)  # sum over possible inputs
#       Q0 <- lapply(1:nrow(temp_combos_current_0), function(s) {
#         each_row <- temp_combos_current_0[s, ]
#         if (length(ancesters_needed) == 0) {
#           loc_tasks <- left_join(each_row, 
#                                  temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
#           
#           # for all non-A, non-0 variables, calculate the variable by rule
#           # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
#           ifelse_vec(length(loc_Z_needed) == 0, 
#                      data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_RLY_needed])
#                      ) %>% pmap_dbl(prod), 
#                      data.frame(self$observed_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_Z_needed]), 
#                                 self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                          temp_node_names[loc_RLY_needed])
#                      ) %>% pmap_dbl(prod)
#           )
#           
#         } else {
#           combos_to_search <- data.frame(obs_data %>% select(ancesters_needed), 
#                                          each_row)
#           # ZW todo: for multivariate nodes
#           names(combos_to_search) <- names(temp_combos)
#           loc_tasks <- left_join(combos_to_search, 
#                                  temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
#           sapply(1:length(loc_tasks), function(i) {
#             # for all non-A, non-0 variables, calculate the variable by rule
#             # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
#             loc_Z_needed <- loc_Z[loc_Z > ind_var]  # only product children variables
#             loc_RLY_needed <- loc_RLY[loc_RLY > ind_var]
#             ifelse_vec(length(loc_Z_needed) == 0, 
#                        prod(list(self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                           temp_node_names[loc_RLY_needed])[i]
#                        ) %>% unlist), 
#                        prod(list(self$observed_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
#                                                                           temp_node_names[loc_Z_needed])[i], 
#                                  self$observed_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
#                                                                           temp_node_names[loc_RLY_needed])[i]
#                        ) %>% unlist)
#             )
#             
#           })
#         }
#       }) %>% pmap_dbl(sum)  # sum over possible inputs
#       list_deltaQ[[ind_var]] <- Q1 - Q0 
#     }
#   }
# }
# 
# list_D <- list()
# for (ind_var in 1:length(list_H)) {
#   if(!is.null(list_H[[ind_var]])) {
#     # ZW todo: for discretized variables
#     current_ind <- (obs_data[[tmle_task$npsem[[ind_var]]$variables]] == 1)*1
#     
#     if (ind_var %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else 
#       temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
#     temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)
#     list_D[[ind_var]] <- (current_ind - temp_p) * 
#       list_H[[ind_var]] * (list_deltaQ[[ind_var]])
#   }
# }
# list_D[[1]] <- pmap_dbl(list_prods, sum) - psi
# 
# vec_D <- list_D %>% compact %>% pmap_dbl(sum)
# IC <- vec_D

result <- list(psi = psi, IC = 0)

return(result)