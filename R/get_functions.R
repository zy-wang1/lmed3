get_obs_Q <- function(tmle_task, obs_data, list_H, 
                      intervention_variables, intervention_levels_treat, intervention_levels_control, 
                      list_all_predicted_lkd, 
                      lt) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  list_Q <- list()
  for (loc_node in 1:length(list_H)) {
    if(!is.null(list_H[[loc_node]])) {
      # the Q integral at the previous variable; current inserted as lt; remove the lt setting if current Q is wanted
      # for the order of ARZL, at each RZL node, integrate out all children and set current = lt
      # note that no A probs will be involved in the product
      # debug: set all A probs as 1, because of Q definition
      
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      loc_current_var <- which(obs_variable_names == tmle_task$npsem[[loc_node]]$variables)
      all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ), 
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables), 
                                           rule_values = c(1, lt, intervention_levels_treat))
      all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:(loc_current_var-1) ),
                                           rule_variables = c(last(obs_variable_names), obs_variable_names[loc_current_var], intervention_variables), 
                                           rule_values = c(1, lt, intervention_levels_control))
      
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- obs_data[1:(loc_current_var-1)] %>% unique
      library_output <- data.frame(unique_input, output = 
                                     map_dbl(1:nrow(unique_input), function(which_row) {
                                       # probs in the integrals, A=1 or A=0 is inserted
                                       temp_all_comb_0 <- data.frame(unique_input[which_row, ], all_possible_RZLY_0)
                                       temp_all_comb_1 <- data.frame(unique_input[which_row, ], all_possible_RZLY_1)
                                       for (i in 1:length(intervention_variables_loc)) {
                                         temp_all_comb_0[, intervention_variables_loc[i]] <- intervention_levels_control[i]
                                         temp_all_comb_1[, intervention_variables_loc[i]] <- intervention_levels_treat[i]
                                       }
                                       # for all non-A, non-0 variables, calculate the variable by rule
                                       # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                       loc_Z_needed <- loc_Z[loc_Z > loc_node]  # only product children variables
                                       temp_list_0 <- lapply(loc_Z, 
                                                             function(each_t) {
                                                               left_join(temp_all_comb_0, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       loc_RLY_needed <- loc_RLY[loc_RLY > loc_node]
                                       temp_list_1 <- lapply(loc_RLY_needed, 
                                                             function(each_t) {
                                                               left_join(temp_all_comb_1, list_all_predicted_lkd[[each_t]])$output
                                                             })
                                       temp_list <- c(temp_list_0, temp_list_1)
                                       pmap_dbl(temp_list, prod) %>% sum %>% return
                                     })
      )
      list_Q[[loc_node]] <- left_join(obs_data[1:(loc_current_var-1)], library_output)$output
    }
  }
}




get_obs_H <- function(tmle_task, obs_data, current_likelihood, 
                      cf_task_treatment, cf_task_control, 
                      intervention_variables, intervention_levels_treat, intervention_levels_control
) {
  obs_variable_names <- colnames(obs_data)
  temp_node_names <- names(tmle_task$npsem)
  loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
  loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
  intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
  
  # get a list of corresponding H covariates; ordered by nodes, not variables
  list_H <- list()
  # calculate RLY nodes
  for (temp_ind in loc_RLY) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    # this is the At indicators for H_RLY; now
    A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] ==
      intervention_levels_treat[tmle_task$npsem[[last(loc_A_needed)]]$variables]
    # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
    # these A probs will be taken as product
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])) %>% pmap_dbl(prod)  # this is the likelihood of being 1
    part_Z <- lapply(loc_Z_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k]) / 
        current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    if(length(part_Z) == 0) part_Z <- 1
    
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
  }
  # calculate Z nodes
  for (temp_ind in loc_Z) {
    loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
    A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] ==
      intervention_levels_control[tmle_task$npsem[[last(loc_A_needed)]]$variables]
    part_A <- lapply(loc_A_needed, function(k) current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])) %>% pmap_dbl(prod)
    part_RLY <- lapply(loc_RLY_needed, function(k) {
      current_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k]) / 
        current_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
    }) %>% pmap_dbl(prod)
    list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
  }
  return(list_H)
}
