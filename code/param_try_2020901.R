# self$observed_likelihood -> initial_likelihood
# self$intervention_list_treatment
# self$ -> tmle_params$




if (is.null(tmle_task)) {
  tmle_task <- initial_likelihood$training_task
}


intervention_nodes <- union(names(tmle_params$intervention_list_treatment), names(tmle_params$intervention_list_control))

# # clever_covariates happen here (for this param) only, but this is repeated computation
# HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
# 
# 
# # todo: make sure we support updating these params
# pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

# todo: extend for stochastic
cf_task_treatment <- tmle_params$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- tmle_params$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]

Y <- tmle_task$get_tmle_node(tmle_params$outcome_node)

# EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
# EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
# EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)

# ZW todo: store all possible inputs and outputs as df, and use left_join to call for outputs
# R and L(t!=0) nodes: treated
# Z nodes: control
# the function bellow from Likelihood class might be used to get all wanted interventions
# one group of interventions for treated, one group of interventions for controls; each group contains all combos

# self$observed_likelihood is just a likelihood; it has get_possible_counterfactuals method

if (is.null(tmle_task)) {
  tmle_task <- initial_likelihood$training_task
}

# all not A, not t=0 nodes
temp_node_names <- names(tmle_task$npsem)
loc_A <- grep("A", temp_node_names)
if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
nodes_to_combo <- temp_node_names[ -c(loc_A, which(!if_not_0)) ]
temp_combos <- initial_likelihood$get_possible_counterfactuals(nodes = nodes_to_combo)
combos_treat <- combos_control <- data.frame(matrix(0, nrow(temp_combos), length(loc_A)), 
                                             temp_combos)  # first tau cols are treat/control
names(combos_treat)[1:length(loc_A)] <- names(combos_control)[1:length(loc_A)] <- temp_node_names[loc_A]
# ZW todo: stochastic interventions
combos_treat[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) tmle_params$intervention_list_treatment[[k]]$value %>% as.character %>% as.numeric)
combos_control[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) tmle_params$intervention_list_control[[k]]$value %>% as.character %>% as.numeric)
combo_node_names <- temp_node_names[if_not_0]
combos_treat <- combos_treat[ combo_node_names ]
combos_control <- combos_control[ combo_node_names ]
# only difference between these two groups of combos is in intervention nodes

# get the list of LF_static interventions
temp_list <- list()
temp_list[[1]] <- lapply(1:nrow(combos_treat), function(row_id) {
  temp_row <- combos_treat[row_id, ]
  temp_treat <- lapply(1:length(temp_row), function(k) {
    define_lf(LF_static, combo_node_names[k], value = temp_row[k])
  })
  names(temp_treat) <- combo_node_names
  return(temp_treat)
})
temp_list[[2]] <- lapply(1:nrow(combos_control), function(row_id) {
  temp_row <- combos_control[row_id, ]
  temp_control <- lapply(1:length(temp_row), function(k) {
    define_lf(LF_static, combo_node_names[k], value = temp_row[k])
  })
  names(temp_control) <- combo_node_names
  return(temp_control)
})
# temp_list[[1]]  # each slot of it is a intervention list, behaves like treatment or control

# list of CF_Likelihood objects
cf_list <- list()
cf_list[[1]] <- lapply(temp_list[[1]], function(s) {
  CF_Likelihood$new(initial_likelihood, s)
})
cf_list[[2]] <- lapply(temp_list[[2]], function(s) {
  CF_Likelihood$new(initial_likelihood, s)
})
# ZW todo: save combo fittings
private$.cf_likelihood_combo_list <- cf_list

cf_task_list <- list()
cf_task_list[[1]] <- lapply(cf_list[[1]], function(s) {
  s$enumerate_cf_tasks(tmle_task)[[1]]
})
cf_task_list[[2]] <- lapply(cf_list[[2]], function(s) {
  s$enumerate_cf_tasks(tmle_task)[[1]]
})



cf_list <- tmle_params$cf_likelihood_combo_list

loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))

# the list of prob products to be summed
list_prods <- lapply(1:nrow(combos_treat), function(s) {
  if(last(combos_treat[s, ]) == 0) {
    # for Y_tau = 0, return 0
    return(0)
  } else {
    current_combo <- combos_treat[s,]
    current_combo <- current_combo[-grep("A", names(current_combo))]
    # task_1 <- cf_list[[1]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=1 for RLY nodes
    # task_2 <- cf_list[[2]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=0 for Z nodes
    # take Z with control, take RLY with treat
    temp_return <- data.frame(initial_likelihood$get_likelihoods(cf_task_list[[1]][[s]], temp_node_names[loc_RLY], fold_number), 
                              initial_likelihood$get_likelihoods(cf_task_list[[2]][[s]], temp_node_names[loc_Z], fold_number)
    )
    temp_return <- temp_return[names(current_combo)]  # match for next step
    # decide p or 1-p
    temp_return <- lapply(1:length(current_combo), function(k) if (current_combo[k] == 1) temp_return[k] else 1 - temp_return[k] ) %>% data.frame
    apply(temp_return, 1, prod)
  }
})

psi <- pmap_dbl(list_prods, sum) %>% mean

# ZW todo: get true IC
obs_data <- tmle_task$data
# vars_A <- lapply(loc_A, function(s) tmle_task$npsem[[s]]$variables) %>% unlist  # if raw variable name is inquired from obs_data

# get a list of corresponding H covariates; ordered by nodes, not variables
list_H <- list()
# calculate RLY nodes
for (temp_ind in loc_RLY) {
  loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
  loc_Z_needed <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
  # this is the At indicators for H_RLY; now 
  A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == 
    tmle_params$intervention_list_treatment[[ temp_node_names[last(loc_A_needed)] ]]$value
  # A_ind <- obs_data[[temp_node_names[last(loc_A_needed)]]]  # using variable names (rather than node names) to inquire obs_data
  # these A probs will be taken as product
  part_A <- lapply(loc_A_needed, function(k) initial_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])) %>% pmap_dbl(prod)
  part_Z <- ifelse_vec(length(loc_Z_needed) == 0, rep(1, length(part_A)), 
                       lapply(loc_Z_needed, function(k) {
                         initial_likelihood$get_likelihoods(cf_task_control, temp_node_names[k]) / 
                           initial_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k])
                       }) %>% pmap_dbl(prod))
  list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_Z, 0) %>% as.vector
}
# calculate Z nodes
for (temp_ind in loc_Z) {
  loc_A_needed <- loc_A[loc_A < temp_ind]  # all needed A nodes
  loc_RLY_needed <- loc_RLY[loc_RLY < temp_ind]
  A_ind <- obs_data[[tmle_task$npsem[[last(loc_A_needed)]]$variables]] == 
    tmle_params$intervention_list_control[[ temp_node_names[last(loc_A_needed)] ]]$value
  part_A <- lapply(loc_A_needed, function(k) initial_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])) %>% pmap_dbl(prod)
  part_RLY <- lapply(loc_RLY_needed, function(k) {
    initial_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[k]) / 
      initial_likelihood$get_likelihoods(cf_task_control, temp_node_names[k])
  }) %>% pmap_dbl(prod)
  list_H[[temp_ind]] <- ifelse(A_ind, 1/part_A*part_RLY, 0) %>% as.vector
}


# get a list of needed deltaQ
list_deltaQ <- list()
for (ind_var in 1:length(list_H)) {
  if(!is.null(list_H[[ind_var]])) {
    
    loc_Z_needed <- loc_Z[loc_Z > ind_var]  # only product children variables
    loc_RLY_needed <- loc_RLY[loc_RLY > ind_var]
    
    # deltaQ part in EIC
    # for the order of ARZLY, at each RZLY node, integrate out all children
    # note that no A probs will be involved in the product
    loc_in_nodes_to_combo <- which(nodes_to_combo == temp_node_names[ind_var])
    # Y_tau = 1, current Lt = 1 or 0
    if ( (loc_in_nodes_to_combo + 1) > (length(nodes_to_combo) - 1) ) {
      temp_combos_current_1 <- data.frame(1, 
                                          1)
      colnames(temp_combos_current_1)[c(1, ncol(temp_combos_current_1))] <- c(nodes_to_combo[loc_in_nodes_to_combo], 
                                                                              last(nodes_to_combo)
      )
      temp_combos_current_0 <- temp_combos_current_1
      temp_combos_current_0[, 1] <- 0
    } else {
      temp_combos_current_1 <- initial_likelihood$get_possible_counterfactuals(nodes = nodes_to_combo[(loc_in_nodes_to_combo + 1) : (length(nodes_to_combo) - 1)])
      temp_combos_current_1 <- data.frame(1, 
                                          temp_combos_current_1, 
                                          1)
      colnames(temp_combos_current_1)[c(1, ncol(temp_combos_current_1))] <- c(nodes_to_combo[loc_in_nodes_to_combo], 
                                                                              last(nodes_to_combo)
      )
      temp_combos_current_0 <- temp_combos_current_1
      temp_combos_current_0[, 1] <- 0
    }
    
    if (ind_var == length(list_H)) {
      list_deltaQ[[ind_var]] <- rep(1, length(list_H[[ind_var]]))
    } else {
      # For Q1, use current Lt = 1 (_1) for sum
      ancesters_needed <- lapply(nodes_to_combo[ifelse_vec(loc_in_nodes_to_combo > 1, 1:(loc_in_nodes_to_combo - 1), 0)], function(temp_name) tmle_task$npsem[[ temp_name ]]$variables) %>% unlist
      Q1 <- lapply(1:nrow(temp_combos_current_1), function(s) {
        each_row <- temp_combos_current_1[s, ]
        if (length(ancesters_needed) == 0) {
          loc_tasks <- left_join(each_row, 
                                 temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
          
          # for all non-A, non-0 variables, calculate the variable by rule
          # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
          ifelse_vec(length(loc_Z_needed) == 0, 
                     data.frame(initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_RLY_needed])), 
                     data.frame(initial_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_Z_needed]), 
                                initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_RLY_needed])
                     )  
          ) %>% pmap_dbl(prod)
        } else {
          combos_to_search <- data.frame(obs_data %>% select(ancesters_needed), 
                                         each_row)
          # ZW todo: for multivariate nodes
          names(combos_to_search) <- names(temp_combos)
          loc_tasks <- left_join(combos_to_search, 
                                 temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
          
          sapply(1:length(loc_tasks), function(i) {
            # for all non-A, non-0 variables, calculate the variable by rule
            # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
            
            ifelse_vec(length(loc_Z_needed) == 0, 
                       prod(list(initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_RLY_needed])[i]
                       ) %>% unlist), 
                       prod(list(initial_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_Z_needed])[i], 
                                 initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_RLY_needed])[i]
                       ) %>% unlist)
            )
            
            
          })
        }
      }) %>% pmap_dbl(sum)  # sum over possible inputs
      Q0 <- lapply(1:nrow(temp_combos_current_0), function(s) {
        each_row <- temp_combos_current_0[s, ]
        if (length(ancesters_needed) == 0) {
          loc_tasks <- left_join(each_row, 
                                 temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
          
          # for all non-A, non-0 variables, calculate the variable by rule
          # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
          ifelse_vec(length(loc_Z_needed) == 0, 
                     data.frame(initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_RLY_needed])
                     ) %>% pmap_dbl(prod), 
                     data.frame(initial_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_Z_needed]), 
                                initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                   temp_node_names[loc_RLY_needed])
                     ) %>% pmap_dbl(prod)
          )
          
        } else {
          combos_to_search <- data.frame(obs_data %>% select(ancesters_needed), 
                                         each_row)
          # ZW todo: for multivariate nodes
          names(combos_to_search) <- names(temp_combos)
          loc_tasks <- left_join(combos_to_search, 
                                 temp_combos %>% mutate(which_row = 1:nrow(temp_combos)))$which_row
          sapply(1:length(loc_tasks), function(i) {
            # for all non-A, non-0 variables, calculate the variable by rule
            # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
            loc_Z_needed <- loc_Z[loc_Z > ind_var]  # only product children variables
            loc_RLY_needed <- loc_RLY[loc_RLY > ind_var]
            ifelse_vec(length(loc_Z_needed) == 0, 
                       prod(list(initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_RLY_needed])[i]
                       ) %>% unlist), 
                       prod(list(initial_likelihood$get_likelihoods(cf_task_list[[2]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_Z_needed])[i], 
                                 initial_likelihood$get_likelihoods(cf_task_list[[1]][[loc_tasks[i]]], 
                                                                    temp_node_names[loc_RLY_needed])[i]
                       ) %>% unlist)
            )
            
          })
        }
      }) %>% pmap_dbl(sum)  # sum over possible inputs
      list_deltaQ[[ind_var]] <- Q1 - Q0 
    }
  }
}


list_D <- list()
for (ind_var in 1:length(list_H)) {
  if(!is.null(list_H[[ind_var]])) {
    # ZW todo: for discretized variables
    current_ind <- (obs_data[[tmle_task$npsem[[ind_var]]$variables]] == 1)*1
    
    if (ind_var %in% loc_Z) temp_p <- initial_likelihood$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else 
      temp_p <- initial_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
    temp_p <- ifelse(current_ind == 1, temp_p, 1 - temp_p)
    list_D[[ind_var]] <- (current_ind - temp_p) * 
      list_H[[ind_var]] * (list_deltaQ[[ind_var]])
  }
}
list_D[[1]] <- pmap_dbl(list_prods, sum) - psi

vec_D <- list_D %>% compact %>% pmap_dbl(sum)



IC <- vec_D

result <- list(psi = psi, IC = IC)
