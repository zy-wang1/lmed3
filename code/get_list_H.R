get_list_H <- function(all_observed_1, all_observed_0, data_wide, variables, tau, list_predicted_probs) {
  list_H <- list()
  
  # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
  for (temp_ind in which_variable_drop(variables, to_drop_variable = c("A", "Z"), to_drop_time = 0, to_drop = paste0("Y_", tau))) {
    loc_A <- which_variable_take(variables, "A")
    loc_A <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z <- which_variable_take(variables, "Z")
    loc_Z <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    
    # all predicted probs at these locs are needed
    # needed inputs are in either all_observed_1 or all_observed_0
    
    # these A probs will be taken as product
    part_A <- lapply(loc_A, function(k) {
      left_join(all_observed_1, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    # these ratios Z probs will be taken as product
    part_Z <- ifelse_vec(length(loc_Z) == 0, rep(1, length(part_A)), 
                         lapply(loc_Z, function(k) {
                           left_join(all_observed_0, list_predicted_probs[[k]])$output / 
                             left_join(all_observed_1, list_predicted_probs[[k]])$output
                         }) %>% pmap_dbl(prod))
    
    list_H[[temp_ind]] <- ifelse(data_wide[last(loc_A)] == 1, 1/part_A*part_Z, 0) %>% as.vector
  }
  # Z nodes
  for (temp_ind in which_variable_take(variables, "Z")) {
    loc_A <- which_variable_take(variables, "A")
    loc_A <- loc_A[loc_A < temp_ind]
    # all the non-A, non-Z, not t=0, nodes before temp_ind; Y_tau won't be involved
    loc_LR <- which_variable_drop(variables, c("A", "Z"), 0)
    loc_LR <- loc_LR[loc_LR < temp_ind]
    
    # these A probs will be taken as product
    part_A <- lapply(loc_A, function(k) {
      left_join(all_observed_0, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    # these ratios Z probs will be taken as product
    part_LR <- lapply(loc_LR, function(k) {
      left_join(all_observed_1, list_predicted_probs[[k]])$output / 
        left_join(all_observed_0, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    list_H[[temp_ind]] <- ifelse(data_wide[last(loc_A)] == 0, 1/part_A*part_LR, 0) %>% as.vector
  }
  return(list_H)
}


# get Q for the previous variable
get_list_Q <- function(data_wide, variables, tau, list_predicted_probs, list_H, lt = 1) {
  list_Q <- list()
  # data_1 and data_0 are observed values with A=1 and A=0; stil need to impute these in probabilities
  # lt = 1 for now
  
  # data_1 <- data_0 <- data_wide
  # data_1[, which_variable_take(variables, "A")] <- 1
  # data_0[, which_variable_take(variables, "A")] <- 0
  
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      # the Q integral at the previous variable; current inserted as lt; remove the lt setting if current Q is wanted
      # for the order of ARZL, at each RZL node, integrate out all children and set current = lt
      # note that no A probs will be involved in the product
      # debug: set all A probs as 1, because of Q definition
      
      # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
      all_possible_rlz_1 <- expand_values(variables, to_drop = c(1:(ind_var-1) ), 
                                          A = 1, 
                                          rule_variables = c(last(variables), variables[ind_var]), rule_values = c(1, lt))
      all_possible_rlz_0 <- expand_values(variables, to_drop = c(1:(ind_var-1) ), 
                                          A = 0, 
                                          rule_variables = c(last(variables), variables[ind_var]), rule_values = c(1, lt))
      
      # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
      unique_input <- data_wide[1:(ind_var-1)] %>% unique
      library_output <- data.frame(unique_input, output = 
                                 map_dbl(1:nrow(unique_input), function(which_row) {
                                   # probs in the integrals, A=1 or A=0 is inserted
                                   temp_all_comb_0 <- cbind(unique_input[which_row, ], all_possible_rlz_0)
                                   temp_all_comb_1 <- cbind(unique_input[which_row, ], all_possible_rlz_1)
                                   temp_all_comb_0[which_variable_take(variables, "A")] <- 0
                                   temp_all_comb_1[which_variable_take(variables, "A")] <- 1
                                   
                                   # for all non-A, non-0 variables, calculate the variable by rule
                                   # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                   loc_Z <- which_variable_take(variables, to_take_variable = "Z")
                                   loc_Z <- loc_Z[loc_Z > ind_var]  # only product children variables
                                   temp_list_0 <- lapply(loc_Z, 
                                                         function(each_t) {
                                                           left_join(temp_all_comb_0, list_predicted_probs[[each_t]])$output
                                                         })
                                   loc_other <- which_variable_drop(variables, c("A", "Z"), 0)
                                   loc_other <- loc_other[loc_other > ind_var]
                                   temp_list_1 <- lapply(loc_other, 
                                                         function(each_t) {
                                                           left_join(temp_all_comb_1, list_predicted_probs[[each_t]])$output
                                                         })
                                   temp_list <- c(temp_list_0, temp_list_1)
                                   pmap_dbl(temp_list, prod) %>% sum %>% return
                                 })
      )
      list_Q[[ind_var]] <- left_join(data_wide[1:(ind_var-1)], library_output)$output
    }
  }
  
  
  return(list_Q)
}




# get all Q's for the current variable; taking previous variables as inputs
# debug it is still the next non-A Q that is needed; for example at L_t, R_{t+1} is needed, so current value is not integrated out
get_Q_current <- function(data_temp, variables, tau, list_predicted_probs) {
  # data_temp is now only a part of it; inputs from list of probs
  ind_var <- ncol(data_temp)
  
  # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1;
  all_possible_rlz_1 <- expand_values(variables, to_drop = c(1:(ind_var) ), 
                                      A = 1, 
                                      rule_variables = c(last(variables)), rule_values = c(1))
  all_possible_rlz_0 <- expand_values(variables, to_drop = c(1:(ind_var) ), 
                                      A = 0, 
                                      rule_variables = c(last(variables)), rule_values = c(1))
  
  # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
  unique_input <- data_temp[1:(ind_var)] %>% unique
  library_output <- data.frame(unique_input, output = 
                                 map_dbl(1:nrow(unique_input), function(which_row) {
                                   # probs in the integrals, A=1 or A=0 is inserted
                                   temp_all_comb_0 <- cbind(unique_input[which_row, ], all_possible_rlz_0)
                                   temp_all_comb_1 <- cbind(unique_input[which_row, ], all_possible_rlz_1)
                                   temp_all_comb_0[which_variable_take(variables, "A")] <- 0
                                   temp_all_comb_1[which_variable_take(variables, "A")] <- 1
                                   
                                   # for all non-A, non-0 variables, calculate the variable by rule
                                   # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                                   loc_Z <- which_variable_take(variables, to_take_variable = "Z")
                                   loc_Z <- loc_Z[loc_Z > ind_var]  # only integrate previous variables; debug: not including current
                                   temp_list_0 <- lapply(loc_Z, 
                                                         function(each_t) {
                                                           left_join(temp_all_comb_0, list_predicted_probs[[each_t]])$output
                                                         })
                                   loc_other <- which_variable_drop(variables, c("A", "Z"), 0)
                                   loc_other <- loc_other[loc_other > ind_var]
                                   temp_list_1 <- lapply(loc_other, 
                                                         function(each_t) {
                                                           left_join(temp_all_comb_1, list_predicted_probs[[each_t]])$output
                                                         })
                                   temp_list <- c(temp_list_0, temp_list_1)
                                   pmap_dbl(temp_list, prod) %>% sum %>% return
                                 })
  )
  left_join(data_temp, library_output)$output %>% return
}





get_H_current <- function(data_temp, variables, tau, list_predicted_probs) {
  # impute current data_temp
  
  temp_ind <- ncol(data_temp)
  
  all_observed_1 <- all_observed_0 <- data_temp
  loc_A <- which_variable_take(variables, "A")
  loc_A <- loc_A[loc_A <= temp_ind]
  all_observed_1[, loc_A] <- 1
  all_observed_0[, loc_A] <- 0
  
  # R and L (including Y before t = tau) (or not Z not A not t=0, not Y_tau) can be calculated in the same way
  if (temp_ind %in% which_variable_drop(variables, to_drop_variable = c("A", "Z"), to_drop_time = 0, to_drop = paste0("Y_", tau))) {
    loc_A <- which_variable_take(variables, "A")
    loc_A <- loc_A[loc_A < temp_ind]  # all needed A nodes
    loc_Z <- which_variable_take(variables, "Z")
    loc_Z <- loc_Z[loc_Z < temp_ind]  # all needed Z nodes
    
    # all predicted probs at these locs are needed
    # needed inputs are in either all_observed_1 or all_observed_0
    
    # these A probs will be taken as product
    part_A <- lapply(loc_A, function(k) {
      left_join(all_observed_1, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    # these ratios Z probs will be taken as product
    part_Z <- ifelse_vec(length(loc_Z) == 0, rep(1, length(part_A)), 
                         lapply(loc_Z, function(k) {
                           left_join(all_observed_0, list_predicted_probs[[k]])$output / 
                             left_join(all_observed_1, list_predicted_probs[[k]])$output
                         }) %>% pmap_dbl(prod))
    
    return(ifelse(data_temp[last(loc_A)] == 1, 1/part_A*part_Z, 0) %>% as.vector)
  }
  # Z nodes
  if (temp_ind %in% which_variable_take(variables, "Z")) {
    loc_A <- which_variable_take(variables, "A")
    loc_A <- loc_A[loc_A < temp_ind]
    # all the non-A, non-Z, not t=0, nodes before temp_ind; Y_tau won't be involved
    loc_LR <- which_variable_drop(variables, c("A", "Z"), 0)
    loc_LR <- loc_LR[loc_LR < temp_ind]
    
    # these A probs will be taken as product
    part_A <- lapply(loc_A, function(k) {
      left_join(all_observed_0, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    # these ratios Z probs will be taken as product
    part_LR <- lapply(loc_LR, function(k) {
      left_join(all_observed_1, list_predicted_probs[[k]])$output / 
        left_join(all_observed_0, list_predicted_probs[[k]])$output
    }) %>% pmap_dbl(prod)
    
    return(ifelse(data_temp[last(loc_A)] == 0, 1/part_A*part_LR, 0))
  }
}
