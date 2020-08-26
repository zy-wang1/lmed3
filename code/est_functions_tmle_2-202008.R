# one-step MLE grid search, jointly update, (1 + de D)p
est_density_tmle_2 <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  # process simulated data
  data_wide <- data.frame(data_sim)
  variables <- colnames(data_wide)
  tau <- length(data_sim) - 1  # the first slot in data_sim is L_0
  
  # initial fit; only for t!=0 nodes
  list_density <- lapply(1:length(variables), function(temp_loc) {
    if (temp_loc > ncol(data_sim[[1]])) {
      glm(current~., family = binomial, 
          data = data.frame(current = data_wide[[temp_loc]], data_wide[1:(temp_loc - 1)]))
    }
  })
  
  # get initial predicted probs
  list_predicted_probs <- lapply(1:length(variables), function(temp_loc) {
    if (!is.null(list_density[[temp_loc]])) {
      temp_input <- expand_values(variables = variables[1:temp_loc])
      temp_output <- predict.glm(list_density[[temp_loc]], newdata = temp_input, type = "response")
      temp_output <- ifelse(temp_input[[temp_loc]] == 1, temp_output, 1 - temp_output)  # when current value is 0, takes the 1 - predicted
      # need to update for discretized values
      data.frame(temp_input, output = temp_output) %>% return
    }
  })
  
  # generate all possible 0, 1 valued RZL(Y); set A to 0 first; set y_tau to always 1, because y_tau=0 makes whole product 0
  all_possible_rlz_1 <- expand_values(variables, to_drop = c(1:length(data_sim[[1]]) ), 
                                      A = 1, 
                                      rule_variables = last(variables), rule_values = 1)
  all_possible_rlz_0 <- expand_values(variables, to_drop = c(1:length(data_sim[[1]]) ), 
                                      A = 0, 
                                      rule_variables = last(variables), rule_values = 1)
  
  # getting H's  
  all_observed_1 <- all_observed_0 <- data_wide
  all_observed_1[, which_variable_take(variables, "A")] <- 1
  all_observed_0[, which_variable_take(variables, "A")] <- 0
  # when A=1 and A=0 data tables are inserted the same as the first n-1 cols in list_predicted_probs, H will be saved as a function in the form of a list
  list_H <- get_list_H(all_observed_1, all_observed_0, data_wide, variables, tau, list_predicted_probs)
  # need to update where list_H is not NULL
  
  # the imputed Q_{R t+1} for Lt nodes; lt is the inserted Lt value
  list_Q_1 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 1)
  list_Q_0 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 0)
  
  
  # targeting step; get 1 epsilon for each non-null list_H variable
  # try p_update = (1 + D_current * 0.01)p_current and D_update = D(p_update)
  
  list_predicted_probs_update <- list_predicted_probs_current <- list_predicted_probs
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      temp_current <- list_predicted_probs_current[[ind_var]]  # this is the full input output p^0 at current variable
      
      # where E Lt needs to be updated
      Xt_input <- temp_current[ind_var]
      
      loc_to_update <- Xt_input == 1
      df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
      df_to_update_0 <- df_to_update
      df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0
      
      Q_current_1 <- get_Q_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs_current)  # Q with lt = 1
      Q_current_0 <- get_Q_current(df_to_update_0[1:ind_var], variables, tau, list_predicted_probs_current)  # Q with lt = 1
      Q_current <- Q_current_1 - Q_current_0
      
      H_current <- get_H_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs_current)
      df_to_update$output <- (1 + 0.01 * H_current * Q_current * (1 - temp_current$output[loc_to_update])) * temp_current$output[loc_to_update]
        # expit(logit(df_to_update$output) + list_e[[ind_var]] * H_current * Q_current)
      raw_probs <- left_join(temp_current[1:(ind_var-1)], df_to_update)$output
      list_predicted_probs_update[[ind_var]]$output <- ifelse(loc_to_update, raw_probs, 1-raw_probs)
    }
  }

  # get list of observed updated D parts; for convergence
  list_D <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      # for non-Z variables, use A=1 data; otherwise use A=0 data
      if(substr(variables[ind_var], 1, 1) == "Z") {
        temp_data <- all_observed_0
      } else {
        temp_data <- all_observed_1
      }
      # calculate D's except for D_{L_0}; 
      temp_data[[ind_var]] <- 1  # update prob of being 1
      # for each Z or RLY node/variable, call predicted probs and decide p or 1-p
      temp_p <- left_join(temp_data, list_predicted_probs_update[[ind_var]])$output
      temp_p <- ifelse(temp_data[[ind_var]] == 1, temp_p, 1-temp_p)
      # for now, these are binary variables
      # get updated lists first
      list_H <- get_list_H(all_observed_1, all_observed_0, data_wide, variables, tau, list_predicted_probs_update)
      list_Q_1 <- get_list_Q(data_wide, variables, tau, list_predicted_probs_update, list_H, lt = 1)
      list_Q_0 <- get_list_Q(data_wide, variables, tau, list_predicted_probs_update, list_H, lt = 0)
      
      list_D[[ind_var]] <- (data_wide[[ind_var]] - temp_p) * 
        list_H[[ind_var]] * (list_Q_1[[ind_var]] - list_Q_0[[ind_var]])
    }
  }
  # stop if D(p_update) < sd(D(p_update)) / log(n) 
  vec_D <- list_D %>% compact %>% pmap_dbl(sum)
  counter <- 1
  while (abs(mean(vec_D)) > sd(vec_D)/log(length(vec_D)) & counter <= 10) {
    list_predicted_probs_current <- list_predicted_probs_update
    for (ind_var in 1:length(list_H)) {
      if(!is.null(list_H[[ind_var]])) {
        temp_current <- list_predicted_probs_current[[ind_var]]  # this is the full input output p^0 at current variable
        
        # where E Lt needs to be updated
        Xt_input <- temp_current[ind_var]
        
        loc_to_update <- Xt_input == 1
        df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
        df_to_update_0 <- df_to_update
        df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0
        
        Q_current_1 <- get_Q_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs_current)  # Q with lt = 1
        Q_current_0 <- get_Q_current(df_to_update_0[1:ind_var], variables, tau, list_predicted_probs_current)  # Q with lt = 1
        Q_current <- Q_current_1 - Q_current_0
        
        H_current <- get_H_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs_current)
        df_to_update$output <- (1 + 0.01 * H_current * Q_current * (1 - temp_current$output[loc_to_update])) * temp_current$output[loc_to_update]
        # expit(logit(df_to_update$output) + list_e[[ind_var]] * H_current * Q_current)
        raw_probs <- left_join(temp_current[1:(ind_var-1)], df_to_update)$output
        list_predicted_probs_update[[ind_var]]$output <- ifelse(loc_to_update, raw_probs, 1-raw_probs)
      }
    } 
    for (ind_var in 1:length(list_H)) {
      if(!is.null(list_H[[ind_var]])) {
        # for non-Z variables, use A=1 data; otherwise use A=0 data
        if(substr(variables[ind_var], 1, 1) == "Z") {
          temp_data <- all_observed_0
        } else {
          temp_data <- all_observed_1
        }
        # calculate D's except for D_{L_0}; 
        temp_data[[ind_var]] <- 1  # update prob of being 1
        # for each Z or RLY node/variable, call predicted probs and decide p or 1-p
        temp_p <- left_join(temp_data, list_predicted_probs_update[[ind_var]])$output
        temp_p <- ifelse(temp_data[[ind_var]] == 1, temp_p, 1-temp_p)
        # for now, these are binary variables
        # get updated lists first
        list_H <- get_list_H(all_observed_1, all_observed_0, data_wide, variables, tau, list_predicted_probs_update)
        list_Q_1 <- get_list_Q(data_wide, variables, tau, list_predicted_probs_update, list_H, lt = 1)
        list_Q_0 <- get_list_Q(data_wide, variables, tau, list_predicted_probs_update, list_H, lt = 0)
        
        list_D[[ind_var]] <- (data_wide[[ind_var]] - temp_p) * 
          list_H[[ind_var]] * (list_Q_1[[ind_var]] - list_Q_0[[ind_var]])
      }
    }
    vec_D <- list_D %>% compact %>% pmap_dbl(sum)
    counter <- counter + 1
  }
  
  
  
  
  
  # for each observed L_0 vector, generate all needed combinations, one version for A = 1, one version for A = 0
  unique_L0 <- data_sim[[1]] %>% unique
  library_L0 <- data.frame(unique_L0, output = 
                             map_dbl(1:nrow(unique_L0), function(which_row) {
                               temp_all_comb_0 <- cbind(unique_L0[which_row, ], all_possible_rlz_0)
                               temp_all_comb_1 <- cbind(unique_L0[which_row, ], all_possible_rlz_1)
                               # for all non-A, non-0 variables, calculate the variable by rule
                               # for Z's, use A = 0 values; outputs are predicted probs at each possible comb
                               temp_list_0 <- lapply(which_variable_take(variables, to_take_variable = "Z"), 
                                                     function(each_t) {
                                                       left_join(temp_all_comb_0, list_predicted_probs_update[[each_t]])$output
                                                     })
                               temp_list_1 <- lapply(which_variable_drop(variables, c("A", "Z"), 0), 
                                                     function(each_t) {
                                                       left_join(temp_all_comb_1, list_predicted_probs_update[[each_t]])$output
                                                     })
                               temp_list <- c(temp_list_0, temp_list_1)
                               pmap_dbl(temp_list, prod) %>% sum %>% return
                             })
  )
  
  if (warn == F) options(warn=0) else {}
  # substitution estimator
  left_join(data_sim[[1]], library_L0)$output %>% mean %>% return
  
}

