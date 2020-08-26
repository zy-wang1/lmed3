# separated epsilon, logistic MLE, first step
est_density_tmle_1 <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  data_wide <- data.frame(data_sim)
  variables <- colnames(data_wide)
  tau <- length(data_sim) - 1  # the first slot in data_sim is L_0
  
  list_density <- lapply(1:length(variables), function(temp_loc) {
    if (temp_loc > ncol(data_sim[[1]])) {
      glm(current~., family = binomial, 
          data = data.frame(current = data_wide[[temp_loc]], data_wide[1:(temp_loc - 1)]))
    }
  })
  
  list_predicted_probs <- lapply(1:length(variables), function(temp_loc) {
    if (!is.null(list_density[[temp_loc]])) {
      temp_input <- expand_values(variables = variables[1:temp_loc])
      temp_output <- predict.glm(list_density[[temp_loc]], newdata = temp_input, type = "response")
      temp_output <- ifelse(temp_input[[temp_loc]] == 1, temp_output, 1 - temp_output)  # when current value is 0, takes the 1 - predicted
      # need to update for discretized values
      data.frame(temp_input, output = temp_output) %>% return
    }
  })
  
  # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1
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
  
  list_Q_1 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 1)
  list_Q_0 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 0)
  list_Q_diff <- map2(list_Q_1, list_Q_0, ~.x-.y)
  
  # targeting step; get 1 epsilon for each non-null list_H variable
  list_e <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      # for non-Z variables, take A=1 probs
      if(substr(variables[ind_var], 1, 1) == "Z") {
        temp_data <- all_observed_0
      } else {
        temp_data <- all_observed_1
      }
      
      temp_data[[ind_var]] <- 1  # this is updateing E Lt, so we need to set Lt = 1 when calling predicted probs
      
      temp_model <- glm(current~-1 + x + offset(off), family = binomial,
                        data = data.frame(current = data_wide[[ind_var]],  # later might need to transform it as indicators
                                          x = list_H[[ind_var]] * list_Q_diff[[ind_var]],
                                          off = left_join(temp_data, list_predicted_probs[[ind_var]])$output %>% logit
                        )
      )
      
      list_e[[ind_var]] <- coef(temp_model)[1]
    }
  }
  
  
  list_predicted_probs_updated <- list_predicted_probs
  # update list of probs by expit(logit(E Lt) + e * correct Q w.r.t. current lt )
  for (ind_var in 1:length(list_e)) {
    if(!is.null(list_e[[ind_var]])) {
      temp_current <- list_predicted_probs[[ind_var]]  # this is the full input output p^0 at current variable
      
      # where E Lt needs to be updated
      Xt_input <- temp_current[ind_var]
      
      loc_to_update <- Xt_input == 1
      df_to_update <- temp_current[loc_to_update, ]  # the possible input where L_t = 1
      df_to_update_0 <- df_to_update
      df_to_update_0[ind_var] <- 0  # replacing the above input with Lt = 0
      
      Q_current_1 <- get_Q_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs)  # Q with lt = 1
      Q_current_0 <- get_Q_current(df_to_update_0[1:ind_var], variables, tau, list_predicted_probs)  # Q with lt = 1
      Q_current <- Q_current_1 - Q_current_0
      
      H_current <- get_H_current(df_to_update[1:ind_var], variables, tau, list_predicted_probs)
      df_to_update$output <- expit(logit(df_to_update$output) + list_e[[ind_var]] * H_current * Q_current)
      raw_probs <- left_join(temp_current[1:(ind_var-1)], df_to_update)$output
      list_predicted_probs_updated[[ind_var]]$output <- ifelse(loc_to_update, raw_probs, 1-raw_probs)
    }
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
                                                       left_join(temp_all_comb_0, list_predicted_probs_updated[[each_t]])$output
                                                     })
                               temp_list_1 <- lapply(which_variable_drop(variables, c("A", "Z"), 0), 
                                                     function(each_t) {
                                                       left_join(temp_all_comb_1, list_predicted_probs_updated[[each_t]])$output
                                                     })
                               temp_list <- c(temp_list_0, temp_list_1)
                               pmap_dbl(temp_list, prod) %>% sum %>% return
                             })
  )
  
  if (warn == F) options(warn=0) else {}
  # substitution estimator
  left_join(data_sim[[1]], library_L0)$output %>% mean %>% return
  
}



# updated substitution estimator
est_density_sub <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  data_wide <- data.frame(data_sim)
  variables <- colnames(data_wide)
  
  list_density <- lapply(1:length(variables), function(temp_loc) {
    if (temp_loc > ncol(data_sim[[1]])) {
      glm(current~., family = binomial, 
          data = data.frame(current = data_wide[[temp_loc]], data_wide[1:(temp_loc - 1)]))
    }
  })
  
  list_predicted_probs <- lapply(1:length(variables), function(temp_loc) {
    if (!is.null(list_density[[temp_loc]])) {
      temp_input <- expand_values(variables = variables[1:temp_loc])
      temp_output <- predict.glm(list_density[[temp_loc]], newdata = temp_input, type = "response")
      temp_output <- ifelse(temp_input[[temp_loc]] == 1, temp_output, 1 - temp_output)  # when current value is 0, takes the 1 - predicted
      # need to update for discretized values
      data.frame(temp_input, output = temp_output) %>% return
    }
  })
  
  # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1
  all_possible_rlz_1 <- expand_values(variables, to_drop = c(1:length(data_sim[[1]]) ), 
                                      A = 1, 
                                      rule_variables = last(variables), rule_values = 1)
  all_possible_rlz_0 <- expand_values(variables, to_drop = c(1:length(data_sim[[1]]) ), 
                                      A = 0, 
                                      rule_variables = last(variables), rule_values = 1)
  
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
                                                       left_join(temp_all_comb_0, list_predicted_probs[[each_t]])$output
                                                     })
                               temp_list_1 <- lapply(which_variable_drop(variables, c("A", "Z"), 0), 
                                                     function(each_t) {
                                                       left_join(temp_all_comb_1, list_predicted_probs[[each_t]])$output
                                                     })
                               temp_list <- c(temp_list_0, temp_list_1)
                               pmap_dbl(temp_list, prod) %>% sum %>% return
                             })
  )

  if (warn == F) options(warn=0) else {}
  # substitution estimator
  left_join(data_sim[[1]], library_L0)$output %>% mean %>% return
}


expand_values <- function(variables, to_drop = NULL, values = NULL, rule_variables = NULL, rule_values = NULL, ...) {
  # variables is the vector of variable names
  # drop is either a vector of which variables to drop, or their indices
  # values are the possible values to expand; if null, generate binary values
  # ... are other value rules
  
  # input_list <- list(A = 1)
  
  input_list <- list(...)
  rules_list <- lapply(names(input_list), function(eachName) {
    if (length(grep("_", eachName)) > 0) eachName else variables[grep(eachName, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))]
  })
  if(!is.null(rule_variables)) {
    temp_list <- as.list(rule_values)
    names(temp_list) <- rule_variables
    input_list <- c(input_list, temp_list)
    rules_list <- c(rules_list, rule_variables)
  }
  
  # drop variables that don't want to generate
  # to_drop <- c("A", "L1_0")
  if (is.null(to_drop)) variables_to_generate <- variables else
    if(is.numeric(to_drop)) variables_to_generate <- variables[-to_drop] else 
      if (is.character(to_drop)) {
        ind_to_drop <- map(to_drop, ~which(variables == .x)) %>% compact %>% unlist
        if (length(ind_to_drop) == 0) variables_to_generate <- variables else 
          variables_to_generate <- variables[-ind_to_drop]
      } else {
        # consider how to identify error later
        print("invalid variables to drop")
        variables_to_generate <- variables
      }
  
  all_possible_values <- map(variables_to_generate, function(eachVar) {
    test_rules <- which(map_dbl(rules_list, ~length(grep(eachVar, .x))) != 0)
    if (length(test_rules) == 0) return(0:1) else {
      return(input_list[test_rules] %>% as.numeric)
    }
  }) %>% expand.grid() %>% data.frame
  colnames(all_possible_values) <- variables_to_generate
  
  # to add values 
  # to add variable name rules
  
  return(all_possible_values)
}

# values <- expand_values(variables, A = 1, rule_variables = last(variables), rule_values = 1)
# merge(values[1:5, ], values)

which_variable <- function(variables, target_variable, timepoint) {
  grep(paste0(target_variable, "_", timepoint), variables)
}

# get the orders of variables after droping by name or time
which_variable_drop <- function(variables, to_drop_variable = NULL, to_drop_time = NULL, to_drop = NULL) {
  if (!is.null(to_drop_variable)) {
    temp_drop_variable <- lapply(to_drop_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else 
    temp_drop_variable <- NULL
  if (!is.null(to_drop_time)) {
    temp_drop_time <- lapply(to_drop_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_drop_time <- NULL
  if (!is.null(to_drop)) {
    temp_drop_both <- lapply(to_drop, function(s) grep(s, variables)) %>% unlist
  } else temp_drop_both <- NULL
  temp_drop <- c(temp_drop_variable, temp_drop_time, temp_drop_both)
  if (!is.null(temp_drop)) (1:length(variables))[-temp_drop] %>% return else 1:length(variables)
}

# which_variable_drop(variables, "L")

# take variables by name or time
which_variable_take <- function(variables, to_take_variable = NULL, to_take_time = NULL, logic = "or") {
  if (!is.null(to_take_variable)) {
    temp_take_variable <- lapply(to_take_variable, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][1]))) %>% unlist
  } else 
    temp_take_variable <- NULL
  if (!is.null(to_take_time)) {
    temp_take_time <- lapply(to_take_time, function(s) grep(s, map_chr(variables, ~strsplit(.x, "_")[[1]][2]))) %>% unlist
  }  else temp_take_time <- NULL
  # default is to take all selected, by name or by time
  # and is to take the intersection
  if (logic == "or") temp_take <- c(temp_take_variable, temp_take_time) else 
    if (logic == "and") temp_take <- intersect(temp_take_variable, temp_take_time)
  if (!is.null(temp_take)) sort(unique(temp_take)) else NULL
}

# which_variable_take(variables, to_take_variable = "Z")
# which_variable_take(variables, to_take_time = 0)


# sequential regression
est_seq_reg <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  data_seq_reg <- data.frame(data_sim)
  
  # R_{tau+1}
  temp_loc <- ncol(data_seq_reg)
  temp_pred <- data_seq_reg[[temp_loc]]  
  temp_loc <- temp_loc - 1
  
  # L_{t}, Z_t, R_t to R_1
  while(temp_loc > ncol(data_sim[[1]])) {
    if (substr(colnames(data_seq_reg[temp_loc]) == "A", 1, 1)) {} else {
      temp_regressor <- data_seq_reg[1:(temp_loc - 1)]
      temp_model <- glm(current~., family = binomial(), 
                        data = data.frame(current = temp_pred, temp_regressor))
      temp_input <- temp_regressor
      if (substr(colnames(data_seq_reg)[temp_loc], 1, 1) == "Z") {
        temp_input[which(substr(colnames(temp_input), 1, 1) == "A")] <- 0  # A set to untreated for mediator models
      } else {
        temp_input[which(substr(colnames(temp_input), 1, 1) == "A")] <- 1  # A set to treated
      }
      temp_pred <- predict.glm(temp_model, newdata = temp_input, type = "response")  # update current prediction
      temp_loc <- temp_loc - 1
    }
  }
  if (warn == F) options(warn=0) else {}
  predict.glm(temp_model, newdata = data_sim[[1]], type = "response") %>% mean
}