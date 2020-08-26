# 1-by-1 logistic, one-step, jointly update
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
  
  list_Q_1 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 1)
  list_Q_0 <- get_list_Q(data_wide, variables, tau, list_predicted_probs, list_H, lt = 0)
  
  # targeting step; get 1 epsilon for each non-null list_H variable
  list_D_1 <- list_D_0 <- list()
  for (ind_var in 1:length(list_H)) {
    if(!is.null(list_H[[ind_var]])) {
      # for non-Z variables, take A=1 probs
      if(substr(variables[ind_var], 1, 1) == "Z") {
        temp_data <- all_observed_0
      } else {
        temp_data <- all_observed_1
      }
      
      # for each lt value
      temp_data[[ind_var]] <- 1  # lt
      # D with Lt == lt
      list_D_1[[ind_var]] <- ((data_wide[[ind_var]] == 1)*1 - left_join(temp_data, list_predicted_probs[[ind_var]])$output) * list_H[[ind_var]] * list_Q_1[[ind_var]]
      
      # for each lt value
      temp_data[[ind_var]] <- 0  # lt
      # D with Lt == lt
      list_D_0[[ind_var]] <- ((data_wide[[ind_var]] == 0)*1 - left_join(temp_data, list_predicted_probs[[ind_var]])$output) * list_H[[ind_var]] * list_Q_1[[ind_var]]
    }
  }
  
  list_observed_p <- list()
  
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
  left_join(data_sim[[1]], library_L0)$output
  
  
  
  
  list_predicted_probs_updated <- list_predicted_probs
  # update list of probs by expit(logit(prob) + correct e w.r.t. current lt * correct Q w.r.t. current lt )
  for (ind_var in 1:length(list_e_0)) {
    if(!is.null(list_e_0[[ind_var]])) {
      temp_current <- list_predicted_probs[[ind_var]]
      Q_current <- get_Q_current(temp_current[1:ind_var], variables, tau, list_predicted_probs)
      H_current <- get_H_current(temp_current[1:ind_var], variables, tau, list_predicted_probs)
      list_predicted_probs_updated[[ind_var]]$output <- expit(logit(temp_current$output) + 
                                                                ifelse(temp_current[[ind_var]] == 1, 
                                                                       list_e_1[[ind_var]], 
                                                                       list_e_0[[ind_var]]) * H_current * Q_current)
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

