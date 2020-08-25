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


# density substitution
est_density_sub <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  data_wide <- data.frame(data_sim)
  variables <- colnames(data_wide)
  
  # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1
  all_possible_rlz <- expand.grid(map(variables[-c(1:ncol(data_sim[[1]]), length(variables))], 
                                      ~ifelse_vec(substr(.x, 1, 1) == "A", 0, 0:1)))
  all_possible_rlz <- cbind(all_possible_rlz, 1)
  
  data_0s <- lapply(1:nrow(all_possible_rlz), function(which_row) {
    temp_data <- data_wide
    temp_values <- all_possible_rlz[which_row, ]
    # [substr(variables, 1, 1) != "A"]  # 
    temp_data[-(1:ncol(data_sim[[1]]))] <- temp_values
    temp_data %>% return
  })
  
  data_1s <- lapply(data_0s, function(each_df) {
    each_df[substr(colnames(each_df), 1, 1) == "A"] <- 1
    return(each_df)
  })
  
  list_density <- lapply(1:ncol(data_wide), function(temp_loc) {
    if (temp_loc > ncol(data_sim[[1]])) {
      glm(current~., family = binomial, 
          data = data.frame(current = data_wide[[temp_loc]], data_wide[1:(temp_loc - 1)]))
    }
  })
  
  to_be_summed <- lapply(1:length(data_0s), function(temp_loc) {
    temp_list_out <- map(1:length(list_density), function(.x) {
      current_density <- list_density[[.x]]
      # omit t = 0 and A nodes
      if(is.null(current_density)) NULL else if(substr(variables[.x], 1, 1) == "A") NULL else if(substr(variables[.x], 1, 1) == "Z") 
        # .x is the loc of current variable
        ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                   predict(current_density, newdata = data_0s[[temp_loc]], type = "response"), 
                   1 - predict(current_density, newdata = data_0s[[temp_loc]], type = "response")
                   )
         else 
           ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                      predict(current_density, newdata = data_1s[[temp_loc]], type = "response"), 
                      1 - predict(current_density, newdata = data_1s[[temp_loc]], type = "response")
           )
    })
    # current value can also be 0 or 1
    pmap_dbl(temp_list_out %>% compact, prod)
  })
  if (warn == F) options(warn=0) else {}
  pmap_dbl(to_be_summed, sum) %>% mean
}




# density, one-by-one maximize log likelihood
est_density_onestep <- function(data_sim, warn = F) {
  if (warn == F) options(warn=-1) else {}
  
  data_wide <- data.frame(data_sim)
  variables <- colnames(data_wide)
  
  # generate all possible 0, 1 valued rlz; set A to 0 first; set y_tau to always 1
  all_possible_rlz <- expand.grid(map(variables[-c(1:ncol(data_sim[[1]]), length(variables))], 
                                      ~ifelse_vec(substr(.x, 1, 1) == "A", 0, 0:1)))
  all_possible_rlz <- cbind(all_possible_rlz, 1)
  
  data_0s <- lapply(1:nrow(all_possible_rlz), function(which_row) {
    temp_data <- data_wide
    temp_values <- all_possible_rlz[which_row, ]
    # [substr(variables, 1, 1) != "A"]  # 
    temp_data[-(1:ncol(data_sim[[1]]))] <- temp_values
    temp_data %>% return
  })
  
  data_1s <- lapply(data_0s, function(each_df) {
    each_df[substr(colnames(each_df), 1, 1) == "A"] <- 1
    return(each_df)
  })
  
  list_density <- lapply(1:ncol(data_wide), function(temp_loc) {
    if (temp_loc > ncol(data_sim[[1]])) {
      glm(current~., family = binomial, 
          data = data.frame(current = data_wide[[temp_loc]], data_wide[1:(temp_loc - 1)]))
    }
  })
  
  to_be_summed <- lapply(1:length(data_0s), function(temp_loc) {
    temp_list_out <- map(1:length(list_density), function(.x) {
      current_density <- list_density[[.x]]
      # omit t = 0 and A nodes
      if(is.null(current_density)) NULL else if(substr(variables[.x], 1, 1) == "A") NULL else if(substr(variables[.x], 1, 1) == "Z") 
        # .x is the loc of current variable
        ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                   predict(current_density, newdata = data_0s[[temp_loc]], type = "response"), 
                   1 - predict(current_density, newdata = data_0s[[temp_loc]], type = "response")
        )
      else 
        ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                   predict(current_density, newdata = data_1s[[temp_loc]], type = "response"), 
                   1 - predict(current_density, newdata = data_1s[[temp_loc]], type = "response")
        )
    })
    # current value can also be 0 or 1
    temp_list_out
  })
  
  # to_be_summed %>% length
  # to_be_summed[[1]] %>% length
  # to_be_summed[[2]][[4]] %>% length
  
  lll <- to_be_summed %>% map_depth(.depth = 1, ~pmap_dbl(.x %>% compact, prod))

  
  # pmap_dbl(temp_list_out %>% compact, prod)
  
  if (warn == F) options(warn=0) else {}
  pmap_dbl(lll, sum) %>% mean
  
  
  # to_be_summed %>% vec_depth()  # 128 combos, variables, 400/1000 subjects
  
  # at each combo, at each relevant variable, need to update the prob
  # DV is the observed variable data value
  # off set is just the current list of probs
  # covariates involve last Q with current value inserted
  {
    # list of current Q functions; each Q function can be fed with newdata (as list)
    Q_functions <- lapply(1:length(list_density), function(loc_each_density) {
      # updated Q's, each take a new data_sim
      if(is.null(list_density[loc_each_density])) NULL else {
        function(data_sim) {
          data_wide <- data.frame(data_sim)
          variables <- colnames(data_wide)
          
          if(loc_each_density == length(variables)) {return(data_wide[[loc_each_density]])} else {
            # integrate out future (not including current) variables
            all_possible_rlz <- expand.grid(map(variables[-c(1:(loc_each_density-1), length(variables))], 
                                                ~ifelse_vec(substr(.x, 1, 1) == "A", 0, 0:1)))
            all_possible_rlz <- cbind(all_possible_rlz, 1)
            
            data_0s <- lapply(1:nrow(all_possible_rlz), function(which_row) {
              temp_data <- data_wide
              temp_values <- all_possible_rlz[which_row, ]
              # [substr(variables, 1, 1) != "A"]  # 
              temp_data[-(1:ncol(data_sim[[1]]))] <- temp_values
              temp_data %>% return
            })
            
            data_1s <- lapply(data_0s, function(each_df) {
              each_df[substr(colnames(each_df), 1, 1) == "A"] <- 1
              return(each_df)
            })
            
            to_be_summed <- lapply(1:length(data_0s), function(temp_loc) {
              temp_list_out <- map(1:length(list_density), function(.x) {
                current_density <- list_density[[.x]]
                # omit t = 0 and A nodes
                if(is.null(current_density)) NULL else if(substr(variables[.x], 1, 1) == "A") NULL else if(substr(variables[.x], 1, 1) == "Z") 
                  # .x is the loc of current variable
                  ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                             predict(current_density, newdata = data_0s[[temp_loc]], type = "response"), 
                             1 - predict(current_density, newdata = data_0s[[temp_loc]], type = "response")
                  )
                else 
                  ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                             predict(current_density, newdata = data_1s[[temp_loc]], type = "response"), 
                             1 - predict(current_density, newdata = data_1s[[temp_loc]], type = "response")
                  )
              })
              # current value can also be 0 or 1
              temp_list_out
            })
            lll <- to_be_summed %>% map_depth(.depth = 1, ~pmap_dbl(.x %>% compact, prod))
            pmap_dbl(lll, sum)
          }
          
        }
      }      
    })
  }
  
  # for each non-A non-L0 variable, there should also be a covariate function
  {
    # list of current Q functions; each Q function can be fed with newdata (as list)
    H_functions <- lapply(1:length(list_density), function(loc_each_density) {
      # updated H's, each take a new data_sim
      if(is.null(list_density[loc_each_density])) NULL else {
        function(data_sim) {
          data_wide <- data.frame(data_sim)
          variables <- colnames(data_wide)
          
          # no covariate for A's
          if(substr(variables[loc_each_density], 1, 1) == "A") return(NULL) else {
            if (substr(variables[loc_each_density], 1, 1) == "Z") {
              # need to consider both L and Y
              current_t <- strsplit(variables[loc_each_density], "_")[[1]][2] %>% as.numeric
              # A_t
              # data_wide[paste0("A_", current_t)]
              # two dataset to impute
              data_wide0 <- data_wide
              data_wide0[substr(variables, 1, 1) == "A"] <- 0
              data_wide1 <- data_wide
              data_wide1[substr(variables, 1, 1) == "A"] <- 1
              # product of p_A's; note this is control group probability
              temp_pA <- map(list_density[substr(variables, 1, 1) == "A"], ~1 - predict(.x, newdata = data_wide0, type = "response"))[1:current_t] %>% pmap_dbl(prod)
              # L ratios; only later L's; only to t-1
              temp1 <- map(list_density[substr(variables, 1, 1) == "L" & map_lgl(variables, ~strsplit(.x, "_")[[1]][2] != 0) ], ~predict(.x, newdata = data_wide1, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              temp2 <- map(list_density[substr(variables, 1, 1) == "L" & map_lgl(variables, ~strsplit(.x, "_")[[1]][2] != 0) ], ~predict(.x, newdata = data_wide0, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              temp3 <- map(list_density[substr(variables, 1, 1) == "Y" & map_lgl(variables, ~strsplit(.x, "_")[[1]][2] != 0) ], ~predict(.x, newdata = data_wide1, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              temp4 <- map(list_density[substr(variables, 1, 1) == "Y" & map_lgl(variables, ~strsplit(.x, "_")[[1]][2] != 0) ], ~predict(.x, newdata = data_wide0, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              temp_L_ratio <- temp1/temp2*temp3/temp4
              # R ratios; to t
              temp1 <- map(list_density[substr(variables, 1, 1) == "R"], ~predict(.x, newdata = data_wide1, type = "response"))[1:(current_t)] %>% pmap_dbl(prod)
              temp2 <- map(list_density[substr(variables, 1, 1) == "R"], ~predict(.x, newdata = data_wide0, type = "response"))[1:(current_t)] %>% pmap_dbl(prod)
              temp_R_ratio <- temp1/temp2
              ifelse(data_wide[paste0("A_", current_t)] == 0, 1/temp_pA*temp_L_ratio*temp_R_ratio, 0) %>% return
            } else if (substr(variables[loc_each_density], 1, 1) == "R") {
              current_t <- strsplit(variables[loc_each_density], "_")[[1]][2] %>% as.numeric
              # A_t's
              # data_wide[paste0("A_", current_t)]
              # two dataset to impute
              data_wide0 <- data_wide
              data_wide0[substr(variables, 1, 1) == "A"] <- 0
              data_wide1 <- data_wide
              data_wide1[substr(variables, 1, 1) == "A"] <- 1
              # product pA's
              temp_pA <- map(list_density[substr(variables, 1, 1) == "A"], ~predict(.x, newdata = data_wide1, type = "response"))[1:current_t] %>% pmap_dbl(prod)
              # ratio of pZ's; only to t-1
              temp1 <- map(list_density[substr(variables, 1, 1) == "Z"], ~predict(.x, newdata = data_wide0, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              temp2 <- map(list_density[substr(variables, 1, 1) == "Z"], ~predict(.x, newdata = data_wide1, type = "response"))[1:(current_t-1)] %>% pmap_dbl(prod)
              ifelse(data_wide[paste0("A_", current_t)] == 1, 1/temp_pA*temp1/temp2, 0) %>% return
            } else 
              # if(substr(variables[loc_each_density], 1, 1) == "L")   # L  or Y
                {
              current_t <- strsplit(variables[loc_each_density], "_")[[1]][2] %>% as.numeric
              # A_t's
              # data_wide[paste0("A_", current_t)]
              # two dataset to impute
              data_wide0 <- data_wide
              data_wide0[substr(variables, 1, 1) == "A"] <- 0
              data_wide1 <- data_wide
              data_wide1[substr(variables, 1, 1) == "A"] <- 1
              # product pA's
              temp_pA <- map(list_density[substr(variables, 1, 1) == "A"], ~predict(.x, newdata = data_wide1, type = "response"))[1:current_t] %>% pmap_dbl(prod)
              # ratio of pZ's
              temp1 <- map(list_density[substr(variables, 1, 1) == "Z"], ~predict(.x, newdata = data_wide0, type = "response"))[1:current_t] %>% pmap_dbl(prod)
              temp2 <- map(list_density[substr(variables, 1, 1) == "Z"], ~predict(.x, newdata = data_wide1, type = "response"))[1:current_t] %>% pmap_dbl(prod)
              ifelse(data_wide[paste0("A_", current_t)] == 1, 1/temp_pA*temp1/temp2, 0) %>% return
            }
            
            
          }
        }
      }      
    })
  }
  
  # D functions
  {
    # list of current D
    D_functions <- lapply(1:length(list_density), function(loc_each_density) {
      if(is.null(list_density[loc_each_density])) NULL else {
        function(data_sim) {
          data_wide <- data.frame(data_sim)
          variables <- colnames(data_wide)
          # no D for the last Y, or any A
          if(substr(variables[loc_each_density], 1, 1) == "A" | loc_each_density == length(variables)) {
            return(NULL)
            } else {
            if(substr(variables[loc_each_density], 1, 1) == "Y") {
              (Q_functions[[loc_each_density + 2]](data_sim) - Q_functions[[loc_each_density]](data_sim)) * H_functions[[loc_each_density]](data_sim)
            } else {
              (Q_functions[[loc_each_density + 1]](data_sim) - Q_functions[[loc_each_density]](data_sim)) * H_functions[[loc_each_density]](data_sim)
            }
          }

        }
      }      
    })
  }
  
  full_D <- function(data_sim) {
    data_wide <- data.frame(data_sim)
    variables <- colnames(data_wide)
    
    lapply(1:length(D_functions), function(k){
      if ((substr(variables, 1, 1) != "A" & map_lgl(variables, ~as.numeric(strsplit(.x, "_")[[1]][2]) != "0") & 
          variables != variables[length(variables)])[k]) D_functions[[k]](data_sim) else 0
    })
  }
  
  D_functions[[11]](data_sim)
  
  # detailed_D <- full_D(data_sim)
  # current_D <- pmap_dbl(detailed_D, sum)
  # 
  # epsilon <- 0.001
  
  p_full <- map(list_density %>% compact, ~predict.glm(.x, type = "response", newdata = data_wide)) %>% pmap_dbl(prod)
  
  updated_list_density <- lapply(1:length(list_density), function(k) {
    if(is.null(list_density[[k]])) NULL else if (is.null(D_functions[[k]](data_sim))) {
      function(data_sim, epsilon = 0) {
        data_wide <- data.frame(data_sim)
        predict.glm(list_density[[k]], type = "response", data_wide)
      }
      
    } else 
      function(data_sim, epsilon) {
        data_wide = data.frame(data_sim)
        return((1 + D_functions[[k]](data_sim) * epsilon) * predict.glm(list_density[[k]], type = "response", newdata = data_wide))
      }
  })
  
  search_func <- function(epsilon) {
    p_full <- map(updated_list_density %>% compact, ~.x(data_sim, epsilon)) %>% pmap_dbl(prod)
    sum((p_full))
  }
  
  
  recover_list <- function(x) {
    largest_t <- strsplit(colnames(x)[ncol(x)], "_")[[1]][2]
    list_table <- lapply(0:largest_t, function(t) {
      temp_table <- data_wide[as.numeric(map_chr(colnames(x), ~strsplit(.x, "_")[[1]][2])) == t]
      colnames(temp_table) <- colnames(x)[as.numeric(map_chr(colnames(x), ~strsplit(.x, "_")[[1]][2])) == t]
      temp_table
    })
    return(list_table)
  }

  
  to_be_summed <- lapply(1:length(data_0s), function(temp_loc) {
    temp_list_out <- map(1:length(list_density), function(.x) {
      current_density <- updated_list_density[[.x]]
      # omit t = 0 and A nodes
      if(is.null(current_density)) NULL else if(substr(variables[.x], 1, 1) == "A") NULL else if(substr(variables[.x], 1, 1) == "Z") {
        # .x is the loc of current variable
        ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                   current_density(recover_list(data_0s[[temp_loc]]), 0.001), 
                   1 - current_density(recover_list(data_0s[[temp_loc]]), 0.001)
        )
      } else 
        ifelse_vec(data_0s[[temp_loc]][[.x]][1] == 1, 
                   current_density(recover_list(data_1s[[temp_loc]]), 0.001), 
                   1 - current_density(recover_list(data_1s[[temp_loc]]), 0.001)
        )
    })
    # current value can also be 0 or 1
    temp_list_out
  })
  
  # to_be_summed %>% length
  # to_be_summed[[1]] %>% length
  # to_be_summed[[2]][[4]] %>% length
  
  lll <- to_be_summed %>% map_depth(.depth = 1, ~pmap_dbl(.x %>% compact, prod))
  
  
  # pmap_dbl(temp_list_out %>% compact, prod)
  
  if (warn == F) options(warn=0) else {}
  pmap_dbl(lll, sum) %>% mean
  
  
  
  # get updated densities
  
  
  # update all densities by (1 + eD) * P
  
  
  
}



