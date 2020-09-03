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