#' Helper Functions for Point Treatment
#'
#' Handles the common W (covariates), A (treatment/intervention), Y (outcome) data structure
#'
#' @param data a \code{data.frame}, or \code{data.table} containing data for use in estimation
#' @param node_list a list of character vectors, listing the variables that comprise each node
#' @param variable_types a list of variable types, one for each node. If missing, variable types will be guessed
#' @param tmle_task a \code{\link{tmle3_Task}} as constructed via \code{point_tx_task}
#' @param learner_list a list of sl3 learners, one for A and one for Y to be used for likelihood estimation
#' @param ... extra arguments.
#' @export
#' @rdname point_tx
middle_npsem <- function(node_list, variable_types = NULL) {
  # make tmle_task
  npsem <- c(define_node("L_0", node_list$L_0, variable_type = variable_types$L_0),
             lapply(2:length(node_list), function(k) {
               if (k < length(node_list)) {
                 define_node(names(node_list)[k], 
                             node_list[[k]], 
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]])
               } else {
                 define_node(names(node_list)[k], 
                             node_list[[k]], 
                             names(node_list)[1:(k-1)],
                             variable_type = variable_types[[ names(node_list)[k] ]], 
                             scale = TRUE)
               }
             })
  )
  
  # npsem <- list(
  #   define_node("W", node_list$W, variable_type = variable_types$W),
  #   define_node("A", node_list$A, c("W"), variable_type = variable_types$A),
  #   define_node("Y", node_list$Y, c("A", "W"), variable_type = variable_types$Y, scale = TRUE)
  # )
  
  return(npsem)
}

#' @export
#' @rdname point_tx
middle_task <- function(data, node_list, variable_types = NULL, ...) {
  setDT(data)
  
  npsem <- middle_npsem(node_list, variable_types)
  
  if (!is.null(node_list$id)) {
    tmle_task <- lmed3_Task$new(data, npsem = npsem, id = node_list$id, ...)
  } else {
    tmle_task <- lmed3_Task$new(data, npsem = npsem, ...)
  }
  
  return(tmle_task)
}

#' @export
#' @rdname point_tx
middle_likelihood <- function(tmle_task, learner_list) {
  factor_list <- list()
  
  # covariates
  L_0_factor <- define_lf(LF_emp, "L_0")
  factor_list[[1]] <- L_0_factor
  
  # treatment (bound likelihood away from 0 (and 1 if binary))
  # A_type <- tmle_task$npsem[["A"]]$variable_type
  A_type <- tmle_task$npsem[[ grep("A", names(tmle_task$npsem))[1] ]]$variable_type  # evaluate the first A node
  if (A_type$type == "continous") {
    A_bound <- c(1 / tmle_task$nrow, Inf)
  } else if (A_type$type %in% c("binomial", "categorical")) {
    A_bound <- 0.025
  } else {
    A_bound <- NULL
  }
  
  temp_names <- names(tmle_task$npsem)
  loc_A <- grep("A", temp_names)
  factor_list[loc_A] <- lapply(loc_A, function(k) {
    define_lf(LF_fit, temp_names[k], learner = learner_list[[ temp_names[k] ]], bound = A_bound)
  })
  # A_factor <- define_lf(LF_fit, "A", learner = learner_list[["A"]], bound = A_bound)
  
  # others
  loc_others <- (1:length(temp_names))[-c(grep("A", temp_names), 1)]
  factor_list[loc_others] <- lapply(loc_others, function(k) {
    define_lf(LF_fit, temp_names[k], learner = learner_list[[ temp_names[k] ]], type = "mean")
  })
  
  # # outcome
  # Y_factor <- define_lf(LF_fit, "Y", learner = learner_list[["Y"]], type = "mean")
  # 
  # # construct and train likelihood
  # factor_list <- list(W_factor, A_factor, Y_factor)
  
  likelihood_def <- Likelihood$new(factor_list)
  likelihood <- likelihood_def$train(tmle_task)
  return(likelihood)
}

