#' Average Treatment Effect
#'
#' Parameter definition for the Average Treatment Effect (ATE).
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_param(Param_ATT, observed_likelihood, intervention_list, ..., outcome_node)}
#'
#'   \describe{
#'     \item{\code{observed_likelihood}}{A \code{\link{Likelihood}} corresponding to the observed likelihood
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention.
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention.
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{outcome_node}}{character, the name of the node that should be treated as the outcome
#'     }
#'     }
#'

#' @section Fields:
#' \describe{
#'     \item{\code{cf_likelihood_treatment}}{the counterfactual likelihood for the treatment
#'     }
#'     \item{\code{cf_likelihood_control}}{the counterfactual likelihood for the control
#'     }
#'     \item{\code{intervention_list_treatment}}{A list of objects inheriting from \code{\link{LF_base}}, representing the treatment intervention
#'     }
#'     \item{\code{intervention_list_control}}{A list of objects inheriting from \code{\link{LF_base}}, representing the control intervention
#'     }
#' }
#' @export
Param_middle <- R6Class(
  classname = "Param_middle",
  portable = TRUE,
  class = TRUE,
  inherit = Param_base,
  public = list(
    initialize = function(observed_likelihood, intervention_list_treatment, intervention_list_control, outcome_node = "Y") {
      super$initialize(observed_likelihood, list(), outcome_node = outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full", update = F) {
      # self$observed_likelihood -> initial_likelihood
      # self$ -> tmle_params$
      
      # if update == T, calculate based on current full lkd, and write .list_newH
      # if update == F, cache or calculate based on current full lkd in targeted likelihood
      if (!is.null(private$.list_newH) & update == F) {
        return(private$.list_newH)
      } else {
        if (is.null(tmle_task)) {
          tmle_task <- self$observed_likelihood$training_task
        }
        
        intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
        
        # todo: extend for stochastic
        cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
        cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
        
        Y <- tmle_task$get_tmle_node(self$outcome_node)
        
        # all not A, not t=0 nodes
        temp_node_names <- names(tmle_task$npsem)
        loc_A <- grep("A", temp_node_names)
        loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
        loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
        if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
        
        # get list of all possible predicted lkds
        obs_data <- tmle_task$data
        obs_variable_names <- colnames(obs_data)
        # ZW todo: see if observed_likelihood needs to change to targeted likelihood
        private$.list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd
        if (!is.null(private$.list_all_predicted_lkd)) {
          list_all_predicted_lkd <- private$.list_all_predicted_lkd
        } else {
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
          private$.list_all_predicted_lkd <- list_all_predicted_lkd
        }
        
        intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
        intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
        intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
        intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
        
        
        list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood, 
                            cf_task_treatment, cf_task_control, 
                            intervention_variables, intervention_levels_treat, intervention_levels_control)
        list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H, 
                              intervention_variables, intervention_levels_treat, intervention_levels_control, 
                              list_all_predicted_lkd, 
                              lt = 1)
        list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H, 
                              intervention_variables, intervention_levels_treat, intervention_levels_control, 
                              list_all_predicted_lkd, 
                              lt = 0)
        
        list_newH <- list()
        for (ind_var in 1:length(list_H)) {
          if(!is.null(list_H[[ind_var]])) {
            list_newH[[ind_var]] <- ( list_H[[ind_var]] * (list_Q_1[[ind_var]] - list_Q_0[[ind_var]]) ) %>% as.matrix
          }
        }
        names(list_newH) <- temp_node_names
        
        private$.list_newH <- list_newH
        return(list_newH)
      }
      
      
    },
    estimates = function(tmle_task = NULL, fold_number = "full", update = F) {
      # self$observed_likelihood -> initial_likelihood
      # self$ -> tmle_params[[1]]$
      
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
      
      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      
      # all not A, not t=0 nodes
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
      loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      
      # get list of all possible predicted lkds
      obs_data <- tmle_task$data
      obs_variable_names <- colnames(obs_data)
      
      private$.list_all_predicted_lkd <- self$observed_likelihood$list_all_predicted_lkd
      # only calculate list of lkd here when it is null; otherwise only update it in updater$apply_update_all
      if (!is.null(private$.list_all_predicted_lkd)) {
        list_all_predicted_lkd <- private$.list_all_predicted_lkd
      } else {
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
        private$.list_all_predicted_lkd <- list_all_predicted_lkd
      }
      
      # recalculate if list_D or result is null, or if we force it to update
      # make sure we force update this after each updating step
      # this helps speed up updater$check_convergence
      if (!is.null(private$.list_D) & !is.null(private$.result) & update == F) {
        return(private$.result)
      } else {
        intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
        intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
        intervention_levels_treat <- map_dbl(self$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
        intervention_levels_control <- map_dbl(self$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
        # nodes to integrate out in the target identification
        # only support univaraite node for now; assume treatment level is one
        all_possible_RZLY_1 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                             rule_variables = c(intervention_variables, 
                                                                last(obs_variable_names)), 
                                             rule_values = c(intervention_levels_treat, 1))
        all_possible_RZLY_0 <- expand_values(obs_variable_names, to_drop = c(1:length(tmle_task$npsem[[1]]$variables) ), 
                                             rule_variables = c(intervention_variables, 
                                                                last(obs_variable_names)), 
                                             rule_values = c(intervention_levels_control, 1))
        
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
        vec_est <- left_join(obs_data[, tmle_task$npsem[[1]]$variables, with = F], library_L0)$output
        psi <- mean(vec_est)
        
        # get true IC
        
        # # getting H's  
        # all_observed_1 <- all_observed_0 <- obs_data
        # for (temp_A in intervention_variables) {
        #   all_observed_1 <- all_observed_1 %>% mutate(!!temp_A := 1)
        #   all_observed_0 <- all_observed_0 %>% mutate(!!temp_A := 0)
        # }
        # 
        # list_H <- get_obs_H(tmle_task, obs_data, current_likelihood = self$observed_likelihood, 
        #                     cf_task_treatment, cf_task_control, 
        #                     intervention_variables, intervention_levels_treat, intervention_levels_control)
        # # get a list of needed deltaQ
        # list_Q_1 <- get_obs_Q(tmle_task, obs_data, list_H, 
        #                       intervention_variables, intervention_levels_treat, intervention_levels_control, 
        #                       list_all_predicted_lkd, 
        #                       lt = 1)
        # list_Q_0 <- get_obs_Q(tmle_task, obs_data, list_H, 
        #                       intervention_variables, intervention_levels_treat, intervention_levels_control, 
        #                       list_all_predicted_lkd, 
        #                       lt = 0)
        
        list_newH <- self$clever_covariates(tmle_task, fold_number)
        
        list_D <- list()
        for (ind_var in 1:length(list_newH)) {
          if(!is.null(list_newH[[ind_var]])) {
            # ZW todo: for discretized variables
            current_ind <- (obs_data[[tmle_task$npsem[[ind_var]]$variables]] == 1)*1
            if (ind_var %in% loc_Z) temp_p <- self$observed_likelihood$get_likelihoods(cf_task_control, temp_node_names[ind_var]) else 
              temp_p <- self$observed_likelihood$get_likelihoods(cf_task_treatment, temp_node_names[ind_var])
            temp_p <- ifelse(current_ind == 1, temp_p, 1 - ttemp_p)
            list_D[[ind_var]] <- (current_ind - temp_p) * 
              list_newH[[ind_var]]
            # list_H[[ind_var]] * (list_Q_1[[ind_var]] - list_Q_0[[ind_var]])
          }
        }
        list_D[[1]] <- vec_est - psi
        names(list_D) <- names(list_newH)
        
        vec_D <- list_D %>% compact %>% pmap_dbl(sum)
        IC <- vec_D
        
        result <- list(psi = psi, IC = IC)
        
        # these are cached; unless likelihood is updated, or we force it to update, they shouldn't be changed
        private$.list_D <- list_D
        private$.result <- result
        
        return(result)
      }
    }
  ),
  active = list(
    name = function() {
      param_form <- sprintf("ATE[%s_{%s}-%s_{%s}]", self$outcome_node, self$cf_likelihood_treatment$name, self$outcome_node, self$cf_likelihood_control$name)
      return(param_form)
    },
    cf_likelihood_treatment = function() {
      return(private$.cf_likelihood_treatment)
    },
    cf_likelihood_control = function() {
      return(private$.cf_likelihood_control)
    },
    list_all_predicted_lkd = function() {
      return(private$.list_all_predicted_lkd)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    list_D = function() {
      return(private$.list_D)
    },
    update_nodes = function() {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_update <- temp_node_names[if_not_0 & !((1:length(temp_node_names)) %in% loc_A)]
      # nodes_to_update <- nodes_to_update[-length(nodes_to_update)]
      return(nodes_to_update)
    }, 
    list_newH = function() {
      return(private$.list_newH)
    }
  ),
  private = list(
    .type = "middle",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL, 
    .list_all_predicted_lkd = NULL, 
    .list_newH = NULL,  # the clever covariates
    .list_D = NULL, 
    .result = NULL
  )
)




#' Base Class for Defining Parameters
#'
#' A parameter is a function of the likelihood. Once given a \code{\link{Likelihood}} object, a parameter will a value.
#' These objects also contain information about the efficient influence function (EIF) of a parameter, as well as its clever covariate(s).
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Parameters
#' @keywords data
#'
#' @return \code{Param_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Param_base_extra
#' @export
Param_base <- R6Class(
  classname = "Param_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(observed_likelihood, ..., outcome_node = "Y") {
      private$.observed_likelihood <- observed_likelihood
      private$.outcome_node <- outcome_node
      if (inherits(observed_likelihood, "Targeted_Likelihood")) {
        # register parameter with updater
        observed_likelihood$updater$register_param(self)
      } else if (inherits(observed_likelihood, "Likelihood")) {
        warning("Parameter was passed a non-Targeted Likelihood object so estimates cannot be updated from initial")
      } else {
        stop("Invalid Likelihood class: ", class(observed_likelihood))
      }
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      stop("Param_base is a base class")
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      stop("Param_base is a base class")
    },
    print = function() {
      cat(sprintf("%s: %s\n", class(self)[1], self$name))
    }
  ),
  active = list(
    name = function() {
      return(private$.type)
    },
    type = function() {
      return(private$.type)
    },
    
    observed_likelihood = function() {
      return(private$.observed_likelihood)
    },
    outcome_node = function() {
      return(private$.outcome_node)
    }
  ),
  private = list(
    .type = "undefined",
    .observed_likelihood = NULL,
    .outcome_node = NULL
  )
)


#' Define a Parameter
#'
#' @param Param_class the class of the Parameter. Should inherit from \code{\link{Param_base}}
#' @param ... arguments that define the parameter See the constructor for the specified \code{Parameter}.
#' @family Parameters
#' @export
#
define_param <- function(Param_class, ...) {
  return(Param_class$new(...))
}