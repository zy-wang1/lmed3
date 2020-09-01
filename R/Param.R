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
      super$initialize(observed_likelihood, list(), outcome_node)
      private$.cf_likelihood_treatment <- CF_Likelihood$new(observed_likelihood, intervention_list_treatment)
      private$.cf_likelihood_control <- CF_Likelihood$new(observed_likelihood, intervention_list_control)
    },
    clever_covariates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
      
      pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      
      HA <- (cf_pA_treatment - cf_pA_control) / pA
      
      # ZW todo: get correct covariates
      temp_list <- list(rep(1, length(HA[[1]])))
      names(temp_list) <- outcome_node
      return(temp_list)
    },
    estimates = function(tmle_task = NULL, fold_number = "full") {
      if (is.null(tmle_task)) {
        tmle_task <- self$observed_likelihood$training_task
      }
      
      
      intervention_nodes <- union(names(self$intervention_list_treatment), names(self$intervention_list_control))
      
      # # clever_covariates happen here (for this param) only, but this is repeated computation
      # HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
      # 
      # 
      # # todo: make sure we support updating these params
      # pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      # cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)
      
      # todo: extend for stochastic
      cf_task_treatment <- self$cf_likelihood_treatment$enumerate_cf_tasks(tmle_task)[[1]]
      cf_task_control <- self$cf_likelihood_control$enumerate_cf_tasks(tmle_task)[[1]]
      
      Y <- tmle_task$get_tmle_node(self$outcome_node)
      
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
        tmle_task <- self$observed_likelihood$training_task
      }
      
      # all not A, not t=0 nodes
      temp_node_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_node_names)
      if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
      nodes_to_combo <- temp_node_names[ -c(loc_A, which(!if_not_0)) ]
      temp_combos <- self$observed_likelihood$get_possible_counterfactuals(nodes = nodes_to_combo)
      combos_treat <- combos_control <- data.frame(matrix(0, nrow(temp_combos), length(loc_A)), 
                                                   temp_combos)  # first tau cols are treat/control
      names(combos_treat)[1:length(loc_A)] <- names(combos_control)[1:length(loc_A)] <- temp_node_names[loc_A]
      # ZW todo: stochastic interventions
      combos_treat[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) self$intervention_list_treatment[[k]]$value %>% as.character %>% as.numeric)
      combos_control[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) self$intervention_list_control[[k]]$value %>% as.character %>% as.numeric)
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
        CF_Likelihood$new(self$observed_likelihood, s)
      })
      cf_list[[2]] <- lapply(temp_list[[2]], function(s) {
        CF_Likelihood$new(self$observed_likelihood, s)
      })
      # ZW todo: save combo fittings
      private$.cf_likelihood_combo_list <- cf_list
      
      
      

      
      cf_list <- self$cf_likelihood_combo_list

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
          task_1 <- cf_list[[1]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=1 for RLY nodes
          task_2 <- cf_list[[2]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=0 for Z nodes
          # take Z with control, take RLY with treat
          temp_return <- data.frame(self$observed_likelihood$get_likelihoods(task_1, temp_node_names[loc_RLY], fold_number), 
                                    self$observed_likelihood$get_likelihoods(task_2, temp_node_names[loc_Z], fold_number)
          )
          temp_return <- temp_return[names(current_combo)]  # match for next step
          # decide p or 1-p
          temp_return <- lapply(1:length(current_combo), function(k) if (current_combo[k] == 1) temp_return[k] else 1 - temp_return[k] ) %>% data.frame
          apply(temp_return, 1, prod)
        }
      })
      
      psi <- pmap_dbl(list_prods, sum) %>% mean

      # ZW todo: get true IC
      IC <- 0
      
      result <- list(psi = psi, IC = IC)
      return(result)
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
    cf_likelihood_combo_list = function() {
      return(private$.cf_likelihood_combo_list)
    },
    intervention_list_treatment = function() {
      return(self$cf_likelihood_treatment$intervention_list)
    },
    intervention_list_control = function() {
      return(self$cf_likelihood_control$intervention_list)
    },
    update_nodes = function() {
      return(c(self$outcome_node))
    }
  ),
  private = list(
    .type = "middle",
    .cf_likelihood_treatment = NULL,
    .cf_likelihood_control = NULL, 
    .cf_likelihood_combo_list = NULL
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