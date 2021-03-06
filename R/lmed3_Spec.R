#' Defines a TML Estimator (except for the data)
#'
#' Current limitations: pretty much tailored to \code{Param_TSM}
#'
#' @importFrom R6 R6Class
#'
#' @export
#
lmed3_Spec <- R6Class(
  classname = "lmed3_Spec",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(likelihood_override = NULL,
                          variable_types = NULL, ...) {
      private$.options <- list(
        likelihood_override = likelihood_override,
        variable_types = NULL, ...
      )
    },
    make_tmle_task = function(data, node_list, ...) {
      variable_types <- self$options$variable_types
      tmle_task <- middle_task(data, node_list, variable_types)
      return(tmle_task)
    },
    make_initial_likelihood = function(tmle_task, learner_list = NULL) {
      # produce trained likelihood when likelihood_def provided
      
      if (!is.null(self$options$likelihood_override)) {
        likelihood <- self$options$likelihood_override$train(tmle_task)
      } else {
        likelihood <- middle_likelihood(tmle_task, learner_list)  # see middle_helper
      }
      
      return(likelihood)
    },
    make_updater = function(...) {
      updater <- tmle3_Update$new(...)
    },
    make_targeted_likelihood = function(likelihood, updater) {
      targeted_likelihood <- Targeted_Likelihood$new(likelihood, updater)
      return(targeted_likelihood)
    },
    make_params = function(tmle_task, targeted_likelihood) {
      stop("this is a base class, try tsm_Spec_TSM_all")
    }
  ),
  active = list(
    options = function() {
      return(private$.options)
    }
  ),
  private = list(
    .options = NULL
  )
)