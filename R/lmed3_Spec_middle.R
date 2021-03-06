#' Defines a TML Estimator (except for the data)
#'
#'
#' @importFrom R6 R6Class
#'
#' @export
#
lmed3_Spec_middle <- R6Class(
  classname = "lmed3_Spec_middle",
  portable = TRUE,
  class = TRUE,
  inherit = lmed3_Spec,
  public = list(
    initialize = function(treatment_level, control_level, ...) {
      super$initialize(
        treatment_level = treatment_level,
        control_level = control_level, ...
      )
    },
    make_params = function(tmle_task, likelihood) {
      temp_names <- names(tmle_task$npsem)
      loc_A <- grep("A", temp_names)
      # ZW todo: in future can be dynamic
      treatment_value <- self$options$treatment_level
      control_value <- self$options$control_level
      A_levels <- tmle_task$npsem[[ temp_names[loc_A[1]] ]]$variable_type$levels
      if (!is.null(A_levels)) {
        treatment_value <- factor(treatment_value, levels = A_levels)
        control_value <- factor(control_value, levels = A_levels)
      }
      # list of intervention nodes as LF_static objects
      treatment <- lapply(temp_names[loc_A], function(eachA) {
        define_lf(LF_static, eachA, value = treatment_value)
      })
      control <- lapply(temp_names[loc_A], function(eachA) {
        define_lf(LF_static, eachA, value = control_value)
      })
      names(treatment) <- names(control) <- temp_names[loc_A]
      # treatment <- define_lf(LF_static, "A", value = treatment_value)
      # control <- define_lf(LF_static, "A", value = control_value)
      middle <- Param_middle$new(likelihood, treatment, control, outcome_node = last(temp_names))
      tmle_params <- list(middle)
      return(tmle_params)
    }
  ),
  active = list(),
  private = list()
)

#' All Treatment Specific Means
#'
#' O=(W,A,Y)
#' W=Covariates
#' A=Treatment (binary or categorical)
#' Y=Outcome (binary or bounded continuous)
#' @importFrom sl3 make_learner Lrnr_mean
#' @param treatment_level the level of A that corresponds to treatment
#' @param control_level the level of A that corresponds to a control or reference level
#' @export
lmed_middle <- function(treatment_level, control_level) {
  # TODO: unclear why this has to be in a factory function
  lmed3_Spec_middle$new(treatment_level, control_level)
}
