
#' Likelihood Factor Estimated from Data using sl3.
#'
#' Uses an \code{sl3} learner to estimate a likelihood factor from data.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_fit, name, learner, ..., type = "density")}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{learner}}{An sl3 learner to be used to estimate the factor
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     \item{\code{type}}{character, either "density", for conditional density or, "mean" for conditional mean
#'     }
#'     }
#'
#' @section Fields:
#' \describe{
#'     \item{\code{learner}}{The learner or learner fit object}
#'     }
#'
#' @export
LF_fit <- R6Class(
  classname = "LF_fit",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, learner, ..., type = "density") {
      super$initialize(name, ..., type = type)
      private$.learner <- learner
    },
    delayed_train = function(tmle_task) {
      # just return prefit learner if that's what we have
      # otherwise, make a delayed fit and return that
      if (self$learner$is_trained) {
        return(self$learner)
      }

      outcome_node <- self$name
      
      # fit scaled task for bounded continuous
      learner_task <- tmle_task$get_regression_task(outcome_node, scale = TRUE)
      learner_fit <- delayed_learner_train(self$learner, learner_task)
      return(learner_fit)
    },
    train = function(tmle_task, learner_fit) {
      super$train(tmle_task)
      private$.learner <- learner_fit
    },
    get_mean = function(tmle_task, fold_number) {
      learner_task <- tmle_task$get_regression_task(self$name)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)
      
      # unscale preds (to handle bounded continuous)
      preds_unscaled <- tmle_task$unscale(preds, self$name)
      return(preds_unscaled)
    },
    get_density = function(tmle_task, fold_number) {
      learner_task <- tmle_task$get_regression_task(self$name)
      learner <- self$learner
      preds <- learner$predict_fold(learner_task, fold_number)
      
      outcome_type <- self$learner$training_task$outcome_type
      observed <- outcome_type$format(learner_task$Y)
      if (outcome_type$type == "binomial") {
        likelihood <- ifelse(observed == 1, preds, 1 - preds)
      } else if (outcome_type$type == "categorical") {
        unpacked <- sl3::unpack_predictions(preds)
        index_mat <- cbind(seq_along(observed), observed)
        likelihood <- unpacked[index_mat]
      } else if (outcome_type$type == "continuous") {
        likelihood <- unlist(preds)
      } else {
        stop(sprintf("unsupported outcome_type: %s", outcome_type$type))
      }
      return(likelihood)
    }
  ),
  active = list(
    learner = function() {
      return(private$.learner)
    }
  ),
  private = list(
    .name = NULL,
    .learner = NULL
  )
)




#' Likelihood Factor Estimated using Empirical Distribution
#'
#' Uses the empirical probability distribution (puts mass \eqn{1/n} on each of the observations, or uses weights if specified) to estimate a marginal density.
#' Inherits from \code{\link{LF_base}}; see that page for documentation on likelihood factors in general.
#' Only compatible with marginal likelihoods (no parent nodes). Only compatible with densities (no conditional means).
#' The \code{type} argument will be ignored if specified.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @section Constructor:
#'   \code{define_lf(LF_emp, name, ...)}
#'
#'   \describe{
#'     \item{\code{name}}{character, the name of the factor. Should match a node name in the nodes specified by \code{\link{tmle3_Task}$npsem}
#'     }
#'     \item{\code{...}}{Not currently used.
#'     }
#'     }
#'
#' @export
LF_emp <- R6Class(
  classname = "Lf_emp",
  portable = TRUE,
  class = TRUE,
  inherit = LF_base,
  public = list(
    initialize = function(name, ...) {
      super$initialize(name, ..., type = "density")
      private$.name <- name
    },
    get_mean = function(tmle_task, fold_number = "full") {
      stop("nothing to predict")
    },
    get_density = function(tmle_task, fold_number = "full") {
      weights <- tmle_task$weights
      return(weights / sum(weights))
    }
  ),
  active = list(),
  private = list(
    .name = NULL
  )
)




#' Base Class for Defining Likelihood Factors
#'
#' A Likelihood factor models a conditional density function.
#' The conditioning set is defined as all parent nodes (defined in \code{\link{tmle3_Task}}). In the case of a continuous
#' outcome variable, where a full density isn't needed, this can also model a conditional mean. This is the base class, which
#' is intended to be abstract. See below for a list of possible likelihood factor classes.
#'
#' @importFrom R6 R6Class
#' @importFrom uuid UUIDgenerate
#' @importFrom methods is
#' @family Likelihood objects
#' @keywords data
#'
#' @return \code{LF_base} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template LF_base_extra
#'
#' @export
LF_base <- R6Class(
  classname = "LF_base",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(name, bound = NULL, ..., type = "density", cache = TRUE) {
      private$.name <- name
      private$.type <- type
      private$.bound <- bound
      private$.uuid <- UUIDgenerate(use.time = TRUE)
      private$.cache <- cache
    },
    delayed_train = function(tmle_task) {
      return(list())
    },
    train = function(tmle_task, ...) {
      # get possible values from task if discrete
      tmle_node <- tmle_task$npsem[[self$name]]
      private$.variable_type <- tmle_node$variable_type
      
      # subclasses may do more, like fit sl3 models
    },
    get_density = function(tmle_task, fold_number) {
      stop("density not supported")
    },
    get_mean = function(tmle_task, fold_number) {
      stop("mean not supported")
    },
    get_likelihood = function(tmle_task, fold_number = "full") {
      if (self$type == "mean") {
        values <- self$get_mean(tmle_task, fold_number)
      } else {
        values <- self$get_density(tmle_task, fold_number)
      }
      if (!is.null(self$bound)) {
        values <- bound(values, self$bound)
      }
      
      return(values)
    },
    
    cf_values = function(tmle_task) {
      stop(sprintf("%s is not a valid intervention type", class(self)[1]))
    },
    print = function() {
      cat(sprintf("%s: %s\n", self$name, class(self)[1]))
    }
  ),
  active = list(
    name = function() {
      return(private$.name)
    },
    variable_type = function() {
      return(private$.variable_type)
    },
    type = function() {
      return(private$.type)
    },
    values = function() {
      variable_type <- self$variable_type
      if (!is.null(variable_type)) {
        return(variable_type$levels)
      } else {
        return(NULL)
      }
    },
    uuid = function() {
      return(private$.uuid)
    },
    bound = function() {
      return(private$.bound)
    },
    cache = function() {
      return(private$.cache)
    }
  ),
  private = list(
    .name = NULL,
    .variable_type = c(),
    .memoized_values = list(),
    .type = NULL,
    .uuid = NULL,
    .bound = NULL,
    .cache = TRUE
  )
)

#' Define a Likelihood Factor
#'
#' @param LF_class the class of likelihood factor. Should inherit from \code{\link{LF_base}}
#' @param ... arguments that define the likelihood factor. See the constructor for the specified \code{LF_class}.
#' @family Likelihood objects
#' @export
#
define_lf <- function(LF_class, ...) {
  return(LF_class$new(...))
}




