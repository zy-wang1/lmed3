#' Class for Likelihood
#'
#' This object represents an estimate of the relevant factors of the likelihood estimated from data, or based on \emph{a priori} knowledge where appropriate.
#' That is, it represents some subset of $P_n$. This object inherits from \code{\link[sl3]{Lrnr_base}}, and so shares some properties with \code{sl3} learners.
#' Specifically, to fit a likelihood object to data, one calls \code{likelihood$train(tmle3_task)}.
#' Each likelihood factor is represented by an object inheriting from \code{\link{LF_base}}.
#'
#' @docType class
#'
#' @importFrom R6 R6Class
#' @importFrom sl3 Lrnr_base
#' @importFrom assertthat assert_that is.count is.flag
#' @importFrom delayed bundle_delayed
#' @import data.table
#' @family Likelihood objects
#' @export
#'
#' @keywords data
#'
#' @return \code{Likelihood} object
#'
#' @format \code{\link{R6Class}} object.
#'
#' @template Likelihood_extra
#'
#' @export
Likelihood <- R6Class(
  classname = "Likelihood",
  portable = TRUE,
  class = TRUE,
  inherit = Lrnr_base,
  public = list(
    initialize = function(factor_list, cache = NULL, ...) {
      params <- args_to_list()
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }
      
      factor_names <- sapply(factor_list, `[[`, "name")
      names(factor_list) <- factor_names
      params$factor_list <- factor_list
      if (is.null(cache)) {
        cache <- Likelihood_cache$new()
      }
      private$.cache <- cache
      
      super$initialize(params)
    },
    print = function() {
      lapply(self$factor_list, print)
      invisible(NULL)
    },
    validate_task = function(tmle_task) {
      assert_that(is(tmle_task, "lmed3_Task"))
      
      factor_list <- self$factor_list
      factor_names <- names(factor_list)
      task_nodes <- names(tmle_task$npsem)
      if (!all(factor_names %in% task_nodes)) {
        stop("factor_list and task$npsem must have matching names")
      }
    },
    # get_likelihood = function(tmle_task, node, fold_number = "full") {
    #   likelihood_factor <- self$factor_list[[node]]
    #   # first check for cached values for this task
    #   likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number)
    #   
    #   if (is.null(likelihood_values)) {
    #     # if not, generate new ones
    #     likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number)
    #     self$cache$set_values(likelihood_factor, tmle_task, 0, fold_number, likelihood_values)
    #   }
    #   
    #   return(likelihood_values)
    # },
    # get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full") {
    #   if (is.null(nodes)) {
    #     nodes <- self$nodes
    #   }
    #   
    #   if (length(nodes) > 1) {
    #     all_likelihoods <- lapply(nodes, function(node) {
    #       self$get_likelihood(tmle_task, node, fold_number)
    #     })
    #     likelihood_dt <- as.data.table(all_likelihoods)
    #     setnames(likelihood_dt, nodes)
    #     return(likelihood_dt)
    #   } else {
    #     return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
    #   }
    # },
    get_likelihood = function(tmle_task, node, fold_number = "full") {
      likelihood_factor <- self$factor_list[[node]]
      # first check for cached values for this task
      likelihood_values <- self$cache$get_values(likelihood_factor, tmle_task, fold_number)
      
      if (is.null(likelihood_values)) {
        # if not, generate new ones
        likelihood_values <- likelihood_factor$get_likelihood(tmle_task, fold_number)
        self$cache$set_values(likelihood_factor, tmle_task, 0, fold_number, likelihood_values)
      }
      
      return(likelihood_values)
    },
    get_likelihoods = function(tmle_task, nodes = NULL, fold_number = "full") {
      if (is.null(nodes)) {
        nodes <- self$nodes
      }
      
      if (length(nodes) > 1) {
        all_likelihoods <- lapply(nodes, function(node) {
          self$get_likelihood(tmle_task, node, fold_number)
        })
        likelihood_dt <- as.data.table(all_likelihoods)
        setnames(likelihood_dt, nodes)
        return(likelihood_dt)
      } else {
        return(self$get_likelihood(tmle_task, nodes[[1]], fold_number))
      }
    },
    get_possible_counterfactuals = function(nodes = NULL) {
      
      # get factors for nodes
      factor_list <- self$factor_list
      if (!is.null(nodes)) {
        factor_list <- factor_list[nodes]
      }
      
      all_levels <- lapply(factor_list, function(likelihood_factor) {
        likelihood_factor$variable_type$levels
      })
      all_levels <- all_levels[ !(sapply(all_levels, is.null))]
      level_grid <- expand.grid(all_levels)
      return(level_grid)
    },
    base_train = function(task, pretrain) {
      self$validate_task(task)
      fit_object <- private$.train(task, pretrain)
      new_object <- self$clone() # copy parameters, and whatever else
      new_object$set_train(fit_object, task)
      return(new_object)
    },
    add_factors = function(factor_list) {
      if (inherits(factor_list, "LF_base")) {
        factor_list <- list(factor_list)
      }
      
      factor_names <- sapply(factor_list, `[[`, "name")
      
      # train factors if necessary
      factor_list <- lapply(factor_list, train_lf, self$training_task)
      
      # add factors to list of factors
      private$.params$factor_list[factor_names] <- factor_list
    }
  ),
  active = list(
    factor_list = function() {
      return(self$params$factor_list)
    },
    nodes = function() {
      return(names(self$factor_list))
    },
    cache = function() {
      return(private$.cache)
    }, 
    list_all_predicted_lkd = function() {
      if (!is.null(private$.list_all_predicted_lkd)) {
        return(private$.list_all_predicted_lkd)
      } else {
        # all not A, not t=0 nodes
        tmle_task <- self$training_task
        temp_node_names <- names(tmle_task$npsem)
        loc_A <- grep("A", temp_node_names)
        loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
        loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))
        if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
        
        # get list of all possible predicted lkds
        obs_data <- tmle_task$data
        obs_variable_names <- colnames(obs_data)
        list_all_predicted_lkd <- lapply(1:length(temp_node_names), function(loc_node) {
          if (loc_node > 1) {
            # currently only support univariate node for t>0
            current_variable <- tmle_task$npsem[[loc_node]]$variables
            temp_input <- expand_values(variables = obs_variable_names[1:which(obs_variable_names == current_variable)])  # all possible inputs
            temp_task <- lmed3_Task$new(temp_input, tmle_task$npsem[1:loc_node])
            temp_output <- self$factor_list[[loc_node]]$get_likelihood(temp_task, fold_number = "full")  # corresponding outputs
            data.frame(temp_input, output = temp_output) %>% return
          }
        })
        names(list_all_predicted_lkd) <- temp_node_names
        private$.list_all_predicted_lkd <- list_all_predicted_lkd
        return(list_all_predicted_lkd)  
      }
    }
  ),
  private = list(
    .train_sublearners = function(tmle_task) {
      factor_fits <- lapply(self$factor_list, function(factor) factor$delayed_train(tmle_task))
      result <- bundle_delayed(factor_fits)
      return(result)
    },
    .train = function(tmle_task, factor_fits) {
      factor_list <- self$factor_list
      for (i in seq_along(factor_list)) {
        factor_list[[i]]$train(tmle_task, factor_fits[[i]])
      }
      # TODO: mutating factor list of Lrnr_object instead of returning a fit
      #       which is not what sl3 Lrnrs usually do
      
      return("trained")
    },
    .predict = function(tmle_task) {
      stop("predict method doesn't work for Likelihood. See Likelihood$get_likelihoods for analogous method")
    },
    .chain = function(tmle_task) {
      stop("chain method doesn't work for Likelihood. Currently, no analogous functionality")
    },
    .list_all_predicted_lkd = NULL,
    .cache = NULL
  )
)

#' @param ... Passes all arguments to the constructor. See documentation for the
#'  Constructor below.
#'
#' @rdname Likelihood
#'
#' @export
#
make_Likelihood <- Likelihood$new



#' 
#' #' Cache Likelihood values, update those values
#' #' @docType class
#' #'
#' #' @importFrom R6 R6Class
#' #' @importFrom sl3 Lrnr_base
#' #' @importFrom assertthat assert_that is.count is.flag
#' #' @importFrom delayed bundle_delayed
#' #' @export
#' Likelihood_cache <- R6Class(
#'   classname = "Likelihood_cache",
#'   portable = TRUE,
#'   class = TRUE,
#'   public = list(
#'     initialize = function() {
#'       private$.cache <- new.env()
#'     },
#'     tasks_at_step = function(current_step) {
#'       self$tasks[task_uuids]
#'     },
#'     get_update_step = function(likelihood_factor, tmle_task, fold_number) {
#'       key <- self$key(likelihood_factor, tmle_task, fold_number)
#'       step_key <- sprintf("%s_%s", key, "step")
#'       get0(step_key, self$cache, inherits = FALSE)
#'     },
#'     key = function(likelihood_factor, tmle_task, fold_number) {
#'       key <- sprintf("%s_%s_%s", likelihood_factor$uuid, tmle_task$uuid, fold_number)
#'       return(key)
#'     },
#'     set_values = function(likelihood_factor, tmle_task, update_step = 0, fold_number, values) {
#'       self$cache_task(tmle_task)
#'       
#'       # respect likelihood factors that don't want to cache
#'       if (!likelihood_factor$cache) {
#'         return(0)
#'       }
#'       key <- self$key(likelihood_factor, tmle_task, fold_number)
#'       assign(key, values, self$cache)
#'       
#'       step_key <- sprintf("%s_%s", key, "step")
#'       assign(step_key, update_step, self$cache)
#'       
#'       return(1)
#'     },
#'     get_values = function(likelihood_factor, tmle_task, fold_number) {
#'       # matching_index <- self$find_match(likelihood_factor, tmle_task, fold_number)
#'       key <- self$key(likelihood_factor, tmle_task, fold_number)
#'       values <- get0(key, self$cache, inherits = FALSE)
#'       
#'       return(values)
#'     },
#'     cache_lf = function(likelihood_factor) {
#'       private$.lfs[likelihood_factor$uuid] <- likelihood_factor
#'     },
#'     cache_task = function(task) {
#'       if (!(task$uuid %in% names(private$.tasks))) {
#'         private$.tasks[[task$uuid]] <- task
#'       }
#'     }
#'   ),
#'   active = list(
#'     cache = function() {
#'       return(private$.cache)
#'     },
#'     tasks = function() {
#'       return(private$.tasks)
#'     }
#'   ),
#'   private = list(
#'     .tasks = list(),
#'     .lfs = list(),
#'     .cache = NULL
#'   )
#' )
#' 
#' 
#' 
#' 
