#' Defines an update procedure (submodel+loss function)
#'
#' Current Limitations:
#' loss function and submodel are hard-coded (need to accept arguments for these)
#' @section Constructor:
#'   \code{define_param(maxit, cvtmle, one_dimensional, constrain_step, delta_epsilon, verbose)}
#'
#'   \describe{
#'     \item{\code{maxit}}{The maximum number of update iterations
#'     }
#'     \item{\code{cvtmle}}{If \code{TRUE}, use CV-likelihood values when
#'        calculating updates.
#'     }
#'     \item{\code{one_dimensional}}{If \code{TRUE}, collapse clever covariates
#'        into a one-dimensional clever covariate scaled by the mean of their
#'        EIFs.
#'     }
#'     \item{\code{constrain_step}}{If \code{TRUE}, step size is at most
#'        \code{delta_epsilon} (it can be smaller if a smaller step decreases
#'        the loss more).
#'     }
#'     \item{\code{delta_epsilon}}{The maximum step size allowed if
#'        \code{constrain_step} is \code{TRUE}.
#'     }
#'     \item{\code{convergence_type}}{The convergence criterion to use: (1)
#'        \code{"scaled_var"} corresponds to sqrt(Var(D)/n)/logn (the default)
#'        while (2) \code{"sample_size"} corresponds to 1/n.
#'     }
#'     \item{\code{fluctuation_type}}{Whether to include the auxiliary covariate
#'        for the fluctuation model as a covariate or to treat it as a weight.
#'        Note that the option \code{"weighted"} is incompatible with a
#'        multi-epsilon submodel (\code{one_dimensional = FALSE}).
#'     }
#'     \item{\code{verbose}}{If \code{TRUE}, diagnostic output is generated
#'        about the updating procedure.
#'     }
#'     }
#'
#' @importFrom R6 R6Class
#'
#' @export
#
lmed3_Update <- R6Class(
  classname = "lmed3_Update",
  portable = TRUE,
  class = TRUE,
  public = list(
    initialize = function(maxit = 100, cvtmle = TRUE, one_dimensional = FALSE,
                          constrain_step = FALSE, delta_epsilon = 1e-4,
                          convergence_type = c("scaled_var", "sample_size"),
                          fluctuation_type = c("standard", "weighted"),
                          verbose = FALSE, 
                          d_epsilon = 0.01, 
                          submodel_type = NULL) {
      private$.maxit <- maxit
      private$.cvtmle <- cvtmle
      private$.one_dimensional <- one_dimensional
      private$.constrain_step <- constrain_step
      private$.delta_epsilon <- delta_epsilon
      private$.convergence_type <- match.arg(convergence_type)
      private$.fluctuation_type <- match.arg(fluctuation_type)
      private$.verbose <- verbose
      private$.d_epsilon <- d_epsilon
      private$.submodel_type <- submodel_type
    },
    collapse_covariates = function(estimates, clever_covariates) {
      ED <- ED_from_estimates(estimates)
      EDnormed <- ED / norm(ED, type = "2")
      collapsed_covariate <- clever_covariates %*% EDnormed
      
      return(collapsed_covariate)
    },
    update_step = function(likelihood, tmle_task, fold_number = "full") {
      if (self$submodel_type == "onestep") {
        # update likelihoods
        likelihood$update(new_epsilons = 0, self$step_number, fold_number)  # now full lkd list is updated too
        
        nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
          tmle_param$clever_covariates(tmle_task, fold_number
                                       , update = T
          )  # this updates the covariates; it does not call full list of lkd
        }))
        nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number
                               , update = T
          )  # this updates the D_list and results (est and ICs)
        }))
        
        if (fold_number != "full") {
          # update full fit likelihoods if we haven't already
          likelihood$update(new_epsilons = 0, self$step_number, "full")
        }
        # increment step count
        private$.step_number <- private$.step_number + 1
      } else {
        # get new submodel fit
        all_submodels <- self$generate_submodel_data(
          likelihood, tmle_task,
          fold_number
        )
        new_epsilons <- self$fit_submodels(all_submodels)
        
        # update likelihoods
        likelihood$update(new_epsilons, self$step_number, fold_number)  # now full lkd list is updated too
        
        nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
          tmle_param$clever_covariates(tmle_task, fold_number
                                       , update = T
          )  # this updates the covariates; it does not call full list of lkd
        }))
        nothing <- suppressWarnings(lapply(self$tmle_params, function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number
                               , update = T
          )  # this updates the D_list and results (est and ICs)
        }))
        
        if (fold_number != "full") {
          # update full fit likelihoods if we haven't already
          likelihood$update(new_epsilons, self$step_number, "full")  # now full lkd list is updated too
        }
        # increment step count
        private$.step_number <- private$.step_number + 1
      }
      
      
    },
    generate_submodel_data = function(likelihood, tmle_task,
                                      fold_number = "full") {
      update_nodes <- self$update_nodes
      
      # TODO: support not getting observed for case where we're applying
      #       updates instead of fitting them
      clever_covariates <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$clever_covariates(tmle_task, fold_number
                                     # , update = T
                                     )  # this returns the list of H
      })
      
      observed_values <- lapply(update_nodes, tmle_task$get_tmle_node)
      
      all_submodels <- lapply(update_nodes, function(update_node) {
        node_covariates <- lapply(clever_covariates, `[[`, update_node)
        covariates_dt <- do.call(cbind, node_covariates)  # there might be multiple targets
        # if (self$one_dimensional) {
        #   observed_task <- likelihood$training_task
        #   estimates <- lapply(self$tmle_params, function(tmle_param) {
        #     tmle_param$estimates(observed_task, fold_number)
        #   })
        #   covariates_dt <- self$collapse_covariates(estimates, covariates_dt)
        # }
        
        observed <- tmle_task$get_tmle_node(update_node)  # raw data on this node
        initial <- likelihood$get_likelihood(
          tmle_task, update_node,
          fold_number
        )  # observed (updated or initial) likelihood at this node
        initial <- ifelse(observed == 1, initial, 1 - initial)
        
        # scale observed and predicted values for bounded continuous
        observed <- tmle_task$scale(observed, update_node)
        initial <- tmle_task$scale(initial, update_node)
        
        # protect against qlogis(1)=Inf
        initial <- bound(initial, 0.005)
        
        submodel_data <- list(
          observed = observed,
          H = covariates_dt,
          initial = initial
        )
      })
      
      names(all_submodels) <- update_nodes
      
      return(all_submodels)
    },
    fit_submodel = function(submodel_data) {
      if (self$constrain_step) {
        ncol_H <- ncol(submodel_data$H)
        if (!(is.null(ncol_H) || (ncol_H == 1))) {
          stop(
            "Updater detected `constrain_step=TRUE` but multi-epsilon submodel.\n",
            "Consider setting `collapse_covariates=TRUE`"
          )
        }
        
        risk <- function(epsilon) {
          submodel_estimate <- self$apply_submodel(submodel_data, epsilon)
          loss <- self$loss_function(submodel_estimate, submodel_data$observed)
          mean(loss)
        }
        
        optim_fit <- optim(
          par = list(epsilon = self$delta_epsilon), fn = risk,
          lower = 0, upper = self$delta_epsilon,
          method = "Brent"
        )
        epsilon <- optim_fit$par
        risk_val <- optim_fit$value
        risk_zero <- risk(0)
        
        if (self$verbose) {
          cat(sprintf("risk_change: %e ", risk_val - risk_zero))
        }
      } else {
        if (self$fluctuation_type == "standard") {
          suppressWarnings({
            submodel_fit <- glm(observed ~ H - 1 
                                + offset(off), 
                                data.frame(submodel_data, 
                                           off = submodel_data$initial %>% logit
                                           ),
                                # offset = qlogis(submodel_data$initial),
                                family = binomial()
                                # ,
                                # start = rep(0, ncol(submodel_data$H)
                                            # )
            )
          })
        } else if (self$fluctuation_type == "weighted") {
          if (self$one_dimensional) {
            suppressWarnings({
              submodel_fit <- glm(observed ~ -1, submodel_data,
                                  offset = qlogis(submodel_data$initial),
                                  family = binomial(),
                                  weights = as.numeric(H),
                                  start = rep(0, ncol(submodel_data$H))
              )
            })
          } else {
            warning(
              "Updater detected `fluctuation_type='weighted'` but multi-epsilon submodel.\n",
              "This is incompatible. Proceeding with `fluctuation_type='standard'`."
            )
            suppressWarnings({
              submodel_fit <- glm(observed ~ H - 1, submodel_data,
                                  offset = qlogis(submodel_data$initial),
                                  family = binomial(),
                                  start = rep(0, ncol(submodel_data$H))
              )
            })
          }
        }
        epsilon <- coef(submodel_fit)
        
        # NOTE: this protects against collinear covariates
        # (which we don't care about, we just want an update)
        epsilon[is.na(epsilon)] <- 0
      }
      
      if (self$verbose) {
        cat(sprintf("epsilon: %e ", epsilon))
      }
      
      return(epsilon)
    },
    fit_submodels = function(all_submodels) {
      all_epsilon <- lapply(all_submodels, self$fit_submodel)
      
      names(all_epsilon) <- names(all_submodels)
      private$.epsilons <- c(private$.epsilons, list(all_epsilon))
      
      return(all_epsilon)
    },
    submodel = function(epsilon, initial, H) {
      plogis(qlogis(initial) + H %*% epsilon)
    },
    loss_function = function(estimate, observed) {
      -1 * ifelse(observed == 1, log(estimate), log(1 - estimate))
    },
    apply_submodel = function(submodel_data, epsilon) {
      self$submodel(epsilon, submodel_data$initial, submodel_data$H %>% as.matrix)
    },
    apply_update = function(tmle_task, likelihood, fold_number, all_epsilon) {
      update_nodes <- self$update_nodes
      
      # get submodel data for all nodes
      all_submodel_data <- self$generate_submodel_data(
        likelihood, tmle_task,
        fold_number
      )
      
      # needed to transform between lkd and probs
      observed_values <- lapply(update_nodes, tmle_task$get_tmle_node)
      names(observed_values) <- update_nodes
      
      # apply update to all nodes
      updated_likelihoods <- lapply(update_nodes, function(update_node) {
        submodel_data <- all_submodel_data[[update_node]]
        epsilon <- all_epsilon[[update_node]]
        updated_likelihood <- self$apply_submodel(submodel_data, epsilon)
        
        # we updated the predicted probabilities, or the non-zero lt likelihoods
        # for lt=0 they are just 1-p or 1- sum of other p
        updated_likelihood <- ifelse(observed_values[[update_node]] == 1, updated_likelihood, 1 - updated_likelihood)
        
        # un-scale to handle bounded continuous
        updated_likelihood <- tmle_task$unscale(
          updated_likelihood,
          update_node
        )
      })
      names(updated_likelihoods) <- update_nodes
      
      return(updated_likelihoods)
    },
    apply_update_onestep = function(tmle_task, likelihood, fold_number, d_epsilon) {
      update_nodes <- self$update_nodes
      
      # only need the list of D to get (1 + d_epsilon*D(p_k)) * p_k for each p_k
      list_D <- lapply(self$tmle_params, function(tmle_param) {
        tmle_param$list_D
      })
      
      observed_values <- lapply(update_nodes, tmle_task$get_tmle_node)  # needed to transform lkd to prob
      names(observed_values) <- update_nodes
      
      all_submodels <- lapply(update_nodes, function(update_node) {
        node_D <- lapply(list_D, `[[`, update_node)
        D_dt <- do.call(cbind, node_D)  # there might be multiple targets

        initial <- likelihood$get_likelihood(
          tmle_task, update_node,
          fold_number
        )  # observed (updated or initial) likelihood at this node
        
        # scale observed and predicted values for bounded continuous
        observed <- observed_values[[update_node]]
        # observed <- tmle_task$scale(observed, update_node)
        initial <- tmle_task$scale(initial, update_node)
        
        initial <- ifelse(observed == 1, initial, 1 - initial)
        
        # protect against qlogis(1)=Inf
        initial <- bound(initial, 0.005)
        
        submodel_data <- list(
          observed = observed,
          D = D_dt,
          initial = initial
        )
      })
      names(all_submodels) <- update_nodes
      
      # apply update to all nodes
      updated_likelihoods <- lapply(update_nodes, function(update_node) {
        submodel_data <- all_submodels[[update_node]]
        updated_likelihood <- (1 + submodel_data$D %*% d_epsilon) * submodel_data$initial
        # ZW todo: handle multiple targets
        
        # un-scale to handle bounded continuous
        updated_likelihood <- tmle_task$unscale(
          updated_likelihood,
          update_node
        )
        updated_likelihood <- ifelse(observed_values[[update_node]] == 1, updated_likelihood, 1 - updated_likelihood)
      })
      
      names(updated_likelihoods) <- update_nodes
      
      return(updated_likelihoods)
    },
    # use this to update list of full likelihood: list_all_predicted_lkd
    apply_update_full = function(tmle_task, likelihood, fold_number, all_epsilon) {
      update_nodes <- self$update_nodes
      list_all_predicted_lkd <- likelihood$list_all_predicted_lkd
      tmle_task <- likelihood$training_task
      temp_node_names <- names(tmle_task$npsem)
      obs_data <- tmle_task$data
      obs_variable_names <- names(obs_data)
      
      tmle_params <- self$tmle_params
      intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params$intervention_list_control))
      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      
      # apply update to all nodes
      updated_likelihoods <- lapply(update_nodes, function(update_node) {
        loc_node <- which(temp_node_names == update_node)
        
        # this is where full H list gets updated
        current_newH <- get_current_newH(loc_node, 
                         tmle_task, obs_data,
                         intervention_variables, intervention_levels_treat, intervention_levels_control, 
                         list_all_predicted_lkd
        )
        observed <- list_all_predicted_lkd[[update_node]][[ tmle_task$npsem[[update_node]]$variables ]] %>% as.vector()
        observed <- tmle_task$scale(observed, update_node)
        initial <- list_all_predicted_lkd[[update_node]]$output
        initial <- tmle_task$scale(initial, update_node)
        initial <- bound(initial, 0.005)
        submodel_data <- list(
          observed = observed,
          H = current_newH %>% as.matrix,
          initial = initial
        )
        submodel_data_1 <- list(
          observed = observed[observed == 1],
          H = current_newH %>% as.matrix,
          initial = initial[observed == 1]
        )
        
        epsilon <- all_epsilon[[update_node]]
        updated_likelihood_1 <- self$apply_submodel(submodel_data_1, epsilon)
        updated_likelihood <- submodel_data$initial
        updated_likelihood[observed == 1] <- updated_likelihood_1
        updated_likelihood[observed == 0] <- 1 - updated_likelihood_1  # list of probs are symmetric for now
        
        # un-scale to handle bounded continuous
        updated_likelihood <- tmle_task$unscale(
          updated_likelihood,
          update_node
        )
      })
      names(updated_likelihoods) <- update_nodes
      
      return(updated_likelihoods)
    },
    # use this to update list of full likelihood: list_all_predicted_lkd
    apply_update_full_onestep = function(tmle_task, likelihood, fold_number, d_epsilon) {
      update_nodes <- self$update_nodes
      list_all_predicted_lkd <- likelihood$list_all_predicted_lkd
      tmle_task <- likelihood$training_task
      temp_node_names <- names(tmle_task$npsem)
      obs_data <- tmle_task$data
      obs_variable_names <- names(obs_data)
      
      tmle_params <- self$tmle_params
      intervention_nodes <- union(names(tmle_params[[1]]$intervention_list_treatment), names(tmle_params$intervention_list_control))
      intervention_variables <- map_chr(tmle_task$npsem[intervention_nodes], ~.x$variables)
      intervention_variables_loc <- map_dbl(intervention_variables, ~grep(.x, obs_variable_names))
      intervention_levels_treat <- map_dbl(tmle_params[[1]]$intervention_list_treatment, ~.x$value %>% as.character %>% as.numeric)
      intervention_levels_control <- map_dbl(tmle_params[[1]]$intervention_list_control, ~.x$value %>% as.character %>% as.numeric)
      
      # apply update to all nodes
      updated_likelihoods <- lapply(update_nodes, function(update_node) {
        loc_node <- which(temp_node_names == update_node)
        # for lt = 1, calculate (1 + de D) p
        # for lt = 0, give 1 - updated p
        
        # this is where full H list gets updated
        current_newH <- get_current_newH(loc_node, 
                                         tmle_task, obs_data,
                                         intervention_variables, intervention_levels_treat, intervention_levels_control, 
                                         list_all_predicted_lkd
        )
        
        observed <- list_all_predicted_lkd[[update_node]][[tmle_task$npsem[[update_node]]$variables]] %>% as.vector()
        # observed <- tmle_task$scale(observed, update_node)
        initial <- list_all_predicted_lkd[[update_node]]$output
        initial <- tmle_task$scale(initial, update_node)
        initial <- bound(initial, 0.005)
        submodel_data <- list(
          observed = observed,
          H = current_newH %>% as.matrix,
          initial = initial
        )
        submodel_data_1 <- list(
          observed = observed[observed == 1],
          H = current_newH %>% as.matrix,
          initial = initial[observed == 1]
        )
        # for lt = 1, D = (observed - prob) * newH
        D <- (submodel_data_1$observed - submodel_data_1$initial) * submodel_data_1$H
        # p_{k+1} = (1 + D de) pk
        updated_likelihood_1 <- (1 + D * d_epsilon) * submodel_data_1$initial
        updated_likelihood <- submodel_data$initial
        updated_likelihood[observed == 1] <- updated_likelihood_1
        updated_likelihood[observed == 0] <- 1 - updated_likelihood_1  # list of probs are symmetric for now
        
        # un-scale to handle bounded continuous
        updated_likelihood <- tmle_task$unscale(
          updated_likelihood,
          update_node
        )
      })
      names(updated_likelihoods) <- update_nodes
      
      return(updated_likelihoods)
    },
    check_convergence = function(tmle_task, fold_number = "full") {
      estimates <- lapply(
        self$tmle_params,
        function(tmle_param) {
          tmle_param$estimates(tmle_task, fold_number = fold_number)
        }
      )
      
      if (self$convergence_type == "scaled_var") {
        # NOTE: the point of this criterion is to avoid targeting in an overly
        #       aggressive manner, as we simply need check that the following
        #       condition is met |P_n D*| / SE(D*) =< max(1/log(n), 1/10)
        IC <- do.call(cbind, lapply(estimates, `[[`, "IC"))
        se_Dstar <- sqrt(apply(IC, 2, var) 
                         / tmle_task$nrow
                         )
        ED_threshold <- se_Dstar / min(log(tmle_task$nrow), 10)
      } else if (self$convergence_type == "sample_size") {
        ED_threshold <- 1 / tmle_task$nrow
      }
      
      # get |P_n D*| of any number of parameter estimates
      ED <- ED_from_estimates(estimates)
      ED_criterion <- abs(ED)
      
      if (self$verbose) {
        cat(sprintf("max(abs(ED)): %e\n", ED_criterion))
      }
      return(all(ED_criterion <= ED_threshold))
    },
    update = function(likelihood, tmle_task) {
      update_fold <- self$update_fold
      maxit <- private$.maxit
      for (steps in seq_len(maxit)) {
        self$update_step(likelihood, tmle_task, update_fold)
        if (self$check_convergence(tmle_task, update_fold)) {
          break
        }
      }
    },
    register_param = function(new_params) {
      if (inherits(new_params, "Param_base")) {
        new_params <- list(new_params)
      }
      private$.tmle_params <- c(private$.tmle_params, new_params)
      new_update_nodes <- unlist(lapply(new_params, `[[`, "update_nodes"))
      private$.update_nodes <- unique(c(
        private$.update_nodes,
        new_update_nodes
      ))
    }
  ),
  active = list(
    epsilons = function() {
      return(private$.epsilons)
    },
    tmle_params = function(new_params = NULL) {
      if (!is.null(new_params)) {
        if (inherits(new_params, "Param_base")) {
          new_params <- list(new_params)
        }
        private$.tmle_params <- new_params
        private$.update_nodes <- unique(unlist(lapply(
          new_params, `[[`,
          "update_nodes"
        )))
      }
      return(private$.tmle_params)
    },
    update_nodes = function() {
      return(private$.update_nodes)
    },
    update_fold = function() {
      if (self$cvtmle) {
        # use training predictions on validation sets
        update_fold <- "validation"
      } else {
        # use predictions from full fit
        update_fold <- "full"
      }
    },
    step_number = function() {
      return(private$.step_number)
    },
    maxit = function() {
      return(private$.maxit)
    },
    cvtmle = function() {
      return(private$.cvtmle)
    },
    one_dimensional = function() {
      return(private$.one_dimensional)
    },
    constrain_step = function() {
      return(private$.constrain_step)
    },
    delta_epsilon = function() {
      return(private$.delta_epsilon)
    },
    convergence_type = function() {
      return(private$.convergence_type)
    },
    fluctuation_type = function() {
      return(private$.fluctuation_type)
    },
    verbose = function() {
      return(private$.verbose)
    }, 
    submodel_type = function() {
      return(private$.submodel_type)
    }, 
    d_epsilon = function() {
      return(private$.d_epsilon)
    }
  ),
  private = list(
    .epsilons = list(),
    .tmle_params = NULL,
    .update_nodes = NULL,
    .step_number = 0,
    .maxit = 100,
    .cvtmle = NULL,
    .one_dimensional = NULL,
    .constrain_step = NULL,
    .delta_epsilon = NULL,
    .convergence_type = NULL,
    .fluctuation_type = NULL,
    .verbose = FALSE, 
    .submodel_type = NULL, 
    .d_epsilon = 0.01
  )
)