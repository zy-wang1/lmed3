setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# this is for local running; delete it for cluster file
setwd("..")  # to home
home <- getwd()  # specify it if needed

nCores <- 8

library(parallel)  # mclapply
library(abind)  # abind
library(xtable)  # xtable

# used for barebone lmed3 functions
library(data.table)
library(sl3)
library(tlverse)
library(dplyr)  # %>% 
# other dependency
library(R6)  # R6Class
library(uuid)  # UUIDgenerate
library(delayed)  # bundle_delayed
library(assertthat)  # assert_that
library(methods)  # is
library(hal9001)
library(ranger)
library(digest)

source(file.path(home, "code", "basic_functions-202008.R"))
source(file.path(home, "code", "generate_Zheng_data-202008.R"))
source(file.path(home, "code", "get_list_H.R"))
source(file.path(home, "code", "est_functions_tmle_1-202008.R"))
source(file.path(home, "code", "est_functions_tmle_2-202008.R"))
code_list <- list.files("./R", full.names = T)
for (code in code_list) source(code)

n_sim <- nCores * 13

timepoint <- 2
# configs for a barebone lmed3 non-targeting substitution estimator
# ARZLY model
node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1"
                  ,
                  A_2 = "A_2",
                  R_2 = "R_2",
                  Z_2 = "Z_2",
                  L_2 = "L1_2",
                  Y_2 = "Y_2"
)

  
# # only needed when debugging
# n_sim <- 24
# if_misspec <- F
# sample_size <- 500
# data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
# truth <- data_truth[[timepoint + 1]]$Y %>% mean
# truth

for (if_misspec in c(T
                     # , F
                     )) {
  data_truth <- generate_Zheng_data(B = 100000, tau = timepoint, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
  truth <- data_truth[[timepoint + 1]]$Y %>% mean
  truth
  
  for (sample_size in c(
    # 400
    #                     , 
    4000
                        )) {
    {
      start.time <- Sys.time()
      
      # sample_size <- 500
      RNGkind("L'Ecuyer-CMRG")
      set.seed(1234)
          
      
      results <- mclapply(1:n_sim, function(s) {
        data_sim <- generate_Zheng_data(B = sample_size, tau = timepoint, if_LY_misspec = if_misspec)
        temp_seq_reg <- est_seq_reg(data_sim = data_sim)
        temp_density_sub <- est_density_sub(data_sim = data_sim)
        temp_density_tmle <- est_density_tmle_1(data_sim = data_sim)
        temp_density_onestep <- est_density_tmle_2(data_sim = data_sim)
        
        {
          data_wide <- data.frame(data_sim)
          middle_spec <- lmed_middle(
            treatment_level = 1,
            control_level = 0
          )
          tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
          
          # choose base learners
          lrnr_mean <- make_learner(Lrnr_mean)
          # lrnr_glm <- make_learner(Lrnr_glm)
          lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
          lrnr_ranger50 <- make_learner(Lrnr_ranger, num.trees = 50)
          lrnr_hal_simple <- make_learner(Lrnr_hal9001, max_degree = 3, n_folds = 3)
          lrnr_lasso <- make_learner(Lrnr_glmnet) # al  pha default is 1
          lrnr_ridge <- make_learner(Lrnr_glmnet, alpha = 0)
          lrnr_elasticnet <- make_learner(Lrnr_glmnet, alpha = .5)
          learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
            learners = list(
              # lrnr_mean,
              # lrnr_glm, 
              lrnr_glm_fast
              # , 
              # lrnr_ranger50, 
              # lrnr_hal_simple
              # ,
              # lrnr_lasso, lrnr_ridge, lrnr_elasticnet
              )
          ))
          # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
          names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
          
          initial_likelihood <- middle_spec$make_initial_likelihood(
            tmle_task,
            learner_list
          )
          tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood)
          
          nontargeting <- tmle_params[[1]]$estimates(tmle_task)
          temp_lmed3_nontargeting <- nontargeting$psi
          temp_IC <- nontargeting$IC
          
          var_D <- var(temp_IC)
          n <- length(temp_IC)
          se <- sqrt(var_D / n)  
          
          
          # temp_var <- var(temp_IC)/length(temp_IC)
          CI2 <- temp_lmed3_nontargeting + 1.96 * se
          CI1 <- temp_lmed3_nontargeting - 1.96 * se
        }
        
        {
          # data_wide <- data.frame(data_sim)
          # 
          # middle_spec <- lmed_middle(
          #   treatment_level = 1,
          #   control_level = 0
          # )
          # tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
          # 
          # # choose base learners
          # lrnr_mean <- make_learner(Lrnr_mean)
          # lrnr_glm <- make_learner(Lrnr_glm)
          # lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
          # lrnr_hal_simple <- make_learner(Lrnr_hal9001, degrees = 1, n_folds = 2)
          # learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          #   learners = list(
          #     lrnr_mean,
          #     # lrnr_glm, 
          #     lrnr_hal_simple,
          #     lrnr_glm_fast)
          # ))          # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
          # names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
          # 
          # initial_likelihood <- middle_spec$make_initial_likelihood(
          #   tmle_task,
          #   learner_list
          # )
          updater <- lmed3_Update$new(maxit = 1, convergence_type = "scaled_var",
                                      fluctuation_type = "standard"
                                      # , cvtmle = T
                                      )
          targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
          tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
          updater$tmle_params <- tmle_params
          test <- fit_lmed3(tmle_task, targeted_likelihood, tmle_params, updater)
          firsttargeting <- test$estimates[[1]]
          temp_lmed3_est <- firsttargeting$psi
          temp_IC <- firsttargeting$IC
          
          var_D <- var(temp_IC)
          n <- length(temp_IC)
          se <- sqrt(var_D / n)  
          
          
          # temp_var <- var(temp_IC)/length(temp_IC)
          CI2_first <- firsttargeting$psi + 1.96 * se
          CI1_first <- firsttargeting$psi - 1.96 * se
        }
        
        {
          # data_wide <- data.frame(data_sim)
          # 
          # middle_spec <- lmed_middle(
          #   treatment_level = 1,
          #   control_level = 0
          # )
          # tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
          # 
          # # choose base learners
          # lrnr_mean <- make_learner(Lrnr_mean)
          # lrnr_glm <- make_learner(Lrnr_glm)
          # lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
          # lrnr_hal_simple <- make_learner(Lrnr_hal9001, degrees = 1, n_folds = 2)
          # learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          #   learners = list(
          #     lrnr_mean,
          #     # lrnr_glm, 
          #     lrnr_hal_simple,
          #     lrnr_glm_fast)
          # ))
          # # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
          # names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
          # 
          # initial_likelihood <- middle_spec$make_initial_likelihood(
          #   tmle_task,
          #   learner_list
          # )
          updater <- lmed3_Update$new(maxit = 10, convergence_type = "scaled_var",
                                      fluctuation_type = "standard")
          targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
          tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
          updater$tmle_params <- tmle_params
          test <- fit_lmed3(tmle_task, targeted_likelihood, tmle_params, updater)
          tentargeting <- test$estimates[[1]]
          temp_lmed3_est <- tentargeting$psi
          temp_IC <- tentargeting$IC
          
          var_D <- var(temp_IC)
          n <- length(temp_IC)
          se <- sqrt(var_D / n)  
          
          
          # temp_var <- var(temp_IC)/length(temp_IC)
          CI2_ten <- tentargeting$psi + 1.96 * se
          CI1_ten <- tentargeting$psi - 1.96 * se
        }
        
        {
          # data_wide <- data.frame(data_sim)
          # 
          # middle_spec <- lmed_middle(
          #   treatment_level = 1,
          #   control_level = 0
          # )
          # tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
          # 
          # # choose base learners
          # lrnr_mean <- make_learner(Lrnr_mean)
          # lrnr_glm <- make_learner(Lrnr_glm)
          # lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
          # lrnr_hal_simple <- make_learner(Lrnr_hal9001, degrees = 1, n_folds = 2)
          # learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          #   learners = list(
          #     lrnr_mean,
          #     # lrnr_glm, 
          #     lrnr_hal_simple,
          #     lrnr_glm_fast)
          # ))
          # # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
          # names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
          # 
          # initial_likelihood <- middle_spec$make_initial_likelihood(
          #   tmle_task,
          #   learner_list
          # )
          updater <- lmed3_Update$new(maxit = 100, convergence_type = "scaled_var",
                                      fluctuation_type = "standard", submodel_type = "onestep", d_epsilon = 0.01)
          targeted_likelihood <- Targeted_Likelihood$new(initial_likelihood, updater)
          tmle_params <- middle_spec$make_params(tmle_task, targeted_likelihood)
          updater$tmle_params <- tmle_params
          test <- fit_lmed3(tmle_task, targeted_likelihood, tmle_params, updater)
          onestep_targeting <- test$estimates[[1]]
          temp_lmed3_est <- onestep_targeting$psi
          temp_IC <- onestep_targeting$IC
          
          var_D <- var(temp_IC)
          n <- length(temp_IC)
          se <- sqrt(var_D / n)  
          
          
          # temp_var <- var(temp_IC)/length(temp_IC)
          CI2_onestep <- onestep_targeting$psi + 1.96 * se
          CI1_onestep <- onestep_targeting$psi - 1.96 * se
        }
        
        
        return(list(est = c(temp_seq_reg, 
                            temp_density_sub,
                            temp_density_tmle, 
                            temp_density_onestep,
                            nontargeting$psi, 
                            firsttargeting$psi, 
                            tentargeting$psi, 
                            onestep_targeting$psi
        ), 
        ci = c(rep(NA, 8), 
               CI1, CI2, 
               CI1_first, CI2_first, 
               CI1_ten, CI2_ten, 
               CI1_onestep, CI2_onestep)
        ))
      }, mc.cores = nCores)
      estimations <- lapply(results, function(x) x[[1]]) %>% abind(along = 0)
      cis <- lapply(results, function(x) x[[2]]) %>% abind(along = 0)
      
      end.time <- Sys.time()
      time.taken <- end.time - start.time
    }
    
    time.taken
    rowmax <- apply(estimations, 1, max)
    rowmax_noNA <- apply(estimations, 1, function(x) max(x, na.rm = T))
    ifNA <- is.na(rowmax)
    ifLarge <- rowmax_noNA > 10
    
    
    estimations <- estimations[!(ifNA & ifLarge), ]
    report <- data.frame(Bias = apply(estimations, 2, function(s) mean(s, na.rm = T)) - truth, 
                         lapply(1:ncol(estimations), function(which_col) c(mean((estimations[, which_col] - truth)^2, na.rm = T), 
                                                                           sd(estimations[, which_col], na.rm = T))) %>% abind(along = 0)
    )
    
    names(report)[2:3] <- c("MSE", "SD")
    # report <- report[-((nrow(report) - 1):nrow(report)), ]
    rownames(report) <- c("Non-targeted Sequential Regression", 
                          "Non-targeted Density",
                          "First-step Logistic MLE, Density",
                          "One-step MLE, Density",
                          "lmed3, Non-targeted", 
                          "lmed3, First-step Logistic", 
                          "lmed3, Iterative Logistic maxit 10", 
                          "lmed3, one-step common de 0.01, maxit = 100"
    )
    report <- report[, c(2, 1, 3)]
    
    report <- data.frame(report, 
                         coverage = sapply(1:nrow(report), function(k) {
                           mean(cis[, 2*k-1] < truth & cis[, 2*k] > truth)
                         })
    )
    report
    results
    ifNA
    ifLarge
    
    report %>% xtable(type = "latex", caption = paste0("Sample size ", sample_size, "; iteration: ", n_sim, "; run time: ", round(time.taken, 2), " ", units(time.taken),
                                                       "; NA and Large:  ",   sum(ifNA), " and ", sum(ifLarge)), digits = 6) %>% print(caption.placement = "top"
                                                                                                                                            ,
                                                                                                                                            file = paste0("./temp/", sample_size, "_LY_misspec_", ifelse(if_misspec, "misspecified", "correct"), "_sl_20200916_simplesl_oldseD.tex")
                                                       )
  }
}

