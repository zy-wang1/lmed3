setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# this is for local running; delete it for cluster file
setwd("..")  # to home
home <- getwd()  # specify it if needed

library(parallel)  # mclapply
library(abind)  # abind
library(xtable)  # xtable
nCores <- 8
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


source(file.path(home, "code", "basic_functions-202008.R"))
source(file.path(home, "code", "generate_Zheng_data-202008.R"))
source(file.path(home, "code", "get_list_H.R"))
source(file.path(home, "code", "est_functions_tmle_1-202008.R"))
source(file.path(home, "code", "est_functions_tmle_2-202008.R"))
code_list <- list.files("./R", full.names = T)
lapply(code_list, source)



data_truth <- generate_Zheng_data(B = 100000, tau = 2, seed = 202008, setAM = c(1, 0))
# data_truth <- generate_Zheng_data(B = 100000, tau = 2, seed = 202008, setAM = c(1, 1))
# data_truth <- generate_Zheng_data(B = 100000, tau = 2, seed = 202008, setAM = c(0, 0))
# data_truth %>% length
truth <- data_truth[[3]]$Y %>% mean

# checking generated data
glm(data_truth[[2]]$R ~ data_truth[[2]]$A + data_truth[[1]]$L1 + data_truth[[1]]$L2, family = binomial)
glm(data_truth[[2]]$Z ~ data_truth[[1]]$L2 + data_truth[[2]]$R, family = binomial)
glm(data_truth[[2]]$L ~ data_truth[[1]]$L2 + data_truth[[2]]$Z, family = binomial)
glm(data_truth[[2]]$Y ~ data_truth[[1]]$L2 + data_truth[[2]]$R + data_truth[[2]]$L1 + data_truth[[2]]$Z + data_truth[[2]]$A:data_truth[[2]]$Z, family = binomial)

glm(data_truth[[3]]$R ~ data_truth[[2]]$L1 + data_truth[[2]]$R, family = binomial)
glm(data_truth[[3]]$Z ~ data_truth[[1]]$L2 + data_truth[[3]]$R, family = binomial)
glm(data_truth[[3]]$L ~ data_truth[[1]]$L2 + data_truth[[2]]$Z, family = binomial)
glm(data_truth[[3]]$Y ~ data_truth[[1]]$L2 + data_truth[[3]]$R + data_truth[[3]]$L1 + data_truth[[3]]$Z + data_truth[[3]]$A:data_truth[[3]]$Z + data_truth[[2]]$R, family = binomial)

# configs for a barebone lmed3 non-targeting substitution estimator
# ARZLY model
node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1", 
                  A_2 = "A_2", 
                  R_2 = "R_2", 
                  Z_2 = "Z_2", 
                  L_2 = "L1_2", 
                  Y_2 = "Y_2" 
)




n_sim <- 200
# sample_size <- 50
for (sample_size in c(50
                      , 100, 400
)) {
  {
    start.time <- Sys.time()
    
    # sample_size <- 100
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
    
    results <- mclapply(1:n_sim, function(s) {
      data_sim <- generate_Zheng_data(B = sample_size, tau = 2)
      temp_seq_reg <- est_seq_reg(data_sim = data_sim)
      temp_density_sub <- est_density_sub(data_sim = data_sim)
      temp_density_tmle <- est_density_tmle_1(data_sim = data_sim)
      # temp_onestep <- est_density_onestep(data_sim = data_sim)
      
      {
        data_wide <- data.frame(data_sim)
        middle_spec <- lmed_middle(
          treatment_level = 1,
          control_level = 0
        )
        tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
        
        # choose base learners
        lrnr_mean <- make_learner(Lrnr_mean)
        lrnr_glm <- make_learner(Lrnr_glm)
        lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
        learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          learners = list(
            lrnr_mean,
            # lrnr_glm, 
            lrnr_glm_fast)
        ))
        # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
        names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
        
        initial_likelihood <- middle_spec$make_initial_likelihood(
          tmle_task,
          learner_list
        )
        
        tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood)[[1]]
        temp_lmed3_nontargeting <- tmle_params$estimates(tmle_task)$psi
      }
      
      return(c(temp_seq_reg, temp_density_sub, temp_density_tmle, temp_lmed3_nontargeting
      ))
    }, mc.cores = nCores)
    results <- results %>% abind(along = 0)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
  }
  
  time.taken
  report <- data.frame(Bias = apply(results, 2, function(s) mean(s, na.rm = T)) - truth, 
                       lapply(1:ncol(results), function(which_col) c(mean((results[, which_col] - truth)^2, na.rm = T), sd(results[, which_col], na.rm = T))) %>% abind(along = 0)
  )
  names(report)[2:3] <- c("MSE", "SD")
  rownames(report) <- c("Non-targeted Sequential Regression", "Non-targeted Density", "First-step Logistic MLE, Density", 
                        "Non-targeted lmed3 functions"
  )
  report <- report[, c(2, 1, 3)]
  report
  
  report %>% xtable(type = "latex", caption = paste0("Sample size ", sample_size, "; run time: ", round(time.taken, 2), " ", units(time.taken)), digits = 6) %>% print(caption.placement = "top",
                                                                                                                                                                       file = paste0("./temp/", sample_size, ".tex")
  )
}





# first one-step simulation
n_sim <- 8
# sample_size <- 100
for (sample_size in c(50
                      , 100, 400
)) {
  {
    start.time <- Sys.time()
    
    # sample_size <- 100
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
    
    results <- mclapply(1:n_sim, function(s) {
      data_sim <- generate_Zheng_data(B = sample_size, tau = 2)
      temp_seq_reg <- est_seq_reg(data_sim = data_sim)
      # temp_density_sub <- est_density_sub(data_sim = data_sim)
      # temp_density_tmle <- est_density_tmle_1(data_sim = data_sim)
      temp_onestep <- est_density_tmle_2(data_sim = data_sim)
      
      
      return(c(temp_seq_reg, 
               # temp_density_sub, temp_density_tmle, temp_lmed3_nontargeting
               temp_onestep
      ))
    }, mc.cores = nCores)
    results <- results %>% abind(along = 0)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
  }
  
  time.taken
  report <- data.frame(Bias = apply(results, 2, function(s) mean(s, na.rm = T)) - truth, 
                       lapply(1:ncol(results), function(which_col) c(mean((results[, which_col] - truth)^2, na.rm = T), sd(results[, which_col], na.rm = T))) %>% abind(along = 0)
  )
  names(report)[2:3] <- c("MSE", "SD")
  rownames(report) <- c("Non-targeted Sequential Regression", 
                        # "Non-targeted Density", "First-step Logistic MLE, Density", 
                        # "Non-targeted lmed3 functions"
                        "One-step MLE grid search"
  )
  report <- report[, c(2, 1, 3)]
  report
  
  report %>% xtable(type = "latex", caption = paste0("Sample size ", sample_size, "; run time: ", round(time.taken, 2), " ", units(time.taken)), digits = 6) %>% print(caption.placement = "top",
                                                                                                                                                                       file = paste0("./temp/", sample_size, "_onestep.tex")
  )
}





# LY model misspec

node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1", 
                  A_2 = "A_2", 
                  R_2 = "R_2", 
                  Z_2 = "Z_2", 
                  L_2 = "L1_2", 
                  Y_2 = "Y_2" 
)

data_truth <- generate_Zheng_data(B = 100000, tau = 2, seed = 202008, setAM = c(1, 0), if_LY_misspec = T)
truth <- data_truth[[3]]$Y %>% mean
truth

n_sim <- 48
# sample_size <- 50
for (sample_size in c(
  # 50
  # ,
  100, 
  200
  # , 400
  # 1000, 4000
)) {
  {
    start.time <- Sys.time()
    
    # sample_size <- 100
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
    
    results <- mclapply(1:n_sim, function(s) {
      data_sim <- generate_Zheng_data(B = sample_size, tau = 2, if_LY_misspec = T)
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
        lrnr_glm <- make_learner(Lrnr_glm)
        lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
        learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          learners = list(
            lrnr_mean,
            # lrnr_glm, 
            lrnr_glm_fast)
        ))
        # learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)  # simplest learner
        names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates
        
        initial_likelihood <- middle_spec$make_initial_likelihood(
          tmle_task,
          learner_list
        )
        
        tmle_params <- middle_spec$make_params(tmle_task, initial_likelihood)[[1]]
        nontargeting <- tmle_params$estimates(tmle_task)
        temp_lmed3_nontargeting <- nontargeting$psi
        temp_IC <- nontargeting$IC
        temp_var <- var(temp_IC)/length(temp_IC)
        CI2 <- temp_lmed3_nontargeting + 1.96 * sqrt(temp_var)
        CI1 <- temp_lmed3_nontargeting - 1.96 * sqrt(temp_var)
      }
      
      return(c(temp_seq_reg, temp_density_sub, temp_density_tmle, temp_density_onestep, temp_lmed3_nontargeting, CI1, CI2
      ))
    }, mc.cores = nCores)
    results <- results %>% abind(along = 0)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
  }
  
  time.taken
  rowmax <- apply(results, 1, max)
  rowmax_noNA <- apply(results, 1, function(x) max(x, na.rm = T))
  ifNA <- is.na(rowmax)
  ifLarge <- rowmax_noNA > 10
  
  
  # results <- results[!(ifNA & ifLarge), ]
  report <- data.frame(Bias = apply(results, 2, function(s) mean(s, na.rm = T)) - truth, 
                       lapply(1:ncol(results), function(which_col) c(mean((results[, which_col] - truth)^2, na.rm = T), 
                                                                     sd(results[, which_col], na.rm = T))) %>% abind(along = 0)
  )
  report
  
  names(report)[2:3] <- c("MSE", "SD")
  report <- report[1:5, ]
  rownames(report) <- c("Non-targeted Sequential Regression", "Non-targeted Density", "First-step Logistic MLE, Density", 
                        "One-step MLE, Density", 
                        "lmed3, Non-targeted"
  )
  report <- report[, c(2, 1, 3)]
  report <- data.frame(report, coverage = c(NA, NA, NA, NA, 
                                            mean(results[, 6] < truth & results[, 7] > truth)))
  
  
  report %>% xtable(type = "latex", caption = paste0("Sample size ", sample_size, "; run time: ", round(time.taken, 2), " ", units(time.taken), 
                                                     "; # of NA and Large:  ",   sum(ifNA), " and ", sum(ifLarge)), digits = 6) %>% print(caption.placement = "top",
                                                                                                                                          file = paste0("./temp/", sample_size, "LY_misspec_20200902.tex")
                                                     )
}










# test lmed3 functions

if_misspec <- F

node_list <- list(L_0 = c("L1_0", "L2_0"), 
                  A_1 = "A_1",
                  R_1 = "R_1",
                  Z_1 = "Z_1", 
                  L_1 = "L1_1", 
                  Y_1 = "Y_1", 
                  A_2 = "A_2", 
                  R_2 = "R_2", 
                  Z_2 = "Z_2", 
                  L_2 = "L1_2", 
                  Y_2 = "Y_2" 
)

data_truth <- generate_Zheng_data(B = 100000, tau = 2, seed = 202008, setAM = c(1, 0), if_LY_misspec = if_misspec)
truth <- data_truth[[3]]$Y %>% mean
truth

n_sim <- 8
sample_size <- 100
# for (sample_size in c(
#   # 50
#   # ,
#   100, 
#   200
#   # , 400
#   # 1000, 4000
# )) {
  {
    start.time <- Sys.time()
    
    # sample_size <- 100
    RNGkind("L'Ecuyer-CMRG")
    set.seed(123)
    
    results <- mclapply(1:n_sim, function(s) {
      data_sim <- generate_Zheng_data(B = sample_size, tau = 2, if_LY_misspec = if_misspec)
      temp_seq_reg <- est_seq_reg(data_sim = data_sim)
      # temp_density_sub <- est_density_sub(data_sim = data_sim)
      # temp_density_tmle <- est_density_tmle_1(data_sim = data_sim)
      # temp_density_onestep <- est_density_tmle_2(data_sim = data_sim)
      
      {
        data_wide <- data.frame(data_sim)
        middle_spec <- lmed_middle(
          treatment_level = 1,
          control_level = 0
        )
        tmle_task <- middle_spec$make_tmle_task(data_wide, node_list)
        
        # choose base learners
        lrnr_mean <- make_learner(Lrnr_mean)
        lrnr_glm <- make_learner(Lrnr_glm)
        lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
        learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
          learners = list(
            lrnr_mean,
            # lrnr_glm, 
            lrnr_glm_fast)
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
      
      return(c(temp_seq_reg, 
               # temp_density_sub, 
               # temp_density_tmle, temp_density_onestep, 
               temp_lmed3_nontargeting, CI1, CI2
      ))
    }, mc.cores = nCores)
    results <- results %>% abind(along = 0)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
  }
  
  time.taken
  rowmax <- apply(results, 1, max)
  rowmax_noNA <- apply(results, 1, function(x) max(x, na.rm = T))
  ifNA <- is.na(rowmax)
  ifLarge <- rowmax_noNA > 10
  
  
  # results <- results[!(ifNA & ifLarge), ]
  report <- data.frame(Bias = apply(results, 2, function(s) mean(s, na.rm = T)) - truth, 
                       lapply(1:ncol(results), function(which_col) c(mean((results[, which_col] - truth)^2, na.rm = T), 
                                                                     sd(results[, which_col], na.rm = T))) %>% abind(along = 0)
  )
  
  names(report)[2:3] <- c("MSE", "SD")
  report <- report[-((nrow(report) - 1):nrow(report)), ]
  rownames(report) <- c("Non-targeted Sequential Regression", 
                        # "Non-targeted Density", 
                        # "First-step Logistic MLE, Density",
                        # "One-step MLE, Density",
                        "lmed3, Non-targeted"
  )
  report <- report[, c(2, 1, 3)]
  report <- data.frame(report, coverage = c(NA, 
                                            # NA,
                                            # NA, NA, 
                                            mean(results[, ncol(results) - 1] < truth & results[, ncol(results)] > truth)))
  report
  results
  ifNA
  ifLarge
  
  # report %>% xtable(type = "latex", caption = paste0("Sample size ", sample_size, "; run time: ", round(time.taken, 2), " ", units(time.taken), 
  #                                                    "; # of NA and Large:  ",   sum(ifNA), " and ", sum(ifLarge)), digits = 6) %>% print(caption.placement = "top",
  #                                                                                                                                         file = paste0("./temp/", sample_size, "LY_misspec_20200902.tex")
  #                                                    )
# }
