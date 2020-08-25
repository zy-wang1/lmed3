setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# this is for local running; delete it for cluster file
setwd("..")  # to home
home <- getwd()  # specify it if needed

library(parallel)  # mclapply
library(abind)  # abind
library(xtable)  # xtable
nCores <- 8

source(file.path(home, "code", "basic_functions-202008.R"))
source(file.path(home, "code", "generate_Zheng_data-202008.R"))
source(file.path(home, "code", "get_list_H.R"))
source(file.path(home, "code", "est_functions_tmle_1-202008.R"))


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

{
start.time <- Sys.time()

n_sim <- 200
sample_size <- 400
RNGkind("L'Ecuyer-CMRG")
set.seed(123)

results <- mclapply(1:n_sim, function(s) {
  data_sim <- generate_Zheng_data(B = sample_size, tau = 2)
  temp_seq_reg <- est_seq_reg(data_sim = data_sim)
  temp_density_sub <- est_density_sub(data_sim = data_sim)
  temp_density_tmle <- est_density_tmle_1(data_sim = data_sim)
  # temp_onestep <- est_density_onestep(data_sim = data_sim)
  return(c(temp_seq_reg, temp_density_sub, temp_density_tmle
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
rownames(report) <- c("Non-targeted Sequential Regression", "Non-targeted Density", "First-step Logistic MLE, Density")
report <- report[, c(2, 1, 3)]

# report50 <- report
# report100 <- report
report400 <- report

report400 %>% xtable(type = "latex", caption = "Sample size 400.", digits = 6) %>% print(caption.placement = "top")
report100 %>% xtable(type = "latex", caption = "Sample size 100.", digits = 6) %>% print(caption.placement = "top")
report50 %>% xtable(type = "latex", caption = "Sample size 50.", digits = 6) %>% print(caption.placement = "top")


# report50 <- report
# report100 <- report
# report500 <- report
report500 %>% xtable(type = "latex", caption = "Sample Size 500", digits = 6) %>% print(caption.placement = "top")
report100 %>% xtable(type = "latex", caption = "Sample Size 100", digits = 6) %>% print(caption.placement = "top")
report50 %>% xtable(type = "latex", caption = "Sample Size 50", digits = 6) %>% print(caption.placement = "top")

results %>% head


