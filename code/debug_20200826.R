tmle_task <- middle_spec$make_tmle_task(data_wide, node_list, variable_types = "binomial")

# choose base learners
lrnr_mean <- make_learner(Lrnr_mean)
lrnr_glm <- make_learner(Lrnr_glm)
lrnr_glm_fast <- Lrnr_glm_fast$new(outcome_type = "binomial")
# learner_list <- lapply(1:length(tmle_task$npsem), function(s) Lrnr_sl$new(
#   learners = list(
#     lrnr_mean,
#                   lrnr_glm)
# ))
learner_list <- lapply(1:length(tmle_task$npsem), function(s) lrnr_glm_fast)
names(learner_list) <- names(tmle_task$npsem)  # the first will be ignored; empirical dist. will be used for covariates

initial_likelihood <- middle_spec$make_initial_likelihood(
  tmle_task,
  learner_list
)

temp_names <- names(tmle_task$npsem)
loc_A <- grep("A", temp_names)
# ZW todo: in future can be dynamic
treatment_value <- 1
control_value <- 0
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


intervention_nodes <- union(names(treatment), names(control))

# # clever_covariates happen here (for this param) only, but this is repeated computation
# HA <- self$clever_covariates(tmle_task, fold_number)[[self$outcome_node]]
# 
# 
# # todo: make sure we support updating these params
# pA <- self$observed_likelihood$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_treatment <- self$cf_likelihood_treatment$get_likelihoods(tmle_task, intervention_nodes, fold_number)
# cf_pA_control <- self$cf_likelihood_control$get_likelihoods(tmle_task, intervention_nodes, fold_number)

cf_lkd_trt <- CF_Likelihood$new(initial_likelihood, treatment)
cf_lkd_ctrl <- CF_Likelihood$new(initial_likelihood, control)

# todo: extend for stochastic
cf_task_treatment <- cf_lkd_trt$enumerate_cf_tasks(tmle_task)[[1]]
cf_task_control <- cf_lkd_ctrl$enumerate_cf_tasks(tmle_task)[[1]]

Y <- tmle_task$get_tmle_node("Y_2")

# EY <- self$observed_likelihood$get_likelihood(tmle_task, self$outcome_node, fold_number)
# EY1 <- self$observed_likelihood$get_likelihood(cf_task_treatment, self$outcome_node, fold_number)
# EY0 <- self$observed_likelihood$get_likelihood(cf_task_control, self$outcome_node, fold_number)

# ZW todo: store all possible inputs and outputs as df, and use left_join to call for outputs
# R and L(t!=0) nodes: treated
# Z nodes: control
# the function bellow from Likelihood class might be used to get all wanted interventions
# one group of interventions for treated, one group of interventions for controls; each group contains all combos

# self$observed_likelihood is just a likelihood; it has get_possible_counterfactuals method

# all not A, not t=0 nodes
temp_node_names <- names(tmle_task$npsem)
loc_A <- grep("A", temp_node_names)
if_not_0 <- sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][2] != 0)
nodes_to_combo <- temp_node_names[ -c(loc_A, which(!if_not_0)) ]
temp_combos <- initial_likelihood$get_possible_counterfactuals(nodes = nodes_to_combo)
combos_treat <- combos_control <- data.frame(matrix(0, nrow(temp_combos), length(loc_A)), 
                                             temp_combos)  # first tau cols are treat/control
names(combos_treat)[1:length(loc_A)] <- names(combos_control)[1:length(loc_A)] <- temp_node_names[loc_A]
# ZW todo: stochastic interventions
combos_treat[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) treatment[[k]]$value %>% as.character %>% as.numeric)
combos_control[1:length(loc_A)] <- lapply(1:length(loc_A), function(k) control[[k]]$value %>% as.character %>% as.numeric)
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


test <- CF_Likelihood$new(initial_likelihood, temp_list[[1]][[2]])
task1 <- test$enumerate_cf_tasks(tmle_task)[[1]]
task1$data
initial_likelihood$get_likelihoods(task1, nodes = "A_1", fold_number = "full")
initial_likelihood$factor_list$L_1$type

# list of CF_Likelihood objects
cf_list <- list()
cf_list[[1]] <- lapply(temp_list[[1]], function(s) {
  CF_Likelihood$new(initial_likelihood, s)
})
cf_list[[2]] <- lapply(temp_list[[2]], function(s) {
  CF_Likelihood$new(initial_likelihood, s)
})
# ZW todo: save combo fittings
# private$.cf_likelihood_combo_list <- cf_list

loc_Z <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] == "Z"))
loc_RLY <- which(sapply(temp_node_names, function(s) strsplit(s, "_")[[1]][1] %in% c("R", "L", "Y") & strsplit(s, "_")[[1]][2] != 0))

fold_number <- "full"
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
    temp_return <- data.frame(initial_likelihood$get_likelihoods(task_1, temp_node_names[loc_RLY], fold_number), 
                              initial_likelihood$get_likelihoods(task_2, temp_node_names[loc_Z], fold_number)
    )
    temp_return <- temp_return[names(current_combo)]
    # decide p or 1-p
    temp_return <- lapply(1:length(current_combo), function(k) if (current_combo[k] == 1) temp_return[k] else 1 - temp_return[k] ) %>% data.frame
    apply(temp_return, 1, prod)
  }
})
s <- 129
current_combo <- combos_treat[s,]
current_combo <- current_combo[-grep("A", names(current_combo))]

task_1 <- cf_list[[1]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=1 for RLY nodes
task_2 <- cf_list[[2]][[s]]$enumerate_cf_tasks(tmle_task)[[1]]  # plug in A=0 for Z nodes
temp_return <- data.frame(initial_likelihood$get_likelihoods(task_1, temp_node_names[loc_RLY], fold_number), 
                          initial_likelihood$get_likelihoods(task_2, temp_node_names[loc_Z], fold_number)
)
temp_return

initial_likelihood$get_likelihoods(task_1, temp_node_names[loc_RLY], fold_number)
initial_likelihood$factor_list[["L_1"]]$get_likelihood(task_1)


psi <- pmap_dbl(list_prods, sum) %>% mean

