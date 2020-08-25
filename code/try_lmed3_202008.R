library(data.table)
library(sl3)
library(tlverse)
library(R6)  # R6Class
library(dplyr)  # %>% 



home <- getwd()  # specify it if needed
source(file.path(home, "code", "basic_functions-202008.R"))
source(file.path(home, "code", "generate_Zheng_data-202008.R"))
source("./code/try_lmed3_to_import_202008.R")

sample_size <- 100
data_sim <- generate_Zheng_data(B = sample_size, tau = 2)
data_wide <- data_sim %>% data.frame

# all the imported class
{
  
}

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

