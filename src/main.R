##################################################################################################################
# R code to run pomp model and each of the methods
#
# This code runs the following methods:
#  - Pearson's correlation coefficient
#  - GAM estimated correlation coefficient
#  - Transfer entropy
#  - Granger causality analysis 
#  - convergent cross mapping
# 
# Created by: Sarah Pirikahu 
# Creation date: 28 Feb 2024
##################################################################################################################

# # Code to run through all jobids locally (uncomment and copy-paste in console):
# for (jobid_use in 1:16) {
#   source('src/main.R')
#   detachAllPackages <- function() {
#     basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
#     package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
#     package.list <- setdiff(package.list,basic.packages)
#     if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
#   }
#   detachAllPackages()
#   # source: https://stackoverflow.com/questions/7505547/detach-all-packages-while-working-in-r
# }

# Setup
tic_all <- Sys.time()

#---- load libraries ----#
library(foreach)
library(doParallel)
library(doSNOW)
library(doMC)

print(detectCores())

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size
run_local <- as.logical(Sys.getenv("RUNLOCAL")); print(run_local)

#---- run local or on cluster? ----#
if (is.na(run_local)) {
  run_local <- TRUE
  
  if (exists('jobid_use')) {
    jobid <- jobid_use
  } else {
    jobid <- 1
  }
  print(jobid)
  
}

#---- set parameters and generate synthetic data ----#
source('src/generate_data.R')

#---- set up list to store all results ----#
results <- vector(mode = 'list', length = 7)
names(results) <- c('true_param', 'data', 'cor', 'gam_cor', 'granger', 'transfer_entropy', 'CCM')

results$true_param <- true_params
results$data <- dat

##################################################################################################################

# Run various methods to determine causality

#---- Correlation coefficents ----#

corr_func <- function(data){
  
  # Function to calculate correlation
  # param data: Tibble containing time series of viruses 1 and 2
  # returns: Tibble of Pearson's r with confidence intervals and p-values
  
  cor_raw <- data %>% group_by(.id) %>%
    dplyr::select(.id, V1_obs:V2_obs) %>%
    group_map(~ cor.test(.$V1_obs, .$V2_obs))
  
  temp_res <- bind_cols(cor = lapply(cor_raw, getElement, 'estimate') %>%
                          unlist()) %>%
    bind_cols(CI_lower_95 = lapply(lapply(cor_raw, getElement, 'conf.int'), '[[', 1) %>%
                unlist()) %>%
    bind_cols(CI_upper_95 = lapply(lapply(cor_raw, getElement, 'conf.int'), '[[', 2) %>%
                unlist()) %>%
    bind_cols(p_value = lapply(cor_raw, getElement, 'p.value') %>%
                unlist()) %>%
    mutate(.id = 1:length(cor_raw), .before = cor)
  
  return(temp_res)
}

# Apply correlation function to all simulated datasets and save results:
if (run_local) {
  tic <- Sys.time()
  results$cor <- corr_func(dat)
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'secs'
  print(etime)
}

#---- GAM approach ----#
source('src/methods/gam_cor.R')

if (!run_local) {
  
  # setting up parallelism for the foreach loop
  # registerDoParallel(cl <- makeCluster(50))
  registerDoMC(50)
  
  # apply the GAM correlation approach to each simulated data set and save the results
  tic <- Sys.time()
  res_gam_cor <- foreach(i = 1:n_sim, .packages=c('tidyverse', 'mgcv', 'vars', 'boot')) %dopar% {
    
    dat %>% filter(.id == i) %>% gam_cor()
    
    # # if the dataset was removed because the outbreak died out then skip it
    # if(dim(results$data %>% filter(.id==i))[1]!=0){
    #   results$data %>% filter(.id==i) %>% gam_cor(.)
    # }
    
  }
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)
  
  # compile all results
  results$gam_cor <- do.call(rbind, res_gam_cor) %>%
    mutate(.id = 1:n_sim, .before = cor)
  
}

#---- Granger causality analysis  ----#
# apply granger analysis to each simulated data set and save the results
if (run_local) {
  source('src/methods/granger_analysis.R')
  
  tic <- Sys.time()
  results$granger <- dat %>% group_by(.id) %>% do(granger_func(.))
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)
}

#---- Transfer entropy analysis ----#
if (run_local) {
  
  source('src/methods/transfer_entropy_jidt.R')
  
  tic <- Sys.time()
  
  # lag = 1
  res_te_1 <- dat %>% group_by(.id) %>% do(te_jidt(., lag = '1'))
  # lag = 2
  res_te_2 <- dat %>% group_by(.id) %>% do(te_jidt(., lag = '2'))
  # lag = 4
  res_te_4 <- dat %>% group_by(.id) %>% do(te_jidt(., lag = '4'))
  # lag = 6
  res_te_6 <- dat %>% group_by(.id) %>% do(te_jidt(., lag = '6'))
  
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)
  
  # combine results and store
  results$transfer_entropy <- bind_rows(res_te_1, res_te_2, res_te_4, res_te_6)
  rm(res_te_1, res_te_2, res_te_4, res_te_6)
  
}

#---- Convergent Cross mapping analysis ----#
source('src/methods/CCM.R')

if (!run_local) {
  tic <- Sys.time()
  results$CCM <- dat %>% group_by(.id) %>% do(ccm_func(.))
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'hours'
  print(etime)
}

# save out results
write_rds(results, file=sprintf('results/results_%s_%s.rds', jobid, run_local))

#---- Clean up ----#
toc_all <- Sys.time()
etime <- toc_all - tic_all
units(etime) <- 'mins'
print(etime)

rm(list = ls())

print('Done!')
