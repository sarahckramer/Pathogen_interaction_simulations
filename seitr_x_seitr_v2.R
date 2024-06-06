##################################################################################################################
# R code to run pomp model and each of the methods
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/seitr_x_seitr.cpp
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

# load libraries
library(tidyverse)
library(testthat)
library(pomp)
library(janitor)
library(ggfortify)
library(future) # allows for parallel processing
library(foreach)
library(doParallel)

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 

#--- reading in CSnippets ---# 
# read in the C code for the pomp model 
mod_code <- readLines('seitr_x_seitr.cpp')

# pull out the various components of the C code ready to feed into pomp
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'rsim', 'skel', 'toest', 'fromest')
# initialise list
components_l <- vector(mode = 'list', length = length(components_nm))
names(components_l) <- components_nm
# create list with code components
for (nm in components_nm) {
  components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
    str_flatten(collapse = '\n')
  components_l[[nm]] <- Csnippet(text = components_l[[nm]])
}

#---- setting parameter input values ----# 

# set seed:
set.seed(1234)

# set noise parameters
beta_sd1 <- 0 
beta_sd2 <- 0

# total number of simulated datasets to create for each parameter input 
nsim <- 100

# total number of seasons
tot_weeks <- 625
tot_seasons <- round((tot_weeks/52) - 2)

# initialize time of surges (based on week) from start of season (1 July)
# by drawing from a normal distribution 
n_surge <- round(tot_weeks/52) - 1 # total number of surges
mu_Imloss <- 38 # average surge occuring in mid Oct
sd_Imloss <- 4 # allow variation of 4 weeks

t_si <- rnorm(n=n_surge, mean=mu_Imloss,sd=sd_Imloss)
# correcting t_si to give t based on week of the year rather than 
# week from start of season (note: July 1 is week 26)
t_si <- round(seq(26, 52*n_surge, by=52) + t_si)

# remove all t_si which are less than 2years - allowing this amount
# of time for the system to reach an equilibrium (2 yrs ~ 104 weeks)
t_si <- t_si[-which(t_si <= 104)]

# up dating the number of surges to input into parameter vector
n_surge <- length(t_si)

# initialize the rate of loss of immunity corresponding to each of the 
# surge times 
delta_i <- runif(n=length(t_si), min = 0.01*7, max=0.1*7)

# parameter inputs 
theta_lambda1 <- c(0,0.5,1,2,4)
theta_lambda2 <- c(0,0.5,1,2,4)
delta_1 <- c(1,1/4,1/12)
delta_2 <- c(1,1/4,1/12)

# Get all combinations of the interaction parameters
all_param_comb <- expand.grid(theta_lambda1, theta_lambda2, delta_1, delta_2)
names(all_param_comb) <- c("theta_lambda1", "theta_lambda2", "delta_1", "delta_2")

# remove parameter vectors 
rm(theta_lambda1, theta_lambda2, delta_1, delta_2)

# look at just symmetric interactions 
all_param_comb <- all_param_comb %>% filter(theta_lambda1 == theta_lambda2)

# function to create list of true parameter inputs and simulated data 
# function takes a vector of the interaction parameters 
sim_data <- function(tot_weeks,theta_lambda1,theta_lambda2,delta_1,delta_2,beta_sd1,beta_sd2,n_surge,components_l=components_l,nsim){
  set.seed(2908)
  
  # setting parameters to weekly rates - params listed as daily in 
  # spreadsheet list of model parameters.xlsx
  # note also v1 = influenza; v2 = RSV
  true_params <- c(Ri1=1.1, Ri2=1.7,
                   sigma1=7, sigma2=7/5,
                   gamma1=7/5, gamma2=7/10,
                   delta1=delta_1, delta2=delta_2,
                   w1=1/52, w2=1/52,
                   mu = 0.0002, nu = 0.0002,
                   rho1 = 0.002, rho2 = 0.002,
                   theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                   A1=0.2, phi1=26,
                   A2=0.2, phi2=26,
                   beta_sd1=beta_sd1, beta_sd2=beta_sd2, 
                   N=3700000,
                   E01=0.001, E02=0.001,
                   R01=0.4, R02=0.25, R12=0.001, nsurges=n_surge,
                   t_si_=t(t_si), delta_i_=t(delta_i))
  
  #---- Create list to save the parameter sets and results of our different methods ---# 
  
  results <- vector(mode = "list", length = 7)
  results[[1]] <- true_params 
  names(results) <- c("true_param", "data", "cor", "gam_cor", "transfer_entropy", "CCM","granger")
  
  #---- create pomp object ---# 
  po <- pomp(data = data.frame(time = seq(from = 0, to = tot_weeks, by = 1), v1_obs = NA, v2_obs = NA),
             times = "time",
             t0 = 0,
             obsnames = c('v1_obs', 'v2_obs'),
             accumvars = c('v1_T', 'v2_T'),
             statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                            'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                            'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                            'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                            'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                            'v1_T', 'v2_T'),
             paramnames = names(true_params),
             params = true_params,
             partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
             globals = components_l[['globs']],
             dmeasure = components_l[['dmeas']],
             rmeasure = components_l[['rmeas']],
             rprocess = euler(step.fun = components_l[['rsim']], delta.t = 1),
             skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
             rinit = components_l[['rinit']]
  )
  
  # ----simulating data----#
  s1 <- simulate(po, nsim=nsim, times=1:tot_weeks, format="data.frame")
  
  # remove first 2 years where simulation isn't yet at equilibrium 
  s1 <- s1 %>% filter(time > 104) 
  
  # make time into dates based off week number
  s1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(s1$time)
  
  # save results
  results$data <- s1  
  return(results)
}

# generate single true parameter sets and simulate data 
theta_lambda1 <- all_param_comb[jobid,]$theta_lambda1
theta_lambda2 <- all_param_comb[jobid,]$theta_lambda2
delta_1 <- all_param_comb[jobid,]$delta_1
delta_2 <- all_param_comb[jobid,]$delta_2
results <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                    delta_1=delta_1, delta_2=delta_2, beta_sd1=beta_sd1, beta_sd2=beta_sd2,
                    n_surge = n_surge, components_l = components_l, nsim = nsim)
# if outbreak never takes off or dies out then remove it 
temp <- results$data %>% group_by(.id) %>% mutate(sum_v1 = sum(v1_obs), sum_v2=sum(v2_obs))
results$data <- temp %>% filter(sum_v1 !=0 & sum_v2 != 0)

#---- Correlation coefficents --------# 
# function to estimate correlation 
# inputs: v1_obs = time series for v1_obs
#         v2_obs = time series for v2_obs
# corr_func <- function(v1_obs, v2_obs){
#   # calculated correlation coefficent
#   cor_raw <- cor.test(v1_obs, v2_obs)
#   # pull together results in data frame
#   temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
#   names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
#   return(temp_res)
# }
# 
# # apply correlation function to all simulated datasets and save results
# results$cor <- results$data %>% group_by(.id) %>% do((corr_func(.$v1_obs,.$v2_obs)))


#--- GAM approach ---# 
# source("./methods/gam_cor.R")

# setting up parallelism for the foreach loop
# registerDoParallel(cl <- makeCluster(10))
# # apply the GAM correlation approach to each simulated data set and save the results
# res_gam_cor <- foreach(i=1:nsim, .packages=c("tidyverse","mgcv","vars","boot")) %dopar%{
#   # if the dataset was removed because the outbreak died out then skip it
#   if(dim(results$data %>% filter(.id==i))[1]!=0){
#     results$data %>% filter(.id==i) %>% gam_cor(.)
#   }
# }
# results$gam_cor <- do.call(rbind, res_gam_cor)

#----- Transfer entropy analysis ------# 
# source("./methods/transfer_entropy_jidt.R")
# 
# # lag = 1
# results$transfer_entropy <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="1"))
# # lag = 2
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="2"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)
# # lag = 4
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="4"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)
# # lag = 6
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="6"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)

#---- Granger causality analysis  ----# 
# source("./methods/granger_analysis.R")
# 
# # apply granger analysis to each simulated data set and save the results
# results$granger <- results$data %>% group_by(.id) %>% do(granger_func(.))

#------- Convergent Cross mapping analysis -------# 
source("./methods/CCM.R")

results$CCM <- results$data %>% group_by(.id) %>% do(ccm_func(.))

# save out results
save(results, file=sprintf('results_%s_%s_%s_%s_%s.RData',jobid, theta_lambda1,theta_lambda2,delta_1,delta_2))  
