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
#  - likelihood estimation (this method is going to be done separately due to the increased computational power required)
#
# GAM, Grangers, CCM and likelihood methods have been coded up seprately and sourced in due to their increased 
# complexity
# 
# Created by: Sarah Pirikahu 
# Creation date: 22 May 2023
##################################################################################################################

# load libraries
library(tidyverse)
library(testthat)
library(pomp)
library(janitor)
library(ggfortify)
library(vars)
library(RTransferEntropy) 
library(future) # allows for parallel processing


#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 
likelihood <- as.logical(Sys.getenv("LIKELIHOOD")); print(likelihood) # will be TRUE or FALSE
# how many different starting params are we going to run for numerical optimizer for each job -- likelihood approach
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)  
# number of weeks of data to simulate - note: first 2yrs of data discarded to allow mechanistic model to achieve equilibrium
tot_weeks <- as.integer(Sys.getenv("WEEKSSIM")); print(tot_weeks)  
# the amount of demographic stochasticity 
beta_sd1 <- as.integer(Sys.getenv("BETASD1")); print(beta_sd1)  
beta_sd2 <- as.integer(Sys.getenv("BETASD1")); print(beta_sd2)  
# the period for the surrogate generation in CCM 
Tperiod_v1 <- as.integer(Sys.getenv("TPERIODV1")); print(beta_sd1)  
Tperiod_v2 <- as.integer(Sys.getenv("TPERIODV2")); print(beta_sd2)  
# the amount of noise to allow into the CCM surrogates 
alpha_v1 <- as.integer(Sys.getenv("ALPHAV1")); print(alpha_v1)  
alpha_v2 <- as.integer(Sys.getenv("ALPHAV2")); print(alpha_v2)  

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
set.seed(2908)

# total number of seasons
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

# parameter inputs for simulation:
theta_lambda1 <- c(0,0.5,1,2,4)
theta_lambda2 <- c(0,0.5,1,2,4)
delta_1 <- c(1,1/4,1/24)
delta_2 <- c(1,1/4,1/24)

# Get all combinations of the interaction parameters
all_param_comb <- expand.grid(theta_lambda1, theta_lambda2, delta_1, delta_2)
names(all_param_comb) <- c("theta_lambda1", "theta_lambda2", "delta_1", "delta_2")

# for now just keep symmetric interactions 
all_param_comb <- all_param_comb %>% filter(theta_lambda1 == theta_lambda2 & delta_1 == delta_2)
# remove parameter vectors 
rm(theta_lambda1, theta_lambda2, delta_1, delta_2)

# function to create list of true parameter inputs and simulated data 
# function takes a vector of the interaction parameters 
sim_data <- function(tot_weeks,theta_lambda1,theta_lambda2,delta_1,delta_2,beta_sd1,beta_sd2,n_surge,components_l=components_l){
  set.seed(2908)
  
  # setting parameters to weekly rates - params listed as daily in 
  # spreadsheet list of model parameters.xlsx
  # note also v1 = influenza; v2 = RSV
  true_params <- c(Ri1=1.1, Ri2=1.7,
                   sigma1=7, sigma2=7/5,
                   gamma1=7/5, gamma2=7/10,
                   delta1=delta_1, delta2=delta_2,
                   w1=1/78, w2=1/52,
                   mu = 0.0002, nu = 0.0002,
                   rho1 = 0.002, rho2 = 0.002,
                   theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                   A1=0.01, phi1=26,
                   A2=0.2, phi2=20,
                   beta_sd1=beta_sd1, beta_sd2=beta_sd2, 
                   N=3700000,
                   E01=0.001, E02=0.001,
                   R01=0.4, R02=0.25, R12=0.001, nsurges=n_surge,
                   t_si_=t(t_si), delta_i_=t(delta_i))
  
  #---- Create list to save the parameter sets and results of our different methods ---# 
  
  results <- vector(mode = "list", length = 8)
  results[[1]] <- true_params 
  names(results) <- c("true_param", "data", "cor", "gam_cor", "transfer_entropy", "CCM","granger","likelihood")
  
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
  s1 <- simulate(po, times=1:tot_weeks, format="data.frame")
  # deterministic simulation 
  # d1 <- trajectory(po, times=1:tot_weeks, format = "data.frame") %>% dplyr::select(-'.id') %>%
  #   mutate(v1_obs = rbinom(n=length(v1_T),size=round(v1_T), prob=true_params["rho1"]),
  #          v2_obs = rbinom(n=length(v2_T),size=round(v2_T), prob=true_params["rho2"]))
  
  # remove first 2 years where simulation isn't yet at equilibrium 
  s1 <- s1 %>% filter(time > 104) 
  #d1 <- d1 %>% filter(time > 104)
  
  # make time into dates based off week number
  s1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(s1$time)
  #d1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(d1$time)
  
  # save results
  results$data <- s1  
  #results$data <- d1  
  return(results)
}

# generate single true parameter sets and simulate data 
theta_lambda1 <- all_param_comb[jobid,]$theta_lambda1
theta_lambda2 <- all_param_comb[jobid,]$theta_lambda2
delta_1 <- all_param_comb[jobid,]$delta_1
delta_2 <- all_param_comb[jobid,]$delta_2
results <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                    delta_1=delta_1, delta_2=delta_2, beta_sd1=beta_sd1, beta_sd2=beta_sd2,
                    n_surge = n_surge, components_l = components_l)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

# Since the likelihood approach requires a significantly larger amount of compute power will only 
# perform this method if it is specifically asked for at the cmd line generally on an ad hoc basis

# If likelihood is true we don't run the non-likelihood methods at all 
if(likelihood==FALSE){ 
  
  #------------ setup ---------------#
  # Automatically determine the best lag doing several models with lags
  # 1-5 (approximately 1 month) then choose the best lag number based on BIC
  
  # initialising lists to put results in 
  lags <- list()
  lag_v1 <- list()
  lag_v2 <- list()
  
  # determine the number of lags for each simulated dataset
  df <- results$data %>% dplyr::select(v1_obs, v2_obs)
  
  lags <- lapply(df, VARselect, lag.max=5) # lag of approx 1 month
  rm(df)
  # pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1 <- as.numeric(lags$v1_obs$selection[3])
  lag_v2 <- as.numeric(lags$v2_obs$selection[3])
  
  rm(lags)
  
  #---- Correlation coefficents --------# 
  # function to estimate correlation 
  # inputs: v1_obs = time series for v1_obs
  #         v2_obs = time series for v2_obs
  corr_func <- function(v1_obs, v2_obs){
    # calculated correlation coefficent 
    cor_raw <- cor.test(v1_obs, v2_obs)
    # pull together results in data frame 
    temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
    names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
    return(temp_res)
  }
  
  # apply correlation function to all simulated datasets and save results 
  results$cor <- corr_func(v1_obs = results$data$v1_obs, v2_obs = results$data$v2_obs)
  
  
  #--- GAM approach ---# 
  source("./methods/gam_cor.R")
  
  # apply the GAM correlation approach to each simulated data set and save the results
  data <- results$data %>% dplyr::select(time,v1_obs, v2_obs)
  results$gam_cor <- gam_cor(data=data)
  
  #----- Transfer entropy analysis ------# 
  # create function to give transfer entropy results
  # inputs: v1_obs = time series for v1_obs
  #         v2_obs = time series for v2_obs
  #         lag_v1 = total number lag to use with v1_obs time series
  #         lag_v2 = total number lag to use with v2_obs time series
  te_func <- function(v1_obs, v2_obs, lag_v1, lag_v2){
    # Interpreting transfer entropy (note: TE \in [0,1]):
    # If test significant suggests T_{X->Y} > 0 and the uncertainty about 
    # Y is reduced by the addition of X, that is X causes Y.
    
    # Output: provides not the transfer entropy and bias corrected effective transfer entropy  
    # Transfer entropy estimates are biased by small sample sizes. For large sample sizes TE and ETE 
    # will be approximately the same. For a single season the sample size is quite small so we want to 
    # go with ETE... see Behrendt et al. 2019 for more details
    shannon_te <- transfer_entropy(v1_obs, v2_obs, lx = min(lag_v1, lag_v2), ly=min(lag_v1, lag_v2))
    temp_res <- data.frame(coef(shannon_te))
    
    # creating the 95% CIs about ETE - note that in the code for transfer_entropy to calculate the 
    # se they simply look at the sd of the bootstrap samples NOT SE=sd/sqrt(n)
    temp_res$lower95 <- temp_res$ete - 1.96*temp_res$se
    temp_res$upper95 <- temp_res$ete + 1.96*temp_res$se
    return(temp_res)
  }
  
  # apply transfer entropy function to all simulated datasets and save results
  results$transfer_entropy <- te_func(v1_obs = results$data$v1_obs, v2_obs = results$data$v2_obs, 
                                      lag_v1 = lag_v1, lag_v2 = lag_v2) 
  
  #---- Granger causality analysis  ----# 
  source("./methods/granger_analysis.R")
  # apply the granger analysis to each simulated data set and save the results
  data <- results$data %>% dplyr::select(v1_obs, v2_obs)
  results$granger <- granger_func(data = data, lag_v1 = lag_v1, lag_v2 = lag_v2)
  
  #------- Convergent Cross mapping analysis -------# 
  source("./methods/CCM.R")
  
  # apply the CCM approach to each simulated data set 
  data <- results$data %>% dplyr::select(time, v1_obs, v2_obs)
  # run CCM 
  results$CCM <- ccm_func(data = data, Tperiod_v1=Tperiod_v1, Tperiod_v2=Tperiod_v2,
                          alpha_v1 = alpha_v1, alpha_v2 = alpha_v2, tot_weeks=tot_weeks)
  
  #------Finish up by saving out results -------#
  # save out the results: results_jobid_totweeks_noise
  save(results, file=sprintf('results_%s_%s_%s.RData',jobid, tot_weeks, beta_sd1*100))  
  
  # run results extraction to get a spreadsheet with just summary stats for each method
  #source("results_extraction.R")
  
} else {
  
  #----- Likelihood approach -----# 
  # Note this is very computationally intensive and may not always achieve convergence to the MLE
  # so was only run ad hoc - not for the entire simulation study
  
  # max time that we want to allow the optmizer to run for
  maxtime <- 12*60*60  # 12 hrs - want this to run for as long as possible so that I don't have to re-run to get the MLE
  
  # dataset to apply method to 
  data <- results$data %>% dplyr::select(time, v1_obs, v2_obs)
  # true parameters used to create the simulated data set
  true_params <- results$true_param
  
  # run likelihood estimation - these results will not be saved to the overall results
  # list here but instead to individual RDS files
  source("./methods/Likelihood_estimation.R")
  results$likelihood <- lik(data=data, true_params, components_l = components_l, sobol_size, jobid, no_jobs = no_jobs, maxtime)
  
  # save out the results
  # if I do it this way then I will have a list which has input params the simulated data and the 
  # output of the different input parameters for the likelihood estimation in their own results file
  # total output: number of jobs x sobol size results files  
  save(results, file=sprintf('results_%s.RData',jobid)) 
}




