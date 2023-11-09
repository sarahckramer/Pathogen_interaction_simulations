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
library(gridExtra)

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 
likelihood <- as.logical(Sys.getenv("LIKELIHOOD")); print(likelihood) # will be TRUE or FALSE
# how many different starting params are we going to run for numerical optimizer for each job 
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)  

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

# total number of weeks of data we are going to want 
tot_weeks <- 625 # 12 years 
#tot_weeks <- 1145 # 22 years 
#tot_weeks <- 5304 # 102 years 

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

# create a function to specify multiple sets of parameter inputs
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
sim_data <- function(tot_weeks,theta_lambda1,theta_lambda2,delta_1,delta_2,n_surge,components_l=components_l){
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
                   beta_sd1=0, beta_sd2=0, 
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
                    delta_1=delta_1, delta_2=delta_2, n_surge = n_surge, components_l = components_l)


# ---- Plotting simulated data ----#
# changing the surge times to dates
t_si_date <- lubridate::ymd("2012-July-01") + lubridate::weeks(t_si)

# order starting parameters by theta
all_param_comb <- all_param_comb[order(all_param_comb$theta_lambda1),]

# creating multiple plots at once
temp <- vector(mode = "list", length = 15)
plot_list <- vector(mode = "list", length = 15)
attack_plots <- vector(mode = "list", length = 15)
res_all <- NULL
for(i in 1:15){
  theta_lambda1 <- all_param_comb[i,]$theta_lambda1
  theta_lambda2 <- all_param_comb[i,]$theta_lambda2
  delta_1 <- all_param_comb[i,]$delta_1
  delta_2 <- all_param_comb[i,]$delta_2
  temp[[i]] <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2,
                        delta_1=delta_1, delta_2=delta_2, n_surge=n_surge, components_l=components_l)
  data <- temp[[i]]$data

  legend_colors <- c("v1_obs" = "black", "v2_obs" = "blue")
  plot_list[[i]] <- ggplot(aes(x=time_date, y=v1_obs, colour="v1_obs"),data=data) + geom_line() + geom_line(aes(x=time_date, y=v2_obs,colour="v2_obs")) +
    ggtitle(paste("theta_lambda1 and theta_lambda2 =", temp[[i]]$true_param["theta_lambda1"],
                  "AND delta_1 = delta_2 =", temp[[i]]$true_param["delta1"])) + labs(y="observed cases") +
    scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")  +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +  geom_vline(xintercept = t_si_date, linetype="dotted") +
    scale_colour_manual(values=legend_colors) + labs(colour="")


 # also estimate attack rates by year for each plot...... NOT WORKING
  data$season <- c(rep(1:tot_seasons, each=52),tot_seasons+1)
  #data$season <- rep(1:tot_seasons, each=52)
  seasonal_incidence <- data %>% group_by(season) %>%
                          summarise(obs_v1 = sum(v1_obs), obs_v2 = sum(v2_obs),
                                    tot_v1 = sum(v1_T), tot_v2 = sum(v2_T))
  # trying to calculate the attack rate based on observed data
  obs_v1_attack <- seasonal_incidence$obs_v1/3700000 * 100
  obs_v2_attack <- seasonal_incidence$obs_v2/3700000 * 100

  range_obs_v1_att <- range(obs_v1_attack[-length(obs_v1_attack)]) # 0.14 - 0.19
  range_obs_v2_att <- range(obs_v2_attack[-length(obs_v2_attack)]) # 0.18 - 0.20

  # trying to calculate the attack rate based on true number of cases from the model
  tot_v1_attack <- seasonal_incidence$tot_v1/3700000 * 100
  tot_v2_attack <- seasonal_incidence$tot_v2/3700000 * 100
  tot_v1_attack <- tot_v1_attack[-c(length(tot_v1_attack))] 
  tot_v2_attack <- tot_v2_attack[-c(length(tot_v2_attack))] 
  
  plot_dat <- data.frame(cbind(tot_v1_attack = tot_v1_attack, tot_v2_attack = tot_v2_attack))
  attack_plots[[i]] <- ggplot(aes(x=tot_v2_attack,y=tot_v1_attack), data=plot_dat) + geom_point() +
    ggtitle(paste("theta_lambda1 and theta_lambda2 =", temp[[i]]$true_param["theta_lambda1"],
                  "AND delta_1 = delta_2 =", temp[[i]]$true_param["delta1"])) 
  
  range_tot_v1_att <- range(tot_v1_attack[-length(tot_v1_attack)]) # 71 - 95
  range_tot_v2_att <- range(tot_v2_attack[-length(tot_v2_attack)]) # 90 - 97

  res <- cbind(all_param_comb[i,],
               range_obs_v1_att[1],range_obs_v1_att[2],
               range_obs_v2_att[1],range_obs_v2_att[2],
               range_tot_v1_att[1],range_tot_v1_att[2],
               range_tot_v2_att[1],range_tot_v2_att[2])
  res_all <- rbind(res_all, res)

}

# plot simulated timeseries data
grid.arrange(plot_list[[1]],plot_list[[2]],plot_list[[3]],ncol=1)
grid.arrange(plot_list[[4]],plot_list[[5]],plot_list[[6]],ncol=1)
grid.arrange(plot_list[[7]],plot_list[[8]],plot_list[[9]],ncol=1)
grid.arrange(plot_list[[10]],plot_list[[11]],plot_list[[12]],ncol=1)
grid.arrange(plot_list[[13]],plot_list[[14]],plot_list[[15]],ncol=1)

# plot scatter plots of seasonal attack rates
grid.arrange(attack_plots[[1]],attack_plots[[2]],attack_plots[[3]],ncol=1)
grid.arrange(attack_plots[[4]],attack_plots[[5]],attack_plots[[6]],ncol=1)
grid.arrange(attack_plots[[7]],attack_plots[[8]],attack_plots[[9]],ncol=1)
grid.arrange(attack_plots[[10]],attack_plots[[11]],attack_plots[[12]],ncol=1)
grid.arrange(attack_plots[[13]],attack_plots[[14]],attack_plots[[15]],ncol=1)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

# Since the likelihood approach requires a significantly larger amount of compute power will only 
# perform this method if it is specifically asked for at the cmd line

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
  # separated this method out as a bit more code required
  source("./methods/granger_analysis.R")
  # apply the granger analysis to each simulated data set and save the results
  data <- results$data %>% dplyr::select(v1_obs, v2_obs)
  results$granger <- granger_func(data = data, lag_v1 = lag_v1, lag_v2 = lag_v2)
  
  #------- Convergent Cross mapping analysis -------# 
  # separated this method out as a bit more code required
  source("./methods/CCM.R")
  
  # apply the CCM approach to each simulated data set 
  data <- results$data %>% dplyr::select(time, v1_obs, v2_obs)
  results$CCM <- ccm_func(data = data)
  
  # save out the results
  save(results, file=sprintf('results_%s.RData',jobid))
} else {
  
  #----- Likelihood approach -----# 
  
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




