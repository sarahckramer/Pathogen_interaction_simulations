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

#--- reading in CSnippets ---# 
# read in the C code for the pomp model 
mod_code <- readLines('seitr_x_seitr.cpp')

# pull out the various components of the C code ready to feed into pomp
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'rsim', 'skel')
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
tot_weeks <- 365 # 7 yrs

# initialize time of surges (based on week) from start of season (1 July)
# by drawing from a normal distribution 
n_surge <- round(tot_weeks/52) - 1 # total number of surges
mu_Imloss <- 36 # average surge occuring in early Oct
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
delta_i <- runif(n=length(t_si), min = 0.01, max=0.12)

# create a function to specify multiple sets of parameter inputs
theta_lambda1 <- c(0,1,2)
theta_lambda2 <- c(0,1,2)
delta_1 <- c(1,1/2,1/3)
delta_2 <- c(1,1/2,1/3)

# Get all combinations of the interaction parameters
all_param_comb <- expand.grid(theta_lambda1, theta_lambda2, delta_1, delta_2)
names(all_param_comb) <- c("theta_lambda1", "theta_lambda2", "delta_1", "delta_2")

# for now just keep symmetric interactions 
all_param_comb <- all_param_comb %>% filter(theta_lambda1 == theta_lambda2 & delta_1 == delta_2)
# remove parameter vectors 
rm(theta_lambda1, theta_lambda2, delta_1, delta_2)

# function to create list of true parameter inputs and simulated data 
# function takes a vector of the interaction parameters 
sim_data <- function(tot_weeks,theta_lambda1, theta_lambda2, delta_1, delta_2, components_l=components_l){
  set.seed(2908)
  
  # setting parameters to weekly rates - params listed as daily in 
  # spreadsheet list of model parameters.xlsx
  # note also v1 = influenza; v2 = RSV
  true_params <- data.frame(Ri1=1.3, Ri2=1.7,
                          sigma1=7, sigma2=7/5,
                          gamma1=7/5, gamma2=7/10,
                          delta1=delta_1, delta2=delta_2,
                          mu = 0.0002, nu=0.0007, 
                          w1=1/52, w2=1/28,
                          rho1 = 0.002, rho2 = 0.002,
                          theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                          A1=0.2, phi1=26,
                          A2=0.2, phi2=20,
                          beta_sd1=0, beta_sd2=0, 
                          N=3700000,
                          E01=0.0001, E02=0.0001,
                          R01=0.4, R02=0.2, R12=0.001,
                          n_surge = n_surge, t_si=t(t_si), delta_i=t(delta_i))

  # replacing . in names of true params with _
  names(true_params) <- gsub(x = names(true_params), pattern = "\\.", replacement = "_") 

#---- Create list to save the parameter sets and results of our different methods ---# 

    results <- vector(mode = "list", length = 8)
    results[[1]] <- true_params[1,] 
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
    # d1 <- trajectory(po, times=1:364, format = "data.frame") %>% dplyr::select(-'.id') %>% 
    #   mutate(v1_obs = rbinom(n=length(v1_T),size=round(v1_T), prob=true_params$rho1),  
    #          v2_obs = rbinom(n=length(v2_T),size=round(v2_T), prob=true_params$rho2))
    
    # remove first 2 years where simulation isn't yet at equilibrium 
    s1 <- s1 %>% filter(time > 104)
    
    # make time into dates based off week number
    s1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(s1$time)
    #d1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(d1$time)
    
  # save results
    results$data <- s1  
    return(results)
}
 
# generate all true parameter sets and simulate data 
theta_lambda1 <- all_param_comb[jobid,]$theta_lambda1
theta_lambda2 <- all_param_comb[jobid,]$theta_lambda2
delta_1 <- all_param_comb[jobid,]$delta_1
delta_2 <- all_param_comb[jobid,]$delta_2
results <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                      delta_1=delta_1, delta_2=delta_2, components_l = components_l)
  
# ---- Plotting simulated data ----#
# changing the surge times to dates
t_si_date <- lubridate::ymd("2012-July-01") + lubridate::weeks(t_si)
  
# creating multiple plots at once
temp <- vector(mode = "list", length = 3)
plot_list <- vector(mode = "list", length = 3)
for(i in 1:3){
  theta_lambda1 <- all_param_comb[i,]$theta_lambda1
  theta_lambda2 <- all_param_comb[i,]$theta_lambda2
  delta_1 <- all_param_comb[i,]$delta_1
  delta_2 <- all_param_comb[i,]$delta_2
  temp[[i]] <- sim_data(theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2,
                      delta_1=delta_1, delta_2=delta_2)
  data <- temp[[i]]$data
  plot_list[[i]] <- ggplot(aes(x=time_date, y=v1_obs),data=data) + geom_line() + geom_line(aes(x=time_date, y=v2_obs), colour="blue") +
    ggtitle(paste("theta_lambda1 and theta_lambda2 =", temp[[i]]$true_param$theta_lambda1,
                  "AND delta_1 = delta_2 =", temp[[i]]$true_param$delta1)) + labs(y="observed cases") +
    scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y") + ylim(0,400) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +  geom_vline(xintercept = t_si_date, linetype="dotted")

  # # also estimate attack rates by year for each plot...... NOT WORKING
  # data$season <- rep(1:5, each=52)
  # seasonal_incidence <- data %>% group_by(season) %>% summarise(tot_v1 = sum(v1_obs), tot_v2 = sum(v2_obs))
  # start_season <- data %>% group_by(season) %>% summarise(min(time_date))
  # start_season <- data %>% filter(time_date %in% start_season$`min(time_date)`)
  # v1_susceptible <- start_season$X_SS + start_season$X_SE + start_season$X_SI + start_season$X_ST + start_season$X_SR
  # v2_susceptible <- start_season$X_SS + start_season$X_ES + start_season$X_IS + start_season$X_TS + start_season$X_RS
  # seasonal_incidence$tot_v1/v1_susceptible*100
  
}
grid.arrange(grobs=plot_list,ncol=1)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

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
source("granger_analysis.R")
# apply the granger analysis to each simulated data set and save the results
data <- results$data %>% dplyr::select(v1_obs, v2_obs)
results$granger <- granger_func(data = data, lag_v1 = lag_v1, lag_v2 = lag_v2)
  
#------- Convergent Cross mapping analysis -------# 
# separated this method out as a bit more code required
source("CCM.R")
  
# apply the CCM approach to each simulated data set 
data <- results$data %>% dplyr::select(time, v1_obs, v2_obs)
results$CCM <- ccm_func(data = data)


#--- GAM approach ---# 
source("gam_cor.R")

# apply the GAM correlation approach to each simulated data set and save the results
results$gam_cor <- gam_cor(data=data)

# save out the results
save(results, file=sprintf('results_%s.RData',jobid))

#----- Likelihood approach -----# 

# Since the likelihood approach requires a significantly larger amount of compute power will only 
# perform this method if it is specifically asked for at the cmd line
likelihood <- as.logical(Sys.getenv("LIKELIHOOD")); print(likelihood) # will be TRUE or FALSE

if(likelihood==TRUE){
  # how many jobs are there in total we want to run? 
  no_jobs <- dim(all_param_comb)[1]
  # how many different starting params are we going to run for numerical optimizer run for each job (i.e. interaction parameter combos)
  sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size) # probably ~10 
  
  
  # run this bit separately if we want. Must
  true_params <- results$true_param
  
}




