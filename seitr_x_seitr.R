##################################################################################################################
# R code to run pomp model and each of the methods
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/seitr_x_seitr.cpp
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
library(ggpubr) # stat_cor
library(vars)
library(RTransferEntropy) 
library(future) # allows for parallel processing
library(mgcv)
library(latex2exp)
library(gridExtra)

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
# jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 
# no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
# sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size) 
# 
# # determine which number job each original jobid, from the array, corresponds to
# jobid <- (jobid - 1) %% no_jobs + 1; print(jobid) 
# 
# # Get unique identifiers for each repetition 
# unique_ids <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)

#--- reading in CSnippets ---# 
# read in the C code for the pomp model 
mod_code <- readLines('/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/seitr_x_seitr.cpp')

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
# initialize time of surges (based on week) from start of season (1 July)
# by drawing from a normal distribution 
n_surge <- 6
mu_Imloss <- 36 # early Oct
sd_Imloss <- 4

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

theta_lambda1 <- c(0,1,4)
theta_lambda2 <- c(0,1,4)
delta_1 <- 1
delta_2 <- 1

# function to create list of true parameter inputs and simulated data 
# function takes a vector of the interaction parameters 
sim_data <- function(theta_lambda1, theta_lambda2, delta_1, delta_2){
  set.seed(2908)
  # can add expand.grid function here to get all combinations of the 
  # interaction parameters
  

  # setting parameters to weekly rates - params listed as daily in 
  # spreadsheet list of model parameters.xlsx
  # note also v1 = influenza; v2 = RSV
  true_params <- data.frame(Ri1=1.3, Ri2=3.5,
                          sigma1=7, sigma2=7/5,
                          gamma1=7/5, gamma2=7/10,
                          delta1=1, delta2=1,
                          mu = 0.0002, nu=0.0007, 
                          w1=1/52, w2=1/28,
                          rho1 = 0.004, rho2 = 0.003,
                          theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                          A1=0.2, phi1=26,
                          A2=0.4, phi2=25,
                          beta_sd1=0, beta_sd2=0, 
                          N=3700000,
                          E01=0.001, E02=0.0007,
                          R01=0.4, R02=0.2, R12=0.001,
                          n_surge = n_surge, t_si=t(t_si), delta_i=t(delta_i))

  # replacing . in names of true params with _
  names(true_params) <- gsub(x = names(true_params), pattern = "\\.", replacement = "_") 

#---- Create list to save the parameter sets and results of our different methods ---# 

  results <- vector(mode = "list", length = dim(true_params)[1])
  for (j in 1:dim(true_params)[1]){
    results[[j]] <- vector(mode = "list", length = 6)
    results[[j]][[1]] <- true_params[1,] 
    names(results[[j]]) <- c("true_param", "data", "cor", "transfer_entropy", "CCM","granger")
  


#---- create pomp object ---# 
      po <- pomp(data = data.frame(time = seq(from = 0, to = 364, by = 1), v1_obs = NA, v2_obs = NA),
           times = "time",
           t0 = 0,
           obsnames = c('v1_obs', 'v2_obs'),
           accumvars = c('v1_T', 'v2_T'),
           statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                          'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                          'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                          'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                          'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                          'v1_T', 'v2_T', 
                          'w', 'delta', 'lambda_1', 'lambda_2', 'gamma_1',
                          'beta_1', 's_1'),
           paramnames = names(true_params),
           params = true_params[j,],
           globals = components_l[['globs']],
           dmeasure = components_l[['dmeas']],
           rmeasure = components_l[['rmeas']],
           rprocess = euler(step.fun = components_l[['rsim']], delta.t = 1),
           skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
           rinit = components_l[['rinit']]
      )

# ----simulating data----#
    s1 <- simulate(po, times=1:364, format="data.frame")
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
    results[[j]]$data <- s1  
      }
  return(results)
}

# generate all true parameter sets and simulate data 
results <- mapply(sim_data, theta_lambda1, theta_lambda2)

# ---- Plotting simulated data ----#
# changing the surge times to dates
t_si_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(t_si)

# plots of single simulation 
# ggplot(aes(x=time, y=v1_obs),data=d1) + geom_line() + geom_line(aes(x=time, y=v2_obs), colour="blue") + 
#   ggtitle("deterministic")
# ggplot(aes(x=time_date, y=v1_obs),data=s1) + geom_line() + geom_line(aes(x=time_date, y=v2_obs), colour="blue") + 
#   ggtitle("stochastic") + labs(y="observed cases") + 
#   scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y") + 
#   theme(axis.text.x=element_text(angle=60, hjust=1)) +  geom_vline(xintercept = t_si_date, linetype="dotted")

# creating multiple plots at once
plot_list <- list() 
for(i in 1:3){
    data <- results[[i]]$data %>% filter(time > 104)
    plot_list[[i]] <- ggplot(aes(x=time_date, y=v1_obs),data=data) + geom_line() + geom_line(aes(x=time_date, y=v2_obs), colour="blue") + 
    ggtitle(paste("theta_lambda1 and theta_lambda2 =", results[[i]]$true_param$theta_lambda1, 
                  "AND delta_1 = delta_2 =", results[[i]]$true_param$delta1)) + labs(y="observed cases") + 
    scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y") + ylim(0,800) +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +  geom_vline(xintercept = t_si_date, linetype="dotted")
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

# loop over each simulated dataset to determine the number of lags for each
for(i in 1:3){
  df <- results[[i]]$data %>% dplyr::select(v1_obs, v2_obs)
  lags[[i]] <- lapply(df, VARselect, lag.max=5) # lag of approx 1 month
  rm(df)
  # pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1[[i]] <- as.numeric(lags[[i]]$v1_obs$selection[3])
  lag_v2[[i]] <- as.numeric(lags[[i]]$v2_obs$selection[3])
  rm(lags)
}

#---- Correlation coefficents --------# 

# function to estimate correlation 
corr_func <- function(v1_obs, v2_obs){
  # calculated correlation coefficent 
  cor_raw <- cor.test(v1_obs, v2_obs)
  # pull together results in data frame 
  temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
  names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
  return(temp_res)
}

# apply correlation function to all simulated datasets and save results 
for(i in 1:3){
  results[[i]]$cor <- corr_func(v1_obs = results[[i]]$data$v1_obs, v2_obs = results[[i]]$data$v2_obs)
}

#----- Transfer entropy analysis ------# 

# create function to give transfer entropy results
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
for(i in 1:3){
  results[[i]]$transfer_entropy <- te_func(v1_obs = results[[i]]$data$v1_obs, v2_obs = results[[i]]$data$v2_obs, 
                                           lag_v1 = lag_v1[[i]], lag_v2 = lag_v2[[i]]) 
}

#---- Granger causality analysis  ----# 
# separated this method out as a bit more code required
source("granger_analysis.R")

# apply the granger analysis to each simulated data set and save the results
for(i in 1:3){
  data <- results[[i]]$data %>% dplyr::select(v1_obs, v2_obs)
  results[[i]]$granger <- granger_func(data = data, lag_v1 = lag_v1[[i]], lag_v2 = lag_v2[[i]])
}

#------- Convergent Cross mapping analysis -------# 
# separated this method out as a bit more code required
source("CCM.R")

# apply the CCM approach to each simulated data set and save the results
for(i in 1:3){
  data <- results[[i]]$data %>% dplyr::select(time, v1_obs, v2_obs)
  results[[i]]$CCM <- ccm_func(data = data)
}


#----- Likelihood approach -----# 



#---- Wavelets analysis  ----# 

my.wc <- analyze.coherency(d_var, my.pair = c("H1_obs","H2_obs"),
                           loess.span = 0,
                           dt = 1, dj = 1/100,
                           make.pval = TRUE, n.sim = 10)

wc.image(my.wc, n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05, ## default values
         legend.params = list(lab = "cross-wavelet power levels"),
         timelab = "")

wc.image(my.wc, which.image = "wc", color.key = "interval", n.levels = 250,
         siglvl.contour = 0.1, siglvl.arrow = 0.05,
         legend.params = list(lab = "wavelet coherence levels"),
         timelab = "")

wc.phasediff.image(my.wc, which.contour = "wc", use.sAngle = TRUE,
                   n.levels = 250, siglvl = 0.1,
                   legend.params = list(lab = "phase difference levels",
                                        lab.line = 3),
                   timelab = "")


#--- GAM approach ---# 

mvn_mod <- gam(formula = list(v1_obs ~ s(time), v2_obs ~ s(time)),
               family = mvn(d = 2),
               data = s1)
summary(mvn_mod)

corr_mat <- mvn_mod$family$data$R %>%
  crossprod() %>%
  solve() %>%
  cov2cor(); corr_mat
