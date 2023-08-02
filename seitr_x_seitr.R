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

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size) 

# determine which number job each original jobid, from the array, corresponds to
jobid <- (jobid - 1) %% no_jobs + 1; print(jobid) 

# Get unique identifiers for each repetition 
unique_ids <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)

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

# initialize time of surges from start of season (1 July)
# by drawing from a normal distribution
n_surge <- 5
mu_Imloss <- 188
sd_Imloss <- 35
t_si <- rnorm(n=n_surge, mean=mu_Imloss,sd=sd_Imloss)
# correcting t_si to give t based on day of the year rather than 
# day from start of season (note: July 1 is the 182 day)
t_si <- round(seq(182, 365*n_surge, by=365) + t_si)

# remove all t_si which are less than 2years - allowing this amount
# of time for the system to reach an equilibrium 
#t_si <- t_si[-which(t_si <= 730)]
# up dating the number of surges to input into parameter vector
#n_surge <- length(t_si)

# initialize the rate of loss of immunity corresponding to each of the 
# surge times 
delta_i <- runif(n=length(t_si), min = 0.1, max=0.12)

# create dataframe of single set of parameter inputs
# v1 = influenza; v2 = RSV
# true_params <- data.frame(Ri1=1.3, Ri2=4.5,
#                           sigma1=1, sigma2=1/5,
#                           gamma1=1/4, gamma2=1/10,
#                           delta1=1/7, delta2=1/7,
#                           w1=1/365, w2=1/190,
#                           rho1 = 0.2, rho2 = 0.2,
#                           theta_lambda1=1, theta_lambda2=1,
#                           A=1/10, phi=100,
#                           beta_sd1=0.1, beta_sd2=0.1,
#                           N=3000000,
#                           E01=0.001, E02=0.001,
#                           R01=0.4, R02=0.1, R12=0.1)

true_params <- data.frame(Ri1=1.3, Ri2=1,
                          sigma1=1, sigma2=1/5,
                          gamma1=1/4, gamma2=1/30,
                          delta1=1/7, delta2=1/7,
                          mu = 0.0001, nu=0.00003, 
                          w1=1/365, w2=1/190,
                          rho1 = 0.2, rho2 = 0.2,
                          theta_lambda1=1, theta_lambda2=1, 
                          A1=1/10, phi1=100,
                          A2=1/10, phi2=100,
                          beta_sd1=0, beta_sd2=0, 
                          N=3000000,
                          E01=0.001, E02=0.001,
                          R01=0.4, R02=0.2, R12=0.1,
                          n_surge = n_surge, t_si=t(t_si), delta_i=t(delta_i))

# replacing . in names of true params with _
names(true_params) <- gsub(x = names(true_params), pattern = "\\.", replacement = "_") 


#---- Create list to save the parameter sets and results of our different methods ---# 

results <- vector(mode = "list", length = dim(true_params)[1])
for (i in 1:dim(true_params)[1]){
  results[[i]] <- vector(mode = "list", length = 6)
  results[[i]][[1]] <- true_params[1,] 
  names(results[[i]]) <- c("true_param", "data", "cor", "transfer_entropy", "CCM","granger")
}

#---- create pomp object ---# 
po <- pomp(data = data.frame(time = seq(from = 0, to = 1825, by = 1), v1_obs = NA, v2_obs = NA),
           times = "time",
           t0 = 0,
           obsnames = c('v1_obs', 'v2_obs'),
           accumvars = c('v1_T', 'v2_T'),
           statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                          'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                          'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                          'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                          'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                          'v1_T', 'v2_T', 'death', 'fIS0','fIS1','fIS2', 'fIE', 'fII', 'fIT' , 'fIR','fSI', 'fEI', 'fTI', 'fRI'),
           paramnames = names(true_params),
           params = true_params,
           globals = components_l[['globs']],
           dmeasure = components_l[['dmeas']],
           rmeasure = components_l[['rmeas']],
           rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.3),
           skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
           rinit = components_l[['rinit']]
)

# set seed:
set.seed(2908)

# ----simulating data----#
s1 <- simulate(po, times=1:1825, format="data.frame")
d1 <- trajectory(po, times=1:1825, format = "data.frame") %>% dplyr::select(-'.id') %>% 
  mutate(v1_obs = rbinom(n=length(v1_T),size=round(v1_T), prob=true_params$rho1),  
         v2_obs = rbinom(n=length(v2_T),size=round(v2_T), prob=true_params$rho2))

# make time into years
d1$time <- d1$time/365
s1$time <- s1$time/365
# remove first 2 years where simulation isn't yet at equilibrium 
d1 <- d1 %>% filter(time > 2)
s1 <- s1 %>% filter(time > 2)

# save results
#results[[i]]$data <- s1  

# ---- Plotting simulated data ----#

# observation model output 
ggplot(aes(x=time, y=v1_obs),data=d1) + geom_line() + geom_line(aes(x=time, y=v2_obs), colour="blue") + 
  ggtitle("deterministic")
ggplot(aes(x=time, y=v1_obs),data=s1) + geom_line() + geom_line(aes(x=time, y=v2_obs), colour="blue") + 
  ggtitle("stochastic")


# process model output - each compartment over time 
d1_long <- gather(d1, compartment, cases, X_SS:v2_T, factor_key=T) # make data long
ggplot(aes(x=time,y=cases), data=d1_long) + geom_line() + facet_wrap(.~compartment, scales="free") + 
  ggtitle("deterministic")
s1_long <- gather(s1, compartment, cases, X_SS:v2_T, factor_key=T)
ggplot(aes(x=time,y=cases), data=s1_long) + geom_line() + facet_wrap(.~compartment, scales="free") + 
  ggtitle("stochastic")

# over plots of stochastic and deterministic models - can help identify problems
d1_long$type <- 'deterministic'
s1_long$type <- 'stochastic'
# combine the two datasets
s1_long$.id <- NULL
long_comb <- rbind(d1_long, s1_long)
# plot 
ggplot(aes(x=time, y=cases, colour=type),data=long_comb) + geom_line() +
  facet_wrap(.~compartment, scales="free") 


# remove datasets no longer going to use
rm(s1,s1_states,components_l,po,components_nm,mod_code,i,nm, true_params)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################
i=1

# create dataset with just the observed cases
d_var <- d1[,c("v1_obs", "v2_obs")]
# specify total number of weeks of data we have 
N <- d_var[,1] %>% length()
  
# Automatically determine the best lag doing several models with lags
# 1-5 then choose the best lag number  based on AIC
lags <- lapply(d_var, VARselect, lag.max=5) 
# pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
# regardless of whether raw of normalised data used the lag chosen is the same
lag_v1 <- as.numeric(lags$v1_obs$selection[3])
lag_v2 <- as.numeric(lags$v2_obs$selection[3])

#---- Correlation coefficents --------# 
cor_raw <- cor.test(d_var$v1_obs, d_var$v2_obs); cor_raw

temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
results[[i]]$cor <- temp_res

rm(temp_res, cor_raw)

#----- Transfer entropy analysis ------# 

# Interpreting transfer entropy (note: TE \in [0,1]):
# If test significant suggests T_{X->Y} > 0 and the uncertainty about 
# Y is reduced by the addition of X, that is X causes Y.

# Output: provides not the transfer entropy and bias corrected effective transfer entropy  
# Transfer entropy estimates are biased by small sample sizes. For large sample sizes TE and ETE 
# will be approximately the same. For a single season the sample size is quite small so we want to 
# go with ETE... see Behrendt et al. 2019 for more details
shannon_te <- transfer_entropy(d_var$v1_obs, d_var$v2_obs, lx = min(lag_v1, lag_v2), ly=min(lag_v1, lag_v2))
temp_res <- data.frame(coef(shannon_te))

# creating the 95% CIs about ETE - note that in the code for transfer_entropy to calculate the 
# se they simply look at the sd of the bootstrap samples NOT SE=sd/sqrt(n)
temp_res$lower95 <- temp_res$ete - 1.96*temp_res$se
temp_res$upper95 <- temp_res$ete + 1.96*temp_res$se

# add to the results list
results[[i]]$transfer_entropy <- temp_res

rm(temp_res, shannon_te)

#---- Granger causality analysis  ----# 
# separated this method out as a bit more code required
source("granger_analysis.R")

#------- Convergent Cross mapping analysis -------# 
# separated this method out as a bit more code required
source("CCM.R")


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

