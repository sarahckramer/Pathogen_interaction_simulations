##################################################################################################################
# R code to run pomp model and each of the methods
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/seitr_x_seitr.cpp
#
# Created by: Sarah Pirikahu 
# Creation date: 22 May 2023
##################################################################################################################

# set seed:
set.seed(2908)

# load libraries
library(tidyverse)
library(testthat)
library(pomp)
library(janitor)
library(ggfortify)
library(ggpubr) # stat_cor
library(vars)
library(RTransferEntropy) 
library(vars)
library(future) # allows for parallel processing

#---- set up cluster ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size 
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size) 

# determine which number job each original jobid, from the array, corresponds to
jobid <- (jobid - 1) %% no_jobs + 1; print(jobid) 

# Get unique identifiers for each repetition 
unique_ids <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)

#--- Setup pomp model and simulation ---# 
# read in the C code for the pomp model 
mod_code <- readLines('/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/seitr_x_seitr.cpp')

# pull out the various components of the C code ready to feed into pomp
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'rsim')
# initialise list
components_l <- vector(mode = 'list', length = length(components_nm))
names(components_l) <- components_nm
# create list with code components
for (nm in components_nm) {
  components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
    str_flatten(collapse = '\n')
  components_l[[nm]] <- Csnippet(text = components_l[[nm]])
}


# create dataframe of single set of parameter inputs
true_params <- data.frame(Ri1=2, Ri2=3,
                          sigma1=1, sigma2=1 ,
                          gamma1=7/5, gamma2=7/10,
                          delta1=0.6, delta2=0.6,
                          rho1 = 0.5, rho2 = 0.2,
                          theta_lambda1=2, theta_lambda2=2, 
                          A=0, phi=0,
                          beta_sd1=0, beta_sd2=0, 
                          N=1000000,
                          E01=0.01, E02=0.01,
                          R01=0.1, R02=0.1, R12=0.001)

#---- Create list to save the parameter sets and results of our different methods ---# 

results <- vector(mode = "list", length = dim(true_params)[1])
for (i in 1:dim(true_params)[1]){
  results[[i]] <- vector(mode = "list", length = 6)
  results[[i]][[1]] <- true_params[1,] 
  names(results[[i]]) <- c("true_param", "data", "cor", "transfer_entropy", "CCM","granger")
}

#---- create pomp object ---# 
po <- pomp(data = data.frame(time = seq(from = 0, to = 52, by = 1), v1_obs = NA, v2_obs = NA),
           times = "time",
           t0 = 0,
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
           #skeleton = vectorfield(components_l[['skel']]),
           rprocess = euler(step.fun = components_l[['rsim']], delta.t = 1),
           rinit = components_l[['rinit']]
)


# simulating multiple seasons and pulling them together to make a single timeseries
s1 <- simulate(po, times=1:52)

# NOTE: H1_obs and H2_obs are the number of positive tests to each virus
s1_states <- as(s1, "data.frame") 
d1 <- s1_states

# normalise case data
d1$v1_obs_NORM <- (d1$v1_obs - mean(d1$v1_obs))/sd(d1$v1_obs)
d1$v2_obs_NORM <- (d1$v2_obs - mean(d1$v2_obs))/sd(d1$v2_obs)

results[[i]]$data <- d1  


# remove datasets no longer going to use
rm(s1,s1_states,components_l,po,components_nm,mod_code,i,nm, true_params)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

# create dataset with just the observed cases
d_var <- d1[,c("v1_obs", "v1_obs_NORM", "v2_obs", "v2_obs_NORM")]
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

