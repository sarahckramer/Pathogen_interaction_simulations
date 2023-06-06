##################################################################################################################
# R code to run pomp model 
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
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'skel', 'rsim')
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
                          delta1=0.7, delta2=0.6,
                          rho1 = 0.5, rho2 = 0.2,
                          theta_lambda1=1, theta_lambda2=1, 
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
  names(results[[i]]) <- c("true_param", "data", "cor", "granger", "transfer_entropy", "CCM")
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

#--- plotting the data---#
# putting data into correct format to plot 
# original simulated data
d1_plot <- d1 %>%  dplyr::select(time,v1_obs, v2_obs) %>%
  pivot_longer(v1_obs:v2_obs, names_to = 'Vir', values_to = 'Inc')
d1_plot$Inc_percent <- (d1_plot$Inc/1000000) * 100
head(d1_plot)
d1_plot$Vir <- as.factor(d1_plot$Vir)
levels(d1_plot$Vir) <- c("virus 1", "virus 2")
# plot out the data
ggplot(aes(x=time,y=Inc_percent,colour=Vir),data=d1_plot) + geom_line() + 
  theme_classic() + theme(legend.position="bottom", legend.title=element_blank()) + labs(y="Incidence (%)")

# Normalised data
d1_plot <- d1 %>%  dplyr::select(time,v1_obs_NORM, v2_obs_NORM) %>%
  pivot_longer(v1_obs_NORM:v2_obs_NORM, names_to = 'Vir', values_to = 'Inc')
d1_plot$Inc_percent <- (d1_plot$Inc/1000000) * 100
head(d1_plot)
d1_plot$Vir <- as.factor(d1_plot$Vir)
levels(d1_plot$Vir) <- c("virus 1", "virus 2")
# plot out the data
ggplot(aes(x=time,y=Inc_percent,colour=Vir),data=d1_plot) + geom_line() + 
  theme_classic() + theme(legend.position="bottom", legend.title=element_blank()) + labs(y="Incidence (%)")

# plotting H_obs v H (i.e. total positive tests vs total number of infections)
d1_plot2 <- d1 %>%  dplyr::select(time,v1_obs, v2_obs, v1_T, v2_T) %>% 
  pivot_longer(v1_obs:v2_T, names_to = 'Vir', values_to = 'Inc')
head(d1_plot2)
d1_plot2$Inc_percent <- (d1_plot2$Inc/1000000) * 100
head(d1_plot2)
d1_plot2$Vir <- as.factor(d1_plot2$Vir)
# plot 
ggplot(aes(x=time, y=Inc_percent), data=d1_plot2) + geom_line() + facet_wrap(.~Vir) +
  theme_bw()

# remove datasets no longer going to use
rm(s1,s2,s1_states,s2_states,d1_plot2,d1_plot,components_l,po,
   components_nm,mod_code,i,nm)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

# create dataset with just the observed cases
d_var <- d1[,c("v1_obs", "v1_obs_NORM", "v2_obs", "v2_obs_NORM")]

# Automatically determine the best lag doing several models with lags
# 1-10 then choose the best lag number  based on AIC
lags <- lapply(d_var, VARselect, lag.max=5) 
# pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
# regardless of whether raw of normalised data used the lag choosen is the same
lag_h1 <- as.numeric(lags$v1_obs$selection[3])
lag_h2 <- as.numeric(lags$v2_obs$selection[3])

#---- Correlation coefficents --------# 
cor_raw <- cor.test(d_var$v1_obs, d_var$v2_obs); cor_raw

temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
results[[i]]$cor <- temp_res

# plot of interaction 
ggplot(aes(x=v1_obs,y=v2_obs),data=d1) + geom_point() + stat_cor()

rm(temp_res)

#---- Granger causality analysis  ----# 
# separated this method out as a bit more code required
source("granger_analysis.R")

#----- Transfer entropy analysis ------# 

# determining the transfer entropy (note: TE \in [0,1] 
# if significant suggests T_{X->Y} > 0 and the uncertainty about 
# Y is reduced by the addition of X, that is X causes Y.
shannon_te <- transfer_entropy(d_var$H1_obs, d_var$H2_obs)
shannon_te

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

#------- Convergent Cross mapping analysis -------# 
# separated this method out as a bit more code required
source("CCM.R")

#----- Bayesian multivariate autoregression -----# 

library(coda)
library(rjags)


# add year to d1 
d1$year <- 2020
d1$month <- lubridate::month(as.Date(paste(d1$year, d1$week_year,1, sep="-"), "%Y-%U-%u"))

# Neighbourhood model
neighbour_model <- model(data){
  # initalising 
  Nobs <- 10000 # total number of people in the population 
  Nyears <- 1 # total number of years we have data for  
  # calculating expected counts from data
  
  # total number of people positive for each virus 
  Nt_state <- d1 %>% group_by(year) %>% summarise() 
  # total number of people in the population in each year 
  Nt <- 10000 # cause we aren't accounting for demographics here
  
  # Nmv = total number of testes for virus v in month m
  # we are treating this as H1 and H2 calculate each of these for each month
  Nmv <- d1 %>% group_by(month) %>% summarise(H1_m = sum(H1), H2_m = sum(H2))
  
  # pmv  = probability of testing positive for virus v in month m 
  # the total number of those that test positive is H1_obs and H2_obs
  obs_tot_m <- d1 %>% group_by(month) %>% summarise(H1_obs_m = sum(H1_obs), H2_obs_m = sum(H2_obs))
  pmv <- obs_tot_m/Nmv
  
  Expected <- Nt*pmv # expected count for month m, year t, virus v
  
  
  
  for (i in 1:Nobs){
    for (j in 1:Nyears){
      Y[i,j] ~ dpois(RR[i,j]*Expected[i,j])
      log(RR[i,j]) <- alpha[virus[i]] +
        phi.month[Month.Virus.phi[i],j]
    }
  }
  ##Virus overall means
  for (i in 1:Nvirus) {
    alpha[i] ~ dnorm(0,0.0001)
  }
  phi.month[1:(Nvirus*Nmonths),1] ~ dmnorm(mean.month,Omega.month)
  #temporal tends
  for (j in 2:Nyears){
    phi.month[1:(Nvirus*Nmonths),j] ~
      dmnorm(smooth*phi.month[1:(Nvirus*Nmonths),j-
                                1],Omega.month)#Y[1:Nobs,j-1]
  }
  for (i in 1:(Nvirus*Nmonths)){
    smooth[i]<-smooth1[virus_order[i]]
  }
  for (i in 1:Nvirus){
    smooth1[i]~dunif(0,1) #temporal smoothing can be different for each virus
  }
  #Work out precisions matrix Omega
  for (i in 1:Nmonths) {
    for (j in 1:Nmonths) {
      for (k in 1:Nvirus) {
        for (l in 1:Nvirus) {
          Omega.month[(i-1)*Nvirus+k,(j-1)*Nvirus+l]<-
            omega.month.part1[i,j]*Lambda[k,l]
        }
      }
    }
  }
  omega.month.part1<-D.month-(lambda.month*W.month)
  lambda.month ~ dunif(0,1) #seasonal smoothing parameter
  ###This is the modified chol decomposition
  LAMBDA[1,1]<-1/LAMBDA1[1,1]
  for (i in 2:Nvirus){
    LAMBDA[i,i]<-1/LAMBDA1[i,i]
    for (j in 1:(i-1)){
      LAMBDA[j,i]<-0
      LAMBDA[i,j]<-0
    }
  }
  for (i in 1:Nvirus){ #gamma priors for standard deviations 
    LAMBDA1[i,i]~dgamma(1,1)
  }
  GAMMA[1,1]<-1
  for (i in 2:Nvirus){
    GAMMA[i,i]<-1
    for (j in 1:(i-1)){
      GAMMA[j,i]~dnorm(0,1)T(-1,1) #truncated normal priors for correlations parameters
      GAMMA[i,j]<-0
    }
  }
  Lambda.inv<-LAMBDA%*%GAMMA%*%t(GAMMA)%*%LAMBDA ##This is the covariance matrix
  Lambda<-inverse(Lambda.inv) ##This is the precision matrix
}






#Autoregressive model
model {
  for (i in 1:Nobs){
    for (j in 1:Nyears){
      Y[i,j] ~ dpois(RR[i,j]*Expected[i,j])
      log(RR[i,j]) <- alpha[virus[i]] +
        phi.month[Month.Virus.phi[i],j]
    }
  }
  ##Virus overall means
  for (i in 1:Nvirus) {
    alpha[i] ~ dnorm(0,0.0001)
  }
  phi.month[1:(Nvirus*Nmonths),1] ~ dmnorm(mean.month,Omega.month)
  for (j in 2:Nyears){
    phi.month[1:(Nvirus*Nmonths),j] ~
      dmnorm(smooth*phi.month[1:(Nvirus*Nmonths),j-
                                1],Omega.month)#Y[1:Nobs,j-1]
  }
  for (i in 1:(Nvirus*Nmonths))
  {
    smooth[i]<-smooth1[virus_order[i]]
  }
  for (i in 1:Nvirus)
  {
    smooth1[i]~dunif(0,1) 
    #temporal smoothing can be different for each virus
  }
  #Work out precision matrix Omega
  for (i in 1:Nmonths) {
    for (j in 1:Nmonths) {
      for (k in 1:Nvirus) {
        for (l in 1:Nvirus) {
          Omega.month[(i-1)*Nvirus+k,(j-1)*Nvirus+l]<-
            omega.month.part1[i,j]*Lambda[k,l]
        }
      }
    }
  }
  omega.month.part1<-D.month-(lambda.month*W.month)
  omega.month.part.W[1,1]<-0
  for (i in 2:Nmonths){
    omega.month.part.W[i,i]<-0
    for (j in 1:(i-1)){
      omega.month.part.W[j,i]<-pow(rho,STEPS[j,i]) 
      #STEP is a matrix that determines the number of steps between months
      omega.month.part.W[i,j]<-pow(rho,STEPS[i,j])
    }
  }
  for (i in 1:Nmonths)
  {
    for (j in 1:Nmonths)
    {
      W.month[i,j]<-
        omega.month.part.W[i,j]#/sum(omega.month.part.W[i,])
    }
  }
  D.month[1,1]<-sum(omega.month.part.W[1,])
  for (i in 2:Nmonths)
  {
    D.month[i,i]<-sum(omega.month.part.W[i,])
    for (j in 1:(i-1))
    {
      D.month[i,j]<-0
      D.month[j,i]<-0
    }
  }
  rho ~ dunif(0.001,0.9) #autoregressive parameter
  lambda.month ~ dunif(0,1) #seasonal smoothing
  #This is the modified chol decomposition
  LAMBDA[1,1]<-1/LAMBDA1[1,1]
  for (i in 2:Nvirus){
    LAMBDA[i,i]<-1/LAMBDA1[i,i]
    for (j in 1:(i-1)){
      LAMBDA[j,i]<-0
      LAMBDA[i,j]<-0
    }
  }
  for (i in 1:Nvirus) #can use Log normal priors for standard
    deviations
  {
    LAMBDA1[i,i]~dgamma(1,1)
  }
  GAMMA[1,1]<-1
  for (i in 2:Nvirus){
    GAMMA[i,i]<-1
    for (j in 1:(i-1)){
      GAMMA[j,i]~dnorm(0,1)T(-1,1)
      GAMMA[i,j]<-0
    }
  }
  Lambda.inv<-LAMBDA%*%GAMMA%*%t(GAMMA)%*%LAMBDA ##This is the covariance matrix
  Lambda<-inverse(Lambda.inv) ##This is the precision matrix
}



