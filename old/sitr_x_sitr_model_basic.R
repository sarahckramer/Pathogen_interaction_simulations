##################################################################################################################
# R code to run pomp model 
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/sitr_x_sitr_basic.cpp
# The basic model does not incorporate seasonality in the transmission function and does not incorporate 
# under-reporting
#
# Created by: Sarah Pirikahu 
# Creation date: 19 December 
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
library(vars)
library(future) # allows for parallel processing

# read in the c code for the pomp model 
mod_code <- readLines('/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/old/sitr_x_sitr_basic.cpp')

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

# create dataframe of input parameters which will be used to simulate data from our pomp model

# specify ranges of all input variables - **Ask Mattheiu for thoughts on 
# suitable ranges 
# Ri1 <- c(1,2,3)
# Ri2 <- c(1,2,3)
# gamma1 <- c(0.5,1.0,1.5,2.0,2.5)
# gamma2 <- c(0.5,1.0,1.5,2.0,2.5)
# delta1 <- c(0.5,1.0,1.5,2.0,2.5)
# d2 <- c(0.5,1.0,1.5,2.0,2.5)
# theta_lambda1 <- c(0,0.25,0.5,0.75,1)
# theta_lambda2 <- c(0,0.25,0.5,0.75,1)
# rho1 <- c(0,0.05,0.1,0.25,0.5,0.75,1)
# rho2 <- c(0,0.05,0.1,0.25,0.5,0.75,1)
# theta_rho1 <- 1
# theta_rho2 <- 1
# beta_sd1 <- c(0,0.25,0.5,0.75,1)
# beta_sd2 <- c(0,0.25,0.5,0.75,1)
# N <- 10000 # change in future
# I10 <- c(0,0.001,0.005,0.01)
# I20 <- c(0,0.001,0.005,0.01)
# R10 <- c(0,0.1,0.2,0.3)
# R20 <- c(0,0.1,0.2,0.3)
# R120 <- c(0,0.1,0.2,0.3)  

# create a grid with all possible combinations of the above values
#true_params <- expand.grid(Ri1, Ri2, gamma1, gamma2, delta1, d2, theta_lambda1, theta_lambda2, rho1, rho2,
#            theta_rho1, theta_rho2, beta_sd1, beta_sd2, N, I10, I20, R10, R20, R120)
## this amount of parameter combos will not be feasible this alone is 220 billion combos. 
## So are we going to keep the season specific parameters? Remove them? keep some parameters
## fixed whilst others are allowed to vary? 

# create dataframe of single set of parameter inputs
true_params <- data.frame(Ri1=2, Ri2=5,
                           gamma1=7/5, gamma2=7/10,
                           delta1=7/5, d2=1.0,
                           theta_lambda1=1, theta_lambda2=1, 
                           rho1=0.15, rho2=0.5,
                           theta_rho1=1, theta_rho2=1, 
                           beta_sd1=0, beta_sd2=0, 
                           N=10000,
                           I10=0.002, I20=0.002,
                           R10=0, R20=0.02, R120=0)

#---- Create list to save the parameter sets and results of our different methods ---# 

results <- vector(mode = "list", length = dim(true_params)[1])
for (i in 1:dim(true_params)[1]){
  results[[i]] <- vector(mode = "list", length = 6)
  results[[i]][[1]] <- true_params[1,] 
  names(results[[i]]) <- c("true_param", "data", "cor", "granger", "transfer_entropy", "CCM")
}

#---- create pomp object ---# 
po <- pomp(data = data.frame(time = seq(from = 0, to = 52, by = 1), H1_obs = NA, H2_obs = NA),
           times = "time",
           t0 = 0,
           accumvars = c('H1_tot', 'H2_tot', 'H1', 'H2'),
           statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                          'X_SI', 'X_II', 'X_TI', 'X_RI', 
                          'X_ST', 'X_IT', 'X_TT', 'X_RT',
                          'X_SR', 'X_IR', 'X_TR', 'X_RR', 
                          'H1_tot', 'H2_tot', 
                          'H1', 'H2'),
           paramnames = names(true_params),
           params = true_params,
           globals = components_l[['globs']],
           dmeasure = components_l[['dmeas']],
           rmeasure = components_l[['rmeas']],
           skeleton = vectorfield(components_l[['skel']]),
           rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
           rinit = components_l[['rinit']]
)


# simulating multiple seasons and pulling them together to make a single timeseries
s1 <- simulate(po, times=1:26)
s2 <- simulate(po, times=1:26) 

# NOTE: H1_obs and H2_obs are the number of positive tests to each virus
s1_states <- as(s1, "data.frame") 
s2_states <- as(s2, "data.frame")

# combining the consecutive series 
d1 <- rbind(s1_states, s2_states)
dim(d1) 
names(d1)
# create a column for week number 
d1$week <- 1:dim(d1)[1]
head(d1)

# normalise case data
d1$H1_obs_NORM <- (d1$H1_obs - mean(d1$H1_obs))/sd(d1$H1_obs)
d1$H2_obs_NORM <- (d1$H2_obs - mean(d1$H2_obs))/sd(d1$H2_obs)

results[[i]]$data <- d1  

#--- plotting the data---#
# putting data into correct format to plot 
# original simulated data
d1_plot <- d1 %>%  dplyr::select(week,H1_obs, H2_obs) %>%
  pivot_longer(H1_obs:H2_obs, names_to = 'Vir', values_to = 'Inc')
d1_plot$Inc_percent <- (d1_plot$Inc/10000) * 100
head(d1_plot)
d1_plot$Vir <- as.factor(d1_plot$Vir)
levels(d1_plot$Vir) <- c("virus 1", "virus 2")
# plot out the data
ggplot(aes(x=week,y=Inc_percent,colour=Vir),data=d1_plot) + geom_line() + 
  theme_classic() + theme(legend.position="bottom", legend.title=element_blank()) + labs(y="Incidence (%)")

# Normalised data
d1_plot <- d1 %>%  dplyr::select(week,H1_obs_NORM, H2_obs_NORM) %>%
  pivot_longer(H1_obs_NORM:H2_obs_NORM, names_to = 'Vir', values_to = 'Inc')
d1_plot$Inc_percent <- (d1_plot$Inc/10000) * 100
head(d1_plot)
d1_plot$Vir <- as.factor(d1_plot$Vir)
levels(d1_plot$Vir) <- c("virus 1", "virus 2")
# plot out the data
ggplot(aes(x=week,y=Inc_percent,colour=Vir),data=d1_plot) + geom_line() + 
  theme_classic() + theme(legend.position="bottom", legend.title=element_blank()) + labs(y="Incidence (%)")

# plotting H_obs v H (i.e. total positive tests vs total number of infections)
d1_plot2 <- d1 %>%  dplyr::select(week,H1_obs, H2_obs, H1, H2) %>% 
  pivot_longer(H1_obs:H2, names_to = 'Vir', values_to = 'Inc')
head(d1_plot2)
d1_plot2$Inc_percent <- (d1_plot2$Inc/10000) * 100
head(d1_plot2)
d1_plot2$Vir <- as.factor(d1_plot2$Vir)
# plot 
ggplot(aes(x=week, y=Inc_percent), data=d1_plot2) + geom_line() + facet_wrap(.~Vir) +
  theme_bw()

# remove datasets no longer going to use
rm(s1,s2,s1_states,s2_states,d1_plot2,d1_plot,components_l,po,
   components_nm,mod_code,i,nm)

##########################################################
## Start testing each method for estimating interaction ##
##########################################################

# create dataset with just the observed cases
d_var <- d1[,c("H1_obs", "H1_obs_NORM", "H2_obs", "H2_obs_NORM")]

# Automatically determine the best lag to use based on AIC
lags <- lapply(d_var, VARselect) 
# pull out the lag with best AIC. Lower AIC = better 
# regardless of whether raw of normalised data used the lag choosen is the same
lag_h1 <- as.numeric(lags$H1_obs$selection[1])
lag_h2 <- as.numeric(lags$H2_obs$selection[2])

#---- Correlation coefficents --------# 
cor_raw <- cor.test(d_var$H1_obs, d_var$H2_obs); cor_raw
cor_NORM <- cor.test(d_var$H1_obs_NORM, d_var$H2_obs_NORM); cor_NORM

temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
results[[i]]$cor <- temp_res

# plot of interaction 
ggplot(aes(x=H1_obs,y=H2_obs),data=d1) + geom_point() + stat_cor()

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
str(d1)
names(d1)

# determine embedding dimension 
rho_E <- EmbedDimension(dataFrame = d1, columns = "H2", target = "H2",
                        lib = "1 156", pred = "1 156", showPlot = TRUE)
rho_E <- EmbedDimension(dataFrame = d1, columns = "H1", target = "H1",
                        lib = "1 30", pred = "1 30", showPlot = TRUE)
# E = 2 best embedding dimension for H2 and E = 4 based on H1
E <- 4



# test for non-linearity (this runs SMap under the hood)
rho_theta_e3 = PredictNonlinear(dataFrame = d1, columns = "H2",
                                target = "H2", lib = "1 156", pred = "1 156", E = E)
rho_theta_e3 = PredictNonlinear(dataFrame = d1, columns = "H1",
                                target = "H1", lib = "1 156", pred = "1 156", E = E)
# both obviously non-linear

# set up to do ccm for all pairs 
vars = colnames(d_var)
var_pairs = combn(vars, 2) # Combinations of vars, 2 at a time
libSize = paste(NROW(d_var) - E, NROW(d_var) - E, 10, collapse = " ")
ccm_matrix = array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,
                                                                            vars))
# do the ccm for all pairs of variables 
for (i in 1:ncol(var_pairs)) {
  ccm_out = CCM(dataFrame = d1, columns = var_pairs[1, i], 
                target = var_pairs[2,i], libSizes = libSize, Tp = 0,
                E = E, sample = 100)
  outVars = names(ccm_out)
  var_out = unlist(strsplit(outVars[2], ":"))
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 2]
  var_out = unlist(strsplit(outVars[3], ":"))
  ccm_matrix[var_out[2], var_out[1]] = ccm_out[1, 3]
}
ccm_matrix

# creating lagged cross-correlation for comparison 
corr_matrix <- array(NA, dim = c(length(vars), length(vars)), dimnames = list(vars,
                                                                              vars))
for (ccm_from in vars) {
  for (ccm_to in vars[vars != ccm_from]) {
    ccf_out <- ccf(d1[, ccm_from], d1[, ccm_to], type = "correlation",
                   lag.max = 6, plot = FALSE)$acf
    corr_matrix[ccm_from, ccm_to] <- max(abs(ccf_out))
  }
}
corr_matrix

# look at convergence in the cross-map predictability (compare rho as a function of L, library size) 
thrips_xmap_maxT <- CCM(dataFrame = d1, E = E, Tp = 0, columns = "H1",
                        target = "H2", libSizes = "10 20 30 40", sample = 100, showPlot = TRUE)
abline(h = corr_matrix["H1", "H2"], col = "black", lty = 2)


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



