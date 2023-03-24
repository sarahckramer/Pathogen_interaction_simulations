##################################################################################################################
# R code to run pomp model 
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/sitr_x_sitr.cpp
#
# Created by: Sarah Pirikahu 
# Creation date: 19 December 
##################################################################################################################

# load libraries
library(tidyverse)
library(pomp)
library(janitor)
library(ggfortify)
library(ggpubr) # stat_cor
library(tseries) # adf test - stationary 
library(forecast) # ndiff - check if differences in lags needed for stationary 
library(lmtest) # granger causality test
library(vars) # VAR model - that which is underlying granger causality test
library(RTransferEntropy) # transfer entropy approach
library(future) # allows for parallel processing
library(WaveletComp) # wavelet analysis
library(rEDM) # convergent cross mapping analysis 

# read in the c code
# home load
#mod_code <- readLines('/Users/spirikahu/Documents/Max planck/Project 1/Analysis/Simulation/sitr_x_sitr_basic.cpp')
# work load
mod_code <- readLines('/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/sitr_x_sitr_basic.cpp')

# pull out the various components of the C code ready to feed into pomp
# components looking for
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'skel', 'rsim')
# initialise list
components_l <- vector(mode = 'list', length = length(components_nm))
names(components_l) <- components_nm

for (nm in components_nm) {
  components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
    str_flatten(collapse = '\n')
  components_l[[nm]] <- Csnippet(text = components_l[[nm]])
}

# create pomp object 
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
           paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                          'gamma1', 'gamma2', # 1 / average infectious periods
                          # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                          'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                          'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                          'rho1', 'rho2', # probs. infection leads to ILI consultation
                          'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                          'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                          'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                          'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                          'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                          'N', # population size
                          'I10', 'I20', # props. infectious at outbreak start
                          'R10', 'R20', 'R120'), # props. recovered at outbreak start
           params = c(Ri1 = 2, Ri2 = 5,
                      #gamma1 = 7 / 5, gamma2 = 7 / 10, # or 4 for flu?
                      gamma1 = 7/5, gamma2 = 7/10, # or 4 for flu?
                      # delta = 7 / 5,
                      #delta1 = 14, d2 = 0.5,
                      delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                      theta_lambda1 = 1, theta_lambda2 = 1,
                      rho1 = 0.15, rho2 = 0.5,
                      alpha = 0, phi = 0,
                      theta_rho1 = 1, theta_rho2 = 1,
                      eta_temp1 = 0, eta_temp2 = 0,
                      eta_ah1 = 0, eta_ah2 = 0,
                      beta_sd1 = 0, beta_sd2 = 0,
                      N = 10000,
                      #I10 = 0.002, I20 = 0.002,
                      I10 = 0.002, I20 = 0.002,
                      R10 = 0, R20 = 0.02, R120 = 0),
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

#--- plotting the data---#
# putting data into correct format to plot 
d1_plot <- d1 %>%  dplyr::select(week,H1_obs, H2_obs) %>%
  pivot_longer(H1_obs:H2_obs, names_to = 'Vir', values_to = 'Inc')
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

# plot of interaction 
ggplot(aes(x=H1_obs,y=H2_obs),data=d1) + geom_point() + stat_cor(method="spearman") # so much higher than I feel it should be... 
# really don't trust this approach at all due to the lack of independence between observations
# try out a linear model --- still not going to be good but lets see
lm1 <- lm(H1 ~ H2, data=d1)
summary(lm1)
autoplot(lm1)

#---- Granger causality analysis  ----# 
### NOTE RATHER THAN USING grangertest in simulation use VAR and pull the p-values and coefficents from the 
### linear regresssion summary


# how many differences needed 
ndiffs(d1$H1, alpha=0.05, test=c("kpss")) # 0 differences required - suggests the series is potentially stationary 
ndiffs(d1$H2, alpha=0.05, test=c("kpss")) # 0 

## Determine lag automatically based on AIC
# create dataset with just the observed cases
d_var <- d1[,c("H1_obs", "H2_obs")]
# determine the lag
lags <- lapply(d_var, VARselect) 
# pull out the lag with best AIC. Lower AIC = better 
lag_h1 <- as.numeric(lags$H1_obs$selection[1])
lag_h2 <- as.numeric(lags$H2_obs$selection[2])

# checking stationarity of my time series
# Test hypotheses
# H0: there is a unit root - i.e. the time series is not stationary 
# H1: time series is stationary 
adf.test(d1$H1_obs,k=lag_h1) # k = lag number 
adf.test(d1$H2_obs,k=lag_h2) 
# both time series are stationary 

# running test; order = lag numbers
# Test hypotheses
# H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
# H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 
grangertest(H1_obs ~ H2_obs, order = lag_h1, data=d1) 

# check out reverse causality 
grangertest(H2_obs ~ H1_obs, order = lag_h2, data=d1) 

# extract p-values from grangertest
gt1 <- grangertest(H1_obs ~ H2_obs, order = lag_h1,data=d1); gt1
str(gt1)
gt1$`Pr(>F)`[2] # this gives the p-value

gt2 <- grangertest(H2_obs ~ H1_obs, order = lag_h2,data=d1); gt2
str(gt2)
gt2$`Pr(>F)`[2] # this gives the p-value

# determining the causal effect size log(RSS_X/RSS_X&Y)
# run the VAR (vector autoregressive) model (i.e. which has both X and Y)
var1 <- VAR(y=d_var, p=min(lag_h1,lag_h2))
summary(var1) 

# run AR model to run univariate analysis (Note: VAR is only for multivariate
# and if you try to run a univariate analysis like this with VAR it will tell 
# you to use ar)
ar_H1=ar(d_var$H1_obs,order=var1$p,aic=F,method="ols")
ar_H2=ar(d_var$H2_obs,order=var1$p,aic=F,method="ols")

log_21_inter=log(sum((ar_H1$resid)^2,na.rm=T)/sum((var1$varresult$H1$residuals)^2,na.rm=T))
log_12_inter=log(sum((ar_H2$resid)^2,na.rm=T)/sum((var1$varresult$H2$residuals)^2,na.rm=T))

log_21_inter # effect of H2 on H1 
log_12_inter # effect of H1 on H2


#----- Transfer entropy analysis ------# 

set.seed(12345)
shannon_te <- transfer_entropy(d_var$H1, d_var$H2, lx=3, ly=3)
shannon_te <- transfer_entropy(d_var$H1, d_var$H2, lx=2, ly=2)
shannon_te

shannon_te <- transfer_entropy(d_var$H1, d_var$H2, lx=5, ly=5)
shannon_te

# plotting the lagged data of X and Y against each other 
H1_minus1 <- d_var$H1[-1]
H2_minus1 <- d_var$H2[-1]
H1 <- d_var$H1[-length(d_var$H1)]
H2 <- d_var$H2[-length(d_var$H2)]
# pull into dataframe 
df_plot <- data.frame(cbind(H1,H2,H1_minus1, H2_minus1))

# plot 
ggplot(aes(x=H1, y=H2_minus1),data=df_plot) + geom_point()
ggplot(aes(x=H1, y=H2),data=df_plot) + geom_point()
ggplot(aes(x=H2, y=H1_minus1),data=df_plot) + geom_point()


#---- Wavelets analysis  ----# 

my.wc <- analyze.coherency(d_var, my.pair = c("H1","H2"),
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



