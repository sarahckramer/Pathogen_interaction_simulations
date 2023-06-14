################################################################################
#                  Bayesian Multivariate Autoregressive Model
#
# Useful documentation for this method: Mair et al. 2019
# https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007492
# Partial code is provided in the appendix 
#
# Note: there are two different ways of deriving the precision matrix for the MCAR
# model 1. neighbourhood strucure and 2. autoregressive structure
# we implement both and will compare
#
# Created by: Sarah Pirikahu
# Creation date: 14 June 2023
################################################################################

# load packages
library(coda)
library(rjags)



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


