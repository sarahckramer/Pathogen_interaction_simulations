################################################################################
#                       Granger causality analysis        
#
# Useful documentation describing and giving examples:
# 
# To run Granger analysis we need the time series to be covariance and mean 
# stationary
#
# Created by: Sarah Pirikahu
# Creation date: 23 March 2023
################################################################################

# load packages
library(tseries) # adf test - stationary 
library(forecast) # ndiff - check if differences in lags needed for stationary 
library(lmtest) # granger causality test
library(vars) # VAR model - that which is underlying granger causality test
library(boot)

# checking stationary of time series

# ADF Test hypotheses
# H0: there is a unit root - i.e. the time series is not stationary 
# H1: time series is stationary 
adf_v1 <- adf.test(d_var$v1_obs,k=lag_v1); adf_v1 # k = lag number 
adf_v2 <- adf.test(d_var$v2_obs,k=lag_v2); adf_v2

# KPSS Test hypotheses
# H0: time series is stationary 
# H1: non staionary
kpss_v1 <- kpss.test(d_var$v1_obs); kpss_v1
kpss_v2 <- kpss.test(d_var$v2_obs); kpss_v2

### IF NOT STATIONARY .... difference and run tests again? 


# running Granger test and extracting p-values
# Test hypotheses
# H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
# H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 


# specifying the lag to be the minimum of the two series
p <- min(lag_v1,lag_v2)
# results the same regardless of performance on the raw or normalised data
gt1 <- grangertest(v1_obs ~ v2_obs, order = p, data=d_var); gt1
p_gt1 <- gt1$`Pr(>F)`[2] # p-value

gt2 <- grangertest(v2_obs ~ v1_obs, order = p, data=d_var); gt2
p_gt2 <- gt2$`Pr(>F)`[2] # p-value


# determining the causal effect size
# run the VAR (vector autoregressive) model (i.e. which has both X and Y)
var1 <- VAR(y=d_var[,c("v1_obs","v2_obs")], p=p)
summary(var1) 
# acf plot to check autocorr - suggest stationary if acf drops off very quickly
acf(residuals(var1))

# run AR for univariate analysis - models of just the single virus against its 
# own lags (VAR is only for multivariate and if you try to run a univariate analysis
# like this with VAR it will tell you to use ar or arima)
ar_v1=ar(d_var$v1_obs,order=var1$p,aic=F,method="ols")
ar_v2=ar(d_var$v2_obs,order=var1$p,aic=F,method="ols")

# Estimating the effect size using preventable fraction:
# 1-(prediction estimate model with Y only/prediction estimate of model with X and Y)
# interpretation: 
#   0 = no interaction 
# -ve = strong negative interaction 
# +ve = strong positive interaction 


# extracting model fit of simulated data for AR univariate model with Y only
# i.e data - residuals for the model

# v1
v1_obs <- d_var$v1_obs
fit_uni_v1  <- v1_obs[-c(1:var1$p)] - resid(ar_v1)[-c(1:var1$p)] # estimating fit whilst adjusting for lag
expect_length(fit_uni_v1, length(v1_obs) - var1$p)

# v2
v2_obs <- d_var$v2_obs
fit_uni_v2  <- v2_obs[-c(1:var1$p)] - resid(ar_v2)[-c(1:var1$p)] # estimating fit whilst adjusting for lag
expect_length(fit_uni_v2, length(v2_obs) - var1$p)

# model fit of simulated data for VAR model with both X and Y 
# note VAR does both models for v1 ~ v2 and v2 ~ v1 which are held in var1
# v1
fit_multi_v1 <- as.vector(var1$varresult$v1_obs$fitted.values) # NOTE: some fits are -ve
str(fit_multi_v1)
expect_true(length(fit_multi_v1)==length(fit_uni_v1))

# v2
fit_multi_v2 <- as.vector(var1$varresult$v2_obs$fitted.values) 
str(fit_multi_v2)
expect_true(length(fit_multi_v2)==length(fit_uni_v2))

# Estimate prevelence fraction and log residual sum of squares (Barraquand et al. 2021)
# v1
prev_frac_v1 <- 1-(sum(fit_uni_v1)/sum(fit_multi_v1)); prev_frac_v1 # prevalence fraction 
ratio_res_v1 <- log(sum(residuals(var1)[,1]^2)/sum(residuals(ar_v1)^2,na.rm=T)); ratio_res_v1 # log RSS

# v2 
prev_frac_v2 <- 1-(sum(fit_uni_v2)/sum(fit_multi_v2)); prev_frac_v2
ratio_res_v2 <- log(sum(residuals(var1)[,2]^2)/sum(residuals(ar_v2)^2,na.rm=T)); ratio_res_v2 

# Estimating the uncertainty around these statistics # 

# estimating uncertainty from the model by assuming the parameters follow a multivariate 
# normal distribution with vcov from our model 
model_coef <- coef(var1)
model_vcov <- data.frame(vcov(var1))

# the number of variables n 
n <- length(model_coef)


# bootstrapping 
orig_data <- d_var[,c("v1_obs","v2_obs")]
# remove the first p rows of data to account for the lag
lag_orig <- orig_data[-c(1:p),]
# pull out the residuals from my original model 
residuals_orig <- data.frame(residuals(var1))

# creating a function to define the statistics for the tsboot function which will do 
# the resampling in a block wise fashion
boot_func <- function(tseries, orig_data, var_model,p) {
  bootstrap_residuals <- tseries
  bootstrap_data_v1 <- as.vector(orig_data$v1_obs - residuals(var_model)[,"v1_obs"] + bootstrap_residuals$v1_obs)
  bootstrap_data_v2 <- as.vector(orig_data$v2_obs - residuals(var_model)[,"v2_obs"] + bootstrap_residuals$v2_obs)
  bootstrap_data <- data.frame(cbind(v1_obs=bootstrap_data_v1, v2_obs=bootstrap_data_v2))
  
  # perform univariate AR models
  bootstrap_AR_v1 <- ar(bootstrap_data$v1_obs,order=p,aic=F,method="ols")
  bootstrap_AR_v2 <- ar(bootstrap_data$v1_obs,order=p,aic=F,method="ols")
  
  # perform multivariate VAR model 
  bootstrap_VAR <- VAR(y=bootstrap_data, p=p)
  fit_multi_v1 <- as.vector(bootstrap_VAR$varresult$v1_obs$fitted.values)
  fit_multi_v2 <- as.vector(bootstrap_VAR$varresult$v2_obs$fitted.values)
  
  # calculate statistics (prev frac and log RSS)
  prev_frac_v1 <- 1-(sum(fit_uni_v1)/sum(fit_multi_v1)) 
  ratio_res_v1 <- log(sum(residuals(bootstrap_VAR)[,1]^2)/sum(residuals(bootstrap_AR_v1)^2,na.rm=T)) 
  
  prev_frac_v2 <- 1-(sum(fit_uni_v2)/sum(fit_multi_v2))
  ratio_res_v2 <- log(sum(residuals(bootstrap_VAR)[,2]^2)/sum(residuals(bootstrap_AR_v2)^2,na.rm=T)); ratio_res_v2 
  
  # combine results into vector for output
  res <- c(prev_frac_v1, ratio_res_v1, prev_frac_v2, ratio_res_v2)
  return(res)
}

# number of bootstrap replicates
R <- 50
# implement the bootstrap - this approach is doing a block bootstrap with fixed 
# blocks of size 4 weeks 
boot_out <- tsboot(tseries=residuals_orig, statistic=boot_func, R = R, sim="fixed", l=4,
               orig_data = lag_orig, var_model=var1, p=p); boot_out

# check out the bootstrap distributions 
bootstrap_samples <- data.frame(boot_out$t)
names(bootstrap_samples) <- c("prev_frac_v1","ratio_res_v1","prev_frac_v2","ratio_res_v2")
# make data long 
boot_long <- bootstrap_samples %>% gather() 
ggplot(aes(x=value),data=boot_long) + geom_histogram() + facet_grid(.~key) 
# the distributions for prev frac just look wrong 


# calculate the 95% CI for each statistic


