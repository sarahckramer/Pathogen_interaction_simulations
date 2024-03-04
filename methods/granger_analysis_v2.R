################################################################################
#                       Granger causality analysis        
#
# To run Granger analysis we need the time series to be covariance and mean 
# stationary
#
# inputs: data = data with time, v1_obs, v2_obs
#
# Created by: Sarah Pirikahu
# Creation date: 1 March 2024
################################################################################

# load packages
library(tseries) 
library(lmtest) 
library(vars) 
library(tidyverse)

granger_func <- function(data){
  
  # Automatically determine the best lag doing several models with lags
  # 1-12 (approximately 3 month) then choose the best lag number based on BIC
  
  # initialising lists to put results in 
  lags <- list()
  lag_v1 <- list()
  lag_v2 <- list()
  
  # determine the number of lags for each simulated dataset
  df <- data %>% dplyr::select(v1_obs, v2_obs)
  df$.id <- NULL
  
  lags <- lapply(df, VARselect, lag.max=12) # lag of approx 3 month
  rm(df)
  # pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1 <- as.numeric(lags$v1_obs$selection[3])
  lag_v2 <- as.numeric(lags$v2_obs$selection[3])
  
  rm(lags)
  
  #---- checking stationary of time series ---# 
  # simply plan to just report if the series is stationary or not rather 
  # than apply difference as this is difficult to do in practice when lags
  # will differ by simulation run
  
  # ADF Test hypotheses
  # H0: there is a unit root - i.e. the time series is not stationary 
  # H1: time series is stationary 
  adf_v1 <- adf.test(data$v1_obs,k=lag_v1); adf_v1 # k = lag number 
  adf_v2 <- adf.test(data$v2_obs,k=lag_v2); adf_v2
  
  # KPSS Test hypotheses
  # H0: time series is stationary 
  # H1: non stationary
  kpss_v1 <- kpss.test(data$v1_obs); kpss_v1
  kpss_v2 <- kpss.test(data$v2_obs); kpss_v2
  
  # ------running Granger test and extracting p-values------#
  # Test hypotheses
  # H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
  # H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 
  
  # specifying the lag to be the minimum of the two series
  p <- min(lag_v1,lag_v2)
  # results the same regardless of performance on the raw or normalised data
  gt1 <- grangertest(v1_obs ~ v2_obs, order = p, data=data); gt1
  p_gt1 <- gt1$`Pr(>F)`[2] # p-value
  
  gt2 <- grangertest(v2_obs ~ v1_obs, order = p, data=data); gt2
  p_gt2 <- gt2$`Pr(>F)`[2] # p-value
  
  
  # ----determining the causal effect size ----# 
  # run the VAR (vector autoregressive) model (i.e. which has both X and Y)
  # this is essentially what the Grangers test does under the hood
  var1 <- VAR(y=data[,c("v1_obs","v2_obs")], p=p)
  summary(var1) 
  
  # run AR for univariate analysis - models of just the single virus against its 
  # own lags (VAR is only for multivariate and if you try to run a univariate analysis
  # like this with VAR it will tell you to use ar or arima)
  ar_v1=ar(data$v1_obs,order=var1$p,aic=F,method="ols")
  ar_v2=ar(data$v2_obs,order=var1$p,aic=F,method="ols")
  
  # Estimating the effect size using log RSS (Barraquand et al. 2021)
  # log RSS = log(RSS_multivariate/RSS_univariate) interpretation: 
  # 0 = no interaction
  # -ve = strong negative interaction 
  # +ve = strong positive interaction 
  
  # v1 (v1 ~ v2)
  logRSS_v1 <- log(sum(residuals(ar_v1)^2,na.rm=T)/sum(residuals(var1)[,"v1_obs"]^2)); logRSS_v1 
  # v2 (v2 ~ v1)
  logRSS_v2 <- log(sum(residuals(ar_v2)^2,na.rm=T)/sum(residuals(var1)[,"v2_obs"]^2)); logRSS_v2 
  
  # output results 
  orig_est <- data.frame(direction=c("v2 -> v1", "v1 -> v2"),logRSS = c(logRSS_v1, logRSS_v2))
  temp_res <- data.frame(cbind(orig_est,  
                               adf_p = c(adf_v1$p.value, adf_v2$p.value),
                               kpss_p = c(kpss_v1$p.value, kpss_v2$p.value)))
  return(temp_res)
}
  
  
  