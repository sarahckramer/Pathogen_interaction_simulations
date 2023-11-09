################################################################################
#                       Granger causality analysis        
#
# To run Granger analysis we need the time series to be covariance and mean 
# stationary
#
# inputs: data = data with time, v1_obs, v2_obs
#         lag_v1 = number of lags to use for v1_obs time series
#         lag_v2 = number of lags to use for v2_obs time series 
#
# Created by: Sarah Pirikahu
# Creation date: 23 March 2023
################################################################################

granger_func <- function(data, lag_v1, lag_v2){
  # load packages
  library(tseries) 
  library(lmtest) 
  library(vars) 
  library(boot)
  
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
  
  #----- Estimating the uncertainty around these statistics ----# 
  # pulling out original data, residuals and estimates 
  orig_data <- data[,c("v1_obs","v2_obs")]
  residuals_orig <- data.frame(residuals(var1))
  orig_est <- c(logRSS_v1, logRSS_v2)

  ##---- Block bootstrapping -----# 
  
  # creating a function to define the statistics for the tsboot function which will do 
  # the resampling in a block wise fashion
  boot_func <- function(tseries, orig_data, var_model,p) {
    bootstrap_residuals <- tseries
    # estimating new simulated data by taking original data - residuals + new bootstrapped residual
    bootstrap_data_v1 <- as.vector(orig_data$v1_obs[-c(1:p)] - residuals(var_model)[,"v1_obs"] + bootstrap_residuals$v1_obs)
    bootstrap_data_v2 <- as.vector(orig_data$v2_obs[-c(1:p)] - residuals(var_model)[,"v2_obs"] + bootstrap_residuals$v2_obs)
    bootstrap_data <- data.frame(cbind(v1_obs=bootstrap_data_v1, v2_obs=bootstrap_data_v2))
    
    # perform univariate AR models
    bootstrap_AR_v1 <- ar(bootstrap_data$v1_obs,order=p,aic=F,method="ols")
    bootstrap_AR_v2 <- ar(bootstrap_data$v2_obs,order=p,aic=F,method="ols")
    
    # perform multivariate VAR model 
    bootstrap_VAR <- VAR(y=bootstrap_data, p=p)
    
    # calculate log RS 
    logRSS_v1 <- log(sum(residuals(bootstrap_AR_v1)^2,na.rm=T)/sum(residuals(bootstrap_VAR)[,1]^2)) 
    logRSS_v2 <- log(sum(residuals(bootstrap_AR_v2)^2,na.rm=T)/sum(residuals(bootstrap_VAR)[,2]^2)) 
    
    # combine results into vector for output
    res <- c(logRSS_v1, logRSS_v2)
    return(res)
  }
 
  # generate bootstrapped replicates in blocks of 4 weeks
  R <- 300 # number of replicates
  boot_out <- tsboot(tseries=residuals_orig, statistic=boot_func, R = R, sim="fixed", l=4,
                 orig_data = orig_data, var_model=var1, p=p); boot_out
  
  # check out the bootstrap distributions 
  bootstrap_samples <- data.frame(boot_out$t)
  names(bootstrap_samples) <- c("logRSS_v1_x_v2","logRSS_v2_x_v1")
  # make data long 
  #boot_long <- bootstrap_samples %>% tidyr::gather() 
  #ggplot(aes(x=value),data=boot_long) + geom_histogram() + facet_grid(.~key, scales="free_x") 
  
  # calculate the 95% CI for each statistic
  # standard approach assuming normality of the sampling dist
  # v1 x v2
  CI_lower95_v1_x_v2 <- logRSS_v1 - 1.96*sd(bootstrap_samples$logRSS_v1_x_v2)
  CI_upper95_v1_x_v2 <- logRSS_v1 + 1.96*sd(bootstrap_samples$logRSS_v1_x_v2)
  # v2 x v1
  CI_lower95_v2_x_v1 <- logRSS_v2 - 1.96*sd(bootstrap_samples$logRSS_v2_x_v1)
  CI_upper95_v2_x_v1 <- logRSS_v2 + 1.96*sd(bootstrap_samples$logRSS_v2_x_v1)
  
  
  # percentile bootstrap 
  # v1 x v2
  CIperc_lower95_v1_x_v2 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.025))[1]
  CIperc_upper95_v1_x_v2 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.975))[1]
  # v2 x v1
  CIperc_lower95_v2_x_v1 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.025))[2]
  CIperc_upper95_v2_x_v1 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.975))[2]
  

  # output results 
  temp_res <- data.frame(cbind(original_estimate = orig_est,  
               blockbootMean = apply(boot_out$t,2,mean),
               blockboot_CI_lower95_v1_x_v2 = c(CI_lower95_v1_x_v2,CI_lower95_v2_x_v1),
               blockboot_CI_upper95_v1_x_v2 = c(CI_upper95_v1_x_v2,CI_upper95_v2_x_v1),
               blockboot_CIperc_lower95_v1_x_v2 = c(CIperc_lower95_v1_x_v2,CIperc_lower95_v2_x_v1),
               blockboot_CIperc_upper95_v1_x_v2 = c(CIperc_upper95_v1_x_v2,CIperc_upper95_v2_x_v1),
               granger_p = c(p_gt1, p_gt2),
               adf_p = c(adf_v1$p.value, adf_v2$p.value),
               kpss_p = c(kpss_v1$p.value, kpss_v2$p.value)))
  temp_res <- data.frame(temp_res, row.names=NULL)
  
  # create a list of outputs which includes the results and the bootstrap 
  # distributions
  res_list <- list(summary = temp_res, block_bootstraps = boot_out$t)
  return(res_list)
}


             