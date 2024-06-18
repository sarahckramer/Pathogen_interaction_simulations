################################################################################
#                       Granger causality analysis        
#
# To run Granger analysis we need the time series to be covariance and mean 
# stationary
#
# inputs: data = data with time, v1_obs, v2_obs
#
# Created by: Sarah Pirikahu
# Creation date: 23 March 2023
################################################################################

# load packages
library(tseries)
library(lmtest) 
library(vars)
library(VARtests)
library(boot)

granger_func <- function(data){
  
  # automatically determine the best lag doing several models with lags
  # 1-12 (approximately 3 month) then choose the best lag number based on BIC
  
  # initialise lists to put results in 
  # lags <- list()
  # lag_v1 <- list()
  # lag_v2 <- list()
  
  # determine the number of lags for each simulated dataset
  df <- data %>% dplyr::select(V1_obs, V2_obs)
  # df$.id <- NULL
  
  lags <- lapply(df, VARselect, lag.max = 15) # lag of approx 3 month
  rm(df)
  
  # pull out the lag with lowest BIC (not BIC is labeled SV)
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1 <- as.numeric(lags$V1_obs$selection[3])
  lag_v2 <- as.numeric(lags$V2_obs$selection[3])
  
  rm(lags)
  
  #---- checking stationary of time series ---# 
  # simply plan to just report if the series is stationary or not rather 
  # than apply difference as this is difficult to do in practice when lags
  # will differ by simulation run
  
  # ADF Test hypotheses
  # H0: there is a unit root - i.e. the time series is not stationary 
  # H1: time series is stationary 
  adf_v1 <- adf.test(data$V1_obs, k = lag_v1); #adf_v1 # k = lag number 
  adf_v2 <- adf.test(data$V2_obs, k = lag_v2); #adf_v2
  
  # KPSS Test hypotheses
  # H0: time series is stationary 
  # H1: non stationary
  kpss_v1 <- kpss.test(data$V1_obs); #kpss_v1
  kpss_v2 <- kpss.test(data$V2_obs); #kpss_v2
  
  # ---- running Granger test and extracting p-values ----#
  # Test hypotheses
  # H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
  # H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 
  
  # estimate seasonal component
  data <- data %>%
    dplyr::select(time, V1_obs, V2_obs) %>%
    mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52 * (time - 26)))
  
  # m1 <- lm(V1_obs ~ lag(V1_obs, 1) + lag(V1_obs, 2) + lag(V1_obs, 3) + V2_obs , data = data)
  # m2 <- lm(V1_obs ~ lag(V1_obs, 1) + lag(V1_obs, 2) + lag(V1_obs, 3) , data = data)
  # anova(m1, m2)
  
  # specifying the lag to be the minimum of the two series
  p <- min(lag_v1, lag_v2)
  
  #---- determine causal effect size ----#
  # run the VAR (vector autoregressive) model (i.e. which has both X and Y)
  # this is essentially what the Granger test does under the hood
  var1 <- VAR(y = data[, c('V1_obs', 'V2_obs')], p = p)
  var1_confound <- VAR(y = data[, c('V1_obs', 'V2_obs')], p = p, exogen = data[, 'seasonal_component'])
  
  # # pull out effect of seasonal component
  # summary(var1_confound)$varresult$V1_obs$coefficients['exo1', 1]
  # summary(var1_confound)$varresult$V2_obs$coefficients['exo1', 1]
  # summary(var1_confound)$varresult$V1_obs$coefficients['exo1', 'Pr(>|t|)']
  # summary(var1_confound)$varresult$V2_obs$coefficients['exo1', 'Pr(>|t|)']
  
  # run AR for univariate analysis - models of just the single virus against its 
  # own lags (VAR is only for multivariate and if you try to run a univariate analysis
  # like this with VAR it will tell you to use ar or arima)
  ar_v1 = ar(data$V1_obs, order = var1$p, aic = F, method="ols")
  ar_v2 = ar(data$V2_obs, order = var1$p, aic = F, method="ols")
  # ar_v1_alt <- VARfit(data$V1_obs, p = var1$p)
  # ar_v2_alt <- VARfit(data$V2_obs, p = var1$p)
  
  ar_v1_confound <- VARfit(data$V1_obs, p = var1$p, exogen = data$seasonal_component)
  ar_v2_confound <- VARfit(data$V2_obs, p = var1$p, exogen = data$seasonal_component)
  
  # Estimating the effect size using log RSS (Barraquand et al. 2021)
  # log RSS = log(RSS_univariate/RSS_multivariate) interpretation: 
  # 0 = no interaction
  # -ve = strong negative interaction 
  # +ve = strong positive interaction 
  logRSS_v1 <- log(sum(ar_v1$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'V1_obs'] ** 2));# logRSS_v1 
  logRSS_v2 <- log(sum(ar_v2$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'V2_obs'] ** 2));# logRSS_v2
  
  logRSS_v1_confound <- log(sum(ar_v1_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var1_confound)[, 'V1_obs'] ** 2))
  logRSS_v2_confound <- log(sum(ar_v2_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var1_confound)[, 'V2_obs'] ** 2))
  
  # get p-values
  # results the same regardless of performance on the raw or normalised data
  p_gt1_wald <- grangertest(V1_obs ~ V2_obs, order = p, data = data)$`Pr(>F)`[2]
  p_gt2_wald <- grangertest(V2_obs ~ V1_obs, order = p, data = data)$`Pr(>F)`[2]
  
  p_gt1_ftest <- causality(var1, cause = 'V1_obs')$Granger$p.value
  p_gt2_ftest <- causality(var1, cause = 'V2_obs')$Granger$p.value
  
  p_gt1_ftest_confound <- causality(var1_confound, cause = 'V1_obs')$Granger$p.value
  p_gt2_ftest_confound <- causality(var1_confound, cause = 'V2_obs')$Granger$p.value
  
  #---- estimating the uncertainty around these statistics ----# 
  # get residuals
  residuals_orig <- data.frame(residuals(var1))
  residuals_orig_confound <- data.frame(residuals(var1_confound))
  
  ## ---- block bootstrapping ----# 
  
  # creating a function to define the statistics for the tsboot function which will do 
  # the resampling in a block wise fashion
  boot_func <- function(tseries, orig_data, var_model, p, exogen = NA) {
    
    bootstrap_residuals <- tseries
    
    # estimate new simulated data by taking original data - residuals + new bootstrapped residual
    bootstrap_data_v1 <- as.vector(orig_data$V1_obs[-c(1:p)] - residuals(var_model)[,'V1_obs'] + bootstrap_residuals$V1_obs)
    bootstrap_data_v2 <- as.vector(orig_data$V2_obs[-c(1:p)] - residuals(var_model)[,'V2_obs'] + bootstrap_residuals$V2_obs)
    bootstrap_data <- data.frame(cbind(V1_obs = bootstrap_data_v1, V2_obs = bootstrap_data_v2))
    
    if (is.na(exogen)) {
      
      # perform univariate AR models
      bootstrap_AR_v1 <- ar(bootstrap_data$V1_obs, order = p, aic = F, method = "ols")
      bootstrap_AR_v2 <- ar(bootstrap_data$V2_obs, order = p, aic = F, method = "ols")
      
      # perform multivariate VAR model 
      bootstrap_VAR <- VAR(y = bootstrap_data, p = p)
      
    } else {
      
      # get seasonal component
      bootstrap_data <- bootstrap_data %>%
        mutate(time = (min(orig_data$time) + p):max(orig_data$time)) %>%
        mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52 * (time - 26))) %>%
        dplyr::select(-time)
      
      # perform univariate AR models
      bootstrap_AR_v1 <- VARfit(bootstrap_data$V1_obs, p = p, exogen = bootstrap_data$seasonal_component)
      bootstrap_AR_v2 <- VARfit(bootstrap_data$V2_obs, p = p, exogen = bootstrap_data$seasonal_component)
      
      # perform multivariate VAR model 
      bootstrap_VAR <- VAR(y = bootstrap_data[, c('V1_obs', 'V2_obs')], p = p, exogen = bootstrap_data[, 'seasonal_component'])
      
    }
    
    # calculate log RS 
    logRSS_v1 <- log(sum(bootstrap_AR_v1$resid ** 2, na.rm=TRUE) / sum(residuals(bootstrap_VAR)[, 1] ** 2)) 
    logRSS_v2 <- log(sum(bootstrap_AR_v2$resid ** 2, na.rm=TRUE) / sum(residuals(bootstrap_VAR)[, 2] ** 2))
    
    # combine results into vector for output
    res <- c(logRSS_v1, logRSS_v2)
    return(res)
    
  }
  
  # generate bootstrapped replicates in blocks of 4 weeks
  boot_out <- tsboot(tseries = residuals_orig, statistic = boot_func, R = 100, sim = 'fixed', l = 4,
                     orig_data = data, var_model = var1, p = p, exogen = NA)
  boot_out_confound <- tsboot(tseries = residuals_orig, statistic = boot_func, R = 100, sim = 'fixed', l = 4,
                              orig_data = data, var_model = var1, p = p, exogen = 'seasonal')
  
  # check out the bootstrap distributions
  bootstrap_samples <- boot_out$t %>%
    as_tibble() %>%
    rename('logRSS_v1_x_v2' = 'V1',
           'logRSS_v2_x_v1' = 'V2')
  bootstrap_samples_confound <- boot_out_confound$t %>%
    as_tibble() %>%
    rename('logRSS_v1_x_v2' = 'V1',
           'logRSS_v2_x_v1' = 'V2')
  
  # calculate the 95% CI for each statistic
  # standard approach assuming normality of the sampling dist
  # v1 x v2
  CI_lower95_v1_x_v2 <- logRSS_v1 - 1.96 * sd(bootstrap_samples$logRSS_v1_x_v2)
  CI_upper95_v1_x_v2 <- logRSS_v1 + 1.96 * sd(bootstrap_samples$logRSS_v1_x_v2)
  # v2 x v1
  CI_lower95_v2_x_v1 <- logRSS_v2 - 1.96 * sd(bootstrap_samples$logRSS_v2_x_v1)
  CI_upper95_v2_x_v1 <- logRSS_v2 + 1.96 * sd(bootstrap_samples$logRSS_v2_x_v1)
  
  # v1 x v2
  CI_lower95_v1_x_v2_confound <- logRSS_v1_confound - 1.96 * sd(bootstrap_samples_confound$logRSS_v1_x_v2)
  CI_upper95_v1_x_v2_confound <- logRSS_v1_confound + 1.96 * sd(bootstrap_samples_confound$logRSS_v1_x_v2)
  # v2 x v1
  CI_lower95_v2_x_v1_confound <- logRSS_v2_confound - 1.96 * sd(bootstrap_samples_confound$logRSS_v2_x_v1)
  CI_upper95_v2_x_v1_confound <- logRSS_v2_confound + 1.96 * sd(bootstrap_samples_confound$logRSS_v2_x_v1)
  
  # output results
  res <- data.frame(cbind(logRSS = c(logRSS_v1, logRSS_v2, logRSS_v1_confound, logRSS_v2_confound),
                          direction = rep(c('v2 -> v1', 'v1 -> v2'), 2),
                          confounding = rep(c('none', 'none', 'seasonal', 'seasonal')),
                          CI_lower = c(CI_lower95_v1_x_v2, CI_lower95_v2_x_v1, CI_lower95_v1_x_v2_confound, CI_lower95_v2_x_v1_confound),
                          CI_upper = c(CI_upper95_v1_x_v2, CI_upper95_v2_x_v1, CI_upper95_v1_x_v2_confound, CI_upper95_v2_x_v1_confound),
                          granger_p = c(p_gt1_wald, p_gt2_wald, NA, NA),
                          ftest_p = c(p_gt1_ftest, p_gt2_ftest, p_gt1_ftest_confound, p_gt2_ftest_confound),
                          adf_p = c(adf_v1$p.value, adf_v2$p.value, NA, NA),
                          kpss_p = c(kpss_v1$p.value, kpss_v2$p.value, NA, NA))) %>%
    as_tibble() %>%
    mutate(logRSS = as.numeric(logRSS),
           CI_lower = as.numeric(CI_lower),
           CI_upper = as.numeric(CI_upper),
           granger_p = as.numeric(granger_p),
           adf_p = as.numeric(adf_p),
           kpss_p = as.numeric(kpss_p))
  return(res)
  
}
