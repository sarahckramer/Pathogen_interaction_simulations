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
  
  # determine the number of lags for each simulated dataset
  df <- data %>% dplyr::select(V1_obs, V2_obs)
  
  lags <- lapply(df, VARselect, lag.max = 15) # lag of approx 3 month
  rm(df)
  
  # pull out the lag with lowest BIC
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1 <- as.numeric(lags$V1_obs$selection[3])
  lag_v2 <- as.numeric(lags$V2_obs$selection[3])
  
  rm(lags)
  
  #---- checking stationary of time series ---# 
  
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
    mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52.25 * (time - 26)))
  
  # specifying the lag to be the minimum of the two series
  p <- min(lag_v1, lag_v2)
  
  #---- determine causal effect size ----#
  # run the VAR (vector autoregressive) model (i.e. which has both X and Y)
  # this is essentially what the Granger test does under the hood
  var1 <- VAR(y = data[, c('V1_obs', 'V2_obs')], p = p)
  var1_confound <- VAR(y = data[, c('V1_obs', 'V2_obs')], p = p, exogen = data[, 'seasonal_component'])
  
  # run AR for univariate analysis - models of just the single virus against its 
  # own lags (VAR is only for multivariate and if you try to run a univariate analysis
  # like this with VAR it will tell you to use ar or arima)
  ar_v1 = ar(data$V1_obs, order = var1$p, aic = F, method="ols")
  ar_v2 = ar(data$V2_obs, order = var1$p, aic = F, method="ols")
  
  sink(file = nullfile())
  ar_v1_confound <- VARfit(data$V1_obs, p = var1$p, exogen = data$seasonal_component)
  ar_v2_confound <- VARfit(data$V2_obs, p = var1$p, exogen = data$seasonal_component)
  sink()
  
  # Estimating the effect size using log RSS (Barraquand et al. 2021)
  # log RSS = log(RSS_univariate/RSS_multivariate)
  logRSS_v1 <- log(sum(ar_v1$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'V1_obs'] ** 2));# logRSS_v1 
  logRSS_v2 <- log(sum(ar_v2$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'V2_obs'] ** 2));# logRSS_v2
  
  logRSS_v1_confound <- log(sum(ar_v1_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var1_confound)[, 'V1_obs'] ** 2))
  logRSS_v2_confound <- log(sum(ar_v2_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var1_confound)[, 'V2_obs'] ** 2))
  
  # get p-values
  # results the same regardless of performance on the raw or normalised data
  p_gt1_wald <- grangertest(V1_obs ~ V2_obs, order = p, data = data)$`Pr(>F)`[2]
  p_gt2_wald <- grangertest(V2_obs ~ V1_obs, order = p, data = data)$`Pr(>F)`[2]
  
  p_gt1_ftest <- causality(var1, cause = 'V2_obs')$Granger$p.value
  p_gt2_ftest <- causality(var1, cause = 'V1_obs')$Granger$p.value
  
  p_gt1_ftest_confound <- causality(var1_confound, cause = 'V2_obs')$Granger$p.value
  p_gt2_ftest_confound <- causality(var1_confound, cause = 'V1_obs')$Granger$p.value
  
  # output results
  res <- data.frame(cbind(logRSS = c(logRSS_v1, logRSS_v2, logRSS_v1_confound, logRSS_v2_confound),
                          direction = rep(c('v2 -> v1', 'v1 -> v2'), 2),
                          confounding = rep(c('none', 'none', 'seasonal', 'seasonal')),
                          granger_p = c(p_gt1_wald, p_gt2_wald, NA, NA),
                          ftest_p = c(p_gt1_ftest, p_gt2_ftest, p_gt1_ftest_confound, p_gt2_ftest_confound),
                          adf_p = c(adf_v1$p.value, adf_v2$p.value, NA, NA),
                          kpss_p = c(kpss_v1$p.value, kpss_v2$p.value, NA, NA))) %>%
    as_tibble() %>%
    mutate(logRSS = as.numeric(logRSS),
           granger_p = as.numeric(granger_p),
           ftest_p = as.numeric(ftest_p),
           adf_p = as.numeric(adf_p),
           kpss_p = as.numeric(kpss_p))
  
  return(res)
  
}
