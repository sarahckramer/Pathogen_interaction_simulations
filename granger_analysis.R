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
library(meboot)

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

# Decided I will just record if stationarity is met or not rather than differencing both series
# Ask Mathheiu... 


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
#acf(residuals(var1))

# run AR for univariate analysis - models of just the single virus against its 
# own lags (VAR is only for multivariate and if you try to run a univariate analysis
# like this with VAR it will tell you to use ar or arima)
ar_v1=ar(d_var$v1_obs,order=var1$p,aic=F,method="ols")
ar_v2=ar(d_var$v2_obs,order=var1$p,aic=F,method="ols")

# Estimating the effect size using preventable fraction and log RSS

# preventable fraction: 1-(prediction estimate model with Y only/prediction estimate of model with X and Y)
# interpretation: 
#   0 = no interaction 
# -ve = strong negative interaction 
# +ve = strong positive interaction 

# log RSS = log(RSS_multivariate/RSS_univariate) interpretation: ** CHECK
# 0 = no interaction
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
# bootstrapping 
orig_data <- d_var[,c("v1_obs","v2_obs")]
# pull out the residuals from my original model 
residuals_orig <- data.frame(residuals(var1))

##---- Maximum entropy bootstrapping ----# 
# number of replicates
R <- 300

# resampling the residuals
meboot_v1out <- meboot(x=residuals_orig$v1_obs, reps=R, trim=0.10, elaps=TRUE)
meboot_v2out <- meboot(x=residuals_orig$v2_obs, reps=R, trim=0.10, elaps=TRUE)

# estimating the statistics
me_temp_res <- NULL
for(j in 1:R){
  meboot_data_v1 <- as.vector(orig_data$v1_obs[-c(1:p)] - residuals(var1)[,"v1_obs"] + meboot_v1out$ensemble[,j])
  meboot_data_v2 <- as.vector(orig_data$v2_obs[-c(1:p)] - residuals(var1)[,"v2_obs"] + meboot_v1out$ensemble[,j])
  meboot_data <- data.frame(cbind(v1_obs=meboot_data_v1, v2_obs=meboot_data_v2))
  
  # perform univariate AR models
  meboot_AR_v1 <- ar(meboot_data$v1_obs,order=p,aic=F,method="ols")
  meboot_AR_v2 <- ar(meboot_data$v2_obs,order=p,aic=F,method="ols")
  mefit_uni_v1  <- orig_data$v1_obs[-c(1:p)] - resid(meboot_AR_v1) 
  mefit_uni_v2  <- orig_data$v2_obs[-c(1:p)] - resid(meboot_AR_v2)
  
  # perform multivariate VAR model 
  meboot_VAR <- VAR(y=meboot_data, p=p)
  mefit_multi_v1 <- as.vector(meboot_VAR$varresult$v1_obs$fitted.values)
  mefit_multi_v2 <- as.vector(meboot_VAR$varresult$v2_obs$fitted.values)
  
  # calculate statistics (prev frac and log RSS)
  me_prev_frac_v1 <- 1-(sum(mefit_uni_v1, na.rm=T)/sum(mefit_multi_v1)) 
  me_ratio_res_v1 <- log(sum(residuals(meboot_VAR)[,1]^2)/sum(residuals(meboot_AR_v1)^2,na.rm=T)) 
  
  me_prev_frac_v2 <- 1-(sum(mefit_uni_v2, na.rm=T)/sum(mefit_multi_v2))
  me_ratio_res_v2 <- log(sum(residuals(meboot_VAR)[,2]^2)/sum(residuals(meboot_AR_v2)^2,na.rm=T)) 
  # collect the results
  res_out <- c(me_prev_frac_v1,me_ratio_res_v1,me_prev_frac_v2, me_ratio_res_v2)
  me_temp_res <- rbind(me_temp_res, res_out)
}
me_temp_res <- data.frame(me_temp_res, row.names=NULL)
names(me_temp_res) <- c("prev_frac_v1","ratio_res_v1","prev_frac_v2", "ratio_res_v2")

# plot bootstrap replicates 
# make data long 
meboot_long <- me_temp_res %>% gather() 
ggplot(aes(x=value),data=meboot_long) + geom_histogram() + facet_grid(.~key, scales="free_x") 

# estimating ME 95% CIs 
sd_meboot <- apply(me_temp_res, 2,sd)
# original estimates
orig_est <- c(prev_frac_v1, ratio_res_v1, prev_frac_v2, ratio_res_v2)
# standard
me_CI_lower95 <- orig_est - 1.96*sd_meboot
me_CI_upper95 <- orig_est + 1.96*sd_meboot
# percentile boot
me_percCI_lower95 <- apply(me_temp_res, 2, quantile, prob=0.025)
me_percCI_upper95 <- apply(me_temp_res, 2, quantile, prob=0.975)

# rm unused variables now
rm(meboot_VAR, mefit_multi_v1, mefit_multi_v2, meboot_AR_v1, meboot_AR_v2, meboot_data,
   meboot_data_v1, meboot_data_v2,mefit_uni_v1, mefit_uni_v2, prev_frac_v1, prev_frac_v2, 
   ratio_res_v1, ratio_res_v2, meboot_long)

##---- Block bootstrapping -----# 

# creating a function to define the statistics for the tsboot function which will do 
# the resampling in a block wise fashion
boot_func <- function(tseries, orig_data, var_model,p) {
  bootstrap_residuals <- tseries
  bootstrap_data_v1 <- as.vector(orig_data$v1_obs[-c(1:p)] - residuals(var_model)[,"v1_obs"] + bootstrap_residuals$v1_obs)
  bootstrap_data_v2 <- as.vector(orig_data$v2_obs[-c(1:p)] - residuals(var_model)[,"v2_obs"] + bootstrap_residuals$v2_obs)
  bootstrap_data <- data.frame(cbind(v1_obs=bootstrap_data_v1, v2_obs=bootstrap_data_v2))
  
  # perform univariate AR models
  bootstrap_AR_v1 <- ar(bootstrap_data$v1_obs,order=p,aic=F,method="ols")
  bootstrap_AR_v2 <- ar(bootstrap_data$v1_obs,order=p,aic=F,method="ols")
  fit_uni_v1  <- orig_data$v1_obs[-c(1:p)] - resid(bootstrap_AR_v1)
  fit_uni_v2  <- orig_data$v2_obs[-c(1:p)] - resid(bootstrap_AR_v2) 
  
  # perform multivariate VAR model 
  bootstrap_VAR <- VAR(y=bootstrap_data, p=p)
  fit_multi_v1 <- as.vector(bootstrap_VAR$varresult$v1_obs$fitted.values)
  fit_multi_v2 <- as.vector(bootstrap_VAR$varresult$v2_obs$fitted.values)
  
  # calculate statistics (prev frac and log RSS)
  prev_frac_v1 <- 1-(sum(fit_uni_v1,na.rm=T)/sum(fit_multi_v1)) 
  ratio_res_v1 <- log(sum(residuals(bootstrap_VAR)[,1]^2)/sum(residuals(bootstrap_AR_v1)^2,na.rm=T)) 
  
  prev_frac_v2 <- 1-(sum(fit_uni_v2,na.rm=T)/sum(fit_multi_v2))
  ratio_res_v2 <- log(sum(residuals(bootstrap_VAR)[,2]^2)/sum(residuals(bootstrap_AR_v2)^2,na.rm=T)) 
  
  # combine results into vector for output
  res <- c(prev_frac_v1, ratio_res_v1, prev_frac_v2, ratio_res_v2)
  return(res)
}

# implement the bootstrap - this approach is doing a block bootstrap with fixed 
# blocks of size 4 weeks 
boot_out <- tsboot(tseries=residuals_orig, statistic=boot_func, R = R, sim="fixed", l=4,
               orig_data = orig_data, var_model=var1, p=p); boot_out

# check out the bootstrap distributions 
bootstrap_samples <- data.frame(boot_out$t)
names(bootstrap_samples) <- c("prev_frac_v1","ratio_res_v1","prev_frac_v2","ratio_res_v2")
# make data long 
boot_long <- bootstrap_samples %>% tidyr::gather() 
ggplot(aes(x=value),data=boot_long) + geom_histogram() + facet_grid(.~key, scales="free_x") 

# calculate the 95% CI for each statistic
# standard approach assuming normality of the sampling dist
CI_lower95 <- summary(boot_out)$original - 1.96*summary(boot_out)$bootSE
CI_upper95 <- summary(boot_out)$original + 1.96*summary(boot_out)$bootSE

# percentile bootstrap 
CIperc_lower95 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.025))
CIperc_upper95 <- as_vector(apply(bootstrap_samples, 2, quantile, probs = 0.975))


# output results 
temp_res <- data.frame(cbind(original_estimate = summary(boot_out)$original,  
             blockbootMean = apply(boot_out$t,2,mean),
             blockbootBias = summary(boot_out)$bootBias, # difference between mean estimate and the original
             blockbootSE = summary(boot_out)$bootSE,
             blockboot_CI_lower95 = CI_lower95,
             blockboot_CI_upper95 = CI_upper95,
             blockboot_CIperc_lower95 = CIperc_lower95,
             blockboot_CIperc_upper95 = CIperc_upper95,
             mebootMean = apply(me_temp_res,2,mean),
             mebootSE = sd_meboot,
             me_CI_lower95 = me_CI_lower95,
             me_CI_upper95 = me_CI_upper95,
             me_percCI_lower95 = me_percCI_lower95,
             me_percCI_upper95 = me_percCI_upper95),
             granger_p = c(rep(p_gt1,2), rep(p_gt2,2)),
             adf_p = c(rep(adf_v1$p.value,2), rep(adf_v2$p.value,2)),
             kpss_p = c(rep(kpss_v1$p.value,2), rep(kpss_v2$p.value,2)))
temp_res$statistic <- rownames(temp_res)
temp_res <- data.frame(temp_res, row.names=NULL)

# create a list of outputs which includes the results and the bootstrap 
# distributions
res_list <- list(summary = temp_res, block_bootstraps = boot_out$t, me_bootstraps = me_temp_res)

results[[i]]$granger <- res_list


rm(adf_v1, adf_v2, kpss_v1, kpss_v2, temp_res, gt1, gt2, p_gt1, p_gt2, orig_data,
   boot_out, meboot_v1out, meboot_v2out, temp_res, me_temp_res, res_list, meboot_long,
   var1, ar_v1, ar_v2, residuals_orig, CI_lower95, CI_upper95, CIperc_lower95,
   CIperc_upper95)
             