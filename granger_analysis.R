################################################################################
#                       Granger causality analysis        
#
# Useful documentation describing and giving examples:
# 
#
# Created by: Sarah Pirikahu
# Creation date: 23 March 2023
################################################################################

# load packages
library(tseries) # adf test - stationary 
library(forecast) # ndiff - check if differences in lags needed for stationary 
library(lmtest) # granger causality test
library(vars) # VAR model - that which is underlying granger causality test


# testing to determine if any differencing is required
ndiffs(d_var$H1_obs, alpha=0.05, test=c("kpss")) # 0 differences suggests the series is potentially stationary 
ndiffs(d_var$H2_obs, alpha=0.05, test=c("kpss"))  

# checking stationarity of time series
# Test hypotheses
# H0: there is a unit root - i.e. the time series is not stationary 
# H1: time series is stationary 
adf.test(d_var$H1_obs,k=lag_h1) # k = lag number 
adf.test(d_var$H2_obs,k=lag_h2) 
# test gives same result regardless of if the normalised or raw data is used
# Note: stationarity isn't always achieved by lag selection as above.... but it seems like the 
# best approach we have for automating this 

# running Granger test and extracting p-values
# Test hypotheses
# H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
# H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 

# results the same regardless of performance on the raw or normalised data
gt1 <- grangertest(H1_obs ~ H2_obs, order = min(lag_h1,lag_h2),data=d_var); gt1
p_gt1 <- gt1$`Pr(>F)`[2] # p-value

gt2 <- grangertest(H2_obs ~ H1_obs, order = min(lag_h1,lag_h2),data=d_var); gt2
p_gt2 <- gt2$`Pr(>F)`[2] # p-value\
## note if multiple lags used may need to consider a Bonferroni adjustment....

# determining the causal effect size
# run the VAR (vector autoregressive) model (i.e. which has both X and Y)
var1 <- VAR(y=d_var[,c("H1_obs","H2_obs")], p=min(lag_h1,lag_h2))
summary(var1) 
# acf plot to check autocorr - suggets stationarity if acf drops off very quickly
acf(residuals(var1))


var1_NORM <- VAR(y=d_var[,c("H1_obs_NORM","H2_obs_NORM")], p=min(lag_h1,lag_h2))
summary(var1_NORM) # estimates differ for normalised data, but t- and p-values the same as expected
# acf plot
acf(residuals(var1_NORM))

# run AR for univariate analysis - models of just the single virus against its 
# own lags (VAR is only for multivariate and if you try to run a univariate analysis
# like this with VAR it will tell you to use ar or arima)
ar_H1=ar(d_var$H1_obs,order=var1$p,aic=F,method="ols")
ar_H1_NORM=ar(d_var$H1_obs_NORM,order=var1$p,aic=F,method="ols")

ar_H2=ar(d_var$H2_obs,order=var1$p,aic=F,method="ols")
ar_H2_NORM=ar(d_var$H2_obs_NORM,order=var1$p,aic=F,method="ols")

# Barraquand et al. 2021 used the log ratio of the residuals as a point estimate for 
# Granger analysis. Try this out here - note result the same regardless of whether raw or 
# normalised data is used. 
# Limitation: we can determine whether there is an interaction and its strength but not
# its direction.
log_21_inter=log(sum((ar_H1$resid)^2,na.rm=T)/sum((var1$varresult$H1$residuals)^2,na.rm=T))
log_12_inter=log(sum((ar_H2$resid)^2,na.rm=T)/sum((var1$varresult$H2$residuals)^2,na.rm=T))

log_21_inter # effect of H2 on H1 
log_12_inter # effect of H1 on H2

# Estimating the effect size using preventable fraction:
# 1-(prediction estimate model with Y only/prediction estimate of model with X and Y)
# interpretation: 
#   0 = no interaction 
# -ve = strong negative interaction 
# +ve = strong positive interaction 

# ---H1 models ---#

# model fit of simulated data for univariate model with Y only
# i.e data - residuals for the model
H1_obs <- d_var$H1_obs
fit_uni_H1  <- H1_obs[-c(1:var1$p)] - resid(ar_H1)[-c(1:var1$p)]
expect_length(fit_uni_H1, length(H1_obs) - var1$p)

H1_obs_NORM <- d_var$H1_obs_NORM
fit_uni_H1_NORM  <- H1_obs_NORM[-c(1:var1$p)] - resid(ar_H1_NORM)[-c(1:var1$p)]
expect_length(fit_uni_H1_NORM, length(H1_obs_NORM) - var1$p)

# model fit of simulated data for VAR model with both X and Y 
# pull out H1 fitted values from var1 
fit_multi_H1 <- as.vector(var1$varresult$H1_obs$fitted.values)
str(fit_multi_H1)
expect_true(length(fit_multi_H1)==length(fit_uni_H1))

fit_multi_H1_NORM <- as.vector(var1_NORM$varresult$H1_obs_NORM$fitted.values)
str(fit_multi_H1_NORM)
expect_true(length(fit_multi_H1_NORM)==length(fit_uni_H1_NORM))

# plot of data v fit
ts.plot(H1_obs)
points(c(rep(NA, var1$p),fit_uni_H1), type = "l", col = 2, lty = 2)
points(c(rep(NA, var1$p),fit_multi_H1), type = "l", col = 3, lty = 2)

# Estimate preventable fraction
ratio <- fit_uni_H1/fit_multi_H1
prev_frac <- 1-ratio
avg_prev_frac <- mean(prev_frac)
sd_prev_frac <- sd(prev_frac)

ratio_NORM <- fit_uni_H1_NORM/fit_multi_H1_NORM
prev_frac_NORM <- 1-ratio_NORM
avg_prev_frac_NORM <- mean(prev_frac_NORM)
sd_prev_frac_NORM <- sd(prev_frac_NORM)


# --- H2 model ---# 
# model fit Y only 
H2_obs <- d_var$H2_obs
fit_uni_H2  <- H2_obs[-c(1:var1$p)] - resid(ar_H2)[-c(1:var1$p)]
expect_length(fit_uni_H2, length(H2_obs) - var1$p)

# model fit X and Y 
fit_multi_H2 <- as.vector(var1$varresult$H2_obs$fitted.values)
str(fit_multi_H2)
expect_true(length(fit_multi_H2)==length(fit_uni_H2))

# plot of data v fit
ts.plot(H2_obs)
points(c(rep(NA, var1$p),fit_uni_H2), type = "l", col = 2, lty = 2)
points(c(rep(NA, var1$p),fit_multi_H2), type = "l", col = 3, lty = 2)

# estimate preventable fraction 
ratio_H2 <- fit_uni_H2/fit_multi_H2
prev_frac_H2 <- 1-ratio_H2
avg_prev_frac_H2 <- mean(prev_frac_H2)
sd_prev_frac_H2 <- sd(prev_frac_H2)

