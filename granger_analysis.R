##########################################
#       Granger causality analysis        
#
# Created by: Sarah Pirikahu
# Creation date: 23 March 2023
##########################################

# testing to determine if any differencing is required
ndiffs(d1$H1, alpha=0.05, test=c("kpss")) # 0 differences required - suggests the series is potentially stationary 
ndiffs(d1$H2, alpha=0.05, test=c("kpss")) # 0 

# create dataset with just the observed cases
d_var <- d1[,c("H1_obs", "H2_obs")]

# Automatically determine the best lag to use based on AIC
lags <- lapply(d_var, VARselect) 
# pull out the lag with best AIC. Lower AIC = better 
lag_h1 <- as.numeric(lags$H1_obs$selection[1])
lag_h2 <- as.numeric(lags$H2_obs$selection[2])

# checking stationarity of time series
# Test hypotheses
# H0: there is a unit root - i.e. the time series is not stationary 
# H1: time series is stationary 
adf.test(d1$H1_obs,k=lag_h1) # k = lag number 
adf.test(d1$H2_obs,k=lag_h2) 

# running Granger test and extracting p-values
# Test hypotheses
# H0: b1 = b2 = ... = bk = 0 the lags of x provide no additional information about y beyond the lags of y 
# H1: there exists 1 <= i <= k so that bi \neq 0 at least one x lag provides additional information 

gt1 <- grangertest(H1_obs ~ H2_obs, order = lag_h1,data=d1); gt1
p_gt1 <- gt1$`Pr(>F)`[2] # p-value

gt2 <- grangertest(H2_obs ~ H1_obs, order = lag_h2,data=d1); gt2
p_gt2 <- gt2$`Pr(>F)`[2] # p-value

# determining the causal effect size
# run the VAR (vector autoregressive) model (i.e. which has both X and Y)
var1 <- VAR(y=d_var, p=min(lag_h1,lag_h2))
summary(var1) 

# run AR for univariate analysis (Note: VAR is only for multivariate
# and if you try to run a univariate analysis like this with VAR it will tell 
# you to use ar or arima)
ar_H1=ar(d_var$H1_obs,order=var1$p,aic=F,method="ols")
ar_H2=ar(d_var$H2_obs,order=var1$p,aic=F,method="ols")




log_21_inter=log(sum((ar_H1$resid)^2,na.rm=T)/sum((var1$varresult$H1$residuals)^2,na.rm=T))
log_12_inter=log(sum((ar_H2$resid)^2,na.rm=T)/sum((var1$varresult$H2$residuals)^2,na.rm=T))

log_21_inter # effect of H2 on H1 
log_12_inter # effect of H1 on H2