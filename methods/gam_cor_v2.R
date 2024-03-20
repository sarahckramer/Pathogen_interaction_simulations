################################################################################
#                       Generalised additive model (GAM)         
#
# The idea is to run the GAM to get the joint covariance matrix and then from
# that covariance matrix calculate the correlation matrix. Therefore we can 
# calculate the correlation whilst taking into account the autocorrelation within
# our data unlike simple Pearson's correlation. In v2 of this function we simply 
# calculate the point estimate without the confidence interval.... will in future 
# include lags
#
# input: data = dataset with time, v1_obs, v2_obs
#
# Created by: Sarah Pirikahu
# Creation date: 25 Aug 2023
################################################################################

# load packages
library(mgcv)
library(vars)

gam_cor <- function(data){ 
  # run gam model
  mvn_mod <- gam(formula = list(v1_obs ~ s(time), v2_obs ~ s(time)),
                 family = mvn(d = 2), # multivariate normal distribution of dimension 2
                 data = data)
  # extracting model residuals
  orig_resid <- residuals(mvn_mod)
  
  # pull out the covariance matrix and then calculate the correlation matrix. 
  # output: 2 x 2 symmetric matrix
  corr_mat <- mvn_mod$family$data$R %>%
    crossprod() %>%
    solve() %>%
    cov2cor() 
  
    cor <- data.frame(cor=corr_mat[2,1])
    return(cor)
}