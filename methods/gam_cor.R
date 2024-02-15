################################################################################
#                       Generalised additive model (GAM)         
#
# The idea is to run the GAM to get the joint covariance matrix and then from
# that covariance matrix calculate the correlation matrix. Therefore we can 
# calculate the correlation whilst taking into account the autocorrelation within
# our data unlike simple Pearson's correlation
#
# input: data = dataset with time, v1_obs, v2_obs
#
# Created by: Sarah Pirikahu
# Creation date: 25 Aug 2023
################################################################################

# load packages
library(mgcv)
library(boot)

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
  
  # estimating confidence interval for elements of correlation matrix
  # using block bootstrapping of the residuals
  R <- 500 # number of bootstrap replicates to do 
  boot_res <- NULL # initialising vector to save results to
  
  # creating function to run in the tsboot function 
  boot_func <- function(tseries, orig_data=data, var_model=mvn_mod) {
    bootstrap_residuals <- data.frame(tseries)
    names(bootstrap_residuals) <- c("v1_obs", "v2_obs")
    bootstrap_data_v1 <- as.vector(orig_data$v1_obs - residuals(var_model)[,1] + bootstrap_residuals$v1_obs)
    bootstrap_data_v2 <- as.vector(orig_data$v2_obs - residuals(var_model)[,2] + bootstrap_residuals$v2_obs)
   
    boot_data <- data.frame(cbind(orig_data$time, bootstrap_data_v1, bootstrap_data_v2))
    names(boot_data) <- c("time", "v1_obs", "v2_obs")
    
    boot_mvn_mod <- gam(formula = list(v1_obs ~ s(time), v2_obs ~ s(time)),
                     family = mvn(d = 2), # multivariate normal distribution of dimension 2
                     data = boot_data)
    
    boot_corr_mat <- boot_mvn_mod$family$data$R %>%
        crossprod() %>%
        solve() %>%
        cov2cor() 
  
    cor <- boot_corr_mat[2,1]
    boot_res <- c(boot_res, cor)  
    return(boot_res)
  }
  
  # do the block resampling with 4 week size blocks
  boot_out <- tsboot(tseries=orig_resid, statistic=boot_func, R = R, sim="fixed", l=4,
                     orig_data = data, var_model=mvn_mod); boot_out
  
  # check out the correlation bootstrap distribution 
  boot_corr <- data.frame(boot_out$t)
  
  # get 95% CIs 
  # standard approach assuming normality of the sampling dist
  CI_lower95 <-  corr_mat[2,1] - 1.96*sd(boot_corr$boot_out.t)
  CI_upper95 <-  corr_mat[2,1] + 1.96*sd(boot_corr$boot_out.t)
  
  # summarise results to output
  res <- data.frame(cbind(cor = corr_mat[2,1],
                          CI_lower95 = CI_lower95, CI_upper95 = CI_upper95), row.names = "")
  return(res)
}