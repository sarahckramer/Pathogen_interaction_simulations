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
  # run gam model no confounding
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
  
  # run gam model with seasonal confounding 
  # estimate the shared seasonal component 
  omega <- (2 * pi)/52
  A1 <- 0.2
  phi1 <- 26
  data$seasonal_component <- 1 + A1 * cos(omega * (data$time - phi1))
  # run model
  mvn_mod_confound <- gam(formula = list(v1_obs ~ s(time) + s(seasonal_component), v2_obs ~ s(time) +s(seasonal_component)),
                 family = mvn(d = 2), # multivariate normal distribution of dimension 2
                 data = data)
  # extracting model residuals
  orig_resid_confound <- residuals(mvn_mod_confound)
  
  # pull out the covariance matrix and then calculate the correlation matrix. 
  # output: 2 x 2 symmetric matrix
  corr_mat_confound <- mvn_mod_confound$family$data$R %>%
    crossprod() %>%
    solve() %>%
    cov2cor() 
  
  # estimating confidence interval for elements of correlation matrix
  # using block bootstrapping of the residuals
  R <- 100 # number of bootstrap replicates to do 
  boot_res <- NULL # initialising vector to save results to
  
  # creating function to run in the tsboot function 
  boot_func <- function(tseries, orig_data=data, var_model=mvn_mod, var_model_confound=mvn_mod_confound) {

    # generating new simulated data not accounting for confounding 
    bootstrap_residuals <- data.frame(tseries[,c(1:2)])
    names(bootstrap_residuals) <- c("v1_obs", "v2_obs")
    bootstrap_data_v1 <- as.vector(orig_data$v1_obs - residuals(var_model)[,1] + bootstrap_residuals$v1_obs)
    bootstrap_data_v2 <- as.vector(orig_data$v2_obs - residuals(var_model)[,2] + bootstrap_residuals$v2_obs)
   
    # generating new simulated data accounting for confounding 
    bootstrap_residuals_confound <- data.frame(tseries[,c(3:4)])
    names(bootstrap_residuals_confound) <- c("v1_obs", "v2_obs")
    bootstrap_data_v1_confound <- as.vector(orig_data$v1_obs - residuals(var_model_confound)[,1] + bootstrap_residuals_confound$v1_obs)
    bootstrap_data_v2_confound <- as.vector(orig_data$v2_obs - residuals(var_model_confound)[,2] + bootstrap_residuals_confound$v2_obs)
    
    # making the dataframes  
    boot_data <- data.frame(cbind(orig_data$time, bootstrap_data_v1, bootstrap_data_v2))
    names(boot_data) <- c("time", "v1_obs", "v2_obs")
    
    # confounding dataframe
    boot_data_confound <- data.frame(cbind(orig_data$time, bootstrap_data_v1_confound, bootstrap_data_v2_confound))
    names(boot_data_confound) <- c("time", "v1_obs", "v2_obs")
    # calculating seasonal component
    omega <- (2 * pi)/52
    A1 <- 0.2
    phi1 <- 26
    boot_data_confound$seasonal_component <- 1 + A1 * cos(omega * (data$time - phi1))
    
    # model not accounting for seasonal confounding
    boot_mvn_mod <- gam(formula = list(v1_obs ~ s(time), v2_obs ~ s(time)),
                     family = mvn(d = 2), # multivariate normal distribution of dimension 2
                     data = boot_data)

    boot_corr_mat <- boot_mvn_mod$family$data$R %>%
        crossprod() %>%
        solve() %>%
        cov2cor() 
  
    # model accounting for seasonal confounding
    boot_mvn_mod_confound <- gam(formula = list(v1_obs ~ s(time) + s(seasonal_component), v2_obs ~ s(time) + s(seasonal_component)),
                        family = mvn(d = 2), # multivariate normal distribution of dimension 2
                        data = boot_data_confound)
    
    boot_corr_mat_confound <- boot_mvn_mod_confound$family$data$R %>%
      crossprod() %>%
      solve() %>%
      cov2cor() 
    
    cor <- boot_corr_mat[2,1]
    cor_confound <- boot_corr_mat_confound[2,1]
    #cor_both <- data.frame(cor = cor, cor_confound = cor_confound)
    #boot_res <- rbind(boot_res, cor_both)  
    boot_res <- c(cor, cor_confound)
    return(boot_res)
  }
  
  # do the block resampling with 4 week size blocks
  boot_out <- tsboot(tseries=cbind(orig_resid,orig_resid_confound), statistic=boot_func, R = R,
                     sim="fixed", l=4, orig_data = data, var_model=mvn_mod,  var_model_confound=mvn_mod_confound); boot_out
  
  # check out the correlation bootstrap distribution 
  boot_corr <- data.frame(boot_out$t)
  names(boot_corr) <- c("cor", "cor_confound")
  
  # get 95% CIs 
  # standard approach assuming normality of the sampling dist
  
  # no confounding 
  CI_lower95 <-  corr_mat[2,1] - 1.96*sd(boot_corr$cor)
  CI_upper95 <-  corr_mat[2,1] + 1.96*sd(boot_corr$cor)
  # confounding
  CI_lower95_confound <-  corr_mat_confound[2,1] - 1.96*sd(boot_corr$cor_confound)
  CI_upper95_confound <-  corr_mat_confound[2,1] + 1.96*sd(boot_corr$cor_confound)
  
  
  # summarise results to output
  res <- data.frame(cbind(cor = corr_mat[2,1], CI_lower95 = CI_lower95, CI_upper95 = CI_upper95,
                          cor_confound = corr_mat_confound[2,1], CI_lower95_confound = CI_lower95_confound, 
                          CI_upper95_confound = CI_upper95_confound), row.names = "")
  return(res)
}