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
  
  # GAM w/o confounding:
  
  # run gam model
  mvn_mod <- gam(formula = list(V1_obs ~ s(time), V2_obs ~ s(time)),
                 family = mvn(d = 2), # multivariate normal distribution of dimension 2
                 data = data)
  
  # extract model residuals
  orig_resid <- residuals(mvn_mod, type = 'response')
  
  # pull out the covariance matrix and then calculate the correlation matrix. 
  # output: 2 x 2 symmetric matrix
  corr_mat <- mvn_mod$family$data$R %>%
    crossprod() %>%
    solve() %>%
    cov2cor()
  
  # GAM w/ seasonal confounding:
  
  # calculate the shared seasonal component 
  data <- data %>%
    mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52.25 * (data$time - 26)))
  
  # run gam model
  mvn_mod_confound <- gam(formula = list(V1_obs ~ s(time) + s(seasonal_component), V2_obs ~ s(time) +s(seasonal_component)),
                          family = mvn(d = 2), # multivariate normal distribution of dimension 2
                          data = data)
  
  # extract model residuals
  orig_resid_confound <- residuals(mvn_mod_confound, type = 'response')
  
  # pull out the covariance matrix and then calculate the correlation matrix. 
  # output: 2 x 2 symmetric matrix
  corr_mat_confound <- mvn_mod_confound$family$data$R %>%
    crossprod() %>%
    solve() %>%
    cov2cor() 
  
  # estimate confidence interval for elements of correlation matrix
  # using block bootstrapping of the residuals
  R <- 100 # number of bootstrap replicates to do 
  boot_res <- NULL # initialising vector to save results to
  
  # creating function to run in the tsboot function 
  boot_func <- function(tseries, orig_data, var_model, var_model_confound) {
    
    # generating new simulated data not accounting for confounding
    bootstrap_residuals <- data.frame(tseries[, c(1:2)])
    names(bootstrap_residuals) <- c("V1_obs", "V2_obs")
    bootstrap_data_v1 <- as.vector(orig_data$V1_obs - residuals(var_model, type = 'response')[,1] + bootstrap_residuals$V1_obs)
    bootstrap_data_v2 <- as.vector(orig_data$V2_obs - residuals(var_model, type = 'response')[,2] + bootstrap_residuals$V2_obs)

    # generating new simulated data accounting for confounding
    bootstrap_residuals_confound <- data.frame(tseries[, c(3:4)])
    names(bootstrap_residuals_confound) <- c("V1_obs", "V2_obs")
    bootstrap_data_v1_confound <- as.vector(orig_data$V1_obs - residuals(var_model_confound, type = 'response')[,1] + bootstrap_residuals_confound$V1_obs)
    bootstrap_data_v2_confound <- as.vector(orig_data$V2_obs - residuals(var_model_confound, type = 'response')[,2] + bootstrap_residuals_confound$V2_obs)
    
    # making the dataframes
    boot_data <- orig_data %>%
      dplyr::select(time) %>%
      mutate(V1_obs = bootstrap_data_v1,
             V2_obs = bootstrap_data_v2)
    
    # confounding dataframe
    boot_data_confound <- orig_data %>%
      dplyr::select(time) %>%
      mutate(V1_obs = bootstrap_data_v1_confound,
             V2_obs = bootstrap_data_v2_confound)
    
    # calculating seasonal component
    boot_data_confound <- boot_data_confound %>%
      mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52.25 * (data$time - 26)))
    
    # model not accounting for seasonal confounding
    boot_mvn_mod <- gam(formula = list(V1_obs ~ s(time), V2_obs ~ s(time)),
                        family = mvn(d = 2), # multivariate normal distribution of dimension 2
                        data = boot_data)
    
    boot_corr_mat <- boot_mvn_mod$family$data$R %>%
      crossprod() %>%
      solve() %>%
      cov2cor() 
    
    # model accounting for seasonal confounding
    boot_mvn_mod_confound <- gam(formula = list(V1_obs ~ s(time) + s(seasonal_component), V2_obs ~ s(time) + s(seasonal_component)),
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
  boot_out <- tsboot(tseries = cbind(orig_resid, orig_resid_confound), statistic = boot_func, R = R, sim = "fixed",
                     l = 10, orig_data = data, var_model = mvn_mod,  var_model_confound = mvn_mod_confound)
  
  # check out the correlation bootstrap distribution 
  boot_corr <- data.frame(boot_out$t)
  names(boot_corr) <- c("cor", "cor_confound")
  
  # get 95% CIs 
  # standard approach assuming normality of the sampling dist
  
  # no confounding 
  CI_lower95 <-  corr_mat[2, 1] - 1.96 * sd(boot_corr$cor)
  CI_upper95 <-  corr_mat[2, 1] + 1.96 * sd(boot_corr$cor)
  
  # confounding
  CI_lower95_confound <-  corr_mat_confound[2, 1] - 1.96 * sd(boot_corr$cor_confound)
  CI_upper95_confound <-  corr_mat_confound[2, 1] + 1.96 * sd(boot_corr$cor_confound)
  
  # summarise results to output
  res <- data.frame(cbind(cor = corr_mat[2, 1], CI_lower95 = CI_lower95, CI_upper95 = CI_upper95,
                          cor_confound = corr_mat_confound[2, 1], CI_lower95_confound = CI_lower95_confound, 
                          CI_upper95_confound = CI_upper95_confound), row.names = '')
  return(res)
  
}
