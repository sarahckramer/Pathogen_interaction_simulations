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
library(gratia)
library(brms)

gam_cor <- function(data){ 
  
  # GAM w/o confounding:
  
  # calculate time of year:
  data <- data %>%
    mutate(week = week(date))
  
  # run gam model (brms - Bayesian)
  mvn_mod_form <- bf(mvbind(V1_obs, V2_obs) ~ s(week, bs = 'cc', k = 53)) + set_rescor(TRUE)
  mvn_mod <- brm(mvn_mod_form, data = data, chains = 4, cores = 4, warmup = 2000, iter = 3000)#, control = list(max_treedepth = 15, adapt_delta = 0.95))
  corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__V1obs__V2obs')$rescor__V1obs__V2obs
  
  # check for divergent transitions
  # https://github.com/paul-buerkner/brms/issues/97
  n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
    sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
  }) %>% unlist() %>% sum()
  
  # get results and 95% CIs
  res <- data.frame(cbind(cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)[1:13]), n_div = n_divergent), row.names = '')
  
  # return all results
  return(res)
  
}
