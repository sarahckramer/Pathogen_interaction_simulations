# ---------------------------------------------------------------------------------------------------------------------
# Code to run generalized additive models (GAMs)
# ---------------------------------------------------------------------------------------------------------------------

# load packages
library(mgcv)
library(gratia)
library(brms)

gam_cor <- function(data){ 
  
  # GAM w/o confounding:
  
  if (nrow(data) > 0) {
    
    # calculate time of year:
    data <- data %>%
      mutate(week = week(date))
    
    # run gam model (brms - Bayesian)
    mvn_mod_form <- bf(mvbind(V1_obs_ln, V2_obs_ln) ~ s(week, bs = 'cc', k = 53)) + set_rescor(TRUE)
    mvn_mod <- brm(mvn_mod_form, data = data, chains = 4, cores = 4, warmup = 2000, iter = 3000)#, control = list(adapt_delta = 0.9))
    corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__V1obsln__V2obsln')$rescor__V1obsln__V2obsln
    
    # check for divergent transitions
    # https://github.com/paul-buerkner/brms/issues/97
    n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
      sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
    }) %>% unlist() %>% sum()
    
    # get results and 95% CIs
    res <- data.frame(cbind(cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                      row.names = '')
    
  } else {
    res <- data.frame(cbind(cor = NA, cor_median = NA, CI_lower95 = NA, CI_upper95 = NA, n_div = NA),
                      row.names = '')
  }
  
  # return all results
  return(res)
  
}
