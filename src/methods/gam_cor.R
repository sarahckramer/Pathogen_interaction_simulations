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

gam_cor <- function(data){ 
  
  # functions
  posterior_samples_ADAPT <- function(mod, n) {
    # adapated from "posterior_samples" function in package "gratia," which doesn't seem to support multivariate normal distributions yet...
    # https://github.com/gavinsimpson/gratia/blob/main/R/posterior-samples.R
    # ...in combination with "simulate_gam" from the package "mgcvUtils"
    # https://github.com/dill/mgcvUtils/blob/master/R/simulate_gam.R
    # param mod: model fitted using mgcv, from which posterior samples should be drawn
    # param n: number of draws to take from the posterior
    # returns: list of simulated values of V1_obs and V2_obs
    
    # set seed:
    set.seed(8401590)
    
    # get posterior draws of coefficients:
    betas <- post_draws(mod, n, method = 'gaussian', frequentist = FALSE, unconditional = TRUE)
    
    # get indices of necessary parameters for each virus:
    lss_idx <- attr(formula(mod), "lpi")
    lss_loc <- lss_idx[1:2]
    
    # get linear prediction matrix:
    Xp <- predict(mod, newdata = data, type = 'lpmatrix', unconditional = TRUE)
    
    # get simulated means for each virus:
    mu_1 <- Xp[, lss_loc[[1]], drop = FALSE] %*% t(betas[, lss_loc[[1]], drop = FALSE])
    mu_2 <- Xp[, lss_loc[[2]], drop = FALSE] %*% t(betas[, lss_loc[[2]], drop = FALSE])
    
    # reformat as list:
    mu_list <- vector('list', length = n)
    for (i in 1:n) {
      mu_list[[i]] <- cbind(mu_1[, i], mu_2[, i])
    }
    
    # draw from posteriors:
    sims <- lapply(mu_list, function(ix) {
      rmvn(nrow(ix), mu = ix, V = solve(crossprod(mod$family$data$R)))
    })
    
    # format and collapse:
    sims <- lapply(1:length(sims), function(ix) {
      sims[[ix]] %>%
        as_tibble() %>%
        rownames_to_column(var = 'time') %>%
        mutate(.draw = ix, .before = V1)
    }) %>%
      bind_rows()
    
    # # plot simulated "data":
    # par(mfrow = c(2, 1))
    # plot(data$V1_obs, pch = 20)
    # for (i in 1:n) {
    #   lines(sims$V1[sims$.draw == i], col = 'lightblue')
    # }
    # plot(data$V2_obs, pch = 20)
    # for (i in 1:n) {
    #   lines(sims$V2[sims$.draw == i], col = 'lightblue')
    # }
    
    # return formatted posterior samples:
    return(sims)
    
  }
  
  boot_func <- function(boot_dat, mod) {
    
    # how many draws from posterior?:
    n <- boot_dat %>% pull(.draw) %>% unique() %>% length()
    
    # split simulated data into list:
    boot_dat <- boot_dat %>%
      split(.$.draw)
    
    boot_mod <- lapply(boot_dat, function(ix) {
      gam(formula = list(V1 ~ s(week, bs = 'cc', k = 53) + s(time, k = 150), V2 ~ s(week, bs = 'cc', k = 53) + s(time, k = 100)),
          family = mvn(d = 2),
          data = ix,
          method = 'REML')
    })
    
    # pull out the covariance matrix and calculate the correlation matrix:
    cors <- lapply(boot_mod, function(ix) {
      ix <- ix$family$data$R %>%
        crossprod() %>%
        solve() %>%
        cov2cor()
      return(ix[2, 1])
    })
    cors <- cors %>% unlist() %>% unname()
    
    # return bootstrapped correlations:
    return(cors)
    
  }
  
  # GAM w/o confounding:
  
  # calculate time of year:
  data <- data %>%
    mutate(week = week(date))
  
  # run gam model
  mvn_mod <- gam(formula = list(V1_obs ~ s(week, bs = 'cc', k = 53) + s(time, k = 150), V2_obs ~ s(week, bs = 'cc', k = 53) + s(time, k = 100)),
                 family = mvn(d = 2), # multivariate normal distribution of dimension 2
                 data = data,
                 method = 'REML')
  
  # pull out the covariance matrix and then calculate the correlation matrix
  # output: 2 x 2 symmetric matrix
  corr_mat <- mvn_mod$family$data$R %>%
    crossprod() %>%
    solve() %>%
    cov2cor()
  
  # simulate "data" from fit models to perform parametric bootstrap
  sim_dat <- posterior_samples_ADAPT(mvn_mod, n = 100) %>%
    mutate(time = as.numeric(time)) %>%
    left_join(data %>% dplyr::select(time, week),
              by = 'time')
  
  # estimate confidence interval for elements of correlation matrix
  boot_corr <- boot_func(sim_dat, mvn_mod)
  
  # get results and 95% CIs
  res <- data.frame(cbind(cor = corr_mat[2, 1], CI_lower95 = quantile(boot_corr, p = 0.025), CI_upper95 = quantile(boot_corr, p = 0.975)), row.names = '')
  
  # return all results
  return(res)
  
}
