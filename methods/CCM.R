###############################################################
#                 Convergent Cross Mapping       
#
# Useful documentation describing CCM and giving examples:
# https://ha0ye.github.io/rEDM/articles/rEDM.html
# 
# inputs: data = data with time, v1_obs, v2_obs 
#
# Created by: Sarah Pirikahu
# Creation date: 24 March 2023
###############################################################

# load packages
library(rEDM)
library(Kendall)
library(tidyverse)

ccm_func <- function(data){
  
  data <- data %>% dplyr::select(time, V1_obs, V2_obs)
  
  # Determining Embedding dimension (i.e. the number of lags used to build up the shadow manifold)
  # Based on the prediction skill of the model. See rEDM vingette https://ha0ye.github.io/rEDM/articles/rEDM.html 
  
  # specify library set = how much data to fit to
  # choosing to fit to half the data here
  lib_max <- round(nrow(data) / 2)
  lib <- paste0("1 ", lib_max)
  
  # specify pred = which data to predict on 
  # Using the other half of the data to predict on
  pred <- paste0(lib_max + 1, ' ', nrow(data))
  
  # Get E (Embedding dimension - the number nearest neighbours to use for prediction) for v1 
  # EmbedDimension is a wrapper around the simplex function to get out E only  
  E_v1 <- EmbedDimension(dataFrame = data, columns = 'V1_obs', target = 'V1_obs',
                         lib = lib, pred = pred, maxE = 20, showPlot = TRUE)
  E_v1 <- E_v1 %>% slice_max(rho) %>% pull(E) # keep the row with max prediction skill
  
  # Get E for v2
  E_v2 <- EmbedDimension(dataFrame = data, columns = "V2_obs", target = "V2_obs",
                         lib = lib, pred = pred, maxE = 20, showPlot = FALSE)
  E_v2 <- E_v2 %>% slice_max(rho) %>% pull(E)
  
  # determining if any time delay needs considering: i.e. tp parameter
  # data <- data %>% dplyr::select(-time)
  vars <- c('V1_obs', 'V2_obs')
  # generate all combinations of lib_column, target_column, tp
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -8:8) # ~3 months either side
  # remove cases where lib == target
  params <- params[params$lib_column != params$target_column, ]
  
  # want E to be that corresponding to the lib column variable (i.e. the 
  # names of the input data used to create the library)
  params <- params %>%
    mutate(E = if_else(lib_column == 'V1_obs', E_v1, E_v2))
  
  # explore prediction skill over range of tp values 
  # for a single library size set it to max
  lib_size_tp <- nrow(data) - (15 - 1) - (max(params$E) - 1) # total number of weeks of data - max tp - default tau (-1) - max embedding dimension
  
  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = data, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]), 
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))
  
  # pull out optimal Tp
  # note: lib_column xmap target column and E is based on lib_column
  v1xv2 <- output %>% filter(E == E_v1)
  optimal_tp_v1xv2 <- v1xv2[which.max(v1xv2$`V1_obs:V2_obs`),]$tp
  
  v2xv1 <- output %>% filter(E == E_v2)
  optimal_tp_v2xv1 <- v2xv1[which.max(v2xv1$`V2_obs:V1_obs`),]$tp
  
  p1x2 <- ggplot(output %>% filter(target_column == 'V2_obs'), aes(x = tp, y = `V1_obs:V2_obs`)) + geom_line() + theme_classic()
  p2x1 <- ggplot(output %>% filter(target_column == 'V1_obs'), aes(x = tp, y = `V2_obs:V1_obs`)) + geom_line() + theme_classic()
  
  library(gridExtra)
  grid.arrange(p1x2, p2x1, ncol = 1)
  
  #----- run CCM ------#
  
  # determine max library size
  lib_max <- nrow(data) - (max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1)) - 1) - (max(abs(E_v1), abs(E_v2)) - 1)
  
  # run the ccm
  # if wish to get CIs change random_libs = TRUE and add num_samples = xx
  v1_xmap_v2 <- ccm(data, E = E_v1, lib_column = "v1_obs", target_column = "v2_obs", 
                    lib_sizes = seq(50, lib_max, 2), tp=optimal_tp_v1xv2, random_libs=FALSE,
                    replace = TRUE, stats_only=FALSE)
  
  v2_xmap_v1 <- ccm(data, E = E_v2, lib_column = "v2_obs", target_column = "v1_obs", 
                    lib_sizes = seq(50, lib_max, 2), tp=optimal_tp_v2xv1, random_libs = FALSE,
                     replace = TRUE, stats_only=FALSE)
  
  # pull out the mean rho for each library size
  mean_rho_v1_xmap_v2 <- v1_xmap_v2$LibMeans
  mean_rho_v2_xmap_v1 <- v2_xmap_v1$LibMeans
  
  # combine the means for the two
  res <- mean_rho_v1_xmap_v2 %>%
    dplyr::select(1:2) %>%
    inner_join(mean_rho_v2_xmap_v1 %>%
                 dplyr::select(1:2),
               by = 'LibSize') %>%
    rename('mean_v1xv2' = 'V1_obs:V2_obs',
           'mean_v2xv1' = 'V2_obs:V1_obs')
  
  # make data long
  res_long <- res %>%
    pivot_longer(-LibSize, names_to = 'direction', values_to = 'rho') %>%
    mutate(direction = if_else(direction == 'mean_v1xv2', 'v2 -> v1', 'v1 -> v2'))
  
  # pull out all the predictions to get bootstrap CIs for the mean rho for each libsize
  all_predictions_v1 <- v1_xmap_v2$CCM1_PredictStat
  all_predictions_v2 <- v2_xmap_v1$CCM1_PredictStat
  
  # calculate median aswell as lower and upper bounds (2.5, 97.5%) on rho for each lib size
  # because I have now changed this to doing no subsampling rho2.5 -> rho97.5 will all be the same
  intervals_perc_v1 <- all_predictions_v1 %>%
    group_by(LibSize) %>%
    summarise(rho_median = quantile(rho, probs = 0.5),
              rho_ci_lower = quantile(rho, probs = 0.025),
              rho_ci_upper = quantile(rho, probs = 0.975)) %>%
    mutate(direction = 'v2 -> v1')
  
  intervals_perc_v2 <- all_predictions_v2 %>%
    group_by(LibSize) %>%
    summarise(rho_median = quantile(rho, probs = 0.5),
              rho_ci_lower = quantile(rho, probs = 0.025),
              rho_ci_upper = quantile(rho, probs = 0.975)) %>%
    mutate(direction = 'v1 -> v2')
  
  res <- res_long %>%
    left_join(bind_rows(intervals_perc_v1, intervals_perc_v2), by = c('LibSize', 'direction'))
  
  # check convergence using Mann Kendall - significant p-value implies convergence achieved
  MannK_v1_xmap_v2 <- res %>%
    filter(direction == 'v2 -> v1') %>%
    pull(rho) %>%
    MannKendall() %>%
    getElement('sl') %>%
    purrr::map(~ .[1])
  MannK_v2_xmap_v1 <- res %>%
    filter(direction == 'v1 -> v2') %>%
    pull(rho) %>%
    MannKendall() %>%
    getElement('sl') %>%
    purrr::map(~ .[1])
  
  # ------Creating the null hypothesis for comparison with our CCM output-----#
  
  num_surr <- 100 # number of surrogate datasets
  
  # Create seasonal surrogate data to describe the null hypothesis of no causal effect
  # between v1 and v2 whilst accounting for shared seasonality  
  
  # generate surrogates
  surr_v1 <- SurrogateData(data$V1_obs, method = "seasonal", num_surr = num_surr, T_period = 52, alpha = 20)
  surr_v2 <- SurrogateData(data$V2_obs, method = "seasonal", num_surr = num_surr, T_period = 52, alpha = 20) 
  
  # SurrogateData(data$V1_obs, method = "seasonal", num_surr = 10, T_period = 52, alpha=20) %>% matplot(type = 'l')
  # SurrogateData(data$V2_obs, method = "seasonal", num_surr = 10, T_period = 52, alpha=20) %>% matplot(type = 'l')
  
  # turn any negative surrogates into 0 - can't have a negative number of cases
  surr_v1b = apply(surr_v1, 2, function(x) {
    x[x < 0] <- 0
    x
  })
  
  surr_v2 = apply(surr_v2, 2, function(x) {
    x[x < 0] <- 0
    x
  })
  
  # get list of surrogate data to use in ccm
  surr_dat_list_v1xv2 <- vector('list', length = num_surr)
  surr_dat_list_v2xv1 <- vector('list', length = num_surr)
  
  for (i in 1:num_surr) {
    
    surr_dat_list_v1xv2[[i]] <- data %>%
      dplyr::select(time, V1_obs) %>%
      mutate(V2_obs = surr_v2[, i])
    
    surr_dat_list_v2xv1[[i]] <- data %>%
      dplyr::select(time, V2_obs) %>%
      mutate(V1_obs = surr_v1[, i])
    
  }
  
  # get largest library size used for CCM
  lib_max_use <- max(res_long$LibSize)
  
  # run ccm for surrogate data
  registerDoParallel(cl <- makeCluster(5))
  
    ccm_out_v1 <- ccm(v1_data, E = E_v1, lib_column = "v1_obs", target_column = targetCol,
                  lib_sizes = seq(50, lib_max_null, 2), tp=optimal_tp_v1xv2,
                  random_libs = FALSE, replace = TRUE)

    ccm_out_v2 <- ccm(v2_data, E = E_v2, lib_column = "v2_obs", target_column = targetCol,
                     lib_sizes = seq(50, lib_max_null, 2), tp=optimal_tp_v2xv1,
                     random_libs = FALSE, replace = TRUE)
  surr_res <- foreach(i = 1:num_surr, .packages = c('rEDM', 'tidyverse')) %dopar% {
    
    
    temp_means <- ccm_v1xv2$LibMeans %>%
      dplyr::select(1:2) %>%
      rename('rho' = 'V1_obs:V2_obs') %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$LibMeans %>%
                  dplyr::select(1:2) %>%
                  rename('rho' = 'V2_obs:V1_obs') %>%
                  mutate(direction = 'v1 -> v2'))
    
    temp_ci <- ccm_v1xv2$CCM1_PredictStat %>%
      summarise(rho_median = quantile(rho, probs = 0.5),
                rho_ci_lower = quantile(rho, probs = 0.025),
                rho_ci_upper = quantile(rho, probs = 0.975)) %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$CCM1_PredictStat %>%
                  summarise(rho_median = quantile(rho, probs = 0.5),
                            rho_ci_lower = quantile(rho, probs = 0.025),
                            rho_ci_upper = quantile(rho, probs = 0.975)) %>%
                  mutate(direction = 'v1 -> v2'))
    
    res_temp <- temp_means %>%
      left_join(temp_ci, by = 'direction')
    
    return(res_temp)
    
  }
  
  stopCluster(cl)
  
  # combine all results
  surr_res <- bind_rows(surr_res)
  
  # estimate p-value using empirical cumulative distribution
  rho_v1_x_v2_surr <- surr_res %>% filter(direction == 'v2 -> v1') %>% pull(rho)
  rho_v2_x_v1_surr <- surr_res %>% filter(direction == 'v1 -> v2') %>% pull(rho)
  
  rho_v1_x_v2 <- res %>% filter(direction == 'v2 -> v1', LibSize == lib_max_use) %>% pull(rho)
  rho_v2_x_v1 <- res %>% filter(direction == 'v1 -> v2', LibSize == lib_max_use) %>% pull(rho)
  
  p_surr_v1_x_v2 <- 1 - ecdf(rho_v1_x_v2_surr)(rho_v1_x_v2)
  p_surr_v2_x_v1 <- 1 - ecdf(rho_v2_x_v1_surr)(rho_v2_x_v1)
  
  p_surr_v1_x_v2_alt <- (sum(rho_v1_x_v2 < rho_v1_x_v2_surr) + 1) / (length(rho_v1_x_v2_surr) + 1)
  p_surr_v2_x_v1_alt <- (sum(rho_v2_x_v1 < rho_v2_x_v1_surr) + 1) / (length(rho_v2_x_v1_surr) + 1)
  
  # write out results
  res_list <- list(res,
                   surr_res,
                   c(p_surr_v1_x_v2, p_surr_v2_x_v1, p_surr_v1_x_v2_alt, p_surr_v2_x_v1_alt),
                   c(unlist(MannK_v1_xmap_v2), unlist(MannK_v2_xmap_v1)),
                   c(E_v1, E_v2, optimal_tp_v1xv2, optimal_tp_v2xv1))
  return(res_list)
  
}

