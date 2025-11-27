# ---------------------------------------------------------------------------------------------------------------------
# Code to run convergent cross-mapping (CCM)

# Useful documentation describing CCM and giving examples:
# https://ha0ye.github.io/rEDM/articles/rEDM.html
# ---------------------------------------------------------------------------------------------------------------------

# load packages
library(rEDM)
library(Kendall)
library(pracma)

ccm_func <- function(data){
  
  print(unique(data$.id))
  data <- data %>%
    dplyr::select(time, V1_obs_ln, V2_obs_ln)
  
  #---- determine Embedding dimension (i.e. the number of lags used to build up the shadow manifold) ----#
  # based on the prediction skill of the model. See rEDM vingette https://ha0ye.github.io/rEDM/articles/rEDM.html 
  
  # specify library set = how much data to fit to
  # use full data
  lib <- paste('1', nrow(data), sep = ' ')
  
  # specify pred = which data to predict on 
  pred <- paste('1', nrow(data), sep = ' ')
  
  # get E (embedding dimension - the number nearest neighbours to use for prediction) for v1 
  # EmbedDimension is a wrapper around the simplex function to get out E only  
  E_v1 <- EmbedDimension(dataFrame = data, columns = 'V1_obs_ln', target = 'V1_obs_ln',
                         lib = lib, pred = pred, maxE = 2, showPlot = FALSE)
  E_v1 <- E_v1 %>% slice_max(rho) %>% pull(E) # keep the row with max prediction skill
  
  # get E for v2
  E_v2 <- EmbedDimension(dataFrame = data, columns = "V2_obs_ln", target = "V2_obs_ln",
                         lib = lib, pred = pred, maxE = 5, showPlot = FALSE)
  E_v2 <- E_v2 %>% slice_max(rho) %>% pull(E)
  
  #---- determine if any time delay needs considering: i.e. tp parameter ----#
  # consider only negative tp (causation must be past -> future)
  vars <- c('V1_obs_ln', 'V2_obs_ln')
  
  # generate all combinations of lib_column, target_column, tp
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -20:0) %>%
    filter(lib_column != target_column) # remove cases where lib == target
  
  # want E to be that corresponding to the lib column variable (i.e. the names
  # of the input data used to create the library)
  params <- params %>%
    mutate(E = if_else(lib_column == 'V1_obs_ln', E_v1, E_v2))
  
  # set maximum possible library size
  lib_size_tp <- nrow(data) - (20 - 1) - (max(params$E) - 1) # total number of weeks of data - max tp - default tau (-1) - max embedding dimension
  
  # explore prediction skill over a range of tp values
  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = data, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]), 
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))
  
  # pull out optimal Tp
  # note: lib_column xmap target column and E is based on lib_column
  optimal_tp_v1xv2 <- output %>% filter(target_column == 'V2_obs_ln') %>% filter(`V1_obs_ln:V2_obs_ln` == max(`V1_obs_ln:V2_obs_ln`)) %>% pull(tp)
  optimal_tp_v2xv1 <- output %>% filter(target_column == 'V1_obs_ln') %>% filter(`V2_obs_ln:V1_obs_ln` == max(`V2_obs_ln:V1_obs_ln`)) %>% pull(tp)
  
  highest_crossmap_cor_v1xv2 <- output %>% filter(target_column == 'V2_obs_ln') %>% pull(`V1_obs_ln:V2_obs_ln`) %>% max()
  highest_crossmap_cor_v2xv1 <- output %>% filter(target_column == 'V1_obs_ln') %>% pull(`V2_obs_ln:V1_obs_ln`) %>% max()
  
  #----- run CCM ------#
  
  # determine max library size
  lib_max <- nrow(data) - (max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1)) - 1) - (max(abs(E_v1), abs(E_v2)) - 1)
  if (optimal_tp_v1xv2 == 0 & optimal_tp_v2xv1 == 0) lib_max <- lib_max - 1
  
  # run the ccm
  v1_xmap_v2 <- CCM(dataFrame = data, E = E_v1, columns = 'V1_obs_ln', target = 'V2_obs_ln', 
                    libSizes = c(E_v1 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v1xv2, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  v2_xmap_v1 <- CCM(dataFrame = data, E = E_v2, columns = 'V2_obs_ln', target = 'V1_obs_ln', 
                    libSizes = c(E_v2 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v2xv1, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  # pull out the mean rho for each library size
  mean_rho_v1_xmap_v2 <- v1_xmap_v2$LibMeans %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize)) # set smallest library size to some arbitrary but consistent value
  mean_rho_v2_xmap_v1 <- v2_xmap_v1$LibMeans %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  
  # combine the means for the two
  res <- mean_rho_v1_xmap_v2 %>%
    dplyr::select(1:2) %>%
    inner_join(mean_rho_v2_xmap_v1 %>%
                 dplyr::select(1:2),
               by = 'LibSize') %>%
    rename('mean_v1xv2' = 'V1_obs_ln:V2_obs_ln',
           'mean_v2xv1' = 'V2_obs_ln:V1_obs_ln')
  
  # make data long
  res_long <- res %>%
    pivot_longer(-LibSize, names_to = 'direction', values_to = 'rho') %>%
    mutate(direction = if_else(direction == 'mean_v1xv2', 'v2 -> v1', 'v1 -> v2'))
  
  # pull out all the predictions to get bootstrap CIs for the mean rho for each libsize
  all_predictions_v1 <- v1_xmap_v2$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  all_predictions_v2 <- v2_xmap_v1$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  
  # calculate median as well as lower and upper bounds (2.5, 97.5%) on rho for each libsize
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
  
  # check that rhos for max libsize greater than for min libsize:
  corrs_Lmin_v1 <- all_predictions_v1 %>% filter(LibSize == min(LibSize)) %>% pull(rho)
  corrs_Lmax_v1 <- all_predictions_v1 %>% filter(LibSize == max(LibSize)) %>% pull(rho)
  
  corrs_Lmin_v2 <- all_predictions_v2 %>% filter(LibSize == min(LibSize)) %>% pull(rho)
  corrs_Lmax_v2 <- all_predictions_v2 %>% filter(LibSize == max(LibSize)) %>% pull(rho)
  
  inverse_quantile <- function(x, y) {
    # Based on: https://github.com/cobeylab/pyembedding/blob/master/statutils.py
    
    x <- sort(x)
    
    if (x[1] == x[length(x)]) {
      print('Error: Same initial and final value for Lmin')
    } else {
      quantiles <- seq(0, 1, length.out = length(x))
    }
    
    y[y < x[1]] <- x[1]
    y[y > x[length(x)]] <- x[length(x)]
    
    inv_q_y <- interp1(x, quantiles, xi = y)
    
    return(inv_q_y)
    
  }
  
  p_conv_v1xv2 <- 1 - mean(inverse_quantile(corrs_Lmin_v1, corrs_Lmax_v1))
  p_conv_v2xv1 <- 1 - mean(inverse_quantile(corrs_Lmin_v2, corrs_Lmax_v2))
  
  #---- create the null hypothesis for comparison with our CCM output ----#
  
  num_surr <- 500 # number of surrogate datasets
  
  # Create seasonal surrogate data to describe the null hypothesis of no causal effect
  # between v1 and v2 whilst accounting for shared seasonality  
  
  # generate surrogates
  surr_v1 <- SurrogateData(data$V1_obs_ln, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  surr_v2 <- SurrogateData(data$V2_obs_ln, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  
  # get list of surrogate data to use in ccm
  surr_dat_list_v1xv2 <- vector('list', length = num_surr)
  surr_dat_list_v2xv1 <- vector('list', length = num_surr)
  
  for (i in 1:num_surr) {
    
    surr_dat_list_v1xv2[[i]] <- data %>%
      dplyr::select(time, V1_obs_ln) %>%
      mutate(V2_obs_ln = surr_v2[, i])
    
    surr_dat_list_v2xv1[[i]] <- data %>%
      dplyr::select(time, V2_obs_ln) %>%
      mutate(V1_obs_ln = surr_v1[, i])
    
  }
  
  # get largest library size used for CCM
  lib_max_use <- max(res_long$LibSize)
  
  # run ccm for surrogate data
  registerDoMC(50)
  
  surr_res <- foreach(i = 1:num_surr, .packages = c('rEDM', 'tidyverse')) %dopar% {
    
    ccm_v1xv2 <- CCM(dataFrame = surr_dat_list_v1xv2[[i]], E = E_v1, columns = 'V1_obs_ln', target = 'V2_obs_ln', 
                     libSizes = lib_max_use, Tp = optimal_tp_v1xv2, random = TRUE, sample = 100, includeData = TRUE)
    ccm_v2xv1 <- CCM(dataFrame = surr_dat_list_v2xv1[[i]], E = E_v2, columns = 'V2_obs_ln', target = 'V1_obs_ln', 
                     libSizes = lib_max_use, Tp = optimal_tp_v2xv1, random = TRUE, sample = 100, includeData = TRUE)
    
    temp_means <- ccm_v1xv2$LibMeans %>%
      dplyr::select(1:2) %>%
      rename('rho' = 'V1_obs_ln:V2_obs_ln') %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$LibMeans %>%
                  dplyr::select(1:2) %>%
                  rename('rho' = 'V2_obs_ln:V1_obs_ln') %>%
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
  
  # combine all results
  surr_res <- bind_rows(surr_res)
  
  # estimate p-value using empirical cumulative distribution
  rho_v1xv2_surr <- surr_res %>% filter(direction == 'v2 -> v1') %>% pull(rho)
  rho_v2xv1_surr <- surr_res %>% filter(direction == 'v1 -> v2') %>% pull(rho)
  
  rho_v1xv2 <- res %>% filter(direction == 'v2 -> v1', LibSize == lib_max_use) %>% pull(rho)
  rho_v2xv1 <- res %>% filter(direction == 'v1 -> v2', LibSize == lib_max_use) %>% pull(rho)
  
  # as in ha0ye tutorial:
  p_surr_v1xv2 <- (sum(rho_v1xv2 < rho_v1xv2_surr) + 1) / (length(rho_v1xv2_surr) + 1)
  p_surr_v2xv1 <- (sum(rho_v2xv1 < rho_v2xv1_surr) + 1) / (length(rho_v2xv1_surr) + 1)
  
  # write out results
  res <- res %>%
    mutate(data = 'obs')
  
  surr_res <- surr_res %>%
    mutate(data = 'surr')
  
  res <- bind_rows(res, surr_res) %>%
    mutate(E = if_else(direction == 'v2 -> v1', E_v1, E_v2),
           tp_use = if_else(direction == 'v2 -> v1', optimal_tp_v1xv2, optimal_tp_v2xv1),
           max_cmc = if_else(direction == 'v2 -> v1', highest_crossmap_cor_v1xv2, highest_crossmap_cor_v2xv1),
           p_conv = if_else(direction == 'v2 -> v1', p_conv_v1xv2, p_conv_v2xv1),
           p_surr = if_else(direction == 'v2 -> v1', p_surr_v1xv2, p_surr_v2xv1))
  
  return(res)
  
}
