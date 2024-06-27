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
library(gridExtra)

ccm_func <- function(data){
  
  print(unique(data$.id))
  data <- data %>% dplyr::select(time, V1_obs, V2_obs)
  
  #---- determine Embedding dimension (i.e. the number of lags used to build up the shadow manifold) ----#
  # based on the prediction skill of the model. See rEDM vingette https://ha0ye.github.io/rEDM/articles/rEDM.html 
  
  # specify library set = how much data to fit to
  # use first two seasons
  lib <- '1 104'
  
  # specify pred = which data to predict on 
  pred <- paste('105', nrow(data), sep = ' ')
  
  # get E (embedding dimension - the number nearest neighbours to use for prediction) for v1 
  # EmbedDimension is a wrapper around the simplex function to get out E only  
  E_v1 <- EmbedDimension(dataFrame = data, columns = 'V1_obs', target = 'V1_obs',
                         lib = lib, pred = pred, maxE = 20, showPlot = FALSE)
  E_v1 <- E_v1 %>% slice_max(rho) %>% pull(E) # keep the row with max prediction skill
  
  # get E for v2
  E_v2 <- EmbedDimension(dataFrame = data, columns = "V2_obs", target = "V2_obs",
                         lib = lib, pred = pred, maxE = 20, showPlot = FALSE)
  E_v2 <- E_v2 %>% slice_max(rho) %>% pull(E)
  
  
  #---- check whether highest cross-map correlation is positive and for a negative lag ----#
  vars <- c('V1_obs', 'V2_obs')
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -52:52) %>%
    filter(lib_column != target_column) %>%
    mutate(E = if_else(lib_column == 'V1_obs', E_v1, E_v2))
  
  # explore prediction skill over range of tp values 
  # for a single library size set it to max
  lib_size_tp <- nrow(data) - (52 - 1) - (max(params$E) - 1) # total number of weeks of data - max tp - default tau (-1) - max embedding dimension

  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = data, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]),
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))

  optimal_tp_v1xv2_wide <- output %>% filter(target_column == 'V2_obs') %>% filter(`V1_obs:V2_obs` == max(`V1_obs:V2_obs`)) %>% pull(tp)
  optimal_tp_v2xv1_wide <- output %>% filter(target_column == 'V1_obs') %>% filter(`V2_obs:V1_obs` == max(`V2_obs:V1_obs`)) %>% pull(tp)

  # p1x2 <- ggplot(output %>% filter(target_column == 'V2_obs'), aes(x = tp, y = `V1_obs:V2_obs`)) + geom_line() + theme_classic()
  # p2x1 <- ggplot(output %>% filter(target_column == 'V1_obs'), aes(x = tp, y = `V2_obs:V1_obs`)) + geom_line() + theme_classic()
  # grid.arrange(p1x2, p2x1, ncol = 1)
  
  #---- determine if any time delay needs considering: i.e. tp parameter ----#
  # generate all combinations of lib_column, target_column, tp
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -52:0) %>% # alternatively, -20:5
    filter(lib_column != target_column) # remove cases where lib == target
  
  # want E to be that corresponding to the lib column variable (i.e. the names
  # of the input data used to create the library)
  params <- params %>%
    mutate(E = if_else(lib_column == 'V1_obs', E_v1, E_v2))
  
  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = data, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]), 
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))
  
  # pull out optimal Tp
  # note: lib_column xmap target column and E is based on lib_column
  optimal_tp_v1xv2 <- output %>% filter(target_column == 'V2_obs') %>% filter(`V1_obs:V2_obs` == max(`V1_obs:V2_obs`)) %>% pull(tp)
  optimal_tp_v2xv1 <- output %>% filter(target_column == 'V1_obs') %>% filter(`V2_obs:V1_obs` == max(`V2_obs:V1_obs`)) %>% pull(tp)
  
  
  # p1x2 <- ggplot(output %>% filter(target_column == 'V2_obs'), aes(x = tp, y = `V1_obs:V2_obs`)) + geom_line() + theme_classic()
  # p2x1 <- ggplot(output %>% filter(target_column == 'V1_obs'), aes(x = tp, y = `V2_obs:V1_obs`)) + geom_line() + theme_classic()
  # grid.arrange(p1x2, p2x1, ncol = 1)
  
  #----- run CCM ------#
  
  # determine max library size
  lib_max <- nrow(data) - (max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1)) - 1) - (max(abs(E_v1), abs(E_v2)) - 1)
  if (optimal_tp_v1xv2 == 0 & optimal_tp_v2xv1 == 0) lib_max <- lib_max - 1
  
  # run the ccm
  # if wish to get CIs change random_libs = TRUE and add num_samples = xx
  v1_xmap_v2 <- CCM(dataFrame = data, E = E_v1, columns = 'V1_obs', target = 'V2_obs', 
                    libSizes = c(seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v1xv2, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  v2_xmap_v1 <- CCM(dataFrame = data, E = E_v2, columns = 'V2_obs', target = 'V1_obs', 
                    libSizes = c(seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v2xv1, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
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
  
  #---- create the null hypothesis for comparison with our CCM output ----#
  
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
  # registerDoParallel(cl <- makeCluster(50))
  registerDoMC(50)
  
  surr_res <- foreach(i = 1:num_surr, .packages = c('rEDM', 'tidyverse')) %dopar% {
    
    ccm_v1xv2 <- CCM(dataFrame = surr_dat_list_v1xv2[[i]], E = E_v1, columns = 'V1_obs', target = 'V2_obs', 
                     libSizes = lib_max_use,#libSizes = c(seq(20, 100, by = 10), seq(125, lib_max, 25)),
                     Tp = optimal_tp_v1xv2, random = TRUE, sample = 100, includeData = TRUE)
    ccm_v2xv1 <- CCM(dataFrame = surr_dat_list_v2xv1[[i]], E = E_v2, columns = 'V2_obs', target = 'V1_obs', 
                     libSizes = lib_max_use,#c(seq(20, 100, by = 10), seq(125, lib_max, 25)),
                     Tp = optimal_tp_v2xv1, random = TRUE, sample = 100, includeData = TRUE)
    
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
  
  # stopCluster(cl)
  
  # combine all results
  surr_res <- bind_rows(surr_res)
  
  # ggplot() +
  #   geom_violin(data = surr_res, aes(x = direction, y = rho, group = direction), fill = 'gray90') +
  #   geom_point(data = res %>% filter(LibSize == max(LibSize)), aes(x = direction, y = rho), col = 'red', size = 3) +
  #   theme_classic()
  
  # estimate p-value using empirical cumulative distribution
  rho_v1xv2_surr <- surr_res %>% filter(direction == 'v2 -> v1') %>% pull(rho)
  rho_v2xv1_surr <- surr_res %>% filter(direction == 'v1 -> v2') %>% pull(rho)
  
  rho_v1xv2 <- res %>% filter(direction == 'v2 -> v1', LibSize == lib_max_use) %>% pull(rho)
  rho_v2xv1 <- res %>% filter(direction == 'v1 -> v2', LibSize == lib_max_use) %>% pull(rho)
  
  # original:
  p_surr_v1xv2 <- 1 - ecdf(rho_v1xv2_surr)(rho_v1xv2)
  p_surr_v2xv1 <- 1 - ecdf(rho_v2xv1_surr)(rho_v2xv1)
  
  # as in ha0ye tutorial:
  p_surr_v1xv2_alt <- (sum(rho_v1xv2 < rho_v1xv2_surr) + 1) / (length(rho_v1xv2_surr) + 1)
  p_surr_v2xv1_alt <- (sum(rho_v2xv1 < rho_v2xv1_surr) + 1) / (length(rho_v2xv1_surr) + 1)
  
  # write out results
  res <- res %>%
    mutate(data = 'obs') %>%
    mutate(MannK = if_else(direction == 'v2 -> v1', unlist(MannK_v1_xmap_v2), unlist(MannK_v2_xmap_v1)))
  
  surr_res <- surr_res %>%
    mutate(data = 'surr') %>%
    mutate(MannK = NA)
  
  res <- bind_rows(res, surr_res) %>%
    mutate(E = if_else(direction == 'v2 -> v1', E_v1, E_v2),
           tp_use = if_else(direction == 'v2 -> v1', optimal_tp_v1xv2, optimal_tp_v2xv1),
           tp_opt = if_else(direction == 'v2 -> v1', optimal_tp_v1xv2_wide, optimal_tp_v2xv1_wide),
           p_surr = if_else(direction == 'v2 -> v1', p_surr_v1xv2, p_surr_v2xv1),
           p_surr_alt = if_else(direction == 'v2 -> v1', p_surr_v1xv2_alt, p_surr_v2xv1_alt))
  
  return(res)
  
}

