###############################################################
#                 Convergent Cross Mapping       
#
# Useful documentation describing and giving examples:
# https://ha0ye.github.io/rEDM/articles/rEDM.html
#
# Created by: Sarah Pirikahu
# Creation date: 24 March 2023
###############################################################

ccm_func <- function(data){
  # load packages
  library(rEDM) 
  
  # Determining Embedding dimension (i.e. the number of lags used to build up the shadow manifold)
  # Based on the prediction skill of the model. See rEDM vingette https://ha0ye.github.io/rEDM/articles/rEDM.html 
  
  # specify library set = how much data to fit too 
  # choosing to fit to half the data here
  lib_max <- dim(data)[1]/2
  lib <- paste0("1 ", lib_max)
  
  # specify pred = which data to predict on 
  # Using the other half of the data to predict on
  pred <- paste0(lib_max-1, " ", dim(data)[1])
  
  # Get E for v1 
  # EmbedDimension is a wrapper around the simplex function to get E only out 
  # ....I do wonder if this is enough data to accurately determine E? 
  E_v1 <- EmbedDimension(dataFrame = data, columns = "v1_obs", target = "v1_obs",
                         lib = lib, pred = pred, showPlot = TRUE)
  E_v1 <- E_v1 %>% slice_max(rho) 
  E_v1 <- E_v1[,1] 
  
  # Get E for v2
  E_v2 <- EmbedDimension(dataFrame = data, columns = "v2_obs", target = "v2_obs",
                          lib = lib, pred = pred, showPlot = TRUE)
  E_v2 <- E_v2 %>% slice_max(rho) # keep the row with max prediction skill
  E_v2 <- E_v2[,1] # keep just the value of E
  
  # determining if any time delay needs considering: i.e. tp parameter
  data <- data %>% dplyr::select(-time)
  vars <- names(data)
  # generate all combinations of lib_column, target_column, tp
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -5:5) # ~1 months either side
  # remove cases where lib == target
  params <- params[params$lib_column != params$target_column, ]
  
  # want E to be that corresponding to the lib column variable (i.e. the 
  # names of the input data used to create the library)
  params$E <- NA
  params[which(params$lib_column=="v1_obs"),]$E <- E_v1
  params[which(params$lib_column=="v2_obs"),]$E <- E_v2
  
  # explore prediction skill over range of tp values 
  # for a single library size set it to max
  lib_size_tp <- dim(data)[1] - 5 - 1 - max(params$E) # total number of weeks of data - max tp - default tau (-1) - max embedding dimension
  
  output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
    ccm(data, E = params$E[i], lib_sizes = lib_size_tp, 
        random_libs = FALSE, lib_column = as.character(params$lib_column[i]),
        target_column = as.character(params$target_column[i]), 
        tp = params$tp[i], silent = TRUE)
  }))
  
  # pull out optimal Tp 
  # note: lib_column xmap target column and E is based on lib_column
  v1xv2 <- output %>% filter(E==E_v1) 
  optimal_tp_v1xv2 <- v1xv2[which.max(v1xv2$`v1_obs:v2_obs`),]$tp
   
  v2xv1 <- output %>% filter(E==E_v2) 
  optimal_tp_v2xv1 <- v1xv2[which.max(v1xv2$`v2_obs:v1_obs`),]$tp
  

  #----- run CCM ------#
  
  # number of samples to do with the ccm 
  R <- 100
  # run the ccm 
  lib_max <- dim(data)[1] - max(abs(optimal_tp_v1xv2),abs(optimal_tp_v2xv1)) -1 - max(abs(E_v1), abs(E_v2))
  v1_xmap_v2 <- ccm(data, E = E_v1, lib_column = "v1_obs", target_column = "v2_obs", 
                    lib_sizes = seq(10, lib_max, 2), num_samples = R, tp=optimal_tp_v1xv2,
                    random_libs = TRUE, replace = TRUE, stats_only=FALSE)
  
  v2_xmap_v1 <- ccm(data, E = E_v2, lib_column = "v1_obs", target_column = "v2_obs", 
                    lib_sizes = seq(10, lib_max, 2), num_samples = R, tp=optimal_tp_v2xv1,
                    random_libs = TRUE, replace = TRUE, stats_only=FALSE)
  
  # pull out the mean rho for each library size
  mean_rho_v1_xmap_v2 <- v1_xmap_v2$LibMeans 
  mean_rho_v2_xmap_v1 <- v2_xmap_v1$LibMeans
  # combine the means for the two 
  mean_preds <- data.frame(cbind(mean_rho_v1_xmap_v2$LibSize, mean_rho_v1_xmap_v2$`v1_obs:v2_obs`, mean_rho_v2_xmap_v1$`v2_obs:v1_obs`))
  names(mean_preds) <- c("LibSize", "mean_v1_obs:v2_obs", "mean_v2_obs:v1_obs")  
  
  # pull out all the predictions to get bootstrap CIs for the mean rho for each libsize
  all_predictions_v1 <- v1_xmap_v2$CCM1_PredictStat
  all_predictions_v2 <- v2_xmap_v1$CCM1_PredictStat
  
  # calculate median aswell as lower and upper bounds (2.5, 97.5%) on rho for each lib size
  intervals_perc_v1 <- all_predictions_v1 %>% group_by(LibSize) %>%
    summarize(lower95_v1_xmap_v2 = quantile(rho, probs = 0.025), 
              median_v1_xmap_v2 = quantile(rho, probs = 0.5),
              upper95_v1_xmap_v2 = quantile(rho, probs = 0.975))
  
  intervals_perc_v2 <- all_predictions_v2 %>% group_by(LibSize) %>%
    summarize(lower95_v2_xmap_v1 = quantile(rho, probs = 0.025), 
              median_v2_xmap_v1 = quantile(rho, probs = 0.5),
              upper95_v2_xmap_v1 = quantile(rho, probs = 0.975))
  
  # join the interval datasets together 
  interval_perc <- intervals_perc_v1 %>% left_join(intervals_perc_v2, by ="LibSize")
  
  # join the CI intervals with the mean
  res <- mean_preds %>% left_join(interval_perc)
  
  # Creating the null hypothesis for comparison with our CCM output
  # using a randomised shuffle approach 
  num_surr <- 100 # number of surrogate datasets
  surr_v1 <- make_surrogate_data(data$v1_obs, method = "random", num_surr = num_surr)
  surr_v2 <- make_surrogate_data(data$v2_obs, method = "random", num_surr = num_surr)
  #str(surr_v1) # 100 columns for each of the surrogate datasets all with 52 rows
  
  # run ccm for surrogate data
  # estimating max lib value 
  abs_max_tp <- max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1))
  lib_max_null <- dim(surr_v1)[1] - abs_max_tp - max(E_v1,E_v2)
  
  # data.frame to hold CCM rho values
  rho_surr_v1_xmap_v2 <- data.frame(LibSize = seq(10, lib_max_null, 2)) 
  rho_surr_v2_xmap_v1 <- data.frame(LibSize = seq(10, lib_max_null, 2)) 
  
  # creating data frame with observed and surrogate data to be used with ccm 
  v1_data <-  data.frame(cbind(seq(1:length(data$v1_obs)), data$v1_obs, surr_v1))
  names(v1_data) <- c('time', 'v1_obs', paste('T', as.character(seq(1,100)), sep = ''))
  
  v2_data <- as.data.frame(cbind(seq(1:length(data$v2_obs)), data$v2_obs, surr_v2))
  names(v2_data) <- c('time', 'v2_obs', paste('T', as.character(seq(1,100)), sep = ''))
  
  # Cross mapping
  for (j in 1:num_surr) {
    targetCol <- paste('T', j, sep = '' ) # as in v1T_data
    ccm_out_v1 <- ccm(v1_data, E = E_v1, lib_column = "v1_obs", target_column = targetCol,  
                  lib_sizes = seq(10, lib_max_null, 2), num_samples = R, tp=optimal_tp_v1xv2,
                  random_libs = TRUE, replace = TRUE)
    
    ccm_out_v2 <- ccm(v2_data, E = E_v2, lib_column = "v2_obs", target_column = targetCol,  
                     lib_sizes = seq(10, lib_max_null, 2), num_samples = R, tp=optimal_tp_v2xv1,
                     random_libs = TRUE, replace = TRUE)
    
    col_v1 <- paste('v1_obs', ':', targetCol, sep = '')
    col_v2 <- paste('v2_obs', ':', targetCol, sep = '')
    # pulling out the quantiles for the null 
    
    rho_surr_v1_xmap_v2 = cbind(rho_surr_v1_xmap_v2,ccm_out_v1[,col_v1])
    rho_surr_v2_xmap_v1 = cbind(rho_surr_v2_xmap_v1,ccm_out_v2[,col_v2])
  }
  names(rho_surr_v1_xmap_v2) <- c("LibSize", paste0("T", 1:num_surr))
  names(rho_surr_v2_xmap_v1) <- c("LibSize", paste0("T", 1:num_surr))
  
  # finding the lower and upper quantiles for a 95% CI
  # v1
  dim(rho_surr_v1_xmap_v2) # 21 x 101
  intervals_surr_v1_xmap_v2 <- apply(rho_surr_v1_xmap_v2, 1, quantile, probs = c(0.025,0.5, 0.975),  na.rm = TRUE)
  intervals_surr_v1_xmap_v2 <- t(intervals_surr_v1_xmap_v2) # transpose to get intervals as columns
  intervals_surr_v1_xmap_v2 <- cbind(LibSizes = seq(10, lib_max_null, 2), intervals_surr_v1_xmap_v2) # add libSizes to df
  intervals_surr_v1_xmap_v2 <- data.frame(intervals_surr_v1_xmap_v2)
  names(intervals_surr_v1_xmap_v2) <- c("LibSizes", "lower95_v1_xmap_v2", "median","upper95_v1_xmap_v2")
  
  # v2
  intervals_surr_v2_xmap_v1 <- apply(rho_surr_v2_xmap_v1, 1, quantile, probs = c(0.025,0.5, 0.975),  na.rm = TRUE)
  intervals_surr_v2_xmap_v1 <- t(intervals_surr_v2_xmap_v1) # transpose to get intervals as columns
  intervals_surr_v2_xmap_v1 <- cbind(LibSizes = seq(10, lib_max_null, 2), intervals_surr_v2_xmap_v1) # add libSizes to df
  intervals_surr_v2_xmap_v1 <- data.frame(intervals_surr_v2_xmap_v1)
  names(intervals_surr_v2_xmap_v1) <- c("LibSizes", "lower95_v2_xmap_v1", "median","upper95_v2_xmap_v1")
  
  
  #--- plotting---#
  # v1 xmap v2
  p_v1_xmap_v2 <- ggplot(aes(x=LibSize, y=`median_v1_xmap_v2`), data=res) + geom_line() + 
  geom_ribbon(aes(ymin=lower95_v1_xmap_v2, ymax=upper95_v1_xmap_v2,alpha=0.05))  + 
  geom_line(aes(x=LibSizes,y=median), data=intervals_surr_v1_xmap_v2, colour = "blue") +
  geom_ribbon(aes(x=LibSizes,ymin=lower95_v1_xmap_v2, ymax=upper95_v1_xmap_v2,alpha=0.05), data=intervals_surr_v1_xmap_v2, inherit.aes = FALSE, fill = "lightblue")  
  
  # v2 xmap v1
  p_v2_xmap_v1 <- ggplot(aes(x=LibSize, y=`median_v2_xmap_v1`), data=res) + geom_line() + 
    geom_ribbon(aes(ymin=lower95_v2_xmap_v1, ymax=upper95_v2_xmap_v1,alpha=0.05))  + 
    geom_line(aes(x=LibSizes,y=median), data=intervals_surr_v2_xmap_v1, colour = "blue") +
    geom_ribbon(aes(x=LibSizes,ymin=lower95_v2_xmap_v1, ymax=upper95_v2_xmap_v1,alpha=0.05), data=intervals_surr_v2_xmap_v1, inherit.aes = FALSE, fill = "lightblue")  
  
  
  # implementing two sample Kolmogorov-Smirnov test 
  # this test determines if the two series have come from the same distribution
  # H0: both samples come from the same distribution 
  # H1: the samples come from different distributions
  
  # if p-value is significant this suggests a causal effect as the null
  # and model cross mapping skill are significantly different
  ks_out_v1_x_v2 <- ks.test(res$median_v1_xmap_v2, intervals_surr_v1_xmap_v2$median, exact=TRUE)
  ks_out_v2_x_v1 <- ks.test(res$median_v2_xmap_v1, intervals_surr_v2_xmap_v1$median, exact=TRUE)
  
  #---- write out results ---# 
  
  # our point estimate is simply going to be the cross mapping skill for the largest library size as we have
  # based max library size on our sample size
  temp_res <- res[dim(res)[1],]
  # add ks p_values
  temp_res <- cbind(temp_res, ks_p_v1_x_v2=ks_out_v1_x_v2$p.value, ks_p_v2_x_v1=ks_out_v2_x_v1$p.value)
  # create list 
  overall_res_list <- list(summary=temp_res, v1_xmap_v2_preds=all_predictions_v1,
                           v2_xmap_v1_preds = all_predictions_v2, 
                           v1_xmap_v2_null = intervals_surr_v1_xmap_v2,
                           v2_xmap_v1_null = intervals_surr_v2_xmap_v1)
  res_list <- list(summary = overall_res_list, 
                   fig_v1_xmap_v2 = p_v1_xmap_v2, fig_v2_xmap_v1 = p_v1_xmap_v2)
  return(res_list)
}

