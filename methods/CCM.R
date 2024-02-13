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

ccm_func <- function(data, Tperiod_v1, Tperiod_v2, tot_weeks){
  # load packages
  library(rEDM) 
  library(foreach)
  library(doParallel)
  library(Kendall)
  library(tidyverse)
  
  # Determining Embedding dimension (i.e. the number of lags used to build up the shadow manifold)
  # Based on the prediction skill of the model. See rEDM vingette https://ha0ye.github.io/rEDM/articles/rEDM.html 
  
  # specify library set = how much data to fit too 
  # choosing to fit to half the data here
  lib_max <- round(dim(data)[1]/2)
  lib <- paste0("1 ", lib_max)
  
  # specify pred = which data to predict on 
  # Using the other half of the data to predict on
  pred <- paste0(lib_max-1, " ", dim(data)[1])
  
  # Get E (Embedding dimension - the number nearest neighbours to use for prediction) for v1 
  # EmbedDimension is a wrapper around the simplex function to get out E only  
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
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -12:12) # ~3 months either side
  # remove cases where lib == target
  params <- params[params$lib_column != params$target_column, ]
  
  # want E to be that corresponding to the lib column variable (i.e. the 
  # names of the input data used to create the library)
  params$E <- NA
  params[which(params$lib_column=="v1_obs"),]$E <- E_v1
  params[which(params$lib_column=="v2_obs"),]$E <- E_v2
  
  # explore prediction skill over range of tp values 
  # for a single library size set it to max
  lib_size_tp <- dim(data)[1] - 12 - 1 - max(params$E) # total number of weeks of data - max tp - default tau (-1) - max embedding dimension
  
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
  optimal_tp_v2xv1 <- v2xv1[which.max(v2xv1$`v2_obs:v1_obs`),]$tp

  #----- run CCM ------#
  
  # number of samples to do with the ccm 
  R <- 100
  # determine max library size
  lib_max <- dim(data)[1] - max(abs(optimal_tp_v1xv2),abs(optimal_tp_v2xv1)) -1 - max(abs(E_v1), abs(E_v2))
  # run the ccm
  v1_xmap_v2 <- ccm(data, E = E_v1, lib_column = "v1_obs", target_column = "v2_obs", 
                    lib_sizes = seq(50, lib_max, 2), num_samples = R, tp=optimal_tp_v1xv2, random_libs=TRUE,
                    replace = TRUE, stats_only=FALSE)
  
  v2_xmap_v1 <- ccm(data, E = E_v2, lib_column = "v2_obs", target_column = "v1_obs", 
                    lib_sizes = seq(50, lib_max, 2), num_samples = R, tp=optimal_tp_v2xv1, random_libs = TRUE,
                     replace = TRUE, stats_only=FALSE)
  
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
    summarize(rho2.5_v1_xmap_v2 = quantile(rho, probs = 0.025), 
              rho50_v1_xmap_v2 = quantile(rho, probs = 0.5),
              rho97.5_v1_xmap_v2 = quantile(rho, probs = 0.975))
  
  intervals_perc_v2 <- all_predictions_v2 %>% group_by(LibSize) %>%
    summarize(rho2.5_v2_xmap_v1 = quantile(rho, probs = 0.025), 
              rho50_v2_xmap_v1 = quantile(rho, probs = 0.5),
              rho97.5_v2_xmap_v1 = quantile(rho, probs = 0.975))
  
  # join the interval datasets together 
  interval_perc <- intervals_perc_v1 %>% left_join(intervals_perc_v2, by ="LibSize")
  
  # join the CI intervals with the mean
  res <- mean_preds %>% left_join(interval_perc)
  
  # ------Creating the null hypothesis for comparison with our CCM output-----#
  
  num_surr <- 100 # number of surrogate datasets
  
  # Create seasonal surrogate data to describe the null hypothesis of no causal effect
  # between v1 and v2 whilst accounting for shared seasonality  

  # work out the amount of variability to allow into our surrogate generation based on the 
  # standard deviation of seasonal cycles
  
  # adding seasons to the data frame
  tot_seasons <- round((tot_weeks/52) - 2)
  data$season <- c(rep(1:tot_seasons, each=52),tot_seasons+1)
  # calculating sd for each season
  season_summary <- data %>% group_by(season) %>% summarise(mean_v1 = mean(v1_obs), sd_v1 = sd(v1_obs), mean_v2 = mean(v2_obs), sd_v2 = sd(v2_obs))
  # sd can be surprisingly variable...
  alpha_v1 <- 10
  #alpha_v2 <- as.numeric(quantile(season_summary$sd_v2, probs=0.25, na.rm=T))
  alpha_v2 <- 20
    
  # generate surrogates
  surr_v1 <- make_surrogate_data(data$v1_obs, method = "seasonal", num_surr = num_surr, T_period = Tperiod_v1, alpha=alpha_v1)
  surr_v2 <- make_surrogate_data(data$v2_obs, method = "seasonal", num_surr = num_surr, T_period = Tperiod_v2, alpha=alpha_v2) 
  
  # turn any negative surrogates into 0 - can't have a negative number of cases
  surr_v1 = apply(surr_v1, 2, function(x) {
    i = which(x < 0)
    x[i] = 0
    x
  })

  surr_v2 = apply(surr_v2, 2, function(x) {
    i = which(x < 0)
    x[i] = 0
    x
  })
  
  # run ccm for surrogate data
  # estimating max lib value 
  abs_max_tp <- max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1))
  lib_max_null <- dim(surr_v1)[1] - abs_max_tp - max(E_v1,E_v2)
  
  # data.frame to hold CCM rho values
  rho_surr_v1_xmap_v2 <- data.frame(LibSize = seq(50, lib_max_null, 2))
  rho_surr_v2_xmap_v1 <- data.frame(LibSize = seq(50, lib_max_null, 2))
  
  # creating data frame with observed and surrogate data to be used with ccm 
  # note for the ccm of v1 x map v2 we want to have the v2 surrogate data
  # (cause v1 xmap v2 is determining if v2 is having a causal effect on v1)  and vice-versa
  v1_data <-  data.frame(cbind(seq(1:length(data$v1_obs)), data$v1_obs, surr_v2))
  names(v1_data) <- c('time', 'v1_obs', paste('T', as.character(seq(1,num_surr)), sep = ''))

  v2_data <- as.data.frame(cbind(seq(1:length(data$v2_obs)), data$v2_obs, surr_v1))
  names(v2_data) <- c('time', 'v2_obs', paste('T', as.character(seq(1,num_surr)), sep = ''))

  # plotting out surrogates
  # v1 xmap v2 (therefore looking at v2 surrogates)
  v1_long <- v1_data %>% gather(.,surr,value, v1_obs:T100)
  v1_long <- v1_long %>% filter(surr!="v1_obs")
  
  plot_v2_surr <- ggplot(aes(x=time,y=value, colour=surr), data=v1_long) + geom_line() + 
    scale_colour_manual(values=rep("grey",num_surr)) +  theme(legend.position="none") + 
    geom_line(aes(x=time, y=v2_obs), colour="black",data=v2_data) + ylab("v2_obs")
  
  v2_long <- v2_data %>% gather(.,surr,value, v2_obs:T100)
  v2_long <- v2_long %>% filter(surr!="v2_obs")
  
  plot_v1_surr <- ggplot(aes(x=time,y=value, colour=surr), data=v2_long) + geom_line() + 
    scale_colour_manual(values=rep("grey",num_surr)) +  theme(legend.position="none") + 
    geom_line(aes(x=time, y=v1_obs), colour="black",data=v1_data) + ylab("v1_obs")
  
  # number of samples to use in the ccm for the surrogates
  # because we are doing a large number of replicates it is ok here to reduce 
  # the number of samples 
  R <- 10
  
  # Cross mapping
  # setting up parallelism for the foreach loop
  registerDoParallel(cl <- makeCluster(100))
  #registerDoParallel(cl <- makeCluster(10))
  results_foreach <- 
    foreach (j=1:num_surr, .packages=c("rEDM","tidyverse")) %dopar% {
      
    # run ccm  
    targetCol <- paste('T', j, sep = '' ) # as in v1T_data
    ccm_out_v1 <- ccm(v1_data, E = E_v1, lib_column = "v1_obs", target_column = targetCol,
                  lib_sizes = seq(50, lib_max_null, 2), num_samples = R, tp=optimal_tp_v1xv2,
                  random_libs = TRUE, replace = TRUE)

    ccm_out_v2 <- ccm(v2_data, E = E_v2, lib_column = "v2_obs", target_column = targetCol,
                     lib_sizes = seq(50, lib_max_null, 2), num_samples = R, tp=optimal_tp_v2xv1,
                     random_libs = TRUE, replace = TRUE)

    col_v1 <- paste('v1_obs', ':', targetCol, sep = '')
    col_v2 <- paste('v2_obs', ':', targetCol, sep = '')
    
    # pulling out the quantiles for the null 
    rho_surr_v1_xmap_v2 = ccm_out_v1 %>% dplyr::select(LibSize,all_of(col_v1),E,tau,tp,nn)
    rho_surr_v2_xmap_v1 = ccm_out_v2 %>% dplyr::select(LibSize,all_of(col_v2),E,tau,tp,nn)
    
    res_temp <- list(rho_surr_v1_xmap_v2 = rho_surr_v1_xmap_v2, 
                rho_surr_v2_xmap_v1 = rho_surr_v2_xmap_v1)
    res_temp
    }
  # collapsing the list so that I have the correct column of output from each 
  # surrogate CCM run for X->Y and Y->X
   temp_output <- do.call(cbind, lapply(results_foreach, function(inner_list) {
   lapply(inner_list, function(dataframe){
   dataframe[,2]})
  }))
  
   # collapsing list down to a dataframe of the output for each surrogate
   # v1 -> v2 
   rho_surr_v1_xmap_v2 <-  do.call(cbind, temp_output[1,])
   colnames(rho_surr_v1_xmap_v2) <- paste0("T", 1:num_surr)
   # cbind the library size
   rho_surr_v1_xmap_v2 <- cbind(LibSize = seq(50, lib_max_null, 2), rho_surr_v1_xmap_v2)
   
   # v2 -> v1
   rho_surr_v2_xmap_v1 <-  do.call(cbind, temp_output[2,])
   colnames(rho_surr_v2_xmap_v1) <- paste0("T", 1:num_surr)
   # cbind the library size
   rho_surr_v2_xmap_v1 <- cbind(LibSize = seq(50, lib_max_null, 2), rho_surr_v2_xmap_v1)

  # finding the lower and upper quantiles for a 95% CI
  # v1
  intervals_surr_v1_xmap_v2 <- apply(rho_surr_v1_xmap_v2[,2:(num_surr+1)], 1, quantile, probs = c(0.025,0.5, 0.975),  na.rm = TRUE)
  intervals_surr_v1_xmap_v2 <- t(intervals_surr_v1_xmap_v2) # transpose to get intervals as columns
  intervals_surr_v1_xmap_v2 <- cbind(intervals_surr_v1_xmap_v2,mean_rho=apply(rho_surr_v1_xmap_v2[,2:(num_surr+1)], 1, mean, na.rm = TRUE))
  intervals_surr_v1_xmap_v2 <- cbind(intervals_surr_v1_xmap_v2,sd_rho=apply(rho_surr_v1_xmap_v2[,2:(num_surr+1)], 1, sd, na.rm = TRUE))
  intervals_surr_v1_xmap_v2 <- cbind(intervals_surr_v1_xmap_v2,mean_rho=apply(rho_surr_v1_xmap_v2[,2:(num_surr+1)], 1, min, na.rm = TRUE))
  intervals_surr_v1_xmap_v2 <- cbind(intervals_surr_v1_xmap_v2,mean_rho=apply(rho_surr_v1_xmap_v2[,2:(num_surr+1)], 1, max, na.rm = TRUE))
  intervals_surr_v1_xmap_v2 <- cbind(LibSizes = seq(50, lib_max_null, 2), intervals_surr_v1_xmap_v2) # add libSizes to df
  intervals_surr_v1_xmap_v2 <- data.frame(intervals_surr_v1_xmap_v2)
  names(intervals_surr_v1_xmap_v2) <- c("LibSizes", "rho2.5_v1_xmap_v2", "rho50_v1_xmap_v2","rho97.5_v1_xmap_v2", "mean_rho", "sd_rho", "min_rho", "max_rho")

  # v2
  intervals_surr_v2_xmap_v1 <- apply(rho_surr_v2_xmap_v1[,2:(num_surr+1)], 1, quantile, probs = c(0.025,0.5, 0.975),  na.rm = TRUE)
  intervals_surr_v2_xmap_v1 <- t(intervals_surr_v2_xmap_v1) # transpose to get intervals as columns
  intervals_surr_v2_xmap_v1 <- cbind(intervals_surr_v2_xmap_v1,mean_rho=apply(rho_surr_v2_xmap_v1[,2:(num_surr+1)], 1, mean, na.rm = TRUE))
  intervals_surr_v2_xmap_v1 <- cbind(intervals_surr_v2_xmap_v1,mean_rho=apply(rho_surr_v2_xmap_v1[,2:(num_surr+1)], 1, min, na.rm = TRUE))
  intervals_surr_v2_xmap_v1 <- cbind(intervals_surr_v2_xmap_v1,mean_rho=apply(rho_surr_v2_xmap_v1[,2:(num_surr+1)], 1, max, na.rm = TRUE))
  intervals_surr_v2_xmap_v1 <- cbind(LibSizes = seq(50, lib_max_null, 2), intervals_surr_v2_xmap_v1) # add libSizes to df
  intervals_surr_v2_xmap_v1 <- data.frame(intervals_surr_v2_xmap_v1)
  names(intervals_surr_v2_xmap_v1) <- c("LibSizes", "rho2.5_v2_xmap_v1", "rho50_v2_xmap_v1","rho97.5_v2_xmap_v1", "mean_rho", "min_rho", "max_rho")
    
  #--- plotting---#
  # v1 xmap v2 - mean
  p_v1_xmap_v2_mean <- ggplot(aes(x=LibSize, y=`mean_v1_obs:v2_obs`), data=res) + geom_line() +
  geom_ribbon(aes(ymin=rho2.5_v1_xmap_v2, ymax=rho97.5_v1_xmap_v2,alpha=0.05))  +
  geom_line(aes(x=LibSizes,y=mean_rho), data=intervals_surr_v1_xmap_v2, colour = "blue") +
  geom_ribbon(aes(x=LibSizes,ymin=rho2.5_v1_xmap_v2, ymax=rho97.5_v1_xmap_v2,alpha=0.05), data=intervals_surr_v1_xmap_v2, inherit.aes = FALSE, fill = "lightblue")
  
  # v2 xmap v1 - mean
  p_v2_xmap_v1_mean <- ggplot(aes(x=LibSize, y=`mean_v2_obs:v1_obs`), data=res) + geom_line() +
  geom_ribbon(aes(ymin=rho2.5_v2_xmap_v1, ymax=rho97.5_v2_xmap_v1,alpha=0.05))  +
  geom_line(aes(x=LibSizes,y=mean_rho), data=intervals_surr_v2_xmap_v1, colour = "blue") +
  geom_ribbon(aes(x=LibSizes,ymin=rho2.5_v2_xmap_v1, ymax=rho97.5_v2_xmap_v1,alpha=0.05), data=intervals_surr_v2_xmap_v1, inherit.aes = FALSE, fill = "lightblue")
  
  
  # estimating p-value using empirical cumulative distribution 
  p_surr_v1_xmap_v2 <- 1 - ecdf(as.numeric(rho_surr_v1_xmap_v2[nrow(res),-1]))(res[dim(res)[1],"mean_v1_obs:v2_obs"])
  p_surr_v2_xmap_v1 <- 1 - ecdf(as.numeric(rho_surr_v2_xmap_v1[nrow(res),-1]))(res[dim(res)[1],"mean_v2_obs:v1_obs"])
   
  # Checking convergence using Mann Kendall - significant p-value implies converegence acheived
  MannK_v1_xmap_v2_data <- MannKendall(res$`mean_v1_obs:v2_obs`)$sl[1] # pulling out p-value only
  MannK_v1_xmap_v2_null <- MannKendall(intervals_surr_v1_xmap_v2$mean_rho)$sl[1]
  
  MannK_v2_xmap_v1_data <- MannKendall(res$`mean_v2_obs:v1_obs`)$sl[1] 
  MannK_v2_xmap_v1_null <- MannKendall(intervals_surr_v2_xmap_v1$mean_rho)$sl[1]
  
  #---- write out results ---# 
  
  # pull out point estimates for the max library size that can be achieved by all possible simulations (i.e. some 
  # simulations will have results for larger library sizes but for the results to be comparable across my simulation 
  # runs we need to be looking at the results for the same library size)
  theoretic_min_lib <- dim(data)[1]-12-1-10 # data size - max abs tp - tau - max abs E
  temp_res <- res %>% filter(LibSize==theoretic_min_lib)
  # add seasonal surrogate test p-values and mann kendall p-values
  temp_res <- cbind(temp_res,ecdf_p_v1_x_v2 = p_surr_v1_xmap_v2, ecdf_p_v2_x_v1 = p_surr_v2_xmap_v1,
                    MannK_v1_xmap_v2_data = MannK_v1_xmap_v2_data, MannK_v1_xmap_v2_null = MannK_v1_xmap_v2_null,
                    MannK_v2_xmap_v1_data = MannK_v2_xmap_v1_data, MannK_v2_xmap_v1_null = MannK_v2_xmap_v1_null,
                    E_v1 = E_v1, E_v2 = E_v2, optimal_tp_v1xv2= optimal_tp_v1xv2, optimal_tp_v2xv1=optimal_tp_v2xv1)
  # create list 
  res_list <- list(summary = temp_res, 
                   fig_v1_xmap_v2_mean = p_v1_xmap_v2_mean,
                   fig_v2_xmap_v1_mean = p_v2_xmap_v1_mean,
                   plot_v1_surr = plot_v1_surr,
                   plot_v2_surr = plot_v2_surr)
  return(res_list)
}

