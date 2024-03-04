###############################################################
#                 Convergent Cross Mapping       
#
# Useful documentation describing CCM and giving examples:
# https://ha0ye.github.io/rEDM/articles/rEDM.html
# 
# inputs: data = data with time, v1_obs, v2_obs 
#
# Created by: Sarah Pirikahu
# Creation date: 1 March 2024
###############################################################

# load packages
library(rEDM) 
library(Kendall)
library(tidyverse)

ccm_func <- function(data){
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
  data <- data %>% dplyr::select(v1_obs,v2_obs)
  data$.id <- NULL
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
  #R <- 100
  R <- 1
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
  res <- mean_preds 

  # make data long
  res_long <- gather(res, direction, rho, `mean_v1_obs:v2_obs`:`mean_v2_obs:v1_obs`, factor_key=TRUE)
  res_long$direction <- gsub("mean_v1_obs:v2_obs","v2 -> v1", res_long$direction)
  res_long$direction <- gsub("mean_v2_obs:v1_obs","v1 -> v2", res_long$direction)

  #--- plotting---#
  # plot of both directions cross mapping skill together 
  # cols <- c("v2 -> v1"="grey","v1 -> v2"="lightblue")
  # fig_rho <- ggplot(aes(x=LibSize, y=`mean_v1_obs:v2_obs`), data=res) + geom_line() +
  #   geom_ribbon(aes(ymin=rho2.5_v1_xmap_v2, ymax=rho97.5_v1_xmap_v2, fill="v2 -> v1"),alpha=0.6) +
  #   geom_line(aes(x=LibSize, y=`mean_v2_obs:v1_obs`), data=res, colour="blue") +
  #   geom_ribbon(aes(ymin=rho2.5_v2_xmap_v1, ymax=rho97.5_v2_xmap_v1, fill="v1 -> v2"), alpha=0.6) +
  #   ylab("cross mapping skill (rho)") +
  #   scale_fill_manual(name="",values=cols)
  # 
  # Checking convergence using Mann Kendall - significant p-value implies converegence acheived
  MannK_v1_xmap_v2_data <- MannKendall(res$`mean_v1_obs:v2_obs`)$sl[1] # pulling out p-value only
  MannK_v2_xmap_v1_data <- MannKendall(res$`mean_v2_obs:v1_obs`)$sl[1] 
 
  #---- write out results ---# 
  
  # pull out point estimates for the max library size that can be achieved by all possible simulations (i.e. some 
  # simulations will have results for larger library sizes but for the results to be comparable across my simulation 
  # runs we need to be looking at the results for the same library size)
  theoretic_min_lib <- dim(data)[1]-12-1-10 # data size - max abs tp - tau - max abs E
  # need to make the theoretical minimum library size even as I only go up in library sizes by 2
  theoretic_min_lib <- 2*round(theoretic_min_lib/2) 
  temp_res <- res_long %>% filter(LibSize==theoretic_min_lib)
  # add seasonal surrogate test p-values and mann kendall p-values
  temp_res <- cbind(temp_res,
                    MannK_v1_xmap_v2_data = MannK_v1_xmap_v2_data, 
                    MannK_v2_xmap_v1_data = MannK_v2_xmap_v1_data, 
                    E_v1 = E_v1, E_v2 = E_v2, 
                    optimal_tp_v1xv2= optimal_tp_v1xv2, optimal_tp_v2xv1=optimal_tp_v2xv1)
  # create list 
 # res_list <- list(summary = temp_res, 
 #                  fig = fig_rho)
  return(temp_res)
}

