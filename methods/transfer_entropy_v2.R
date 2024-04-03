###############################################################
#                 Transfer Entropy       
#
# Created by: Sarah Pirikahu
# Creation date: 28 Feb 2024
###############################################################

# packages
library(RTransferEntropy) 
library(vars) 

te_func <- function(data){

  v1_obs <- data$v1_obs
  v2_obs <- data$v2_obs
  
  # initialising lists to put results in 
  lags <- list()
  lag_v1 <- list()
  lag_v2 <- list()
  
  # determine the number of lags for each simulated dataset
  df <- data.frame(v1_obs=v1_obs, v2_obs=v2_obs)
  
  lags <- lapply(df, VARselect, lag.max=12) # lag of approx 3 month
  rm(df)
  # pull out the lag with best BIC. Lower BIC = better (not BIC is labeled SV)
  # regardless of whether raw of normalised data used the lag chosen is the same
  lag_v1 <- as.numeric(lags$v1_obs$selection[3])
  lag_v2 <- as.numeric(lags$v2_obs$selection[3])
  
  rm(lags)
  # Interpreting transfer entropy (note: TE \in [0,1]):
  # If test significant suggests T_{X->Y} > 0 and the uncertainty about 
  # Y is reduced by the addition of X, that is X causes Y.
  
  # Output: provides not the transfer entropy and bias corrected effective transfer entropy  
  # Transfer entropy estimates are biased by small sample sizes. For large sample sizes TE and ETE 
  # will be approximately the same. For a single season the sample size is quite small so we want to 
  # go with ETE... see Behrendt et al. 2019 for more details
  shannon_te <- transfer_entropy(v1_obs, v2_obs, lx = min(lag_v1, lag_v2), ly=min(lag_v1, lag_v2), bins=10)
  temp_res <- data.frame(coef(shannon_te))
  temp_res <- rownames_to_column(temp_res, "direction")
  # adding lag to output
  temp_res$lag <- min(lag_v1, lag_v2)
  # providing better naming for direction
  temp_res$direction <- gsub("X->Y","v1 -> v2", temp_res$direction)
  temp_res$direction <- gsub("Y->X","v2 -> v1", temp_res$direction)
  
  
  return(temp_res)
}

