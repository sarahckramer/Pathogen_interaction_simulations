################################################################################
#                     Transfer Entropy analysis        
#
# Useful documentation describing and giving examples:
# https://cran.rstudio.com/web/packages/RTransferEntropy/vignettes/transfer-entropy.html
# and useful video to describe transfer entropy in a simple way:
# https://www.youtube.com/watch?v=focIC0v5Rds
# 
# Transfer entropy has the disadvantage that it does not tell us about the sign 
# of the interaction 
# 
# Created by: Sarah Pirikahu
# Creation date: 19 April 2023
################################################################################

# load packages
library(RTransferEntropy) 
library(vars)
library(future) # allows for parallel processing


# pull out the lag with best AIC. Lower AIC = better 
# regardless of whether raw of normalised data used the lag choosen is the same
lag <- min(lag_h1, lag_h2)

# determining the transfer entropy (note: TE \in [-1,1]. TE = 0 means)
shannon_te <- transfer_entropy(d_var$v1_obs, d_var$v2_obs, lx=lag, ly=lag)
shannon_te
coef(shannon_te)




transfer_entropy(d_var$v1_obs, d_var$v2_obs, nboot = 100, seed = 123,
                 quantiles=c(2.5,97.5))

# pulling out coefficents 
res <- data.frame(shannon_te$coef)


