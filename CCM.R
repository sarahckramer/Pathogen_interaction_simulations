###############################################################
#                 Convergent Cross Mapping       
#
# Useful documentation describing and giving examples:
# https://ha0ye.github.io/rEDM/articles/rEDM.html
#
# Created by: Sarah Pirikahu
# Creation date: 24 March 2023
###############################################################

# load packages
library(rEDM) 

# create dataset with just the observed cases
d_var <- d1[,c("H1_obs", "H2_obs")]

# Embedding dimension (i.e. the number of lags used to build up the shadow manifold)
# determined based on the prediction skill of the model.The rEDM vingette 
# https://ha0ye.github.io/rEDM/articles/rEDM.html uses simplex function. 
# EmbedDimension is a wrapper around the simplex function to get E only out 

# specify library set = how much data to fit too
lib_max <- dim(d1)[1]/2
lib <- paste0("1 ", lib_max)

# specify pred = which data to predict on
pred <- paste0(lib_max-1, " ", dim(d1)[1])

# Get E for H2
# ....I do wonder if this is enough data to accurately determine E? 
E_h2 <- EmbedDimension(dataFrame = d1, columns = "H2", target = "H2",
                        lib = lib, pred = pred, showPlot = TRUE)
E_h2 <- E_h2 %>% slice_max(rho) # keep the row with max prediction skill
E_h2 <- E_h2[,1] # keep just the value of E

# Get E for H1 
E_h1 <- EmbedDimension(dataFrame = d1, columns = "H1", target = "H1",
                       lib = lib, pred = pred, showPlot = TRUE)
E_h1 <- E_h1 %>% slice_max(rho) 
E_h1 <- E_h1[,1] 

# determining if any time delay needs considering: i.e. Tp parameter
vars <- names(d_var)
# generate all combinations of lib_column, target_column, tp
params <- expand.grid(lib_column = vars, target_column = vars, tp = -12:12) # 3 months either side
# remove cases where lib == target
params <- params[params$lib_column != params$target_column, ]

# want E to be that corresponding to the lib column variable (i.e. the 
# names of the input data used to create the library)
params$E <- NA
params[which(params$lib_column=="H1_obs"),]$E <- E_h1
params[which(params$lib_column=="H2_obs"),]$E <- E_h2


# explore prediction skill over range of tp values 
output <- do.call(rbind, lapply(seq_len(NROW(params)), function(i) {
  ccm(d_var, E = params$E[i], lib_sizes = 40, 
      random_libs = FALSE, lib_column = as.character(params$lib_column[i]),
      target_column = as.character(params$target_column[i]), 
      tp = params$tp[i], silent = TRUE)
}))

# pull out optimal Tp 
optimal_tp_H1xH2 <- output[which.max(output$`H1_obs:H2_obs`),]$tp
optimal_tp_H2xH1 <- output[which.max(output$`H2_obs:H1_obs`),]$tp

# make data long
plot_data <- gather(output, direction, measurement, `H2_obs:H1_obs`:`H1_obs:H2_obs`,
                    factor_key=TRUE)
# plot
ggplot(plot_data, aes(x = tp, y = measurement, color = direction)) + 
  geom_line() + theme_bw()


#----- run ccm for H2 on H1 (i.e. using H1 data) ------#
H1_xmap_H2 <- ccm(d_var, E = E_h1, lib_column = "H1_obs", target_column = "H2_obs", 
                  lib_sizes = seq(10, 50, 2), num_samples = 100, tp=optimal_tp_H1xH2,
                  random_libs = TRUE, replace = TRUE, stats_only=FALSE)
str(H1_xmap_H2, max.level = 1)
# pull out the mean rho for each library size
mean_preds <- H1_xmap_H2$LibMeans
# pull out all the predictions to get bootstrap CIs for the mean rho for each libsize
all_predictions <- H1_xmap_H2$CCM1_PredictStat
# calculate lower and upper bounds on rho for each lib size
intervals <- all_predictions %>% group_by(LibSize) %>%
  summarize(lower95_H1_x_H2 = quantile(rho, probs = 0.025), 
            upper95_H1_x_H2 = quantile(rho, probs = 0.975))
# join the CI intervals with the mean
res <- mean_preds %>% left_join(intervals)

# Creating the null hypothesis for comparison with our CCM output

# using a randomised shuffle approach 
num_surr <- 100 # number of surrogate datasets
surr_H1 <- make_surrogate_data(d_var$H1_obs,method = "random", num_surr = num_surr)
#str(surr_H1) # 100 columns for each of the surrogate datasets all with 52 rows

# run ccm for surrogate data
# data.frame to hold CCM rho values between Thrips abundance and variable
rho_surr <- data.frame(LibSize = seq(10, 50, 2))

# data.frames with time, Thrips, and 1000 surrogate climate variables for CCM()
H1_data = as.data.frame(cbind(seq(1:length(d_var$H1_obs)), d_var$H1_obs, surr_H1))
names(H1_data) = c('time', 'H1_obs', paste( 'T', as.character(seq(1,100)), sep = ''))

# Cross mapping
for (i in 1:num_surr) {
  targetCol = paste('T', i, sep = '' ) # as in H1T_data
  
  ccm_out = CCM( dataFrame = H1_data, E = E_h1, Tp = optimal_tp_H1xH2,
                 columns = 'H1_obs', target = targetCol,
                 libSizes = seq(10, 50, 2), sample = 1)
  
  col = paste('H1_obs', ':', targetCol, sep = '')
  
  rho_surr = cbind(rho_surr,ccm_out[,col])
}
names(rho_surr) <- c("LibSize", paste0("T", 1:num_surr))

# finding the lower and upper quantiles for a 95% CI
dim(rho_surr) # 21 x 101
intervals_surr <- apply(rho_surr, 1, quantile, probs = c(0.025, 0.975),  na.rm = TRUE)
intervals_surr <- t(intervals_surr)
intervals_surr <- cbind(LibSizes = seq(10, 50, 2), intervals_surr)
intervals_surr <- data.frame(intervals_surr)
names(intervals_surr) <- c("LibSizes", "lower95_H1_x_H2", "upper95_H1_x_H2")

surr_median <- apply(rho_surr, 1, median) ## mean doesn't look right when doing the plot so used median... 
# not sure if this is the correct thing, maybe I should be doing the CIs using the standard approach 1.96*SE?
res_surr <- cbind(intervals_surr, median=surr_median)

# using a seasonal surrogate approach
## need to chat to Mattheiu about this one as normally you would use a seasonal variable to create this one. 

# plot H1 on H2
ggplot(aes(x=LibSize, y=`H1_obs:H2_obs`), data=res) + geom_line() + 
geom_ribbon(aes(ymin=lower95_H1_x_H2, ymax=upper95_H1_x_H2,alpha=0.05))  + 
geom_line(aes(x=LibSizes,y=median), data=res_surr, colour = "blue") +
geom_ribbon(aes(x=LibSizes,ymin=lower95_H1_x_H2, ymax=upper95_H1_x_H2,alpha=0.05), data=res_surr, inherit.aes = FALSE, fill = "lightblue")  


#----- run ccm for H1 on H2 (i.e. using H2 data) ------#


