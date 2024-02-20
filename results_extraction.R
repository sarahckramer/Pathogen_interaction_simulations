################################################################################
#                       Results extraction 
#
# In this script we extract the results of the simulation study
# excluding any ad hoc likelihood computation
# 
# Created by: Sarah Pirikahu
# Creation date: 20 Oct 2023
################################################################################

# load packages
library(tidyverse)
library(openxlsx) 
library(gridExtra)

#------- all methods except likelihood extraction ------# 

# reading in all results 
file_names = list.files(pattern = "*.RData", recursive = F)

# load all the data sets in 
for (i in 1:length(file_names)) {
  load(file_names[i], verbose = TRUE)
  assign(paste0("results", i), results)
  rm(results)
}

# naming the true parameters we want to extract into our results tables
params <- c("delta1", "delta2", "theta_lambda1", "theta_lambda2", "Ri1", "Ri2",
            "sigma1", "sigma2", "gamma1", "gamma2", "w1", "w2", "rho1", "rho2",
            "A1", "phi1", "A2", "phi2", "beta_sd1", "beta_sd2", "E01", "E02",
            "R01", "R02", "R12")

# extract true parameters from each list and put into a data frame
list_names <- paste0("results", 1:length(file_names))
first_elements <- purrr::map(list_names, ~ get(.x)$true_param)
results_df <- bind_rows(first_elements)
results_df <- results_df[,params]

#---- extract Spearmann correlation coefficients ----#
cor_res <- purrr::map(list_names, ~ get(.x)$cor)
cor_res_df <- bind_rows(cor_res)
names(cor_res_df) <- c("cor", "cor_CI_2.5", "cor_CI_97.5", "cor_pvalue")
results_df <- cbind(results_df,cor_res_df)

#---- extract GAM correlations ----# 
gam_res <- purrr::map(list_names, ~ get(.x)$gam_cor)
gam_res_df <- data.frame(bind_rows(gam_res), row.names = NULL)
gam_res_df <- gam_res_df %>% dplyr::select(cor, CI_lower95, CI_upper95)
names(gam_res_df) <- c("gam_cor", "gam_CI_2.5", "gam_CI_97.5")
results_df <- cbind(results_df,gam_res_df)

#---- extract transfer entropy ----# 
te_res_v1_x_v2 <- purrr::map(list_names, ~ get(.x)$transfer_entropy[1,])
te_res_v2_x_v1 <- purrr::map(list_names, ~ get(.x)$transfer_entropy[2,])

te_res_v1_x_v2_df <- data.frame(bind_rows(te_res_v1_x_v2), row.names = NULL)
te_res_v2_x_v1_df <- data.frame(bind_rows(te_res_v2_x_v1), row.names = NULL)

# renaming and just keeping important columns 
te_res_v1_x_v2_df <- te_res_v1_x_v2_df %>% dplyr::select(-se)
te_res_v2_x_v1_df <- te_res_v2_x_v1_df %>% dplyr::select(-se)

names(te_res_v1_x_v2_df) <- c("TE_v1_x_v2_te","TE_v1_x_v2_ete","TE_v1_x_v2_pvalue","TE_v1_x_v2_CI2.5", "TE_v1_x_v2_CI97.5")
names(te_res_v2_x_v1_df) <- c("TE_v2_x_v1_te","TE_v2_x_v1_ete","TE_v2_x_v1_pvalue","TE_v2_x_v1_CI2.5", "TE_v2_x_v1_CI97.5")

## Note:: v1 x v2 here means v1 -> v2 unlike in CCM how it is backwards
results_df <- cbind(results_df, te_res_v1_x_v2_df,te_res_v2_x_v1_df)

#---- extract granger ----# 
granger_res_v1_xmap_v2 <- purrr::map(list_names, ~ get(.x)$granger$summary[1,])
granger_res_v2_xmap_v1 <- purrr::map(list_names, ~ get(.x)$granger$summary[2,])

granger_res_v1_xmap_v2 <- data.frame(bind_rows(granger_res_v1_xmap_v2), row.names = NULL)
names(granger_res_v1_xmap_v2) <- c("granger_est_v1_x_v2", "granger_bkbootmean_v1_x_v2",
                                   "granger_blockboot_CI_2.5_v1_x_v2","granger_blockboot_CI_97.5_v1_x_v2",
                                   "granger_blockboot_percCI_2.5_v1_x_v2","granger_blockboot_percCI_97.5_v1_x_v2",
                                   "granger_p_v1_x_v2", "granger_adf_p_v1_x_v2", "granger_kpss_p_v1_x_v2")
granger_res_v2_xmap_v1 <- data.frame(bind_rows(granger_res_v2_xmap_v1), row.names = NULL)
names(granger_res_v2_xmap_v1) <- c("granger_est_v2_x_v1", "granger_bkbootmean_v2_x_v1",
                                   "granger_blockboot_CI_2.5_v2_x_v1","granger_blockboot_CI_97.5_v2_x_v1",
                                   "granger_blockboot_percCI_2.5_v2_x_v1","granger_blockboot_percCI_97.5_v2_x_v1",
                                   "granger_p_v2_x_v1", "granger_adf_p_v2_x_v1", "granger_kpss_p_v2_x_v1")

results_df <- cbind(results_df, granger_res_v1_xmap_v2, granger_res_v2_xmap_v1)

#------- extract convergent cross mapping ----# 
ccm_res <- purrr::map(list_names, ~ get(.x)$CCM$summary[1,])

ccm_res <- data.frame(bind_rows(ccm_res), row.names = NULL)
names(ccm_res) <- c("ccm_libsize", "ccm_mean_rho_v1_x_v2", "ccm_mean_rho_v2_x_v1", 
                    "ccm_rho2.5_v1_x_v2", "ccm_rho50_v1_x_v2", "ccm_rho97.5_v1_x_v2",
                    "ccm_rho2.5_v2_x_v1", "ccm_rho50_v2_x_v1", "ccm_rho97.5_v2_x_v1",
                    "ccm_p_v1_x_v2", "ccm_p_v2_x_v1", 
                    "ccm_mannken_v1_x_v2_data","ccm_mannken_v1_x_v2_null",
                    "ccm_mannken_v2_x_v1_data","ccm_mannken_v2_x_v1_null",
                    "ccm_E_v1", "ccm_E_v2", "ccm_tp_v1_x_v2", "ccm_tp_v2_x_v1")


results_df <- cbind(results_df, ccm_res)

# order results by delta
results_df <- results_df[order(results_df$delta1, decreasing=T),]
# output as csv
#write.xlsx(results_df, file="results_symmetric_5yrs.xlsx")

##################################################
#-------- Extract cross map skill plots----------#
##################################################

for(i in 1:15){
  # pulling out jobid name since results name differs from jobid due to ordering 
  jobid <- gsub("results_(\\d+)_.*\\.RData", "\\1", file_names)
  # getting file name sorted 
  res_names_v1xv2 <- paste0("jobid", jobid, "_v1xv2_mean.png")
  p1_v1xv2_name <- res_names_v1xv2[i]
  
  # printing v2 -> v1
  png(p1_v1xv2_name, width=800, height=400)
  p1 <- get(paste0("results", i))$CCM$fig_v1_xmap_v2_mean
  print(p1)
  dev.off()
  
  # printing v1 -> v2
  res_names_v2xv1 <- paste0("jobid", jobid, "_v2xv1_mean.png")
  p2_v2xv1_name <- res_names_v2xv1[i]
  png(p2_v2xv1_name, width=800, height=400)
  p2 <- get(paste0("results", i))$CCM$fig_v2_xmap_v1_mean
  print(p2)
  dev.off()
}


##################################################
#------- plotting out the observed data ---------#
##################################################

# getting the timing of the influenza surges
set.seed(2908)
# total number of seasons
tot_seasons <- round((tot_weeks/52) - 2)

# initialize time of surges (based on week) from start of season (1 July)
# by drawing from a normal distribution 
n_surge <- round(tot_weeks/52) - 1 # total number of surges
mu_Imloss <- 38 # average surge occuring in mid Oct
sd_Imloss <- 4 # allow variation of 4 weeks

t_si <- rnorm(n=n_surge, mean=mu_Imloss,sd=sd_Imloss)
# correcting t_si to give t based on week of the year rather than 
# week from start of season (note: July 1 is week 26)
t_si <- round(seq(26, 52*n_surge, by=52) + t_si)

# remove all t_si which are less than 2years - allowing this amount
# of time for the system to reach an equilibrium (2 yrs ~ 104 weeks)
t_si <- t_si[-which(t_si <= 104)]
# changing the surge times to dates
t_si_date <- lubridate::ymd("2012-July-01") + lubridate::weeks(t_si)

# parameter inputs:
#theta_lambda1 <- c(0,0.5,1,2,4)
#theta_lambda2 <- c(0,0.5,1,2,4)
theta_lambda1 <- c(0,1,4)
theta_lambda2 <- c(0,1,4)
delta_1 <- c(1,1/4,1/24)
delta_2 <- c(1,1/4,1/24)

# Get all combinations of the interaction parameters
all_param_comb <- expand.grid(theta_lambda1, theta_lambda2, delta_1, delta_2)
names(all_param_comb) <- c("theta_lambda1", "theta_lambda2", "delta_1", "delta_2")

# order starting parameters by theta
all_param_comb <- all_param_comb[order(all_param_comb$theta_lambda1),]

# for now just keep symmetric interactions 
all_param_comb <- all_param_comb %>% filter(theta_lambda1 == theta_lambda2 & delta_1 == delta_2)

# specifying all global variables
beta_sd1 <- 0
beta_sd2 <- 0
tot_weeks <- 625

# plotting 
temp <- vector(mode = "list", length = dim(all_param_comb)[1])
plot_list <- vector(mode = "list", length = dim(all_param_comb)[1])
attack_plots <- vector(mode = "list", length = dim(all_param_comb)[1])
res_all <- NULL
for(i in 1:dim(all_param_comb)[1]){
  
  # if using data after simulation 
  #data <- get(paste0("results", i))$data
  # if getting plots before simulation 
  theta_lambda1 <- all_param_comb[i,]$theta_lambda1
  theta_lambda2 <- all_param_comb[i,]$theta_lambda2
  delta_1 <- all_param_comb[i,]$delta_1
  delta_2 <- all_param_comb[i,]$delta_2
  temp[[i]] <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2,
                          delta_1=delta_1, delta_2=delta_2, beta_sd1=beta_sd1, beta_sd2=beta_sd2,
                          n_surge=n_surge, components_l=components_l)
  data <- temp[[i]]$data
  
  data <- data %>% dplyr::select(time_date,v1_obs,v2_obs,v1_T,v2_T)

  legend_colors <- c("v1_obs" = "black", "v2_obs" = "blue")
  plot_list[[i]] <- ggplot(aes(x=time_date, y=v1_obs, colour="v1_obs"),data=data) + geom_line() + geom_line(aes(x=time_date, y=v2_obs,colour="v2_obs")) +
    ggtitle(paste("theta_lambda1 =", temp[[i]]$true_param["theta_lambda1"], 
                  "AND theta_lambda2 =",temp[[i]]$true_param["theta_lambda2"],
                  "AND delta_1 =", temp[[i]]$true_param["delta1"],
                  "AND delta_2 =", temp[[i]]$true_param["delta2"]))  + labs(y="observed cases") +
    scale_x_date(date_breaks = "3 month", date_labels =  "%b %Y")  +
    theme(axis.text.x=element_text(angle=60, hjust=1)) +  geom_vline(xintercept = t_si_date, linetype="dotted") +
    scale_colour_manual(values=legend_colors) + labs(colour="")


 # also estimate attack rates by year for each plot...... NOT WORKING
  data$season <- c(rep(1:tot_seasons, each=52),tot_seasons+1)
  #data$season <- c(rep(1:tot_seasons, each=52)) # for 100 years
  
  #data$season <- rep(1:tot_seasons, each=52)
  seasonal_incidence <- data %>% group_by(season) %>%
                          summarise(obs_v1 = sum(v1_obs), obs_v2 = sum(v2_obs),
                                    tot_v1 = sum(v1_T), tot_v2 = sum(v2_T))
  # trying to calculate the attack rate based on observed data
  obs_v1_attack <- seasonal_incidence$obs_v1/3700000 * 100
  obs_v2_attack <- seasonal_incidence$obs_v2/3700000 * 100

  range_obs_v1_att <- range(obs_v1_attack[-length(obs_v1_attack)]) # 0.14 - 0.19
  range_obs_v2_att <- range(obs_v2_attack[-length(obs_v2_attack)]) # 0.18 - 0.20

  # trying to calculate the attack rate based on true number of cases from the model
  tot_v1_attack <- seasonal_incidence$tot_v1/3700000 * 100
  tot_v2_attack <- seasonal_incidence$tot_v2/3700000 * 100
  tot_v1_attack <- tot_v1_attack[-c(length(tot_v1_attack))]
  tot_v2_attack <- tot_v2_attack[-c(length(tot_v2_attack))]

  plot_dat <- data.frame(cbind(tot_v1_attack = tot_v1_attack, tot_v2_attack = tot_v2_attack))
  attack_plots[[i]] <- ggplot(aes(x=tot_v2_attack,y=tot_v1_attack), data=plot_dat) + geom_point() +
    ggtitle(paste("theta_lambda1 =", temp[[i]]$true_param["theta_lambda1"], 
                  "AND theta_lambda2 =",temp[[i]]$true_param["theta_lambda2"],
                  "AND delta_1 =", temp[[i]]$true_param["delta1"],
                   "AND delta_2 =", temp[[i]]$true_param["delta1"]))

  range_tot_v1_att <- range(tot_v1_attack[-length(tot_v1_attack)]) # 71 - 95
  range_tot_v2_att <- range(tot_v2_attack[-length(tot_v2_attack)]) # 90 - 97

  res <- cbind(all_param_comb[i,],
               range_obs_v1_att[1],range_obs_v1_att[2],
               range_obs_v2_att[1],range_obs_v2_att[2],
               range_tot_v1_att[1],range_tot_v1_att[2],
               range_tot_v2_att[1],range_tot_v2_att[2])
  res_all <- rbind(res_all, res)

}

# plot simulated timeseries data
for(i in seq(from=1,to=dim(all_param_comb)[1], by=3)){
  grid.arrange(plot_list[[i]],plot_list[[i+1]],plot_list[[i+2]],ncol=1)
}


# plot scatter plots of seasonal attack rates
for(i in seq(from=1,to=dim(all_param_comb)[1], by=3)){
  grid.arrange(attack_plots[[i]],attack_plots[[i+1]],attack_plots[[i+2]],ncol=1)
}

