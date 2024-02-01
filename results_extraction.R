################################################################################
#                       Results extraction 
#
# In this script we extract the results of the simulation study for likelihood 
# and all other methods seprately - as the likelihood is more computationally 
# expensive it is run on its own, whereas all other methods are performed 
# together
# 
# Created by: Sarah Pirikahu
# Creation date: 20 Oct 2023
################################################################################

# load packages
library(tidyverse)

#------- all methods except likelihood extraction ------# 

# reading in all results for 10 years of data
setwd("~/Documents/Simulation_estimating_interaction/20years")
# pick up all the results file names
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
list_names <- paste0("results", 1:5) # change this to length(file_names) when have all results
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

results_df <- cbind(results_df, te_res_v1_x_v2_df,te_res_v2_x_v1_df)

#---- extract granger ----# 
granger_res_v1_xmap_v2 <- purrr::map(list_names, ~ get(.x)$granger$summary[1,])
granger_res_v2_xmap_v1 <- purrr::map(list_names, ~ get(.x)$granger$summary[2,])

granger_res_v1_xmap_v2 <- data.frame(bind_rows(granger_res_v1_xmap_v2), row.names = NULL)
names(granger_res_v1_xmap_v2) <- c("granger_est_v1_x_v2", "granger_bkbootmean_v1_x_v2",
                                   "granger_blockboot_CI_2.5_v1_x_v2","granger_blockboot_CI_97.5_v1_x_v2",
                                   "granger_blockboot_CIperc_2.5_v1_x_v2", "granger_blockboot_CIperc_97.5_v1_x_v2",
                                   "granger_p_v1_x_v2", "granger_adf_p_v1_x_v2", "granger_kpss_p_v1_x_v2")
granger_res_v2_xmap_v1 <- data.frame(bind_rows(granger_res_v2_xmap_v1), row.names = NULL)
names(granger_res_v2_xmap_v1) <- c("granger_est_v2_x_v1", "granger_bkbootmean_v2_x_v1",
                                   "granger_blockboot_CI_2.5_v2_x_v1","granger_blockboot_CI_97.5_v2_x_v1",
                                   "granger_blockboot_CIperc_2.5_v2_x_v1", "granger_blockboot_CIperc_97.5_v2_x_v1",
                                   "granger_p_v2_x_v1", "granger_adf_p_v2_x_v1", "granger_kpss_p_v2_x_v1")

results_df <- cbind(results_df, granger_res_v1_xmap_v2, granger_res_v2_xmap_v1)

#------- extract convergent cross mapping ----# 








#------- likelihood results extraction ----# 
str(results, max.level=1)

res <- NULL
for(i in 1:500){
  res <- cbind(res, results$likelihood[[i]]$estpars)
  
}

res <- t(res)
res <- data.frame(res)
# make data long so I can plot it

res_long <- gather(res, parameter, estimate, Ri1:w2, factor_key=TRUE)
ggplot(aes(x=estimate), data=res_long) + geom_histogram() + facet_wrap(.~parameter, scales="free")

true_ests <- results$true_param %>% dplyr::select(Ri1,Ri2,E01,E02,R01,R02,R12,rho1,rho2,A1,phi1,A2,
                                                  phi2,delta1,delta2,theta_lambda1,theta_lambda2,w1,w2)

res_long
res_summary <- res_long %>% group_by(parameter) %>% summarise(mean=mean(estimate), median=median(estimate),
                                                              sd=sd(estimate)) 

res_summary <- cbind(res_summary, t(true_ests))
names(res_summary)[5] <- "true"

test <- res_long %>% filter(parameter=="delta1" & estimate < 100)
summary(test)
test <- res_long %>% filter(parameter=="delta2" & estimate < 100)
summary(test)

ggplot(aes(x=estimate),data=test) + geom_histogram(binwidth=0.5) + geom_vline(xintercept=1, colour="blue")


ggplot(aes(x=estimate), data=res_long) + geom_histogram() +
  facet_wrap(.~parameter, scales="free") +
  geom_vline(data = res_summary, mapping = aes(xintercept = true), colour="blue") 
