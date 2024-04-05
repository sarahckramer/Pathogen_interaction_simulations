################################################################################
#                       Results extraction 
#
# In this script we extract the results of the simulation study
# 
# Created by: Sarah Pirikahu
# Creation date: 13 March 2024
################################################################################


# load packages
library(tidyverse)
library(openxlsx) 
library(gridExtra)

# reading in all results 
file_names = list.files(pattern = "*.RData", recursive = F)

results_all <- vector(mode = "list", length = length(file_names))

# load all the data sets in 
for (i in 1:length(file_names)) {
  load(file_names[i], verbose = TRUE)
  results_all[[i]] <- results
  rm(results)
}

results <- results_all

#######################################################
#-------- Extract and collapse results ---------------#
#######################################################

# function to pull out certain element from each list of lists
get_elements <- function(x, element) {
  if(is.list(x))
  {
    if(element %in% names(x)) x[[element]]
    else lapply(x, get_elements, element = element)
  }
}

# Add true parameters to each results data frames
results <- lapply(results, function(x) {
  x$cor <- cbind(x$cor,data.frame(t(x$true_param)))
  x$gam_cor <- cbind(x$gam_cor,data.frame(t(x$true_param)))
  x$transfer_entropy <- cbind(x$transfer_entropy,data.frame(t(x$true_param)))
  x$granger <- cbind(x$granger,data.frame(t(x$true_param)))
  #x$CCM <- cbind(x$CCM,data.frame(t(x$true_param)))
  return(x)
})


#---- extract Spearmann correlation coefficients ----#
cor_res <- bind_rows(get_elements(results, "cor"))

# since correlations aren't directional we are going to claim if there is an interaction 
# in either direction that there is a true interaction present
# e.g. if theta_lambda_i for i \in {1,2} != 1 then interaction present

# creating column for true interaction Y/N
cor_res$true_int <- NA
cor_res[cor_res$theta_lambda1 != 1 | cor_res$theta_lambda2 != 1, ]$true_int <- "Y"
cor_res[cor_res$theta_lambda1 == 1 & cor_res$theta_lambda2 == 1,]$true_int <- "N"

# create column for estimated interaction Y/N 
# if p value is significant then there is cor significantly differnt than 0 suggesting 
# interaction 
cor_res$int_est <- NA
cor_res[cor_res$p_value < 0.05,]$int_est <- "Y"
cor_res[cor_res$p_value >= 0.05,]$int_est <- "N"

# confusion matrix
t_cor <- table(cor_res$true_int, cor_res$int_est, useNA="ifany")
# overall accuracy
(t_cor[1,1] + t_cor[2,2])/sum(t_cor)*100 # 91%

# plotting
ggplot(aes(x=true_int, y = cor, colour=int_est),data=cor_res) + geom_point()


# check out misclassified
cor_missclassified <- cor_res %>% filter(true_int != int_est) 
cor_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
cor_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
cor_missclassified %>% group_by(delta1,delta2) %>% tally() 


#---- extract GAM correlations ----# 
gam_res <-bind_rows(get_elements(results, "gam_cor"))
gam_res <- data.frame(gam_res, row.names = NULL)

# creating column for true interaction Y/N
gam_res$true_int <- NA
gam_res[gam_res$theta_lambda1 != 1 | gam_res$theta_lambda2 != 1, ]$true_int <- "Y"
gam_res[gam_res$theta_lambda1 == 1 & gam_res$theta_lambda2 == 1,]$true_int <- "N"

# create column for estimated interaction Y/N 
# is 0 in the interval
temp <- between(rep(0,dim(gam_res)[1]), gam_res$CI_lower95, gam_res$CI_upper95)
gam_res$int_est <- ifelse(temp ==FALSE, "Y", "N")

# confusion matrix
t_gam <- table(gam_res$true_int, gam_res$int_est, useNA="ifany")
# overall accuracy 
(t_gam[1,1] + t_gam[2,2])/sum(t_gam)*100 # 96%

# plotting
ggplot(aes(x=true_int, y=cor, colour=int_est), data=gam_res) + geom_point()

# check out misclassified
gam_cor_missclassified <- gam_res %>% filter(true_int != int_est) 
gam_cor_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
gam_cor_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
gam_cor_missclassified %>% group_by(delta1,delta2) %>% tally() 

#---- extract transfer entropy----# 
trans_res <- bind_rows(get_elements(results, "transfer_entropy"))

# create interaction Y/N outcome variable 
trans_res$true_int <- NA
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 != 1,]$true_int <- "Y"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 != 1,]$true_int <- "Y"
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 == 1,]$true_int <- "N"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 == 1,]$true_int <- "N"

# create column for estimated interaction Y/N
trans_res$int_est <- NA
trans_res[trans_res$p.value < 0.05,]$int_est <- "Y"
trans_res[trans_res$p.value >= 0.05,]$int_est <- "N"

# confusion matrix
t_te <- table(trans_res$true_int, trans_res$int_est,trans_res$direction)
# overall accuracy 
# v1 -> v2
(t_te[1,1,1] + t_te[2,2,1])/sum(t_te[,,1])*100 # 33%
# v2 -> v1
(t_te[1,1,2] + t_te[2,2,2])/sum(t_te[,,2])*100 # 22%

# plotting
ggplot(aes(x=true_int, y=ete, colour=int_est), data=trans_res) + geom_point() + facet_grid(.~direction)

# check out misclassified
te_missclassified <- trans_res %>% filter(true_int != int_est) 
te_missclassified %>% group_by(direction) %>% tally() %>% mutate(precent=n/sum(n)*100)
te_missclassified %>% group_by(theta_lambda1,theta_lambda2,delta1,delta2) %>% tally() 

#---- extract granger ----# 
granger_res <- bind_rows(get_elements(results, "granger"))
granger_res <- rename(granger_res,c('direction'='original_estimate.direction', 
                                    'logRSS'='original_estimate.logRSS'))

# create interaction Y/N outcome variable 
granger_res$true_int <- NA
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 != 1,]$true_int <- "Y"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 != 1,]$true_int <- "Y"
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 == 1,]$true_int <- "N"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 == 1,]$true_int <- "N"
granger_res$true_int <- as.factor(granger_res$true_int)

# true signed outcome
granger_res$true_int_sign <- NA
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 > 1,]$true_int_sign <- "+ve"
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 < 1,]$true_int_sign <- "-ve"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 > 1,]$true_int_sign <- "+ve"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 < 1,]$true_int_sign <- "-ve"
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 == 1,]$true_int_sign <- "N"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 == 1,]$true_int_sign <- "N"
granger_res$true_int_sign <- as.factor(granger_res$true_int_sign)
# check 
table(granger_res$true_int, granger_res$true_int_sign, useNA="ifany") # correct assignment

# create column for estimated interaction 

# Y/N outcome
temp <- between(rep(0,dim(granger_res)[1]), granger_res$blockboot_CI_lower95_v1_x_v2, granger_res$blockboot_CI_upper95_v1_x_v2)
granger_res$int_est <- ifelse(temp == FALSE, "Y", "N")

# signed outcome 
granger_res$int_est_sign <- NA
granger_res[temp==F & granger_res$logRSS > 0 ,]$int_est_sign <- "+ve"
granger_res[temp==F & granger_res$logRSS < 0 ,]$int_est_sign <- "-ve"
granger_res[temp==T,]$int_est_sign <- "N"

# confusion matrices
# Y/N outcome
t_granger <- table(granger_res$true_int, granger_res$int_est, granger_res$direction)
# overall accuracy 
# v1 -> v2
(t_granger[1,1,1] + t_granger[2,2,1])/sum(t_granger[,,1])*100 # 70%
# v2 -> v1
(t_granger[1,1,2] + t_granger[2,2,2])/sum(t_granger[,,2])*100 # 73%

# plotting 
ggplot(aes(x=true_int,y=logRSS,colour=int_est), data=granger_res) + geom_point() + facet_grid(.~direction)

# signed outcome 
table(granger_res$true_int_sign, granger_res$int_est_sign, granger_res$direction)

# plotting 
ggplot(aes(x=true_int_sign,y=logRSS,colour=int_est_sign), data=granger_res) + geom_point() + facet_grid(.~direction)

# checking out missclassifications
granger_missclassified <- granger_res %>% filter(true_int != int_est) 
granger_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
granger_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
granger_missclassified %>% group_by(delta1,delta2) %>% tally() 

