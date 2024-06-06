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
library(hrbrthemes)

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
  x$data <- cbind(x$data,data.frame(t(x$true_param)))
  #x$cor <- cbind(x$cor,data.frame(t(x$true_param)))
  #x$gam_cor <- cbind(x$gam_cor,data.frame(t(x$true_param)))
  #x$transfer_entropy <- cbind(x$transfer_entropy,data.frame(t(x$true_param)))
  #x$granger <- cbind(x$granger,data.frame(t(x$true_param)))
  #x$CCM <- cbind(x$CCM,data.frame(t(x$true_param)))
  return(x)
})


#--------- plotting out some datasets------------#

data_full <- data.frame()

#for(i in 1:225){
for(i in 1:45){
  print(i)
  data <- results[[i]]$data #%>% filter(.id==5) 
  data$delta1_and_delta2 <- with(data, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
  data$delta1_and_delta2 <- as.factor(data$delta1_and_delta2)
  
  data$theta_lambda1_and_theta_lambda2 <- with(data, paste0("theta_lambda1 = ",theta_lambda1," & ","theta_lambda2 = ",theta_lambda2))
  data <- data %>% dplyr::select(time_date, v1_obs, v2_obs, theta_lambda1, theta_lambda2, delta1, 
                          delta2, delta1_and_delta2, theta_lambda1_and_theta_lambda2)
    data_full <- rbind(data_full, data)
} 

data_full %>% group_by(theta_lambda1_and_theta_lambda2) %>% tally()
data_full %>% group_by(delta1_and_delta2) %>% tally()

# making labels more sensible 
data_full$delta1_and_delta2 <- as.factor(data_full$delta1_and_delta2)
levels(data_full$delta1_and_delta2)  <- c("d1=1 & d2=1",
                                          "d1=1/4 & d2=1",
                                          "d1=1/12 & d2=1",
                                          "d1=1 & d2=1/4",
                                          "d1=1/4 & d2=1/4",
                                          "d1=1/12 & d2=1/4",
                                          "d1=1 & d2=1/12",
                                          "d1=1/4 & d2=1/12",
                                          "d1=1/12 & d2=1/12")

# plot
# symmetric interactions
data_full %>% filter(theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 0 & theta_lambda2 = 0"|
                     theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 1 & theta_lambda2 = 1"|  
                     theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 4 & theta_lambda2 = 4") %>% 
  ggplot(aes(x=time_date, y=v1_obs), data=.) + geom_line() +
  geom_line(aes(x=time_date, y=v2_obs),colour="blue") +
  facet_grid(delta1_and_delta2~theta_lambda1_and_theta_lambda2)

# asymmetric interactions
data_full %>% filter(theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 0 & theta_lambda2 = 4"|
                       theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 0 & theta_lambda2 = 1"|  
                       theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 4 & theta_lambda2 = 1"|
                       theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 1 & theta_lambda2 = 0"|  
                       theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 1 & theta_lambda2 = 4"|
                       theta_lambda1_and_theta_lambda2 == "theta_lambda1 = 4 & theta_lambda2 = 0") %>% 
  ggplot(aes(x=time_date, y=v1_obs), data=.) + geom_line() +
  geom_line(aes(x=time_date, y=v2_obs),colour="blue") +
  facet_grid(delta1_and_delta2~theta_lambda1_and_theta_lambda2)


######################################################
#---- extract Spearmann correlation coefficients ----#
######################################################

# extract all data from list and create dataframe of just Pearsons results
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

# overall proportion correct for each parameter set
# making parameters factors
cor_res$theta_lambda1 <- as.factor(cor_res$theta_lambda1)
cor_res$theta_lambda2 <- as.factor(cor_res$theta_lambda2)
cor_res$delta1 <- as.factor(cor_res$delta1)
cor_res$delta2 <- as.factor(cor_res$delta2)
cor_res$int_est <- as.factor(cor_res$int_est)
cor_res$true_int <- as.factor(cor_res$true_int)

# proportion correct for each paramter set
tot_dat <- cor_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, .drop=FALSE) %>% 
  tally() # don't always have 100 datasets for each paraeter combo
tot_dat <- rename(tot_dat, denom = n)

# numerator 
num <- cor_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, true_int, int_est, .drop=FALSE) %>% 
  summarise(n = n(), .groups = "drop")
# denominator
overall_cor <- num %>% left_join(tot_dat)
overall_cor$perc <- overall_cor$n / overall_cor$denom * 100

# correct 
correct_cor <- overall_cor %>% filter(true_int == int_est)
overall_correct_cor <- correct_cor %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2) %>% 
  summarise(tot_correct = sum(n)) 
overall_correct_cor <- overall_correct_cor %>% left_join(tot_dat)
overall_correct_cor$percent <- overall_correct_cor$tot_correct/overall_correct_cor$denom *100

# heatmaps 

# delta symmetric 
p1 <- overall_correct %>% filter(delta1== 1 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1")

p2 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/4")

p3 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile()  + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/24")

grid.arrange(p1,p2,p3, ncol=1) # don't really need heat map for this - the conclusion is it 
# always assigns significant interaction therefore it only gets no interaction wrong

# delta asymmetric
p4 <- overall_correct %>% filter(delta1== 1 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/24")

p5 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4,  delta2 = 1/24")

p6 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24, delta2 =1")

p7 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24,  delta2 = 1/4")

p8 <- overall_correct %>% filter(delta1== 0.25 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4, delta2 = 1")

p9 <- overall_correct %>% filter(delta1== 1 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/4")

grid.arrange(p4,p5,p6,p7,p8,p9, ncol=1)

# create column with combo of delta 1 and delta 2 to try multi facet grid to put all 
# the above plots together 
overall_correct$delta1_and_delta2 <- with(overall_correct, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
overall_correct$delta1_and_delta2 <- as.factor(overall_correct$delta1_and_delta2)
levels(overall_correct$delta1_and_delta2) <- c("d1=1/24 & d2=1/24", 
                                               "d1=1/24 & d2=1/4",
                                               "d1=1/24 & d2=1",
                                               "d1=1/4 & d2=1/24",
                                               "d1=1/4 & d2=1/4",
                                               "d1=1/4 & d2=1",
                                               "d1=1 & d2=1/24",
                                               "d1=1 & d2=1/4",
                                               "d1=1 & d2=1")

ggplot(data=overall_correct, aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") + facet_grid(delta1_and_delta2~.)


# check out misclassified
cor_missclassified <- cor_res %>% filter(true_int != int_est) 
cor_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
cor_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
cor_missclassified %>% group_by(delta1,delta2) %>% tally() 

####################################
#---- extract GAM correlations ----#
####################################

# extract all data from list and create dataframe of just GAM results
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

# overall proportion correct for each parameter set
# making parameters factors
gam_res$theta_lambda1 <- as.factor(gam_res$theta_lambda1)
gam_res$theta_lambda2 <- as.factor(gam_res$theta_lambda2)
gam_res$delta1 <- as.factor(gam_res$delta1)
gam_res$delta2 <- as.factor(gam_res$delta2)
gam_res$int_est <- as.factor(gam_res$int_est)
gam_res$true_int <- as.factor(gam_res$true_int)

# proportion correct for each paramter set
tot_dat <- gam_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, .drop=FALSE) %>% 
  tally() # don't always have 100 datasets for each paraeter combo
tot_dat <- rename(tot_dat, denom = n)

# numerator 
num <- gam_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, true_int, int_est, .drop=FALSE) %>% 
  summarise(n = n(), .groups = "drop")
# denominator
overall_gam <- num %>% left_join(tot_dat)
overall_gam$perc <- overall_gam$n / overall_gam$denom * 100

# correct 
correct_cor <- overall_gam %>% filter(true_int == int_est)
overall_correct_gam <- correct_cor %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2) %>% 
  summarise(tot_correct = sum(n)) 
overall_correct_gam <- overall_correct_gam %>% left_join(tot_dat)
overall_correct_gam$percent <- overall_correct_gam$tot_correct/overall_correct_gam$denom *100


# heatmaps 

# delta symmetric 
p1 <- overall_correct %>% filter(delta1== 1 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1")

p2 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/4")

p3 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile()  + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/24")

grid.arrange(p1,p2,p3, ncol=1)

# delta asymmetric
p4 <- overall_correct %>% filter(delta1== 1 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/24")

p5 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4,  delta2 = 1/24")

p6 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24, delta2 =1")

p7 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24,  delta2 = 1/4")

p8 <- overall_correct %>% filter(delta1== 0.25 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4, delta2 = 1")

p9 <- overall_correct %>% filter(delta1== 1 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/4")

grid.arrange(p4,p5,p6,p7,p8,p9, ncol=1)


# check out misclassified
gam_cor_missclassified <- gam_res %>% filter(true_int != int_est) 
gam_cor_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
gam_cor_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
gam_cor_missclassified %>% group_by(delta1,delta2) %>% tally() 

###################################
#---- extract transfer entropy----#
###################################

# extract all data from list and create dataframe of just TE results
trans_res <- bind_rows(get_elements(results, "transfer_entropy"))

# create interaction Y/N outcome variable 
trans_res$true_int <- NA
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 != 1,]$true_int <- "Y"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 != 1,]$true_int <- "Y"
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 == 1,]$true_int <- "N"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 == 1,]$true_int <- "N"

# create column for estimated interaction Y/N
trans_res$int_est <- NA
trans_res[trans_res$p_value < 0.05,]$int_est <- "Y"
trans_res[trans_res$p_value >= 0.05,]$int_est <- "N"

# choosing the best lag for each parameter combo
# have lags 1, 2, 4, 6 to choose from 
# the lag which gives the greatest transfer entropy is considered best
max_te_lag <- trans_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, id, direction) %>% 
  summarise(max_te = max(te))
# join the max lag with the overall dataset
trans_res <- trans_res %>% left_join(max_te_lag)
# remove the lags which don't give max transfer entropy
trans_res <- trans_res %>% filter(te == max_te)

# confusion matrix
t_te <- table(trans_res$true_int, trans_res$int_est,trans_res$direction)
# overall accuracy 
# v1 -> v2
(t_te[1,1,1] + t_te[2,2,1])/sum(t_te[,,1])*100 # 65.2%
# v2 -> v1
(t_te[1,1,2] + t_te[2,2,2])/sum(t_te[,,2])*100 # 58.4%

# making parameters factors
trans_res$theta_lambda1 <- as.factor(trans_res$theta_lambda1)
trans_res$theta_lambda2 <- as.factor(trans_res$theta_lambda2)
trans_res$delta1 <- as.factor(trans_res$delta1)
trans_res$delta2 <- as.factor(trans_res$delta2)
trans_res$int_est <- as.factor(trans_res$int_est)
trans_res$true_int <- as.factor(trans_res$true_int)

# confusion matrices for each paramter set
tot_dat <- trans_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, .drop=FALSE) %>% 
  tally() # don't always have 100 datasets for each paraeter combo
tot_dat <- rename(tot_dat, denom = n)

# numerator 
num <- trans_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, true_int, int_est, .drop=FALSE) %>% 
  summarise(n = n(), .groups = "drop")
# denominator
overall_te <- num %>% left_join(tot_dat)
overall_te$perc <- overall_te$n / overall_te$denom * 100

# overall total correct for each parameter combo
correct_te <- overall_te %>% filter(true_int == int_est)
overall_correct_te <- correct_te %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction) %>% 
  summarise(tot_correct = sum(n)) 
overall_correct_te <- overall_correct_te %>% left_join(tot_dat)
overall_correct_te$percent <- overall_correct_te$tot_correct/overall_correct_te$denom *100

# Heatmap 
# delta symmetric 
p1 <- overall_correct %>% filter(delta1== 1 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1")

p2 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/4")

p3 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/24")

grid.arrange(p1,p2,p3, ncol=1)

# delta asymmetric
p4 <- overall_correct %>% filter(delta1== 1 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/24")

p5 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4,  delta2 = 1/24")

p6 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24, delta2 =1")

p7 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24,  delta2 = 1/4")

p8 <- overall_correct %>% filter(delta1== 0.25 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4, delta2 = 1")

p9 <- overall_correct %>% filter(delta1== 1 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/4")

grid.arrange(p4,p5,p6,p7,p8,p9, ncol=1)


# create column with combo of delta 1 and delta 2 to try multi facet grid to put all 
# the above plots together 
overall_correct$delta1_and_delta2 <- with(overall_correct, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
overall_correct$delta1_and_delta2 <- as.factor(overall_correct$delta1_and_delta2)
levels(overall_correct$delta1_and_delta2) <- c("d1=1/24 & d2=1/24", 
                                               "d1=1/24 & d2=1/4",
                                               "d1=1/24 & da2=1",
                                               "d1=1/4 & d2=1/24",
                                               "d1=1/4 & d2=1/4",
                                               "d1=1/4 & d2=1",
                                               "d1=1 & d2=1/24",
                                               "d1=1 & d2=1/4",
                                               "d1=1 & d2=1")

ggplot(data=overall_correct, aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") + facet_grid(delta1_and_delta2~direction)


# plotting
ggplot(aes(x=true_int, y=te, colour=int_est), data=trans_res) + geom_point() + facet_grid(.~direction)

# check out misclassified
te_missclassified <- trans_res %>% filter(true_int != int_est) 
te_missclassified %>% group_by(direction) %>% tally() %>% mutate(precent=n/sum(n)*100)
te_missclassified %>% group_by(theta_lambda1,theta_lambda2,delta1,delta2) %>% tally() 


###########################
#---- extract granger ----#
###########################

# extract all data from list and create dataframe of just Granger results
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

# create column for estimated interaction 

# Y/N outcome
granger_res$int_est <- ifelse(granger_res$granger_p <= 0.05, "Y", "N")

# confusion matrices
# Y/N outcome
t_granger <- table(granger_res$true_int, granger_res$int_est, granger_res$direction)
# overall accuracy 
# v1 -> v2
(t_granger[1,1,1] + t_granger[2,2,1])/sum(t_granger[,,1])*100 # 58%
# v2 -> v1
(t_granger[1,1,2] + t_granger[2,2,2])/sum(t_granger[,,2])*100 # 62%

# plotting 
ggplot(aes(x=true_int,y=logRSS,colour=int_est), data=granger_res) + geom_point() + facet_grid(.~direction)


# overall proportion correct for each parameter set
# making parameters factors
granger_res$theta_lambda1 <- as.factor(granger_res$theta_lambda1)
granger_res$theta_lambda2 <- as.factor(granger_res$theta_lambda2)
granger_res$delta1 <- as.factor(granger_res$delta1)
granger_res$delta2 <- as.factor(granger_res$delta2)
granger_res$int_est <- as.factor(granger_res$int_est)
granger_res$true_int <- as.factor(granger_res$true_int)

# proportion correct for each paramter set
tot_dat <- granger_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, .drop=FALSE) %>% 
  tally() # don't always have 100 datasets for each paraeter combo
tot_dat <- rename(tot_dat, denom = n)

# numerator 
num <- granger_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, true_int, int_est, .drop=FALSE) %>% 
  summarise(n = n(), .groups = "drop")
# denominator
overall_granger <- num %>% left_join(tot_dat)
overall_granger$perc <- overall_granger$n / overall_granger$denom * 100

# correct 
correct_ga <- overall_granger %>% filter(true_int == int_est)
overall_correct_ga <- correct_ga %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction) %>% 
  summarise(tot_correct = sum(n)) 
overall_correct_ga <- overall_correct_ga %>% left_join(tot_dat)
overall_correct_ga$percent <- overall_correct_ga$tot_correct/overall_correct_ga$denom *100

# Heatmap 
# delta symmetric 
p1 <- overall_correct %>% filter(delta1== 1 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1")

p2 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/4")

p3 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = delta2 = 1/24")

grid.arrange(p1,p2,p3, ncol=1)

# delta asymmetric
p4 <- overall_correct %>% filter(delta1== 1 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/24")

p5 <- overall_correct %>% filter(delta1== 0.25 & delta2==0.0416666666666667) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4,  delta2 = 1/24")

p6 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24, delta2 =1")

p7 <- overall_correct %>% filter(delta1== 0.0416666666666667 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/24,  delta2 = 1/4")

p8 <- overall_correct %>% filter(delta1== 0.25 & delta2==1) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1/4, delta2 = 1")

p9 <- overall_correct %>% filter(delta1== 1 & delta2==0.25) %>%
  ggplot(data=., aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + facet_grid(.~direction) + scale_fill_distiller(palette = "RdPu") +
  ggtitle("delta1 = 1, delta2 = 1/4")

grid.arrange(p4,p5,p6,p7,p8,p9, ncol=1)


# create column with combo of delta 1 and delta 2 to try multi facet grid to put all 
# the above plots together 
overall_correct$delta1_and_delta2 <- with(overall_correct, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
overall_correct$delta1_and_delta2 <- as.factor(overall_correct$delta1_and_delta2)
levels(overall_correct$delta1_and_delta2) <- c("d1=1/24 & d2=1/24", 
                                               "d1=1/24 & d2=1/4",
                                               "d1=1/24 & da2=1",
                                               "d1=1/4 & d2=1/24",
                                               "d1=1/4 & d2=1/4",
                                               "d1=1/4 & d2=1",
                                               "d1=1 & d2=1/24",
                                               "d1=1 & d2=1/4",
                                               "d1=1 & d2=1")

ggplot(data=overall_correct, aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") + facet_grid(delta1_and_delta2~direction)

# checking out missclassifications
granger_missclassified <- granger_res %>% filter(true_int != int_est) 
granger_missclassified %>% group_by(theta_lambda1,theta_lambda2, delta1,delta2) %>% tally() 
granger_missclassified %>% group_by(theta_lambda1,theta_lambda2) %>% tally() 
granger_missclassified %>% group_by(delta1,delta2) %>% tally() 

###########################
#---- extract CCM---------#
###########################

# extract all data from list and create dataframe of just CCM results
ccm_res <- bind_rows(get_elements(results, "CCM"))

# create true interaction Y/N outcome variable 
ccm_res$true_int <- NA
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$theta_lambda1 != 1,]$true_int <- "Y"
ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$theta_lambda2 != 1,]$true_int <- "Y"
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$theta_lambda1 == 1,]$true_int <- "N"
ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$theta_lambda2 == 1,]$true_int <- "N"
ccm_res$true_int <- as.factor(ccm_res$true_int)

# creating a column to specify significance 
ccm_res$int_est <- NA
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$ecdf_p_v2_x_v1 <= 0.05,]$int_est <- "Y"
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$ecdf_p_v2_x_v1 > 0.05,]$int_est <- "N"

ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$ecdf_p_v1_x_v2 <= 0.05,]$int_est <- "Y"
ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$ecdf_p_v1_x_v2 > 0.05,]$int_est <- "N"
# overall accuracy 
t_ccm <- table(ccm_res$true_int ,ccm_res$int_est, ccm_res$direction)
# v1 -> v2
(t_ccm[1,1,1] + t_ccm[2,2,1])/sum(t_ccm[,,1])*100 # 72%
# v2 -> v1
(t_ccm[1,1,2] + t_ccm[2,2,2])/sum(t_ccm[,,2])*100 # 61%

# checking to see if there are any surrogate rho values that are extremely different 
# from rho for the data... if this is the case it could suggest the surrogates are 
# not particularly representative of the null distribution taking into consideration 
# the underlying seasonality 
parameter_combos <- ccm_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2) %>% tally()

# Will need to check surrogates if the rho_data is 
# high and rho_surr are all close to 0.
for(i in 1:20){
  param <- merge(ccm_res, parameter_combos[i,1:4])

  print(ggplot(aes(x=rho, y=.id),data=param) + geom_point() + 
    geom_point(aes(x=surr_rho, y=.id, colour=int_est), shape=2) + 
    geom_errorbar(aes(xmin = surr_rho_2.5, xmax = surr_rho_97.5,colour=int_est), linetype="dashed") + 
    facet_grid(.~direction) +
    ggtitle(paste0("theta_lambda1=", param$theta_lambda1,
                   " theta_lambda2=", param$theta_lambda2, 
                   " delta1=", param$delta1,
                   " delta2=", param$delta2)) + 
    theme(plot.title = element_text(size = 10)) )
}
# after checking all the plots the only situations where it looks
# like we could have problems is when one or both of the delta's 
# are 1/24

# overall proportion correct for each parameter set
# making parameters factors
ccm_res$theta_lambda1 <- as.factor(ccm_res$theta_lambda1)
ccm_res$theta_lambda2 <- as.factor(ccm_res$theta_lambda2)
ccm_res$delta1 <- as.factor(ccm_res$delta1)
ccm_res$delta2 <- as.factor(ccm_res$delta2)
ccm_res$int_est <- as.factor(ccm_res$int_est)
ccm_res$true_int <- as.factor(ccm_res$true_int)

# proportion correct for each paramter set
tot_dat <- ccm_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, .drop=FALSE) %>% 
  tally() # don't always have 100 datasets for each paraeter combo
tot_dat <- rename(tot_dat, denom = n)

# numerator 
num <- ccm_res %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction, true_int, int_est, .drop=FALSE) %>% 
  summarise(n = n(), .groups = "drop")
# denominator
overall <- num %>% left_join(tot_dat)
overall$perc <- overall$n / overall$denom * 100

# correct 
correct_ccm <- overall %>% filter(true_int == int_est)
overall_correct <- correct_ccm %>% group_by(theta_lambda1, theta_lambda2, delta1, delta2, direction) %>% 
  summarise(tot_correct = sum(n)) 
overall_correct <- overall_correct %>% left_join(tot_dat)
overall_correct$percent <- overall_correct$tot_correct/overall_correct$denom *100


# create column with combo of delta 1 and delta 2 to try multi facet grid to put all 
# the above plots together 
overall_correct$delta1_and_delta2 <- with(overall_correct, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
overall_correct$delta1_and_delta2 <- as.factor(overall_correct$delta1_and_delta2)
levels(overall_correct$delta1_and_delta2) <- c("d1=1/24 & d2=1/24", 
                                               "d1=1/24 & d2=1/4",
                                               "d1=1/24 & d2=1",
                                               "d1=1/4 & d2=1/24",
                                               "d1=1/4 & d2=1/4",
                                               "d1=1/4 & d2=1",
                                               "d1=1 & d2=1/24",
                                               "d1=1 & d2=1/4",
                                               "d1=1 & d2=1")

ggplot(data=overall_correct, aes(x=theta_lambda1, y=theta_lambda2, fill= percent)) + 
  geom_tile() + scale_fill_distiller(palette = "RdPu") + facet_grid(delta1_and_delta2~direction)



###############################################################################
# Plotting all methods together
###############################################################################

# combining overall data for all the methods 

# create direction coloumn for GAM so that I can combine the dataframes 
overall_gam$direction <- NA 
overall_cor$direction <- NA 
# add method to each dataset
overall_cor$method <- "correlation"
overall_gam$method <- "gam"
overall_te$method <- "transfer entropy"
overall_granger$method <- "granger"

# combine all datasets
all <- rbind(overall_cor, overall_gam, overall_granger, overall_te)
dim(all) # 5400

# keep just the symmetric results
all_symmetric <- all %>% filter(theta_lambda1 == theta_lambda2)
dim(all_symmetric) # 1080


# pulling together data frames that provide over % correct for each parameter combination
# add method names
overall_correct_cor$method <- "correlation"
overall_correct_gam$method <- "gam"
overall_correct_te$method <- "transfer entropy"
overall_correct_ga$method <- "granger"

# combine dataframes
overall_correct <- rbind(overall_correct_cor, overall_correct_gam, overall_correct_te, overall_correct_ga)
# keep just symmetric outcomes
overall_correct_symmetric <- overall_correct %>% filter(theta_lambda1 == theta_lambda2)
# relabel deltas
levels(overall_correct_symmetric$delta1) <- c("6 months", "1 month", "1 week")
levels(overall_correct_symmetric$delta2) <- c("6 months", "1 month", "1 week")
# reorder factor levels 
overall_correct_symmetric$delta1 <- factor(overall_correct_symmetric$delta1, levels=c("1 week", "1 month", "6 months"))
overall_correct_symmetric$delta2 <- factor(overall_correct_symmetric$delta2, levels=c("1 week", "1 month", "6 months"))

# plotting 
# line plots overall % correct by interaction parameters
overall_correct_symmetric %>% filter(method=="granger") %>% ggplot(aes(x=delta2, y=percent,colour=theta_lambda1, group=theta_lambda1), data=.) + 
     geom_line(position=position_jitter(w=0, h=0.8))  + facet_grid(delta1~direction) + theme_bw()

library(ggh4x)
overall_correct_symmetric %>% filter(method=="granger" | method =="transfer entropy" | method == "gam") %>% ggplot(aes(x=delta2, y=percent,colour=theta_lambda1, group=theta_lambda1), data=.) + 
  geom_line(position=position_jitter(w=0, h=0.8))  + facet_nested(delta1 ~ method + direction) + theme_bw()

# sensitivity and specificity 

# numerator
num <- all_symmetric %>% group_by(delta1, delta2, true_int, int_est, method,direction) %>% summarise(n=sum(n)) 
# denominator
denom <- all_symmetric %>% group_by(delta1, delta2, method,direction) %>% summarise(n=sum(n)) 
names(denom)[5] <- "denom"
# join denominator to numerator dataframe
to_plot <- num %>% left_join(denom)
# create percentage
to_plot$perc <- to_plot$n/to_plot$denom*100

# columns to represent Se and Sp type I and type II errors 
to_plot$type <- NA
to_plot[to_plot$true_int=="N" & to_plot$int_est=="N" , ]$type <- "sp"
to_plot[to_plot$true_int=="Y" & to_plot$int_est=="Y" , ]$type <- "se"
to_plot[to_plot$true_int=="Y" & to_plot$int_est=="N" , ]$type <- "type II"
to_plot[to_plot$true_int=="N" & to_plot$int_est=="Y" , ]$type <- "type I"
table(to_plot$type, useNA="ifany")

# make data wide
to_plot_wide <- to_plot
# drop variables dont need in the wide dataset
to_plot_wide$int_est <- NULL
to_plot_wide$true_int <- NULL
to_plot_wide$n <- NULL
to_plot_wide$denom <- NULL
# make the data frame wide
to_plot_wide <- spread(to_plot_wide, type, perc)

to_plot_wide$delta1_and_delta2 <- with(to_plot_wide, paste0("delta1 = ",delta1," & ","delta2 = ",delta2))
to_plot_wide$delta1_and_delta2 <- as.factor(to_plot_wide$delta1_and_delta2)
levels(to_plot_wide$delta1_and_delta2) <- c("d1=1/24 & d2=1/24", 
                                             "d1=1/24 & d2=1/4",
                                             "d1=1/24 & d2=1",
                                             "d1=1/4 & d2=1/24",
                                             "d1=1/4 & d2=1/4",
                                             "d1=1/4 & d2=1",
                                             "d1=1 & d2=1/24",
                                             "d1=1 & d2=1/4",
                                             "d1=1 & d2=1")

# plotting
to_plot_wide %>% filter(method=="granger" | method=="transfer entropy") %>% 
  ggplot(aes(x=sp, y=se, colour=method), data=.) + geom_point() +
  facet_grid(delta1_and_delta2 ~ direction)
# table
tab1 <- to_plot_wide %>% filter(method=="granger" | method=="transfer entropy") %>% 
  select(method,direction,se,sp,`type I`, `type II`, delta1_and_delta2)
tab1$delta1 <- NULL
tab1$delta2 <- NULL

to_plot_wide %>% filter(method=="cor" | method=="gam") %>% 
  ggplot(aes(x=sp, y=se, colour=method), data=.) + geom_point() +
  facet_grid(delta1_and_delta2~.)
