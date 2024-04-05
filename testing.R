# when all parameter combos are run serially 


# load libraries
library(tidyverse)
library(testthat)
library(pomp)
library(janitor)
library(ggfortify)
library(future) # allows for parallel processing
library(foreach)
library(doParallel)

#--- reading in CSnippets ---# 
# read in the C code for the pomp model 
mod_code <- readLines('seitr_x_seitr.cpp')

# pull out the various components of the C code ready to feed into pomp
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'rsim', 'skel', 'toest', 'fromest')
# initialise list
components_l <- vector(mode = 'list', length = length(components_nm))
names(components_l) <- components_nm
# create list with code components
for (nm in components_nm) {
  components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
    str_flatten(collapse = '\n')
  components_l[[nm]] <- Csnippet(text = components_l[[nm]])
}

#---- setting parameter input values ----# 

# set seed:
set.seed(1234)

# set noise parameters
beta_sd1 <- 0 
beta_sd2 <- 0

# total number of simulated datasets to create for each parameter input 
nsim <- 1000

# total number of seasons
tot_weeks <- 625
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

# up dating the number of surges to input into parameter vector
n_surge <- length(t_si)

# initialize the rate of loss of immunity corresponding to each of the 
# surge times 
delta_i <- runif(n=length(t_si), min = 0.01*7, max=0.1*7)

# parameter inputs 
theta_lambda1 <- c(0,0.5,1,2,4)
theta_lambda2 <- c(0,0.5,1,2,4)
delta_1 <- c(1,1/4,1/24)
delta_2 <- c(1,1/4,1/24)

# Get all combinations of the interaction parameters
all_param_comb <- expand.grid(theta_lambda1, theta_lambda2, delta_1, delta_2)
names(all_param_comb) <- c("theta_lambda1", "theta_lambda2", "delta_1", "delta_2")

# to try a small number 
#all_param_comb <- all_param_comb[1:4,]
# remove parameter vectors 
rm(theta_lambda1, theta_lambda2, delta_1, delta_2)

# function to create list of true parameter inputs and simulated data 
# function takes a vector of the interaction parameters 
sim_data <- function(tot_weeks,theta_lambda1,theta_lambda2,delta_1,delta_2,beta_sd1,beta_sd2,n_surge,components_l=components_l,nsim){
  set.seed(2908)
  
  # setting parameters to weekly rates - params listed as daily in 
  # spreadsheet list of model parameters.xlsx
  # note also v1 = influenza; v2 = RSV
  true_params <- c(Ri1=1.1, Ri2=1.7,
                   sigma1=7, sigma2=7/5,
                   gamma1=7/5, gamma2=7/10,
                   delta1=delta_1, delta2=delta_2,
                   w1=1/78, w2=1/52,
                   mu = 0.0002, nu = 0.0002,
                   rho1 = 0.002, rho2 = 0.002,
                   theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                   A1=0.01, phi1=26,
                   A2=0.2, phi2=20,
                   beta_sd1=beta_sd1, beta_sd2=beta_sd2, 
                   N=3700000,
                   E01=0.001, E02=0.001,
                   R01=0.4, R02=0.25, R12=0.001, nsurges=n_surge,
                   t_si_=t(t_si), delta_i_=t(delta_i))
  
  #---- Create list to save the parameter sets and results of our different methods ---# 
  
  results <- vector(mode = "list", length = 7)
  results[[1]] <- true_params 
  names(results) <- c("true_param", "data", "cor", "gam_cor", "transfer_entropy", "CCM","granger")
  
  #---- create pomp object ---# 
  po <- pomp(data = data.frame(time = seq(from = 0, to = tot_weeks, by = 1), v1_obs = NA, v2_obs = NA),
             times = "time",
             t0 = 0,
             obsnames = c('v1_obs', 'v2_obs'),
             accumvars = c('v1_T', 'v2_T'),
             statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                            'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                            'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                            'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                            'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                            'v1_T', 'v2_T'),
             paramnames = names(true_params),
             params = true_params,
             partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
             globals = components_l[['globs']],
             dmeasure = components_l[['dmeas']],
             rmeasure = components_l[['rmeas']],
             rprocess = euler(step.fun = components_l[['rsim']], delta.t = 1),
             skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
             rinit = components_l[['rinit']]
  )
  
  # ----simulating data----#
  s1 <- simulate(po, nsim=nsim, times=1:tot_weeks, format="data.frame")
  
  # remove first 2 years where simulation isn't yet at equilibrium 
  s1 <- s1 %>% filter(time > 104) 
  
  # make time into dates based off week number
  s1$time_date <- lubridate::ymd( "2012-July-01" ) + lubridate::weeks(s1$time)
  
  # save results
  results$data <- s1  
  return(results)
}

results <- vector(mode = "list", length = dim(all_param_comb)[1])
# generate single true parameter sets and simulate data 
for(i in 1:dim(all_param_comb)[1]){
  theta_lambda1 <- all_param_comb[i,]$theta_lambda1
  theta_lambda2 <- all_param_comb[i,]$theta_lambda2
  delta_1 <- all_param_comb[i,]$delta_1
  delta_2 <- all_param_comb[i,]$delta_2
  results[[i]] <- sim_data(tot_weeks = tot_weeks, theta_lambda1=theta_lambda1, theta_lambda2=theta_lambda2, 
                      delta_1=delta_1, delta_2=delta_2, beta_sd1=beta_sd1, beta_sd2=beta_sd2,
                      n_surge = n_surge, components_l = components_l, nsim = nsim)
  # if outbreak never takes off or dies out then remove it 
  temp <- results[[i]]$data %>% group_by(.id) %>% mutate(sum_v1 = sum(v1_obs), sum_v2=sum(v2_obs))
  results[[i]]$data <- temp %>% filter(sum_v1 !=0 & sum_v2 != 0)
}

#---- Correlation coefficents --------# 
# function to estimate correlation 
# inputs: v1_obs = time series for v1_obs
#         v2_obs = time series for v2_obs
corr_func <- function(v1_obs, v2_obs){
  # calculated correlation coefficent 
  cor_raw <- cor.test(v1_obs, v2_obs)
  # pull together results in data frame 
  temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
  names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
  return(temp_res)
}

# apply correlation function to all simulated datasets and save results 
for(i in 1:dim(all_param_comb)[1]){
results[[i]]$cor <- results[[i]]$data %>% group_by(.id) %>% do((corr_func(.$v1_obs,.$v2_obs)))
}


#--- GAM approach ---# 
source("./methods/gam_cor.R")

# setting up parallelism for the foreach loop
registerDoParallel(cl <- makeCluster(10))
for(j in 1:dim(all_param_comb)[1]){
  res <- results[[j]]
  # apply the GAM correlation approach to each simulated data set and save the results
  res_gam_cor <- foreach(i=1:nsim, .packages=c("tidyverse","mgcv","vars", "boot")) %dopar%{
    if(dim(res$data %>% filter(.id==i))[1]!=0){
     res$data %>% filter(.id==i) %>% gam_cor(.)
    }
  }
  results[[j]]$gam_cor <- do.call(rbind, res_gam_cor)
  print(j)
}


#----- Transfer entropy analysis ------# 
source("./methods/transfer_entropy_v2.R")

# setting up parallelism for the foreach loop
registerDoParallel(cl <- makeCluster(10))
for(j in 1:dim(all_param_comb)[1]){
  res <- results[[j]]
  # apply transfer entropy to each simulated data set and save the results
  res_te <- foreach(i=1:nsim, .packages=c("tidyverse","RTransferEntropy","vars")) %dopar%{
    if(dim(res$data %>% filter(.id==i))[1]!=0){
    res$data %>% filter(.id==i) %>% te_func(.)
    }
  }
  results[[j]]$transfer_entropy <- do.call(rbind, res_te)
  print(j)
}

#------- Transfer entropy jidt --------# 
source("./methods/transfer_entropy_jidt.R")

for(j in 1:dim(all_param_comb)[1]){
  # apply transfer entropy to each simulated data set and save the results
  # lag = 1
  results[[j]]$transfer_entropy <- results[[j]]$data %>% group_by(.id) %>% do(te_jidt(., lag="1"))
  # lag = 2
  temp <- results[[j]]$data %>% group_by(.id) %>% do(te_jidt(., lag="2"))
  results[[j]]$transfer_entropy <- rbind(results[[j]]$transfer_entropy, temp)
  # lag = 4
  temp <- results[[j]]$data %>% group_by(.id) %>% do(te_jidt(., lag="4"))
  results[[j]]$transfer_entropy <- rbind(results[[j]]$transfer_entropy, temp)
  # lag = 6
  temp <- results[[j]]$data %>% group_by(.id) %>% do(te_jidt(., lag="6"))
  results[[j]]$transfer_entropy <- rbind(results[[j]]$transfer_entropy, temp)
  print(j)
}


#---- Granger causality analysis  ----# 
source("./methods/granger_analysis.R")

for(i in 1:dim(all_param_comb)[1]){
  # apply granger analysis to each simulated data set and save the results
  results[[i]]$granger <- results[[i]]$data %>% group_by(.id) %>% do(granger_func(.))
  print(i)
}


#------- Convergent Cross mapping analysis -------# 
source("./methods/CCM.R")

for(j in 1:dim(all_param_comb)[1]){
  res <- results[[j]]
  # setting up parallelism for the foreach loop
  registerDoParallel(cl <- makeCluster(10))
  # apply CCM to each simulated data set and save the results
  res_ccm <- foreach(i=1:nsim, .packages=c("tidyverse","rEDM","Kendall")) %dopar%{
    res$data %>% filter(.id==i) %>% ccm_func(.)
  } 
  
  results[[j]]$CCM <- do.call(rbind, res_ccm)
  print(j)
}


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
  #x$transfer_entropy <- cbind(x$transfer_entropy,data.frame(t(x$true_param)))
  x$granger <- cbind(x$granger,data.frame(t(x$true_param)))
  #x$CCM <- cbind(x$CCM,data.frame(t(x$true_param)))
  return(x)
})



#---- extract Spearmann correlation coefficients ----#
cor_res <- bind_rows(get_elements(results, "cor"))


#---- extract GAM correlations ----# 
gam_res <-bind_rows(get_elements(results, "gam_cor"))
gam_res <- data.frame(gam_res, row.names = NULL)

# create interaction outcome variables
if(symmetric==TRUE){
  # Y/N interaction variable 
  gam_res$true_int <- NA 
  gam_res[gam_res$theta_lambda1 != 1,]$true_int <- "Y"
  gam_res[gam_res$theta_lambda1 == 1,]$true_int <- "N"
  
  # +ve, no and -ve interaction variable 
  gam_res$true_int_sign <- NA 
  gam_res[gam_res$theta_lambda1 < 1,]$true_int_sign <- "-ve"
  gam_res[gam_res$theta_lambda1 > 1,]$true_int_sign <- "+ve"
  gam_res[gam_res$theta_lambda1 == 1,]$true_int_sign <- "N"
} else{
  # Y/N interaction variable 
  gam_res$int <- NA 
  gam_res[gam_res$theta_lambda1 != 1 | gam_res$theta_lambda2 != 1, ]$int <- "Y"
  gam_res[gam_res$theta_lambda1 == 1 & gam_res$theta_lambda2 == 1,]$int <- "N"
}


#---- extract transfer entropy ----# 
trans_res <- bind_rows(get_elements(results, "transfer_entropy"))
trans_res <- cbind(true_params_undirect,direction=trans_res$direction, te=trans_res$te, ete=trans_res$ete)
# create interaction Y/N outcome variable 
trans_res$int <- NA
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 != 1,]$int <- "Y"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 != 1,]$int <- "Y"
trans_res[trans_res$direction=="v1 -> v2" & trans_res$theta_lambda1 == 1,]$int <- "N"
trans_res[trans_res$direction=="v2 -> v1" & trans_res$theta_lambda2 == 1,]$int <- "N"
# if there are no vaules for the above combos an error will pop up


#---- extract granger ----# 
granger_res <- bind_rows(get_elements(results, "granger"))
names(granger_res)[1:3] <- c("id", "direction", "logRSS")

# create interaction Y/N outcome variable 
granger_res$true_int <- NA
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 != 1,]$true_int <- "Y"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 != 1,]$true_int <- "Y"
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 == 1,]$true_int <- "N"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 == 1,]$true_int <- "N"
granger_res$true_int <- as.factor(granger_res$true_int)

# create interaction +/N/- outcome variable 
granger_res$true_int_sign <- NA
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 < 1,]$true_int_sign <- "-ve"
granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 > 1,]$true_int_sign <- "+ve"

granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 < 1,]$true_int_sign <- "-ve"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 > 1,]$true_int_sign <- "+ve"

granger_res[granger_res$direction=="v1 -> v2" & granger_res$theta_lambda1 == 1,]$true_int_sign <- "N"
granger_res[granger_res$direction=="v2 -> v1" & granger_res$theta_lambda2 == 1,]$true_int_sign <- "N"
granger_res$true_int_sign <- as.factor(granger_res$true_int_sign)

# check out number of interactions
granger_res %>% group_by(true_int) %>% tally()
granger_res %>% group_by(true_int_sign) %>% tally()

#------- extract convergent cross mapping ----# 
ccm_res <- bind_rows(get_elements(results, "CCM"))
ccm_res <- cbind(true_params_undirect,direction=ccm_res$direction, rho=ccm_res$rho)

# create interaction Y/N outcome variable 
ccm_res$int <- NA
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$theta_lambda1 != 1,]$int <- "Y"
ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$theta_lambda2 != 1,]$int <- "Y"
ccm_res[ccm_res$direction=="v1 -> v2" & ccm_res$theta_lambda1 == 1,]$int <- "N"
ccm_res[ccm_res$direction=="v2 -> v1" & ccm_res$theta_lambda2 == 1,]$int <- "N"

# save out each data set so that we can perform prediction separately 


## ---------------- try out SVM --------------#
library(e1071)

###########################################
# trying out SVM for Pearsons correlation #
###########################################

# look at just symmetric datasets to start with 
symmetric_data <- cor_res %>% filter(theta_lambda1 == theta_lambda2 & delta1 == delta2)
dim(symmetric_data) # 1500
View(symmetric_data)

# create interaction outcome variables: For symmetric we are able to say something about the sign of the interaction 
# Y/N interaction variable 
symmetric_data$true_int <- NA 
symmetric_data[symmetric_data$theta_lambda1 != 1,]$true_int <- "Y"
symmetric_data[symmetric_data$theta_lambda1 == 1,]$true_int <- "N"
symmetric_data$true_int <- as.factor(symmetric_data$true_int)  

# +ve, no and -ve interaction variable 
symmetric_data$true_int_sign <- NA 
symmetric_data[symmetric_data$theta_lambda1 < 1,]$true_int_sign <- "-ve"
symmetric_data[symmetric_data$theta_lambda1 > 1,]$true_int_sign <- "+ve"
symmetric_data[symmetric_data$theta_lambda1 == 1,]$true_int_sign <- "N"
symmetric_data$true_int_sign <- as.factor(symmetric_data$true_int_sign)

# splitting the data into training and test set 
library(caTools)
set.seed(2908)
split = sample.split(symmetric_data$true_int_sign, SplitRatio = 0.75)
training_set = subset(symmetric_data, split == TRUE)
test_set = subset(symmetric_data, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(cor,true_int, true_int_sign)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(cor,true_int, true_int_sign)
test_set_s$.id <- NULL

# splits
training_set_s %>% group_by(true_int) %>% tally() # Y = 900, N = 225
test_set_s %>% group_by(true_int) %>% tally() # Y = 300, N = 75

training_set_s %>% group_by(true_int_sign) %>% tally() # Y = 900, N = 225
test_set_s %>% group_by(true_int_sign) %>% tally() # Y = 300, N = 75

# setting up k-fold cross validation 
# i.e. dividing the training set up into k groups and train the model on each of these groups
# the purpose of this step is to choose the appropriate parameter inputs for our overall model which we will run 
# on the entire training dataset
library(caret)
folds = createFolds(training_set_s$int, k = 3)

cv = lapply(folds, function(x) { # start of function
  # in the next two lines we will separate the Training set into it's 10 pieces
  training_fold = training_set_s[-x, ] # training fold =  training set minus (-) it's sub test fold
  test_fold = training_set_s[x, ] # here we describe the test fold individually
 
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial')
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
  
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 86.0 % missclassification happens only when N is assigned as Y

# try adjusting parameters within SVM 

# trying out using some weights to account for the large difference in no interaction in comparison to yes interaction 
#wts <- 1/table(symmetric_data$int); wts
#wts <- 10/table(symmetric_data$int); wts
#wts <- 100/table(symmetric_data$int); wts

cv = lapply(folds, function(x) {
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # apply the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   class.weights = wts)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

#--- accuracy ---#
# wts = 1/n
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 71.6 % # now there is also misclassification of Y that should be N
# wts = 10/n
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 80.5 
# wts = 100/n
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 65.7 % # now getting N correct always but assigning lots of Y to N 
# adding weights doesn't seem to help with the prediction accuracy 


# try adjusting gamma 
cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = 10)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 86.0 gamma = 1
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 80 gamma = 0.1
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 93.3 gamma = 10 same for gamma = 100

#----- try both wts and gamma together -----# 
wts <- 1/table(symmetric_data$int); wts
wts <- 10/table(symmetric_data$int); wts
wts <- 100/table(symmetric_data$int); wts

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   wts = wts,
                   gamma = 10)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 #93.3 same for all wts... so leave it out 

###############
# try on test #
###############
classifier = svm(int ~ cor,
                 data = training_set_s,
                 type = 'C-classification',
                 kernel = 'radial',
                 gamma = 10)
y_pred = predict(classifier, newdata = test_set_s)
cm = table(test_set_s$int, y_pred);cm
accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1]); accuracy * 100 # 92.5%
# Note if we just assigned all the data to Y we would have gotten 900/1125 *100 = 80 % accuracy

# checking exactly where it is getting it wrong more closely 
test_set$prediction_binary <- y_pred
table(test_set$int, test_set$prediction_binary)
# looking at a picture of the correlation via known interaction coloured by prediction made
ggplot(aes(x=int, y=cor, colour=prediction_binary), data=test_set) + geom_point()


# trying to make a confusion matrix just based on estimate v true without doing the 
# machine learning approach 
test_set$int_est <- NA
test_set[test_set$p_value < 0.05,]$int_est <- "Y"
test_set[test_set$p_value >= 0.05,]$int_est <- "N"
table(test_set$int_est, useNA="ifany")
# confusion matrix
table(test_set$int, test_set$int_est, useNA="ifany") # never assigns no interaction in reality gets it wrong 
# 20% of the time yet the SVM method suggests 95.5% accuracy on this data 

# looking at the entire symmetric dataset 
symmetric_data$int_est <- NA
symmetric_data[symmetric_data$p_value < 0.05,]$int_est <- "Y"
symmetric_data[symmetric_data$p_value >= 0.05,]$int_est <- "N"

# confusion matrix
table(symmetric_data$int, symmetric_data$int_est, useNA="ifany") # 20% if the time it is wrong


####------- kmeans approach -----###
# Decide how many clusters to look at
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)

set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(symmetric_data$cor, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

k = 3
set.seed(123)
data <- symmetric_data
km.out <- kmeans(data$cor, centers = k, nstart = 20)

data$cluster_id <- factor(km.out$cluster)
ggplot(data, aes(y=cor, x=int, color = cluster_id)) +
  geom_point() 

ggplot(data, aes(y=cor, x=int_sign, color = cluster_id)) +
  geom_point() 

data %>% group_by(cluster_id) %>% summarise(min=min(cor), max=max(cor))
# base output on cluster id 
data$int_est <- NA
data[data$cluster_id==1,]$int_est <- "Y"
data[data$cluster_id==2,]$int_est <- "N"
data[data$cluster_id==3,]$int_est <- "Y"

data$int_sign_est <- NA
data[data$cluster_id==1,]$int_sign_est <- "+ve"
data[data$cluster_id==2,]$int_sign_est <- "N"
data[data$cluster_id==3,]$int_sign_est <- "-ve"

table(data$int, data$int_est, useNA="ifany")
prop.table(table(data$int, data$int_est, useNA="ifany"))*100
# overall % accuracy 
(100 + 897)/dim(data)[1]*100 # 66.5%

table(data$int_sign, data$int_sign_est, useNA="ifany")
prop.table(table(data$int_sign, data$int_sign_est, useNA="ifany"))*100
# overall % accuracy 
(100 + 497 + 100)/dim(data)[1]*100 # 46.5

# trying kmeans with all data 
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(cor_res$cor, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

k=3
set.seed(123)
# Y/N interaction variable 
all_cor <- cor_res
all_cor$int <- NA 
all_cor[all_cor$theta_lambda1 != 1 | all_cor$theta_lambda2 != 1, ]$int <- "Y"
all_cor[all_cor$theta_lambda1 == 1 & all_cor$theta_lambda2 == 1,]$int <- "N"
all_cor$int <- as.factor(all_cor$int)

km.out <- kmeans(all_cor$cor, centers = k, nstart = 20)

all_cor$cluster_id <- factor(km.out$cluster)
ggplot(all_cor, aes(y=cor, x=int, color = cluster_id)) +
  geom_point() 

all_cor %>% group_by(cluster_id) %>% summarise(min=min(cor), max=max(cor))

all_cor$int_est <- NA
# k= 3
all_cor[all_cor$cluster_id==1,]$int_est <- "N"
all_cor[all_cor$cluster_id==2,]$int_est <- "Y"
all_cor[all_cor$cluster_id==3,]$int_est <- "Y"

# k=2
all_cor[all_cor$cluster_id==1,]$int_est <- "Y"
all_cor[all_cor$cluster_id==2,]$int_est <- "N"

t1 <- table(all_cor$int, all_cor$int_est, useNA="ifany")
prop.table(t1)*100

(t1[1,1] + t1[2,2])/sum(t1)*100 # 66.5% overall accuracy (k=3), 51.5% (k=2)

t2 <- table(all_cor$int, all_cor$int_est_st, useNA="ifany")
prop.table(t2)*100

(t2[1,1] + t2[2,2])/sum(t2)*100 # 91.3% overall accuracy


#----------------------------------------------------------------------------------------------#

# repeat with signed approach
set.seed(2908)
split = sample.split(symmetric_data$int_sign, SplitRatio = 0.75)
training_set = subset(symmetric_data, split == TRUE)
test_set = subset(symmetric_data, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(cor,int_sign)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(cor,int_sign)
test_set_s$.id <- NULL

# splits
training_set_s %>% group_by(int_sign) %>% tally() # -ve = 450, +ve = 450, N = 225
test_set_s %>% group_by(int_sign) %>% tally() # -ve = 150, +ve = 150, N = 75

# setting up k-fold cross validation 
# i.e. dividing the training set up into k groups and train the model on each of these groups
library(caret)
folds = createFolds(training_set_s$int_sign, k = 3)

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int_sign ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial')
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int_sign, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2] + cm[3,3]) / (cm[1,1] + cm[1,2] + cm[1,3] + cm[2,1] + cm[2,2] + cm[2,3] + cm[3,1] + cm[3,2] + cm[3,3])
  return(accuracy)
})

avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 74 %

# try adjusting parameters within SVM 

# trying out gamma params
cv = lapply(folds, function(x) {
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int_sign ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = 0.1)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int_sign, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2] + cm[3,3]) / (cm[1,1] + cm[1,2] + cm[1,3] + cm[2,1] + cm[2,2] + cm[2,3] + cm[3,1] + cm[3,2] + cm[3,3])
  return(accuracy)
})

avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # gamma= 0.1, 57.4; gamma = 1, 74%; gamma = 10, 90.5%; gamma = 100, 90.4%;

# addition of weights
wts <- 1/table(symmetric_data$int); wts
wts <- 10/table(symmetric_data$int); wts
wts <- 100/table(symmetric_data$int); wts

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int_sign ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = 10,
                   wts = wts)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int_sign, y_pred)
  accuracy = (cm[1,1] + cm[2,2] + cm[3,3]) / (cm[1,1] + cm[1,2] + cm[1,3] + cm[2,1] + cm[2,2] + cm[2,3] + cm[3,1] + cm[3,2] + cm[3,3])
  return(accuracy)
})
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # wts= 1/n, 90.5%; wts = 10/n, 90.5%; wts = 100/n, 90.5%; 
# so no improvement by adding the weights so leave them 

# try on test
classifier = svm(int_sign ~ cor,
                 data = training_set_s,
                 type = 'C-classification',
                 kernel = 'radial',
                 gamma = 10)
y_pred = predict(classifier, newdata = test_set_s)
cm = table(test_set_s$int_sign, y_pred);cm
accuracy =  (cm[1,1] + cm[2,2] + cm[3,3]) / (cm[1,1] + cm[1,2] + cm[1,3] + cm[2,1] + cm[2,2] + cm[2,3] + cm[3,1] + cm[3,2] + cm[3,3])
accuracy * 100 # 89.6%

# checking out more closely which ones are incorrect
test_set$prediction_sign <- y_pred
# test_set %>% dplyr::select(.id,theta_lambda1, theta_lambda2, delta1, delta2, int_sign,prediction_sign) %>% 
#   filter(int_sign != prediction_sign) %>% View()
# looking at a picture of the correlation via known interaction coloured by prediction made
ggplot(aes(x=int_sign, y=cor, colour=prediction_sign), data=test_set) + geom_point() + xlab("true interaction")
ggplot(aes(x=int_sign, y=cor, colour=prediction_sign), data=test_set) + geom_point() + facet_grid(.~delta1)

# trying to make a confusion matrix just based on estimate v true without doing the 
# machine learning approach 
test_set$int_est <- NA
test_set[test_set$cor > 0 & test_set$p_value < 0.05,]$int_est <- "+ve"
test_set[test_set$cor < 0 & test_set$p_value < 0.05,]$int_est <- "-ve"
test_set[test_set$p_value >= 0.05,]$int_est <- "N"
table(test_set$int_est, useNA="ifany")
# confusion matrix
table(test_set$int_sign, test_set$int_est) # never assigns no interaction 
# often assigns no or -ve interactions to positive


##--------------------------------------------------------------------------------------------------##

# looking at all parameter combos with the following outcome definition:
# If the interaction is asymmetric we will not be able to talk about sign and will have to make the assumption that if 
# there is an interaction in either direction then there is an overall interaction 

# Y/N interaction variable 
all_cor <- cor_res
all_cor$int <- NA 
all_cor[all_cor$theta_lambda1 != 1 | all_cor$theta_lambda2 != 1, ]$int <- "Y"
all_cor[all_cor$theta_lambda1 == 1 & all_cor$theta_lambda2 == 1,]$int <- "N"
all_cor$int <- as.factor(all_cor$int)

set.seed(2908)
split = sample.split(all_cor$int, SplitRatio = 0.75)
training_set = subset(all_cor, split == TRUE)
test_set = subset(all_cor, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(cor,int)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(cor,int)
test_set_s$.id <- NULL

# splits
training_set_s %>% group_by(int) %>% tally() # Y = 16166, N = 675
test_set_s %>% group_by(int) %>% tally() # Y = 5388, N = 225


# setting up k-fold cross validation 
folds = createFolds(training_set_s$int, k = 10)

# try random forest over SVM 
library(randomForest)

cv_rf = lapply(folds, function(x) {
  training_fold = training_set_s[-x, ]
  test_fold = training_set_s[x, ]
  rf <- randomForest(int ~ cor, data=training_fold,ntree=750, proximity=FALSE)
  cm <- rf$confusion
  accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy_rf <- mean(as.numeric(cv_rf));avg_accuracy_rf*100 # 92.7 %

# trying gradient boosted random forest
library(gbm)
cv_gbm = lapply(folds, function(x) {
  training_fold = training_set_s[-x, ]
  test_fold = training_set_s[x, ]
  training_fold$int <- ifelse(training_fold$int=="Y",1,0)
  test_fold$int <- ifelse(test_fold$int=="Y",1,0)
  gbm_out <- gbm(int ~ cor, data = training_fold,
                      distribution = "bernoulli", n.trees = 5000,
                      interaction.depth = 2)
  prob_Y <- predict(gbm_out, newdata=test_fold, n.trees=100, type="response")
  # Convert probabilities to binary predictions
  predictions <- factor(ifelse(prob_Y >= 0.5, 1, 0),levels=c(0,1))
  cm <- table(test_fold$int, predictions); print(cm)
  accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy_gbm <- mean(as.numeric(cv_gbm));avg_accuracy_gbm*100 # 96 %

# trying SVM
cv_svm = lapply(folds, function(x) {
  training_fold = training_set_s[-x, ]
  test_fold = training_set_s[x, ]
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial')
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred_svm = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred_svm); print(cm)
  accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})

avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 96 %
# if all in training set were just labbeled Y we would have 16166/(16166+675)*100 = 96% accuracy

# trying out weighting
wts <- 1/table(symmetric_data$int); wts
wts <- 10/table(symmetric_data$int); wts
wts <- 100/table(symmetric_data$int); wts
wts[1] <- 1000
wts[2] <- -1

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = 10,
                   wts = wts)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred);print(cm)
  accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 95.6 for all weights

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(int ~ cor,
                   data = training_fold,
                   type = 'C-classification',
                   kernel = 'radial',
                   gamma = 10)
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$int, y_pred);print(cm)
  accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
  return(accuracy)
})
avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 96% Just assigns all as Y always

# try on test
classifier = svm(int ~ cor,
                 data = training_set_s,
                 type = 'C-classification',
                 kernel = 'radial',
                 gamma = 10)
y_pred = predict(classifier, newdata = test_set_s)
test_set_s$prediction <- y_pred
cm = table(test_set_s$int, y_pred);print(cm)
accuracy =   accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[2,2] + cm[1,2] + cm[2,1])
accuracy * 100 # 96%

# plot what is going on 
ggplot(aes(x=int, y = cor, colour=prediction), data=test_set_s) + geom_point()

# trying to make a confusion matrix just based on estimate v true without doing the 
# machine learning approach 
all_cor$int_est_st <- NA
all_cor[all_cor$p_value < 0.05,]$int_est_st <- "Y"
all_cor[all_cor$p_value >= 0.05,]$int_est_st <- "N"

# confusion matrix
table(all_cor$int, all_cor$int_est) # never assigns no interaction correctly and often misclassifies
# interaction as well 


#-----------------------------------------------------------#


############################## 
# confusion matrices gam cor #
############################## 

# look at just symmetric datasets to start with 
symmetric_data <- gam_res %>% filter(theta_lambda1 == theta_lambda2 & delta1 == delta2)
dim(symmetric_data) # 1500

# create interaction outcome variables: For symmetric we are able to say something about the sign of the interaction 
# Y/N interaction variable 
symmetric_data$true_int <- NA 
symmetric_data[symmetric_data$theta_lambda1 != 1,]$true_int <- "Y"
symmetric_data[symmetric_data$theta_lambda1 == 1,]$true_int <- "N"
symmetric_data$true_int <- as.factor(symmetric_data$true_int)  

# +ve, no and -ve interaction variable 
symmetric_data$true_int_sign <- NA 
symmetric_data[symmetric_data$theta_lambda1 < 1,]$true_int_sign <- "-ve"
symmetric_data[symmetric_data$theta_lambda1 > 1,]$true_int_sign <- "+ve"
symmetric_data[symmetric_data$theta_lambda1 == 1,]$true_int_sign <- "N"
symmetric_data$true_int_sign <- as.factor(symmetric_data$true_int_sign)

# creating estimated interaction columns 

# if a CI contains 0 the product of the upper and lower bounds will always be <= 0 
temp <- symmetric_data$CI_lower95*symmetric_data$CI_upper95
symmetric_data$int_est <- ifelse(temp > 0, "Y", "N")
# confusion matrix 
table(symmetric_data$true_int, symmetric_data$int_est, useNA="ifany") # assigns all to Y interaction 
# % accuracy 1200/1500 *100 = 80%

# where output is +ve, -ve N 
symmetric_data$int_est_sign <- NA 
symmetric_data[temp > 0 & symmetric_data$cor > 0,]$int_est_sign <- "+ve"
symmetric_data[temp > 0 & symmetric_data$cor < 0,]$int_est_sign <- "-ve"
symmetric_data[temp <= 0,]$int_est_sign <- "N"

# confusion matrix
table(symmetric_data$true_int_sign,symmetric_data$int_est_sign, useNA="ifany") # mainly assigning +ve interactions, but does assign some -ve
# % accuracy 700/1500 * 100 = 46.7%


###################
# SVM on grangers #
###################


set.seed(2908)
split = sample.split(granger_res$true_int, SplitRatio = 0.75)
training_set = subset(granger_res, split == TRUE)
test_set = subset(granger_res, split == FALSE)

training_set_s <- training_set %>% dplyr::select(logRSS,direction,true_int, true_int_sign)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(logRSS,direction,true_int, true_int_sign)
test_set_s$.id <- NULL


folds = createFolds(training_set_s$true_int, k = 3)

cv = lapply(folds, function(x) { 
  training_fold = training_set_s[-x, ] 
  test_fold = training_set_s[x, ] 
  
  # now apply (train) the classifer on the training_fold
  classifier = svm(true_int ~ logRSS,
                   data = training_fold,
                   type = 'C-classification',
                   gamma = 10,
                   kernel = 'radial')
  # calculate the predictions and cm and we equate the accuracy
  # note we are training on training_fold and testing its accuracy on the test_fold
  y_pred = predict(classifier, newdata = test_fold)
  cm = table(test_fold$true_int, y_pred); print(cm)
  accuracy = (cm[1,1] + cm[2,2]) / (cm[1,1] + cm[1,2] + cm[2,1] + cm[2,2])
  return(accuracy)
})

avg_accuracy <- mean(as.numeric(cv));avg_accuracy*100 # 80%



### 

####------- kmeans approach -----###

# total number of clusters to consider
n_clusters <- 10

# trying kmeans with all data 
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(gam_res$cor, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

k=3
set.seed(123)
# Y/N interaction variable 
all_cor <- gam_res
all_cor$int <- NA 
all_cor[all_cor$theta_lambda1 != 1 | all_cor$theta_lambda2 != 1, ]$int <- "Y"
all_cor[all_cor$theta_lambda1 == 1 & all_cor$theta_lambda2 == 1,]$int <- "N"
all_cor$int <- as.factor(all_cor$int)

km.out <- kmeans(all_cor$cor, centers = k, nstart = 20)

all_cor$cluster_id <- factor(km.out$cluster)
ggplot(all_cor, aes(y=cor, x=int, color = cluster_id)) +
  geom_point() 

all_cor %>% group_by(cluster_id) %>% summarise(min=min(cor), max=max(cor))

all_cor$int_est <- NA
# k= 3
all_cor[all_cor$cluster_id==1,]$int_est <- "Y"
all_cor[all_cor$cluster_id==2,]$int_est <- "Y"
all_cor[all_cor$cluster_id==3,]$int_est <- "N"

t1 <- table(all_cor$int, all_cor$int_est, useNA="ifany")
prop.table(t1)*100

(t1[1,1] + t1[2,2])/sum(t1)*100 # 80% overall accuracy (k=3)


#-----------------------------------------------------------------------##
## transfer entropy  ## 

# total number of clusters to consider
n_clusters <- 10

# trying kmeans with all data 
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(trans_res$ete, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot


k=2
set.seed(123)
all_cor <- trans_res

km.out <- kmeans(all_cor$ete, centers = k, nstart = 20)

all_cor$cluster_id <- factor(km.out$cluster)
ggplot(all_cor, aes(y=ete, x=true_int, color = cluster_id)) +
  geom_point() + facet_grid(.~direction)


all_cor %>% group_by(cluster_id) %>% summarise(min=min(ete), max=max(ete))

all_cor$int_est_kmean <- NA
all_cor[all_cor$cluster_id==1,]$int_est_kmean <- "N"
all_cor[all_cor$cluster_id==2,]$int_est_kmean <- "Y"

t1 <- table(all_cor$true_int, all_cor$int_est_kmean, all_cor$direction, useNA="ifany")

# v1 -> v2
(t1[1,1,1] + t1[2,2,1])/sum(t1[,,1])*100 # 33% overall accuracy (k=2)
# v2 -> v1
(t1[1,1,2] + t1[2,2,2])/sum(t1[,,2])*100 # 23% overall accuracy (k=2) # similar to significant testing 

t2 <- table(all_cor$true_int, all_cor$int_est, all_cor$direction)

# v1 -> v2
(t2[1,1,1] + t2[2,2,1])/sum(t2[,,1])*100 # 33% overall accuracy 
# v2 -> v1
(t2[1,1,2] + t2[2,2,2])/sum(t2[,,2])*100 # 22% 


#-----------------------------------------------------------------------##
## granger analysis  ## 

# total number of clusters to consider
n_clusters <- 10

# trying kmeans with all data 
# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(granger_res$logRSS, centers = i, algorithm = "Lloyd", iter.max=500, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

k=2

set.seed(123)
all_cor <- granger_res

km.out <- kmeans(all_cor$logRSS, centers = k, nstart = 20)

all_cor$cluster_id <- factor(km.out$cluster)
ggplot(all_cor, aes(y=logRSS, x=true_int, color = cluster_id)) +
  geom_point() + facet_grid(.~direction) + scale_colour_manual(values=c("#00BFC4","#F8766D"))

all_cor %>% group_by(cluster_id) %>% summarise(min=min(logRSS), max=max(logRSS))

all_cor$int_est_kmean <- NA
all_cor[all_cor$cluster_id==1,]$int_est_kmean <- "Y"
all_cor[all_cor$cluster_id==2,]$int_est_kmean <- "N"

t1 <- table(all_cor$true_int, all_cor$int_est_kmean, all_cor$direction, useNA="ifany")

# v1 -> v2
(t1[1,1,1] + t1[2,2,1])/sum(t1[,,1])*100 # 38% overall accuracy (k=2)
# v2 -> v1
(t1[1,1,2] + t1[2,2,2])/sum(t1[,,2])*100 # 27% overall accuracy (k=2) # similar to significant testing 

t2 <- table(all_cor$true_int, all_cor$int_est, all_cor$direction)

# v1 -> v2
(t2[1,1,1] + t2[2,2,1])/sum(t2[,,1])*100 # 70% overall accuracy 
# v2 -> v1
(t2[1,1,2] + t2[2,2,2])/sum(t2[,,2])*100 # 73% 

###------------------------------------------------------------------------------------###

#######################
# logsitic regression #
#######################


## pearson correlation ## 
library(caTools)
set.seed(2908)
split = sample.split(symmetric_data$true_int_sign, SplitRatio = 0.75)
training_set = subset(symmetric_data, split == TRUE)
test_set = subset(symmetric_data, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(cor,true_int, true_int_sign, p_value)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(cor,true_int, true_int_sign, p_value)
test_set_s$.id <- NULL

# splits
training_set_s %>% group_by(true_int) %>% tally() # Y = 900, N = 225
test_set_s %>% group_by(true_int) %>% tally() # Y = 300, N = 75

training_set_s %>% group_by(true_int_sign) %>% tally() # +ve = 750, -ve = 750, N = 750
test_set_s %>% group_by(true_int_sign) %>% tally() # +ve = 250, -ve = 250, N = 250

# logistic regression model training set
library(VGAM) 
train_logit <- vglm(true_int_sign ~ cor, data=training_set_s, family="multinomial")
summary(train_logit)


pred <-  predict(train_logit, newdata = test_set_s, type = "response") 
predictions <- apply(pred, 1, which.max)
predictions[which(predictions=="1")] <- "-ve"
predictions[which(predictions=="2")] <- "+ve"
predictions[which(predictions=="3")] <- "N"

table(test_set_s$true_int_sign, predictions)

test_set_s$significance <- NA
test_set_s[test_set_s$p_value <=0.05 & test_set_s$cor > 0,]$significance <-  "+ve" 
test_set_s[test_set_s$p_value <=0.05 & test_set_s$cor < 0,]$significance <-  "-ve" 
test_set_s[test_set_s$p_value > 0.05 ,]$significance <-  "N" 

table(test_set_s$true_int_sign,test_set_s$significance)

ggplot(aes(y=cor,x=true_int_sign, colour=pred), data=test_set_s) + geom_point()


## GAM correlation ## 
symmetric_data <- gam_res %>% filter(theta_lambda1 == theta_lambda2 & delta1 == delta2)
dim(symmetric_data)  # 3000

set.seed(2908)
split = sample.split(symmetric_data$true_int_sign, SplitRatio = 0.75)
training_set = subset(symmetric_data, split == TRUE)
test_set = subset(symmetric_data, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(cor,true_int, true_int_sign, CI_lower95, CI_upper95)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(cor,true_int, true_int_sign, CI_lower95, CI_upper95)
test_set_s$.id <- NULL

# logistic regression model training set
train_logit <- vglm(true_int_sign ~ cor, data=training_set_s, family="multinomial")
summary(train_logit)


pred <-  predict(train_logit, newdata = test_set_s, type = "response") 
predictions <- apply(pred, 1, which.max)
predictions[which(predictions=="1")] <- "-ve"
predictions[which(predictions=="2")] <- "+ve"
predictions[which(predictions=="3")] <- "N"

table(test_set_s$true_int_sign, predictions)

ggplot(aes(y=cor,x=true_int_sign, colour=pred), data=test_set_s) + geom_point()

# significant testing accuracy 
test_set_s$temp <- test_set_s$CI_lower95*test_set_s$CI_upper95
test_set_s$significance <- NA
test_set_s[test_set_s$temp >0 & test_set_s$cor > 0,]$significance <- "+ve"
test_set_s[test_set_s$temp >0 & test_set_s$cor < 0,]$significance <- "-ve"
test_set_s[test_set_s$temp <=0,]$significance <- "N"

table(test_set_s$true_int_sign, test_set_s$significance)


### Linear discriminant analysis ###
library(MASS)
fit <- lda(true_int_sign~cor, data=training_set_s)
summary(fit)
predictions <- predict(fit, test_set_s)$class

table(test_set_s$true_int_sign,predictions)


#---- Granger analysis ------#
symmetric_data <- granger_res
library(caTools)
set.seed(2908)
split = sample.split(symmetric_data$true_int_sign, SplitRatio = 0.75)
training_set = subset(symmetric_data, split == TRUE)
test_set = subset(symmetric_data, split == FALSE)

# keep just the nesscary parts of the training and testing sets
training_set_s <- training_set %>% dplyr::select(logRSS, direction, true_int, true_int_sign, granger_p)
training_set_s$.id <- NULL
test_set_s <- test_set %>% dplyr::select(logRSS, direction, true_int, true_int_sign, granger_p)
test_set_s$.id <- NULL

# multinomial regression model training set
train_logit <- vglm(true_int_sign ~ logRSS, data=training_set_s, family="multinomial", trace=TRUE)
summary(train_logit)

pred <-  predict(train_logit, newdata = test_set_s, type = "response") 
predictions <- apply(pred, 1, which.max)
predictions[which(predictions=="1")] <- "-ve"
predictions[which(predictions=="2")] <- "+ve"
predictions[which(predictions=="3")] <- "N"
test_set_s$predictions <- predictions

table(test_set_s$true_int_sign, predictions, test_set_s$direction)
prop.table(table(test_set_s$true_int_sign, predictions))*100


# looking at the accuracy of significance testing
test_set_s$significance <- NA
test_set_s[test_set_s$granger_p <= 0.05 & test_set_s$logRSS > 0,]$significance <- "+ve"
test_set_s[test_set_s$granger_p <= 0.05 & test_set_s$logRSS < 0,]$significance <- "-ve"
test_set_s[test_set_s$granger_p > 0.05,]$significance <- "N"

table(test_set_s$true_int_sign, test_set_s$significance, test_set_s$direction)

# plotting
ggplot(aes(y=logRSS,x=true_int_sign, colour=predictions), data=test_set_s) + geom_point() + facet_grid(.~direction)
ggplot(aes(y=logRSS,x=true_int_sign, colour=significance), data=test_set_s) + geom_point() + facet_grid(.~direction)
