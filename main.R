##################################################################################################################
# R code to run pomp model and each of the methods
#
# This code runs the following methods:
#  - Pearson's correlation coefficient
#  - GAM estimated correlation coefficient
#  - Transfer entropy
#  - Granger causality analysis 
#  - convergent cross mapping
# 
# Created by: Sarah Pirikahu 
# Creation date: 28 Feb 2024
##################################################################################################################

# Setup

#---- load libraries ----#
library(tidyverse)
library(testthat)
library(pomp)
library(lubridate)
# library(janitor)
# library(ggfortify)
# library(future) # allows for parallel processing
# library(foreach)
# library(doParallel)

#---- read in relevant functions ----#
source('seitr_x_seitr.R')

#---- set up cluster inputs ---# 
# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # based on array size
jobid <- 1

#---- set global parameters ----#
n_sim <- 100 # total number of simulated datasets
tot_weeks <- 626 # number of weeks to simulate
debug_bool <- FALSE

#---- get interaction parameter values ----#
theta_lambda1 <- c(0, 0.25, 0.5, 1, 2, 4)
theta_lambda2 <- c(0, 0.25, 0.5, 1, 2, 4)
delta1 <- c(1, 1/4, 1/13) # 1 week, 1 month, 3 months
delta2 <- c(1, 1/4, 1/13)

int_params <- expand.grid(theta_lambda1, theta_lambda2, delta1, delta2) %>%
  as_tibble() %>%
  rename('theta_lambda1' = 'Var1',
         'theta_lambda2' = 'Var2',
         'delta1' = 'Var3',
         'delta2' = 'Var4') %>%
  filter(theta_lambda1 == theta_lambda2,
         delta1 == delta2)
rm(theta_lambda1, theta_lambda2, delta1, delta2)

#---- generate timing of surges in immunity loss ----#
set.seed(1234)

# tot_seasons <- round((tot_weeks / 52) - 2) # convert to seasons

n_surge <- round(tot_weeks / 52) # number of surges
mu_Imloss <- 38 # average surge occurring in mid Oct
sd_Imloss <- 4 # standard deviation of 4 weeks

t_si <- rnorm(n = n_surge, mean = mu_Imloss, sd = sd_Imloss) # draw from normal dist

t_si <- t_si - 26 + seq(0, 52 * (n_surge - 1), by = 52)
# t_si <- round(t_si + seq(0, 52 * (n_surge - 1), by = 52))
# seq(26, 52 * n_surge, by = 52) + t_si

t_si <- round(t_si) # make whole numbers

t_si <- t_si[-which(t_si <= 104)] # remove first two years to allow system to reach equilibrium
n_surge <- length(t_si)

w_delta_i <- runif(n = length(t_si), min = 0.01 * 7, max = 0.1 * 7) # yearly rate of immunity loss

rm(mu_Imloss, sd_Imloss)

#---- set all true parameter values ----# 
true_int_params <- int_params[jobid, ]

true_params <- c(Ri1 = 1.1, Ri2 = 1.7,
                 sigma1 = 7, sigma2 = 7/5,
                 gamma1 = 7/5, gamma2 = 7/10,
                 w1 = 1/52, w2 = 1/52,
                 mu = 0.0002, nu = 0.0002,
                 rho1 = 0.002, rho2 = 0.002,
                 theta_lambda1 = true_int_params$theta_lambda1,
                 theta_lambda2 = true_int_params$theta_lambda2,
                 delta1 = true_int_params$delta1,
                 delta2 = true_int_params$delta2,
                 A1=0.2, phi1=26,
                 A2=0.2, phi2=26,
                 beta_sd1 = 0, beta_sd2 = 0, 
                 N = 3700000,
                 E01 = 0.001, E02 = 0.001,
                 R01 = 0.4, R02 = 0.25, R012 = 0.001,
                 nsurges = n_surge,
                 t_si_ = t(t_si), w_delta_i_ = t(w_delta_i))

##################################################################################################################

# Create model and synthetic data

#---- create pomp model object ----#
source('seitr_x_seitr.R')
resp_mod <- create_SEITRxSEITR_mod(tot_weeks, true_params, debug_bool = debug_bool)

#---- test pomp model ----#
check_transformations(resp_mod) # check parameter transformations
expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N')))) # check initial conditions
check_correct_N_CONST(resp_mod, unname(coef(resp_mod, 'N'))) # check constant population size
p_indep <- check_independent_dynamics(resp_mod) # check for independent dynamics
if (debug_bool) print(p_indep)

#---- simulate synthetic data ----#
dat <- simulate(resp_mod, nsim = n_sim, format = 'data.frame')

if (debug_bool) {
  resp_mod@data <- dat %>%
    filter(.id == 1) %>%
    select(V1_obs:V2_obs) %>%
    t()
  ll <- logLik(traj_objfun(data = resp_mod)) # check measurement density model
  print(ll)
}

dat <- dat %>%
  filter(time > 104) # remove first 2 years before simulation at equilibrium

dat <- dat %>%
  mutate(date = ymd('2012-July-01') + weeks(time)) # add dates

dat_red <- dat %>% # remove if outbreak never takes off
  group_by(.id) %>%
  mutate(sum_V1 = sum(V1_obs),
         sum_V2 = sum(V2_obs)) %>%
  filter(sum_V1 > 0 & sum_V2 > 0) %>%
  select(-c(sum_V1:sum_V2)) %>%
  ungroup()
expect_true(all.equal(dim(dat), dim(dat_red)))
rm(dat_red)

if (debug_bool) { # check seasonal attack rates
  season_breaks <- dat %>% filter(str_detect(date, '07-0[1-7]')) %>% pull(date) %>% unique()
  season_breaks <- c(season_breaks, '2024-07-01')
  
  dat %>%
    mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
    group_by(season, .id) %>%
    summarise(V1 = sum(V1),
              V2 = sum(V2)) %>%
    mutate(V1 = V1 / 3700000 * 100,
           V2 = V2 / 3700000 * 100) %>%
    summary()
}

if (debug_bool) { # check that surges in immunity happen correctly
  dat %>%
    mutate(S1 = X_SS + X_SE + X_SI + X_ST + X_SR) %>%
    select(time:.id, S1) %>%
    ggplot(aes(x = time, y = S1, group = .id)) + geom_line() + theme_classic()
}

#---- set up list to store all results ----#
results <- vector(mode = 'list', length = 7)
names(results) <- c("true_param", "data", "cor", "gam_cor", "transfer_entropy", "CCM","granger")

results$true_param <- true_params
results$data <- dat

##################################################################################################################

# Run various methods to determine causality



#---- Correlation coefficents --------# 
# function to estimate correlation 
# inputs: v1_obs = time series for v1_obs
#         v2_obs = time series for v2_obs
# corr_func <- function(v1_obs, v2_obs){
#   # calculated correlation coefficent
#   cor_raw <- cor.test(v1_obs, v2_obs)
#   # pull together results in data frame
#   temp_res <- data.frame(cbind(as.numeric(cor_raw$estimate), cor_raw$conf.int[1], cor_raw$conf.int[2]), cor_raw$p.value)
#   names(temp_res) <- c("cor", "CI_lower_95", "CI_upper_95", "p_value")
#   return(temp_res)
# }
# 
# # apply correlation function to all simulated datasets and save results
# results$cor <- results$data %>% group_by(.id) %>% do((corr_func(.$v1_obs,.$v2_obs)))

#--- GAM approach ---# 
# source("./methods/gam_cor.R")

# setting up parallelism for the foreach loop
# registerDoParallel(cl <- makeCluster(10))
# # apply the GAM correlation approach to each simulated data set and save the results
# res_gam_cor <- foreach(i=1:nsim, .packages=c("tidyverse","mgcv","vars","boot")) %dopar%{
#   # if the dataset was removed because the outbreak died out then skip it
#   if(dim(results$data %>% filter(.id==i))[1]!=0){
#     results$data %>% filter(.id==i) %>% gam_cor(.)
#   }
# }
# results$gam_cor <- do.call(rbind, res_gam_cor)

#----- Transfer entropy analysis ------# 
# source("./methods/transfer_entropy_jidt.R")
# 
# # lag = 1
# results$transfer_entropy <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="1"))
# # lag = 2
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="2"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)
# # lag = 4
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="4"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)
# # lag = 6
# temp <- results$data %>% group_by(.id) %>% do(te_jidt(., lag="6"))
# results$transfer_entropy <- rbind(results$transfer_entropy, temp)

#---- Granger causality analysis  ----# 
# source("./methods/granger_analysis.R")
# 
# # apply granger analysis to each simulated data set and save the results
# results$granger <- results$data %>% group_by(.id) %>% do(granger_func(.))

#------- Convergent Cross mapping analysis -------# 
# source("./methods/CCM.R")
# 
# results$CCM <- results$data %>% group_by(.id) %>% do(ccm_func(.))
# 
# # save out results
# save(results, file=sprintf('results_%s_%s_%s_%s_%s.RData',jobid, theta_lambda1,theta_lambda2,delta_1,delta_2))  
