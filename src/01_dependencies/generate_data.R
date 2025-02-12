# ---------------------------------------------------------------------------------------------------------------------
# Code to set relevant model parameters and generate synthetic data
# ---------------------------------------------------------------------------------------------------------------------

# Setup

#--- load libraries ----#
library(tidyverse)
library(testthat)
library(pomp)
library(lubridate)
library(gridExtra)

#---- read in relevant functions ----#
source('src/01_dependencies/seitr_x_seitr.R')

# ---------------------------------------------------------------------------------------------------------------------

# Set all relevant model parameters

#---- set global parameters ----#
n_sim <- 100 # total number of simulated datasets
tot_weeks <- 522 # number of weeks to simulate
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
         delta1 == delta2) %>%
  filter(!(theta_lambda1 == 1 & theta_lambda2 == 1 & delta1 < 1))
rm(theta_lambda1, theta_lambda2, delta1, delta2)

#---- generate timing of surges in immunity loss ----#
set.seed(1234)

n_surge <- round(tot_weeks / 52) # number of surges
mu_Imloss <- 12 # average surge occurring 12 weeks into season
sd_Imloss <- 4 # standard deviation of 4 weeks

t_si_mat <- matrix(nrow = n_surge, ncol = n_sim)
w_delta_i_mat <- matrix(nrow = n_surge, ncol = n_sim)

for (i in 1:n_sim) {
  
  t_si <- rnorm(n = n_surge, mean = mu_Imloss, sd = sd_Imloss) # draw from normal dist
  t_si <- t_si + seq(1, 53 * (n_surge - 1), by = 52)
  t_si <- round(t_si) # make whole numbers
  
  t_si[2:10] <- t_si[2:10] + 1 # account for years with 53 weeks
  t_si[8:10] <- t_si[8:10] + 1 # account for years with 53 weeks
  
  # w_delta_i <- runif(n = length(t_si), min = 0.005 * 7, max = 0.075 * 7) # yearly surge in rate of immunity loss
  w_delta_i <- runif(n = length(t_si), min = 0.05, max = 0.3) # yearly surge in rate of immunity loss
  
  t_si_mat[, i] <- t_si
  w_delta_i_mat[, i] <- w_delta_i
  
}

rm(mu_Imloss, sd_Imloss)

#---- generate range of values for Ri1/Ri2/w2/R02 ----#
r_eff_vals <- sobol_design(lower = setNames(c(1.05, 1.6, 1/(52.25 * 1.0)), c('Ri1', 'Ri2', 'w2')),
                           upper = setNames(c(1.4, 2.0, 1/(52.25 * 0.6)), c('Ri1', 'Ri2', 'w2')),
                           nseq = n_sim)

#---- set all true parameter values ----#
true_int_params <- int_params[jobid, ]

true_params_init <- c(Ri1 = r_eff_vals[1, 1], Ri2 = r_eff_vals[1, 2],
                      sigma1 = 7, sigma2 = 7/5,
                      gamma1 = 7/5, gamma2 = 7/10,
                      w1 = 1/52.25, w2 = 1/52.25,
                      mu = 0.0002, nu = 0.0002,
                      rho1 = 0.15, rho2 = 0.05,
                      theta_lambda1 = true_int_params$theta_lambda1,
                      theta_lambda2 = true_int_params$theta_lambda2,
                      delta1 = true_int_params$delta1,
                      delta2 = true_int_params$delta2,
                      A1=0.20, phi1=26,
                      A2=0.20, phi2=26,
                      k1 = 0.04, k2 = 0.02,
                      beta_sd1 = 0.1 * 0.1, beta_sd2 = 0.05 * 0.1,
                      N = 5000000,
                      E01 = 0.001, E02 = 0.001,
                      R01 = 0.30, R02 = 0.25, R012 = 0.001,
                      nsurges = n_surge,
                      t_si_ = t_si_mat[, 1], w_delta_i_ = w_delta_i_mat[, 1])

true_params <- parmat(true_params_init, nrep = n_sim)

true_params['Ri1', ] <- r_eff_vals[, 1]
true_params['Ri2', ] <- r_eff_vals[, 2]
true_params['w2', ] <- r_eff_vals[, 3]

true_params[str_detect(rownames(true_params), 't_si_'), ] <- t_si_mat
true_params[str_detect(rownames(true_params), 'w_delta_i_'), ] <- w_delta_i_mat

# ---------------------------------------------------------------------------------------------------------------------

# Create model and synthetic data

#---- create pomp model object ----#
resp_mod <- create_SEITRxSEITR_mod(tot_weeks, start_time = -1043, true_params_init, debug_bool = debug_bool)
# 25 years burn-in: -1304; 50 years burn-in: -2609

#---- test pomp model ----#
check_transformations(resp_mod) # check parameter transformations
expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N')))) # check initial conditions
check_correct_N_CONST(resp_mod, unname(coef(resp_mod, 'N'))) # check constant population size
p_indep <- check_independent_dynamics(resp_mod) # check for independent dynamics
if (debug_bool) print(p_indep)

#---- simulate synthetic data ----#
tic <- Sys.time()
dat <- simulate(resp_mod, params = true_params, nsim = 1, format = 'data.frame')
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'mins'
print(etime)

if (debug_bool) {
  resp_mod@data <- dat %>%
    filter(.id == 1) %>%
    select(V1_obs:V2_obs) %>%
    t()
  ll <- logLik(traj_objfun(data = resp_mod)) # check measurement density model
  print(ll)
  rm(ll)
}

# # Get weeks included in each season:
# dat <- dat %>%
#   mutate(date = ymd('2012-June-24') + weeks(time)) # add dates
# 
# season_breaks <- dat %>% filter(str_detect(date, '-07-0[1-7]')) %>% pull(date) %>% unique()
# season_breaks <- c(season_breaks, '2022-07-03')
# 
# dat <- dat %>%
#   as_tibble() %>%
#   mutate(season = cut(date, breaks = season_breaks, labels = 1:30, include.lowest = TRUE)) %>%
#   filter(!is.na(season))
# 
# t <- dat %>% group_by(season) %>% summarise(min_t = min(time), max_t = max(time)) %>% mutate(n = max_t - min_t, wk_12 = min_t + 12)
# write_csv(t, file = 'data/dates_and_seasons_UPDATE.csv')

dat <- dat %>%
  filter(time > 0) # remove burn-in period

dat <- dat %>%
  mutate(date = ymd('2012-June-24') + weeks(time)) # add dates

#---- check for yearly outbreaks ----#
season_breaks <- dat %>% filter(str_detect(date, '-07-0[1-7]')) %>% pull(date) %>% unique()
season_breaks <- c(season_breaks, '2022-07-03')

dat_red <- dat %>%
  mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
  group_by(season, .id) %>%
  summarise(V1 = sum(V1),
            V2 = sum(V2)) %>%
  filter(V1 < 300000)
# dat_red <- dat %>%
#   mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
#   group_by(season, .id) %>%
#   summarise(V1_obs = sum(V1_obs),
#             V2_obs = sum(V2_obs)) %>%
#   filter(V1_obs < 1000)

if (true_int_params$theta_lambda1 == 1.0 & true_int_params$theta_lambda2 == 1.0) {
  expect_true(nrow(dat_red) == 0)
}
rm(dat_red, season_breaks)

# ---------------------------------------------------------------------------------------------------------------------

# If debug_bool==TRUE, get various plots of the data

if (debug_bool) {
  
  #---- visualize outbreaks ----#
  dat %>%
    select(time:.id, V1_obs:V2_obs) %>%
    pivot_longer(V1_obs:V2_obs) %>%
    ggplot(aes(x = time, y = value, group = paste(name, .id))) +
    geom_line() +
    facet_wrap(~ name) +
    theme_classic()
  
  dat %>%
    select(time:.id, V1_obs:V2_obs) %>%
    pivot_longer(V1_obs:V2_obs) %>%
    ggplot(aes(x = time, y = value, group = paste(name, .id), color = name)) +
    geom_line() +
    theme_classic()
  
  dat %>%
    filter(.id %in% 1:12) %>%
    select(time:.id, V1_obs:V2_obs) %>%
    pivot_longer(V1_obs:V2_obs) %>%
    ggplot(aes(x = time, y = value, group = paste(name, .id), color = name)) +
    geom_line() +
    facet_wrap(~ .id) +
    theme_classic()
  
  #---- check that surges in immunity happen correctly ----#
  dat %>%
    filter(.id %in% 1:12) %>%
    mutate(S1 = X_SS + X_SE + X_SI + X_ST + X_SR) %>%
    select(time:.id, S1) %>%
    ggplot(aes(x = time, y = S1 / 5000000, group = .id)) + geom_line() + facet_wrap(~ .id) + theme_classic()
  
  #---- check seasonal attack rates ----#
  season_breaks <- dat %>% filter(str_detect(date, '-07-0[1-7]')) %>% pull(date) %>% unique()
  season_breaks <- c(season_breaks, '2022-07-03')
  
  attack_rates <- dat %>%
    mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
    group_by(season, .id) %>%
    summarise(V1 = sum(V1),
              V2 = sum(V2)) %>%
    mutate(V1 = V1 / 5000000 * 100,
           V2 = V2 / 5000000 * 100)
  print(summary(attack_rates))
  
  #---- check influence of parameter values on attack rates ----#
  params_df <- true_params %>%
    t() %>%
    as_tibble() %>%
    select('Ri1', 'Ri2', 'w2') %>%
    mutate(.id = 1:n_sim)
  
  attack_rates <- attack_rates %>%
    mutate(.id = as.integer(.id)) %>%
    inner_join(params_df, by = '.id')
  
  p_ar1 <- ggplot(data = attack_rates %>%
                    pivot_longer(Ri1:w2),
                  aes(x = value, y = V1)) +
    geom_point() +
    facet_wrap(~ name, scales = 'free_x', nrow = 1) +
    theme_classic()
  p_ar2 <- ggplot(data = attack_rates %>%
                    pivot_longer(Ri1:w2),
                  aes(x = value, y = V2)) +
    geom_point() +
    facet_wrap(~ name, scales = 'free_x', nrow = 1) +
    theme_classic()
  grid.arrange(p_ar1, p_ar2, ncol = 1)
  
  #---- check relative peak timing of flu vs. rsv ----#
  rsv_first <- dat %>%
    mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
    select(time:.id, season, V1_obs:V2_obs) %>%
    group_by(.id, season) %>%
    summarise(pt1 = which.max(V1_obs), pt2 = which.max(V2_obs)) %>%
    mutate(pt_diff = pt2 - pt1 < 0) %>%
    group_by(.id) %>%
    summarise(count = sum(pt_diff)) %>%
    filter(count >= 8) %>%
    pull(.id) %>%
    unique() %>%
    as.numeric()
  
  params_df <- params_df %>%
    mutate(rsv_first = .id %in% rsv_first)
  
  ggplot(data = params_df %>%
           pivot_longer(Ri1:w2),
         aes(x = rsv_first, y = value)) +
    geom_violin(fill = 'gray90') +
    facet_wrap(~ name, scales = 'free_y') +
    theme_classic()
  
}
