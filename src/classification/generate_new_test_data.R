# ---------------------------------------------------------------------------------------------------------------------
# Code to generate synthetic data using a wider range of parameter values, for more rigorous testing of
# machine learning methods
# ---------------------------------------------------------------------------------------------------------------------

# Setup

#--- load libraries ----#
library(tidyverse)
library(testthat)
library(pomp)
library(gridExtra)

#---- read in relevant functions ----#
source('src/seitr_x_seitr.R')

# ---------------------------------------------------------------------------------------------------------------------

# Set all relevant model parameters

#---- set global parameters ----#
n_sim <- 1000 # total number of simulated datasets
tot_weeks <- 626 # number of weeks to simulate
debug_bool <- FALSE

#---- generate timing of surges in immunity loss ----#
set.seed(1234)

n_surge <- round(tot_weeks / 52) # number of surges
mu_Imloss <- 12 # average surge occurring 12 weeks into season
sd_Imloss <- 4 # standard deviation of 4 weeks

t_si_mat <- matrix(nrow = n_surge - 2, ncol = n_sim)
w_delta_i_mat <- matrix(nrow = n_surge - 2, ncol = n_sim)

for (i in 1:n_sim) {
  
  t_si <- rnorm(n = n_surge, mean = mu_Imloss, sd = sd_Imloss) # draw from normal dist
  t_si <- t_si + seq(0, 52 * (n_surge - 1), by = 52)
  t_si <- round(t_si) # make whole numbers
  
  t_si[2:12] <- t_si[2:12] + 1 # account for years with 53 weeks
  t_si[8:12] <- t_si[8:12] + 1
  
  t_si <- t_si[-which(t_si <= 104)] # remove first two years to allow system to reach equilibrium
  w_delta_i <- runif(n = length(t_si), min = 0.005 * 7, max = 0.05 * 7) # yearly surge in rate of immunity loss
  
  t_si_mat[, i] <- t_si
  w_delta_i_mat[, i] <- w_delta_i
  
}

rm(mu_Imloss, sd_Imloss, i, t_si, w_delta_i)

#---- generate range of values for other parameters ----#
parm_vals <- sobol_design(lower = setNames(c(1.05, 1.05, 1/(52.25 * 1.0), log(1.0), 1.0, 0, 11, 0, 0),
                                           c('Ri1', 'Ri2', 'w2', 'theta_lambda', 'delta_recip', 'A', 'phi', 'R01', 'R02')),
                          upper = setNames(c(5.0, 5.0, 1/(52.25 * 0.6), log(10.0), 26.0, 0.5, 41, 0.45, 0.45),
                                           c('Ri1', 'Ri2', 'w2', 'theta_lambda', 'delta_recip', 'A', 'phi', 'R01', 'R02')),
                          nseq = n_sim)
parm_vals[, 'theta_lambda'] <- exp(parm_vals[, 'theta_lambda'])

set.seed(9023)
neg_int <- sample(1:1000, size = 500, replace = FALSE)
parm_vals[neg_int, 'theta_lambda'] <- 1 / parm_vals[neg_int, 'theta_lambda']
set.seed(4891)
null_int <- sample(1:1000, size = 100, replace = FALSE)
parm_vals[null_int, 'theta_lambda'] <- 1.0
rm(neg_int, null_int)

#---- set all true parameter values ----#
true_params_init <- c(Ri1 = 1.2, Ri2 = 1.8,
                      sigma1 = 7, sigma2 = 7/5,
                      gamma1 = 7/5, gamma2 = 7/10,
                      w1 = 1/52.25, w2 = 1/52.25,
                      mu = 0.0002, nu = 0.0002,
                      rho1 = 0.002, rho2 = 0.002,
                      theta_lambda1 = 1.0,
                      theta_lambda2 = 1.0,
                      delta1 = 1/4,
                      delta2 = 1/4,
                      A1=0.20, phi1=26,
                      A2=0.20, phi2=26,
                      k1 = 0.04, k2 = 0.02,
                      beta_sd1 = 0.1 * 0.1, beta_sd2 = 0.05 * 0.1,
                      N = 3700000,
                      E01 = 0.001, E02 = 0.001,
                      R01 = 0.40, R02 = 0.25, R012 = 0.001,
                      nsurges = n_surge - 2,
                      t_si_ = t_si_mat[, 1], w_delta_i_ = w_delta_i_mat[, 1])

true_params <- parmat(true_params_init, nrep = n_sim)

true_params['Ri1', ] <- parm_vals[, 1]
true_params['Ri2', ] <- parm_vals[, 2]
true_params['w2', ] <- parm_vals[, 3]
true_params['theta_lambda1', ] <- parm_vals[, 4]
true_params['theta_lambda2', ] <- parm_vals[, 4]
true_params['delta1', ] <- 1 / parm_vals[, 5]
true_params['delta2', ] <- 1 / parm_vals[, 5]
true_params['A1', ] <- parm_vals[, 6]
true_params['A2', ] <- parm_vals[, 6]
true_params['phi1', ] <- parm_vals[, 7]
true_params['phi2', ] <- parm_vals[, 7]
true_params['R01', ] <- parm_vals[, 8]
true_params['R02', ] <- parm_vals[, 9]

true_params[str_detect(rownames(true_params), 't_si_'), ] <- t_si_mat
true_params[str_detect(rownames(true_params), 'w_delta_i_'), ] <- w_delta_i_mat

# write_rds(true_params, file = 'data/ml_true_params_HIGHR0.rds')
rm(parm_vals, t_si_mat, w_delta_i_mat, n_surge)

# ---------------------------------------------------------------------------------------------------------------------

# Create model and synthetic data

#---- create pomp model object ----#
resp_mod <- create_SEITRxSEITR_mod(tot_weeks, true_params_init, debug_bool = debug_bool)
rm(tot_weeks, true_params_init)

#---- test pomp model ----#
check_transformations(resp_mod) # check parameter transformations
expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N')))) # check initial conditions
check_correct_N_CONST(resp_mod, unname(coef(resp_mod, 'N'))) # check constant population size
p_indep <- check_independent_dynamics(resp_mod) # check for independent dynamics
if (debug_bool) print(p_indep)
rm(p_indep)

#---- simulate synthetic data ----#
tic <- Sys.time()
dat <- bake(file = 'data/ml_test_data_HIGHR0.rds', {
  simulate(resp_mod, params = true_params, nsim = 1, format = 'data.frame')
})
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'secs'
print(etime)
rm(tic, toc, etime, resp_mod)

dat <- dat %>%
  filter(time > 104) # remove first 2 years before simulation at equilibrium

dat <- dat %>%
  mutate(date = ymd('2012-July-01') + weeks(time)) # add dates

#---- check for yearly outbreaks ----#
season_breaks <- dat %>% filter(str_detect(date, '07-0[1-7]')) %>% pull(date) %>% unique()
season_breaks <- c(season_breaks, '2024-07-01')

attack_rates <- dat %>%
  mutate(season = cut(date, breaks = season_breaks, labels = 1:10, include.lowest = TRUE)) %>%
  group_by(season, .id) %>%
  summarise(ar1 = sum(V1) / 3700000,
            ar2 = sum(V2) / 3700000) %>%
  ungroup()

#---- identify ids where there aren't outbreaks every year ----#
ids_to_remove <- attack_rates %>%
  filter(ar1 < 0.05 | ar2 < 0.05) %>%
  pull(.id) %>%
  unique() %>%
  sort()
print(length(ids_to_remove))

#---- remove ids where outbreaks don't occur yearly and recode id numbers ----#
new_ids <- c(1:1000)[!(1:1000 %in% ids_to_remove)] %>% bind_cols(new_id = 1:(1000 - length(ids_to_remove))) %>% rename('.id' = '...1')

dat_new <- dat %>%
  as_tibble() %>%
  filter(!(.id %in% ids_to_remove)) %>%
  mutate(.id = as.numeric(.id)) %>%
  left_join(new_ids, by = '.id') %>%
  mutate(.id = new_id) %>%
  select(-new_id)

true_params_new <- true_params[, new_ids$.id]

rm(new_ids)

#---- checks and plots ----#
if (debug_bool) {
  
  #---- visualize outbreaks ----#
  set.seed(49034)
  dat %>%
    filter(.id %in% sample(c(1:1000)[!(1:1000 %in% ids_to_remove)], size = 25, replace = FALSE)) %>%
    select(time:.id, V1_obs:V2_obs) %>%
    pivot_longer(V1_obs:V2_obs) %>%
    ggplot(aes(x = time, y = value, group = paste(name, .id), color = name)) +
    geom_line() +
    facet_wrap(~ .id) +
    theme_classic()
  
  set.seed(49034)
  dat %>%
    filter(.id %in% sample(c(1:1000)[1:1000 %in% ids_to_remove], size = 25, replace = FALSE)) %>%
    select(time:.id, V1_obs:V2_obs) %>%
    pivot_longer(V1_obs:V2_obs) %>%
    ggplot(aes(x = time, y = value, group = paste(name, .id), color = name)) +
    geom_line() +
    facet_wrap(~ .id) +
    theme_classic()
  
  #---- check influence of parameter values on attack rates ----#
  params_df <- true_params %>%
    t() %>%
    as_tibble() %>%
    select('Ri1', 'Ri2', 'w2', 'theta_lambda1', 'delta1', 'A1', 'phi1', 'R01', 'R02') %>%
    mutate(.id = 1:n_sim)
  
  attack_rates <- attack_rates %>%
    mutate(.id = as.integer(.id)) %>%
    inner_join(params_df, by = '.id')
  
  p_ar1 <- ggplot(data = attack_rates %>%
                    pivot_longer(Ri1:R02) %>%
                    mutate(not_yearly = .id %in% ids_to_remove),
                  aes(x = value, y = ar1, col = not_yearly)) +
    geom_point(size = 0.75) +
    facet_wrap(~ name, scales = 'free_x', nrow = 1) +
    theme_classic() +
    scale_color_manual(values = c('black', 'purple'))
  
  p_ar2 <- ggplot(data = attack_rates %>%
                    pivot_longer(Ri1:R02) %>%
                    mutate(not_yearly = .id %in% ids_to_remove),
                  aes(x = value, y = ar2, col = not_yearly)) +
    geom_point(size = 0.75) +
    facet_wrap(~ name, scales = 'free_x', nrow = 1) +
    theme_classic() +
    scale_color_manual(values = c('black', 'purple'))
  grid.arrange(p_ar1, p_ar2, ncol = 1)
  # not seeing any super clear patterns here; higher A1/A2, Ri1 tends to yield larger outbreaks
  
  attack_rates %>%
    mutate(delta1 = 1 / delta1) %>%
    pivot_longer(Ri1:R02) %>%
    mutate(not_yearly = .id %in% ids_to_remove) %>%
    ggplot(aes(x = not_yearly, y = value, group = not_yearly)) +
    geom_violin(fill = 'gray95') +
    facet_wrap(~ name, scales = 'free_y') +
    theme_classic()
  # simulations that don't have yearly outbreaks are slightly more likely to have negative interactions, higher A1/A2, more extreme phi1/phi2, lower Ri/R0
  
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
           pivot_longer(Ri1:R02),
         aes(x = rsv_first, y = value)) +
    geom_violin(fill = 'gray90') +
    facet_wrap(~ name, scales = 'free_y') +
    theme_classic()
  # only about 120 consistently have RSV first; about 350-550 consistently have flu first
  # rsv first tend to have longer duration of interaction, lower Ri1/higher Ri2, lower phi
  
}

# Clean up:
dat <- dat_new
true_params <- true_params_new

rm(attack_rates, season_breaks, ids_to_remove, true_params_new, dat_new, n_sim, debug_bool)
