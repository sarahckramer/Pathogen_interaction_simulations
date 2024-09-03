# ---------------------------------------------------------------------------------------------------------------------
# Process results of maximum likelihood approach
# 
# Created by: Sarah Kramer
# Creation date: 07 August 2024
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load packages:
library(tidyverse)
library(testthat)
library(viridis)
library(gridExtra)

# Function to load results:
load_and_format_results <- function(res_files, get_top_5_perc = FALSE) {
  
  # # Get list of results files:
  # res_files <- list.files(path = filename, full.names = TRUE)
  
  # Read in all results:
  res_full <- list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Compile:
  res_full <- do.call('c', res_full)
  
  # Remove if fit resulted in an error:
  num_errors <- length(which(res_full == 'error'))
  if (num_errors > 0) {
    res_full <- res_full[-which(res_full == 'error')]
  }
  # print(num_errors)
  
  # Get parameter estimates and log-likelihoods:
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())# %>%
  # bind_cols('niter' = lapply(res_full, getElement, 'niter') %>%
  #             unlist())
  
  expect_true(nrow(pars_df) == (length(res_files) * 50) - num_errors)
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(no_best)
  
  # After first round, take top 10 results to get better starting parameters:
  if (get_top_5_perc) {
    
    if (no_best < 25) {
      no_best <- 25 # 25
    }
    
  }
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
  print(table(pars_top$message))
  
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # return(list(pars_df, pars_top))
  return(pars_top)
  
}

# ---------------------------------------------------------------------------------------------------------------------

# Get true parameter values

# Read in all synthetic data and true parameter values:
data_files <- list.files(path = 'results/', pattern = 'TRUE', full.names = TRUE)

true_params_LIST = data_LIST = vector('list', length= length(data_files))
for (i in 1:length(data_files)) {
  true_params_LIST[[i]] <- read_rds(data_files[i])$true_param
  data_LIST[[i]] <- read_rds(data_files[i])$data
}

# Limit to datasets used for fitting:
set.seed(93859)
ids_to_fit <- sample(1:100, size = 10)

true_params_LIST <- lapply(true_params_LIST, function(ix) {
  ix[, ids_to_fit]
})
data_LIST <- lapply(data_LIST, function(ix) {
  ix %>%
    filter(.id %in% ids_to_fit) %>%
    select(time:.id, date, V1:V2_obs) %>%
    as_tibble()
})

# Reorder:
where_run <- which(!is.na(as.numeric(str_split(data_files, '_')[[1]])))

true_params_LIST <- true_params_LIST[order(as.numeric(unlist(map(str_split(data_files, '_'), where_run))))]
data_LIST <- data_LIST[order(as.numeric(unlist(map(str_split(data_files, '_'), where_run))))]

# Format true params as tibble:
estpars <- c('Ri1', 'Ri2', 'rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'delta2',
             'A1', 'phi1', 'A2', 'phi2', 'k1', 'k2', 'E01', 'E02', 'R01', 'R02', 'R012',
             't_si_1', 't_si_2', 't_si_3', 't_si_4', 't_si_5', 't_si_6',
             't_si_7', 't_si_8', 't_si_9', 't_si_10',
             'w_delta_i_1', 'w_delta_i_2', 'w_delta_i_3', 'w_delta_i_4', 'w_delta_i_5', 'w_delta_i_6',
             'w_delta_i_7', 'w_delta_i_8', 'w_delta_i_9', 'w_delta_i_10')

names(true_params_LIST) <- 1:16
true_params <- lapply(1:length(true_params_LIST), function(ix) {
  true_params_LIST[[ix]][estpars, ] %>%
    t() %>%
    as_tibble() %>%
    rownames_to_column(var = '.id') %>%
    mutate(int_set = ix, .before = .id) %>%
    mutate(.id = ids_to_fit[as.numeric(.id)])
}) %>%
  bind_rows()

# Clean up:
rm(i, data_files, where_run)

# ---------------------------------------------------------------------------------------------------------------------

# Load and process results (round 1)

# Read in results:
res_dir <- 'results/trajectory_matching/round_1/'

# Loop through interaction parameter sets:
ids_to_fit <- ids_to_fit[1:5] # TEMPORARY

res_LIST = vector('list', length = 16)
for (int_set in 1:16) {
  res_list_TEMP = vector('list', length = length(ids_to_fit))
  
  # Loop through synthetic datasets:
  for (i in seq_along(ids_to_fit)) {
    data_id <- ids_to_fit[i]
    pat_temp <- paste(int_set, data_id, sep = '_')
    
    res_list_TEMP[[i]] <- load_and_format_results(list.files(path = res_dir, pattern = pat_temp, full.names = TRUE), get_top_5_perc = TRUE) %>%
      mutate(int_set = int_set, .id = data_id, .before = 'Ri1')
    rm(data_id, pat_temp)
    
  }
  
  res_LIST[[int_set]] <- bind_rows(res_list_TEMP)
  rm(res_list_TEMP)
}
rm(int_set, i)

# Compile all results:
res <- bind_rows(res_LIST)
rm(res_LIST, res_dir)

# Set any unrealistic values to NA:
res$delta1[res$delta1 > 7.0] <- NA
res$delta2[res$delta2 > 7.0] <- NA

res$rho1[res$rho1 == 1.0] <- NA
res$rho2[res$rho2 == 1.0] <- NA

# res <- res %>%
#   mutate(across(contains('w_delta'), ~ if_else(.x == 1.0, NA, .x)))

# Since phi=0 is equivalent to phi=52.25, don't use full range; transform so that we can select from only best-supported range:
phi_ranges_orig <- res %>%
  group_by(int_set, .id) %>%
  summarise(min1 = min(phi1), max1 = max(phi1), min2 = min(phi2), max2 = max(phi2)) %>%
  mutate(range1 = max1 - min1, range2 = max2 - min2) %>%
  select(int_set:.id, range1:range2)

phi_ranges_trans <- res %>%
  mutate(phi1 = if_else(phi1 < 1, phi1 + 52.25, phi1),
         phi2 = if_else(phi2 < 1, phi2 + 52.25, phi2)) %>%
  group_by(int_set, .id) %>%
  summarise(min1 = min(phi1), max1 = max(phi1), min2 = min(phi2), max2 = max(phi2)) %>%
  mutate(range1 = max1 - min1, range2 = max2 - min2) %>%
  select(int_set:.id, range1:range2)

range_reduct <- c(phi_ranges_orig %>%
                    inner_join(phi_ranges_trans, by = c('int_set', '.id')) %>%
                    filter(range1.y < range1.x) %>%
                    nrow(),
                  phi_ranges_orig %>%
                    inner_join(phi_ranges_trans, by = c('int_set', '.id')) %>%
                    filter(range2.y < range2.x) %>%
                    nrow())

if (any(range_reduct > 0)) {
  res <- res %>%
    mutate(phi1 = if_else(phi1 < 1, phi1 + 52.25, phi1),
           phi2 = if_else(phi2 < 1, phi2 + 52.25, phi2))
}

rm(phi_ranges_orig, phi_ranges_trans, range_reduct)

# ---------------------------------------------------------------------------------------------------------------------

# Get start ranges for round 2

# Drop column for log-likelihoods:
res <- res %>%
  select(-loglik)

# List results by int_set/.id:
res_list <- res %>% group_split(int_set, .id)

# Get minimum and maximum start values:
res_list <- lapply(res_list, function(ix) {
  as.data.frame(rbind(summarise(ix %>% group_by(int_set, .id), across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                      summarise(ix %>% group_by(int_set, .id), across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))
})

# Any parameters where minimum and maximum are equal?:
lapply(res_list, function(ix) {
  
  no_range <- c()
  for (i in 3:ncol(ix)) {
    if (identical(ix[1, i], ix[2, i])) {
      no_range <- c(no_range, i)
    }
  }
  
  return(length(no_range))
  
}) %>%
  unlist() %>%
  any(. > 0)

# Reset initial conditions so that sum cannot be >1:
ci_start <- lapply(res_list, function(ix) {
  
  # init_cond_estpars <- c('E01', 'E02', 'R01', 'R02', 'R012')
  
  ix %>%
    mutate(minmax = c('min', 'max')) %>%
    # select(int_set, .id, contains(init_cond_estpars), minmax) %>%
    mutate(sum = E01 + E02 + R01 + R02 + R012)
  
}) %>%
  bind_rows()
rm(res_list)

expect_true(ci_start %>%
              filter(minmax == 'min', sum > 1.0) %>%
              nrow() == 0)

ci_start <- ci_start %>%
  filter(minmax == 'max', sum > 1.0) %>%
  group_by(int_set, .id) %>%
  mutate(across(c('E01', 'E02', 'R01', 'R02', 'R012'),
                ~ .x - ((sum - 0.9999999) * (.x / sum)))) %>%
  bind_rows(ci_start %>%
              filter(minmax == 'min')) %>%
  mutate(minmax = factor(minmax, levels = c('min', 'max'))) %>%
  arrange(minmax, .by_group = TRUE)

# Ensure upper bounds still greater than lower:
expect_true(
  ci_start %>%
    filter(minmax == 'max') %>%
    select(int_set:.id, E01:R012) %>%
    inner_join(ci_start %>%
                 filter(minmax == 'min') %>%
                 select(int_set:.id, E01:R012),
               by = c('int_set', '.id')) %>%
    mutate(E01 = E01.x > E01.y,
           E02 = E02.x > E02.y,
           R01 = R01.x > R01.y,
           R02 = R02.x > R02.y,
           R012 = R012.x > R012.y) %>%
    ungroup() %>%
    select(E01:R012) %>%
    all(TRUE)
)

# Ensure that upper bounds now sum to <1:
expect_equal(ci_start %>%
               mutate(sum = E01 + E02 + R01 + R02 + R012) %>%
               filter(sum > 1.0) %>%
               nrow(), 0)

# Remove 'sum' column:
ci_start <- ci_start %>%
  select(!sum)

# If w_delta_i range greater than c(0, 0.5), reset:
ci_start <- ci_start %>%
  mutate(w_delta_i_1 = if_else(minmax == 'max' & w_delta_i_1 > 0.5, 0.5, w_delta_i_1),
         w_delta_i_2 = if_else(minmax == 'max' & w_delta_i_2 > 0.5, 0.5, w_delta_i_2),
         w_delta_i_3 = if_else(minmax == 'max' & w_delta_i_3 > 0.5, 0.5, w_delta_i_3),
         w_delta_i_4 = if_else(minmax == 'max' & w_delta_i_4 > 0.5, 0.5, w_delta_i_4),
         w_delta_i_5 = if_else(minmax == 'max' & w_delta_i_5 > 0.5, 0.5, w_delta_i_5),
         w_delta_i_6 = if_else(minmax == 'max' & w_delta_i_6 > 0.5, 0.5, w_delta_i_6),
         w_delta_i_7 = if_else(minmax == 'max' & w_delta_i_7 > 0.5, 0.5, w_delta_i_7),
         w_delta_i_8 = if_else(minmax == 'max' & w_delta_i_8 > 0.5, 0.5, w_delta_i_8),
         w_delta_i_9 = if_else(minmax == 'max' & w_delta_i_9 > 0.5, 0.5, w_delta_i_9),
         w_delta_i_10 = if_else(minmax == 'max' & w_delta_i_10 > 0.5, 0.5, w_delta_i_10))

# Write start ranges to file:
write_rds(ci_start, file = 'results/trajectory_matching/round2CI_startvals.rds')
rm(ci_start)

# ---------------------------------------------------------------------------------------------------------------------

# TEMPORARY: Check round 1 accuracy

# Join information on true parameter values:
res <- res %>%
  pivot_longer(Ri1:w_delta_i_10, names_to = 'param') %>%
  inner_join(true_params %>% pivot_longer(Ri1:w_delta_i_10, names_to = 'param', values_to = 'truth'),
             by = c('int_set', '.id', 'param'))

# And also get true interaction parameters for each run:
res <- res %>%
  left_join(true_params %>%
              select(int_set, theta_lambda1, delta1) %>%
              distinct() %>%
              rename('true_theta_lambda' = 'theta_lambda1',
                     'true_delta' = 'delta1') %>%
              mutate(true_delta = 7 / true_delta)
  )

# ---------------------------------------------------------------------------------------------------------------------

# Plot results

# Interaction parameters:
p1 <- ggplot(res %>%
               filter(param == 'theta_lambda1') %>%
               mutate(true_delta = factor(true_delta, levels = c(7, 28, 91)),
                      x_use = case_when(int_set == 7 ~ 2, int_set == 12 ~ 3, int_set == 2 ~ 5, int_set == 8 ~ 6, int_set == 13 ~ 7,
                                        int_set == 3 ~ 9, int_set == 9 ~ 10, int_set == 14 ~ 11, int_set == 4 ~ 13, int_set == 5 ~ 15,
                                        int_set == 10 ~ 16, int_set == 15 ~ 17, int_set == 6 ~ 19, int_set == 11 ~ 20, int_set == 16 ~ 21,
                                        .default = int_set))) +
  geom_boxplot(aes(x = x_use, y = value, fill = true_delta, group = paste(true_delta, x_use))) +
  geom_segment(aes(x = x_use - 0.5, xend = x_use + 0.5, y = truth, yend = truth), lty = 2, linewidth = 0.6) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  facet_wrap(~ .id, ncol = 1) +
  scale_fill_manual(values = viridis(n = 3, option = 'mako')[1:3]) +
  scale_x_continuous(breaks = c(2, 6, 10, 13, 16, 20), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(transform = 'sqrt', breaks = c(0, 0.5, 1, 2, 4, 6, 10)) +
  labs(x = expression('True ' * theta[lambda]), y = 'Fit Value', fill = 'True Duration', title = expression(theta[lambda*1]))

p2 <- ggplot(res %>%
               filter(param == 'delta1') %>%
               mutate(truth_recode = case_when(truth == 1 ~ 1, truth == 7 / 28 ~ 2, truth == 7 / 91 ~ 3),
                      true_theta_lambda = if_else(true_theta_lambda == 0, 0.1, true_theta_lambda))) +
  geom_boxplot(aes(x = truth_recode, y = 7 / value, fill = true_theta_lambda, group = paste(true_theta_lambda, truth_recode))) +
  geom_segment(aes(x = truth_recode - 0.5, xend = truth_recode + 0.5, y = 7 / truth, yend = 7 / truth), lty = 2) +#, linewidth = 1.0) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  facet_wrap(~ .id, ncol = 1) +
  scale_fill_viridis(option = 'mako', trans = 'log', breaks = c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_x_continuous(breaks = 1:3, labels = c(7, 28, 91)) +
  scale_y_continuous(transform = 'sqrt', breaks = c(0, 10, 50, 100, 200)) +
  labs(x = 'True Duration', y = 'Fit Value', fill = expression('True ' * theta[lambda] * '   '), title = expression(delta[1]))

p3 <- ggplot(res %>%
               filter(param == 'theta_lambda2') %>%
               mutate(true_delta = factor(true_delta, levels = c(7, 28, 91)),
                      x_use = case_when(int_set == 7 ~ 2, int_set == 12 ~ 3, int_set == 2 ~ 5, int_set == 8 ~ 6, int_set == 13 ~ 7,
                                        int_set == 3 ~ 9, int_set == 9 ~ 10, int_set == 14 ~ 11, int_set == 4 ~ 13, int_set == 5 ~ 15,
                                        int_set == 10 ~ 16, int_set == 15 ~ 17, int_set == 6 ~ 19, int_set == 11 ~ 20, int_set == 16 ~ 21,
                                        .default = int_set))) +
  geom_boxplot(aes(x = x_use, y = value, fill = true_delta, group = paste(true_delta, x_use))) +
  geom_segment(aes(x = x_use - 0.5, xend = x_use + 0.5, y = truth, yend = truth), lty = 2, linewidth = 0.6) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  facet_wrap(~ .id, ncol = 1) +
  scale_fill_manual(values = viridis(n = 3, option = 'mako')[1:3]) +
  scale_x_continuous(breaks = c(2, 6, 10, 13, 16, 20), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(transform = 'sqrt', breaks = c(0, 0.5, 1, 2, 4, 6, 10)) +
  labs(x = expression('True ' * theta[lambda]), y = 'Fit Value', fill = 'True Duration', title = expression(theta[lambda*2]))

p4 <- ggplot(res %>%
               filter(param == 'delta2') %>%
               mutate(truth_recode = case_when(truth == 1 ~ 1, truth == 7 / 28 ~ 2, truth == 7 / 91 ~ 3),
                      true_theta_lambda = if_else(true_theta_lambda == 0, 0.1, true_theta_lambda))) +
  geom_boxplot(aes(x = truth_recode, y = 7 / value, fill = true_theta_lambda, group = paste(true_theta_lambda, truth_recode))) +
  geom_segment(aes(x = truth_recode - 0.5, xend = truth_recode + 0.5, y = 7 / truth, yend = 7 / truth), lty = 2) +#, linewidth = 1.0) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  facet_wrap(~ .id, ncol = 1) +
  scale_fill_viridis(option = 'mako', trans = 'log', breaks = c(0.1, 0.25, 0.5, 1.0, 2.0, 4.0), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_x_continuous(breaks = 1:3, labels = c(7, 28, 91)) +
  scale_y_continuous(transform = 'sqrt', breaks = c(0, 10, 50, 100, 200, 500)) +
  labs(x = 'True Duration', y = 'Fit Value', fill = expression('True ' * theta[lambda] * '   '), title = expression(delta[2]))

grid.arrange(p1, p2, p3, p4, ncol = 2)

# Other parameters:

ggplot(res %>% filter(param %in% c('Ri1', 'Ri2'))) +
  geom_violin(aes(x = as.character(int_set), y = value, group = int_set)) +
  geom_hline(aes(yintercept = truth), lty = 2) +
  facet_grid(.id ~ param, scales = 'free_y') +
  theme_classic()

ggplot(res %>% filter(param %in% c('E01', 'E02', 'R01', 'R02', 'R012', 'A1', 'A2', 'phi1', 'phi2', 'rho1', 'rho2', 'k1', 'k2'))) +
  geom_violin(aes(x = as.character(.id), y = value, group = paste(int_set, .id), fill = int_set)) +
  geom_hline(aes(yintercept = truth), lty = 2) +
  facet_wrap(~ param, scales = 'free_y', ncol = 2) +
  theme_classic() +
  scale_fill_viridis()

ggplot(res %>% filter(str_detect(param, 't_si'))) +
  geom_violin(aes(x = as.character(int_set), y = value, group = int_set), fill = 'gray90') +
  geom_hline(aes(yintercept = truth), lty = 2) +
  facet_grid(param ~ .id, scales = 'free_y') +
  theme_classic()

ggplot(res %>% filter(str_detect(param, 'w_delta'))) +
  geom_violin(aes(x = as.character(int_set), y = value, group = int_set), fill = 'gray90') +
  geom_hline(aes(yintercept = truth), lty = 2) +
  facet_grid(param ~ .id, scales = 'free_y') +
  theme_classic()
