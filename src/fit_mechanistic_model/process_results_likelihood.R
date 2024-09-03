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



