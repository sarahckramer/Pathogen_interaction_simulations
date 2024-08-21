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
load_and_format_results <- function(res_files) {
  
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
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  
  # Remove where no convergence occurs:
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

# Load and process results

# Read in results:
res_dir <- 'results/trajectory_matching/'

# Loop through interaction parameter sets:
res_LIST = vector('list', length = 16)
for (int_set in 1:16) {
  res_list_TEMP = vector('list', length = length(ids_to_fit))
  
  # Loop through synthetic datasets:
  for (i in seq_along(ids_to_fit)) {
    data_id <- ids_to_fit[i]
    pat_temp <- paste(int_set, data_id, sep = '_')
    
    res_list_TEMP[[i]] <- load_and_format_results(list.files(path = res_dir, pattern = pat_temp, full.names = TRUE)) %>%
      mutate(int_set = int_set, .id = data_id, .before = 'Ri1')
    
  }
  
  res_LIST[[int_set]] <- bind_rows(res_list_TEMP)
  rm(res_list_TEMP)
}
rm(int_set, i)

# Compile all results:
res <- bind_rows(res_LIST)
rm(res_LIST)

# ---------------------------------------------------------------------------------------------------------------------

# Explore mif2 results:
mf <- read_rds('results/particle_filter/res_4_90_1_PARALLEL_mif_TEMP.rds')

m <- lapply(mf, function(ix) {
  ix@params %>% bind_rows() %>% bind_cols(loglik = ix@loglik)
})

estpars <- c('Ri1', 'Ri2', 'rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'delta2',
             'A1', 'phi1', 'A2', 'phi2', 'k1', 'k2', 'beta_sd1', 'beta_sd2',
             't_si_1', 't_si_2', 't_si_3', 't_si_4', 't_si_5', 't_si_6',
             't_si_7', 't_si_8', 't_si_9', 't_si_10',
             'w_delta_i_1', 'w_delta_i_2', 'w_delta_i_3', 'w_delta_i_4', 'w_delta_i_5', 'w_delta_i_6',
             'w_delta_i_7', 'w_delta_i_8', 'w_delta_i_9', 'w_delta_i_10')

m <- bind_rows(m) %>%
  select(all_of(estpars), loglik) %>%
  arrange(desc(loglik))

m_red <- m[1:20, ] %>%
  pivot_longer(-loglik, names_to = 'param')

ggplot(m_red) +
  geom_violin(aes(x = 1, y = value)) +
  # geom_point(aes(x = as.character(.id), y = truth), size = 1.5) +
  facet_wrap(~ param, scales = 'free_y') +
  theme_classic()

trace_df <- lapply(1:length(mf), function(ix) {
  mf[[ix]]@traces %>%
    as.data.frame() %>%
    as_tibble() %>%
    rownames_to_column(var = 'iteration') %>%
    mutate(iteration = as.numeric(iteration),
           .id = ix) %>%
    select(.id, iteration, loglik, all_of(estpars)) %>%
    pivot_longer(-c(iteration, .id, loglik), names_to = 'param')
}) %>%
  bind_rows()

ggplot(trace_df, aes(x = iteration, y = value, group = .id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ param, scales = 'free') +
  theme_classic()

ids_to_plot <- trace_df %>% filter(iteration == 50, loglik > -3500) %>% pull(.id) %>% unique()

ggplot(trace_df %>% filter(.id %in% ids_to_plot),
       aes(x = iteration, y = value, group = .id)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~ param, scales = 'free') +
  theme_classic()
