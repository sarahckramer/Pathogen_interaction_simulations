# ------------------------------------------------------------------------------
# Code to extract and process all results
#
# Created by: Sarah Pirikahu
# Creation date: 13 March 2024
# ------------------------------------------------------------------------------

# Setup

# Load packages
library(tidyverse)
library(gridExtra)
library(testthat)
library(viridis)

# Open pdf to save plots:
date <- format(Sys.Date(), '%d%m%y')
# pdf(file = paste0('results/plots/plot_accuracy_by_method_', date, '_UPDATED.pdf'),
#     width = 16, height = 12)

# ------------------------------------------------------------------------------

# Read in all results

# Get file names:
res_filenames_T <- list.files(path = 'results/', pattern = 'TRUE', full.names = TRUE) # run locally
res_filenames_F <- list.files(path = 'results/', pattern = 'FALSE', full.names = TRUE) # run on cluster

# Read in results:
results_T = results_F = vector('list', length = length(res_filenames_T))

for (i in 1:length(res_filenames_T)) {
  results_T[[i]] <- read_rds(res_filenames_T[i])
  results_F[[i]] <- read_rds(res_filenames_F[i])
}
rm(i)

# Label each results list with run number:
where_run <- which(!is.na(as.numeric(str_split(res_filenames_T, '_')[[1]])))

names(results_T) <- unlist(map(str_split(res_filenames_T, '_'), where_run))
names(results_F) <- unlist(map(str_split(res_filenames_F, '_'), where_run))

rm(res_filenames_T, res_filenames_F, where_run)

# ------------------------------------------------------------------------------

# Extract true parameter values and data

# Get true parameter values:
int_params <- c('theta_lambda1', 'delta1')

res_trueparams <- lapply(1:length(results_T), function(ix) {
  results_T[[ix]]$true_param[int_params, 1]
}) %>%
  bind_rows() %>%
  rename('theta_lambda' = 'theta_lambda1',
         'delta' = 'delta1') %>%
  mutate(run = as.numeric(names(results_T))) %>%
  arrange(run)

res_trueparams_CHECK <- lapply(1:length(results_F), function(ix) {
  results_F[[ix]]$true_param[int_params, 1]
}) %>%
  bind_rows() %>%
  rename('theta_lambda' = 'theta_lambda1',
         'delta' = 'delta1') %>%
  mutate(run = as.numeric(names(results_T))) %>%
  arrange(run)

expect_true(all.equal(res_trueparams, res_trueparams_CHECK))
rm(res_trueparams_CHECK)

# Get data:
data_list = data_list_CHECK = vector('list', length = length(results_T))

for (i in 1:length(results_T)) {
  
  data_list[[i]] <- results_T[[i]]$data %>%
    select(time, date, .id, V1_obs, V2_obs) %>%
    mutate(run = as.numeric(names(results_T)[i]))
  
  data_list_CHECK[[i]] <- results_F[[i]]$data %>%
    select(time, date, .id, V1_obs, V2_obs) %>%
    mutate(run = as.numeric(names(results_F)[i]))
  
}

dat <- bind_rows(data_list) %>%
  arrange(run)
dat_CHECK <- bind_rows(data_list_CHECK) %>%
  arrange(run)
expect_true(all.equal(dat, dat_CHECK))

rm(data_list, data_list_CHECK, i, dat_CHECK, int_params)

# Join true parameter values:
dat <- dat %>%
  as_tibble() %>%
  inner_join(res_trueparams, by = 'run')

# ------------------------------------------------------------------------------

# Visualize some datasets

# Choose five random simulations to plot:
set.seed(490275)
ids_to_plot <- sample(1:100, size = 5)

# Format data for plotting:
dat_plot <- dat %>%
  filter(.id %in% ids_to_plot) %>%
  # filter((theta_lambda == 0 & delta == 0.25) | (theta_lambda == 4 & delta == 0.25) | theta_lambda == 1) %>%
  pivot_longer(V1_obs:V2_obs, names_to = 'virus', values_to = 'obs') %>%
  mutate(virus = if_else(virus == 'V1_obs', 'Influenza', 'RSV'))

# Plot:
p1.1 <- ggplot(data = dat_plot %>%
                 filter((theta_lambda == 0 & delta == 1) | (theta_lambda == 4 & delta == 1) | theta_lambda == 1),
               aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ theta_lambda, ncol = 1) +
  theme_classic() +
  labs(title = 'Duration 1 Week')
p1.2 <- ggplot(data = dat_plot %>%
                 filter((theta_lambda == 0 & delta == 0.25) | (theta_lambda == 4 & delta == 0.25) | theta_lambda == 1),
               aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ theta_lambda, ncol = 1) +
  theme_classic() +
  labs(title = 'Duration 4 Weeks')
p1.3 <- ggplot(data = dat_plot %>%
                 filter((theta_lambda == 0 & delta == 1/13) | (theta_lambda == 4 & delta == 1/13) | theta_lambda == 1),
               aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ theta_lambda, ncol = 1) +
  theme_classic() +
  labs(title = 'Duration 13 Weeks')

p2 <- ggplot(data = dat_plot %>%
               filter((theta_lambda == 0 & delta == 0.25) | (theta_lambda == 0.25 & delta == 0.25) | (theta_lambda == 0.5 & delta == 0.25) | (theta_lambda == 2 & delta == 0.25) | (theta_lambda == 4 & delta == 0.25) | theta_lambda == 1),
             aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ theta_lambda, ncol = 1) +
  theme_classic() +
  labs(title = 'Duration 4 Weeks')

p3.1 <- ggplot(data = dat_plot %>%
                 filter(theta_lambda == 0),
               aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ delta, ncol = 1) +
  theme_classic() +
  labs(title = 'theta_lambda = 0')
p3.2 <- ggplot(data = dat_plot %>%
                 filter(theta_lambda == 4),
               aes(x = date, y = obs, group = paste(virus, .id), color = virus)) +
  geom_line() +
  facet_wrap(~ delta, ncol = 1) +
  theme_classic() +
  labs(title = 'theta_lambda = 4')

print(p1.1)
print(p1.2)
print(p1.3)

# print(p2)

# print(p3.1)
# print(p3.2)

rm(p1.1, p1.2, p1.3, p2, p3.1, p3.2, dat_plot, ids_to_plot)

# ------------------------------------------------------------------------------

# Extract results for all statistical methods

# Correlation coefficients:
res_corr <- lapply(results_T, getElement, 'cor') %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.numeric(run)) %>%
  inner_join(res_trueparams, by = 'run')

# GAMs:
res_gam <- lapply(results_F, getElement, 'gam_cor') %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.numeric(run)) %>%
  inner_join(res_trueparams, by = 'run') %>%
  as_tibble()

# Granger causality:
res_granger <- lapply(results_T, getElement, 'granger') %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.numeric(run)) %>%
  inner_join(res_trueparams, by = 'run')

# Transfer entropy:
res_te <- lapply(results_T, getElement, 'transfer_entropy') %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.numeric(run)) %>%
  inner_join(res_trueparams, by = 'run')

if ('V1 -> V2' %in% unique(res_te$direction)) {
  res_te <- res_te %>%
    mutate(direction = if_else(direction == 'V1 -> V2', 'v1 -> v2', 'v2 -> v1'))
}

# CCM:
res_ccm <- lapply(results_F, getElement, 'CCM') %>%
  bind_rows(.id = 'run') %>%
  mutate(run = as.numeric(run)) %>%
  ungroup() %>%
  inner_join(res_trueparams, by = 'run')

res_ccm_surr <- res_ccm %>%
  filter(data == 'surr')
res_ccm <- res_ccm %>%
  filter(data == 'obs')

# Clean up:
rm(results_T, results_F, res_trueparams)

# ------------------------------------------------------------------------------

# Functions

# Function to calculate accuracy for all values of theta_lambda/delta:
calculate_accuracy_matrix <- function(df) {
  
  all_tl <- sort(unique(df$theta_lambda))
  all_delta <- sort(unique(1 / df$delta))
  
  mat_correct <- matrix(nrow = length(all_tl), ncol = length(all_delta))
  rownames(mat_correct) <- all_tl
  colnames(mat_correct) <- all_delta
  
  for (tl in all_tl) {
    
    if (tl != 1) {
      for (d in all_delta) {
        
        df_temp <- df %>% filter(theta_lambda == tl & delta == 1 / d)
        mat_correct[rownames(mat_correct) == tl, colnames(mat_correct) == d] <- (df_temp %>% filter(int_est == int_true) %>% nrow()) / nrow(df_temp)
        
      }
    } else {
      
      df_temp <- df %>% filter(theta_lambda == tl)
      mat_correct[rownames(mat_correct) == tl, ] <- (df_temp %>% filter(int_est == int_true) %>% nrow()) / nrow(df_temp)
      
    }
    
  }
  
  mat_correct <- mat_correct %>%
    as_tibble(rownames = 'strength') %>%
    pivot_longer(-strength, names_to = 'duration', values_to = 'perc_correct') %>%
    mutate(duration = factor(duration, levels = c(1, 4, 13)))
  return(mat_correct)
  
}

# Function to assess whether higher values of a method's metric are associated with stronger true interaction strenghts:
calculate_assoc_true_strength <- function(df, method, met) {
  
  df <- df %>%
    rename('metric' = all_of(met))
  
  res_temp <- NULL
  
  for (d in unique(df$delta)) {
    
    if (method %in% c('granger', 'te', 'ccm')) {
      
      cor_temp_overall <- df %>%
        mutate(theta_lambda = if_else(theta_lambda < 1, 1 / theta_lambda, theta_lambda)) %>%
        filter(delta == d | theta_lambda == 1) %>%
        # filter(p_value < 0.05) %>%
        cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
      
    } else {
      
      cor_temp_overall <- df %>%
        filter(delta == d | theta_lambda == 1) %>%
        # filter(p_value < 0.05) %>%
        cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
      
    }
    
    cor_temp_neg <- df %>%
      filter(theta_lambda < 1,
             delta == d) %>%
      filter(sig == 'yes') %>%
      cor.test(~ theta_lambda + metric, data = ., method = 'spearman') # lm(data = ., metric ~ theta_lambda)
    cor_temp_pos <- df %>%
      filter(theta_lambda > 1,
             delta == d) %>%
      filter(sig == 'yes') %>%
      cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
    
    # bind_rows(rbind(c(d, cor_temp_neg$estimate, cor_temp_neg$p.value, cor_temp_neg$conf.int[1], cor_temp_neg$conf.int[2]),
    #                 c(d, cor_temp_pos$estimate, cor_temp_pos$p.value, cor_temp_pos$conf.int[1], cor_temp_pos$conf.int[2])) %>%
    #   as_tibble() %>%
    #   mutate(true_int = c('neg', 'pos')))
    
    res_temp <- res_temp %>% bind_rows(rbind(c(d, cor_temp_overall$estimate, cor_temp_overall$p.value),
                                             c(d, cor_temp_neg$estimate, cor_temp_neg$p.value),
                                             c(d, cor_temp_pos$estimate, cor_temp_pos$p.value)) %>%
                                         as_tibble() %>%
                                         mutate(true_int = c('all', 'neg', 'pos')))
    
  }
  
  names(res_temp) <- c('delta', 'rho', 'p_value', 'true_int')
  
  if (method %in% c('granger', 'te', 'ccm')) {
    
    res_temp <- res_temp %>%
      mutate(rho = if_else(true_int == 'neg', -1 * rho, rho))
    
  }
  
  return(res_temp)
  
}

# ------------------------------------------------------------------------------

# Process accuracy of results (correlation coefficients)

# Determine significance/direction of true interaction:
res_corr <- res_corr %>%
  mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
         int_true = if_else(theta_lambda == 1, 'none', int_true))

# Determine significance/direction of detected correlation:
res_corr <- res_corr %>%
  mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
         int_est = if_else(p_value < 0.05, int_est, 'none'))

# Calculate sensitivity/specificity (overall):
print('Sensitivity (Any Interaction) (Overall):')
print((res_corr %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()) / (res_corr %>% filter(int_true != 'none') %>% nrow()))

print('Sensitivity (Correct Direction) (Overall):')
print((res_corr %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_corr %>% filter(int_true == 'neg') %>% nrow()))
print((res_corr %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_corr %>% filter(int_true == 'pos') %>% nrow()))

print('Specificity (Overall):')
print((res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_corr %>% filter(int_true == 'none') %>% nrow()))

# Calculate sensitivity/specificity (by true params):
acc_corr <- calculate_accuracy_matrix(res_corr)

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_corr <- calculate_assoc_true_strength(res_corr %>% mutate(sig = if_else(int_est == 'none', 'no', 'yes')),
                                            method = 'corr', met = 'cor')

# Plot:
p.corr.1 <- ggplot(res_corr %>% mutate(sig = if_else(int_est == 'none', 'no', 'yes'))) +
  geom_violin(aes(x = as.character(theta_lambda), y = cor, group = theta_lambda)) +
  geom_jitter(aes(x = as.character(theta_lambda), y = cor, col = sig)) +
  facet_wrap(~ delta, ncol = 1) +
  theme_classic() +
  labs(title = 'Correlation Coefficients')

p.corr.2 <- ggplot(res_corr %>% mutate(sig = if_else(int_est == 'none', 'no', 'yes')),
                   aes(x = theta_lambda, y = cor, ymin = CI_lower_95, ymax = CI_upper_95, col = sig)) +
  geom_pointrange(position = 'jitter') +
  theme_classic() +
  labs(title = 'Correlation Coefficients')

p.corr.3 <- ggplot(data = acc_corr, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Correlation Coefficients')

print(p.corr.1)
# print(p.corr.2)
# print(p.corr.3)
rm(p.corr.1, p.corr.2, res_corr)
# rm(p.corr.1, p.corr.2, p.corr.3, res_corr, acc_corr, assoc_corr)

# ------------------------------------------------------------------------------

# Process accuracy of results (GAMs)

# Determine significance/direction of true interaction:
res_gam <- res_gam %>%
  mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
         int_true = if_else(theta_lambda == 1, 'none', int_true))

# Determine significance/direction of detected correlation:
res_gam <- res_gam %>%
  mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
         int_est = if_else(CI_lower95 > 0 | CI_upper95 < 0, int_est, 'none')) %>%
  mutate(int_est_confound = if_else(cor_confound > 0, 'pos', 'neg'),
         int_est_confound = if_else(CI_lower95_confound > 0 | CI_upper95_confound < 0, int_est_confound, 'none'))

# Calculate sensitivity/specificity (overall):
print('Sensitivity (Any Interaction) (Overall):')
print((res_gam %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()) / (res_gam %>% filter(int_true != 'none') %>% nrow()))
print((res_gam %>% filter(int_true != 'none' & int_est_confound != 'none') %>% nrow()) / (res_gam %>% filter(int_true != 'none') %>% nrow()))

print('Sensitivity (Correct Direction) (Overall):')
print((res_gam %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_gam %>% filter(int_true == 'neg') %>% nrow()))
print((res_gam %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_gam %>% filter(int_true == 'pos') %>% nrow()))
print((res_gam %>% filter(int_true == 'neg' & int_est_confound == 'neg') %>% nrow()) / (res_gam %>% filter(int_true == 'neg') %>% nrow()))
print((res_gam %>% filter(int_true == 'pos' & int_est_confound == 'pos') %>% nrow()) / (res_gam %>% filter(int_true == 'pos') %>% nrow()))

print('Specificity (Overall):')
print((res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_gam %>% filter(int_true == 'none') %>% nrow()))
print((res_gam %>% filter(int_true == 'none' & int_est_confound == 'none') %>% nrow()) / (res_gam %>% filter(int_true == 'none') %>% nrow()))

# Calculate sensitivity/specificity (by true params):
acc_gam <- calculate_accuracy_matrix(res_gam)
acc_gam_confound <- calculate_accuracy_matrix(res_gam %>%
                                                mutate(int_est = int_est_confound))

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_gam <- calculate_assoc_true_strength(res_gam %>%
                                             mutate(sig = if_else(int_est == 'none', 'no', 'yes')),
                                           method = 'gam', met = 'cor')
assoc_gam_confound <- calculate_assoc_true_strength(res_gam %>%
                                                      mutate(sig = if_else(int_est_confound == 'none', 'no', 'yes')),
                                                    method = 'gam', met = 'cor_confound')

# Plot:
p.gam.1.1 <- ggplot(res_gam %>% mutate(sig = if_else(int_est == 'none', 'no', 'yes'))) +
  geom_violin(aes(x = as.character(theta_lambda), y = cor, group = theta_lambda)) +
  geom_jitter(aes(x = as.character(theta_lambda), y = cor, col = sig)) +
  facet_wrap(~ delta, nrow = 1) +
  theme_classic() +
  labs(title = 'GAMs')
p.gam.1.2 <- ggplot(res_gam %>% mutate(sig = if_else(int_est_confound == 'none', 'no', 'yes'))) +
  geom_violin(aes(x = as.character(theta_lambda), y = cor_confound, group = theta_lambda)) +
  geom_jitter(aes(x = as.character(theta_lambda), y = cor_confound, col = sig)) +
  facet_wrap(~ delta, nrow = 1) +
  theme_classic() +
  labs(title = 'GAMs (Seasonality Controlled)')

p.gam.2.1 <- ggplot(data = acc_gam, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = 'GAMs')
p.gam.2.2 <- ggplot(data = acc_gam_confound, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = 'GAMs (Seasonality Controlled)')

grid.arrange(p.gam.1.1, p.gam.1.2, ncol = 1)
# grid.arrange(p.gam.2.1, p.gam.2.2, nrow = 1)

rm(p.gam.1.1, p.gam.1.2, res_gam)
# rm(p.gam.1.1, p.gam.1.2, p.gam.2.1, p.gam.2.2, res_gam, acc_gam, acc_gam_confound, assoc_gam, assoc_gam_confound)

# ------------------------------------------------------------------------------

# Process accuracy of results (Granger causality)

# Determine significance/direction of true interaction:
res_granger <- res_granger %>%
  # mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
  #        int_true = if_else(theta_lambda == 1, 'none', int_true)) %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'))

# Determine significance/direction of detected correlation:
res_granger <- res_granger %>%
  mutate(int_est = if_else(ftest_p < 0.05, 'interaction', 'none'))

# Check whether any datasets are not stationary:
res_granger %>% filter(adf_p >= 0.05) %>% nrow() %>% print()
res_granger %>% filter(kpss_p < 0.05) %>% nrow() %>% print()

# Get results by direction, and by whether seasonality controlled for:
res_granger_LIST <- vector('list', length = 4)
names(res_granger_LIST) <- c('v1 -> v2 (No confounding)', 'v2 -> v1 (No confounding)', 'v1 -> v2 (Seasonality Controlled)', 'v2 -> v1 (Seasonality Controlled)')

res_granger_LIST[[1]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'none')
res_granger_LIST[[2]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'none')
res_granger_LIST[[3]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'seasonal')
res_granger_LIST[[4]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'seasonal')

# Calculate sensitivity/specificity (overall):
for (i in 1:length(res_granger_LIST)) {
  
  print(names(res_granger_LIST)[i])
  
  print('Sensitivity:')
  print((res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true != 'none') %>% nrow()))
  
  print('Specificity:')
  print((res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()))
  
  print('-----------')
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_granger_LIST <- vector('list', length = 4)
names(acc_granger_LIST) <- names(res_granger_LIST)

for (i in 1:length(acc_granger_LIST)) {
  acc_granger_LIST[[i]] <- calculate_accuracy_matrix(res_granger_LIST[[i]])
}

# Are higher values of logRSS associated with higher true interaction strength?:
assoc_granger_LIST <- vector('list', length = 4)
names(assoc_granger_LIST) <- names(res_granger_LIST)

for (i in 1:length(assoc_granger_LIST)) {
  assoc_granger_LIST[[i]] <- calculate_assoc_true_strength(res_granger_LIST[[i]] %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')),
                                                           method = 'granger', met = 'logRSS')
}
rm(i)

# Plot:
p.granger.1.1 <- res_granger_LIST[[1]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = logRSS, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Granger Causality (', names(res_granger_LIST)[1], ')'))
p.granger.1.2 <- res_granger_LIST[[2]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = logRSS, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Granger Causality (', names(res_granger_LIST)[2], ')'))
p.granger.2.1 <- res_granger_LIST[[3]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = logRSS, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Granger Causality (', names(res_granger_LIST)[3], ')'))
p.granger.2.2 <- res_granger_LIST[[4]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = logRSS, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Granger Causality (', names(res_granger_LIST)[4], ')'))
grid.arrange(p.granger.1.1, p.granger.1.2, p.granger.2.1, p.granger.2.2, ncol = 1)

p.granger.3.1 <- ggplot(data = acc_granger_LIST[[1]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Granger Causality (', names(acc_granger_LIST)[1], ')'))
p.granger.3.2 <- ggplot(data = acc_granger_LIST[[2]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Granger Causality (', names(acc_granger_LIST)[2], ')'))
p.granger.3.3 <- ggplot(data = acc_granger_LIST[[3]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Granger Causality (', names(acc_granger_LIST)[3], ')'))
p.granger.3.4 <- ggplot(data = acc_granger_LIST[[4]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Granger Causality (', names(acc_granger_LIST)[4], ')'))
# grid.arrange(p.granger.3.1, p.granger.3.2, p.granger.3.3, p.granger.3.4, ncol = 2)

rm(p.granger.1.1, p.granger.1.2, p.granger.2.1, p.granger.2.2, res_granger, res_granger_LIST)
# rm(p.granger.1.1, p.granger.1.2, p.granger.2.1, p.granger.2.2, p.granger.3.1, p.granger.3.2, p.granger.3.3, p.granger.3.4,
#    res_granger, res_granger_LIST, acc_granger_LIST, assoc_granger_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (transfer entropy)

# Determine significance/direction of true interaction:
res_te <- res_te %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'))

# Determine significance/direction of detected correlation:
res_te <- res_te %>%
  mutate(int_est = if_else(p_value < 0.05, 'interaction', 'none'))

# Alternatively, base significance on confidence intervals:
res_te <- res_te %>%
  mutate(int_est = if_else(CI_lower > 0, 'interaction', 'none'))
# Lower sensitivity but higher specificity, but only slightly

# Get results by direction and lag:
res_te_LIST <- vector('list', length = 8)
names(res_te_LIST) <- c('v1 -> v2 (lag 1)', 'v1 -> v2 (lag 2)', 'v1 -> v2 (lag 4)', 'v1 -> v2 (lag 6)',
                        'v2 -> v1 (lag 1)', 'v2 -> v1 (lag 2)', 'v2 -> v1 (lag 4)', 'v2 -> v1 (lag 6)')

res_te_LIST[[1]] <- res_te %>% filter(direction == 'v1 -> v2' & lag == '1')
res_te_LIST[[2]] <- res_te %>% filter(direction == 'v1 -> v2' & lag == '2')
res_te_LIST[[3]] <- res_te %>% filter(direction == 'v1 -> v2' & lag == '4')
res_te_LIST[[4]] <- res_te %>% filter(direction == 'v1 -> v2' & lag == '6')
res_te_LIST[[5]] <- res_te %>% filter(direction == 'v2 -> v1' & lag == '1')
res_te_LIST[[6]] <- res_te %>% filter(direction == 'v2 -> v1' & lag == '2')
res_te_LIST[[7]] <- res_te %>% filter(direction == 'v2 -> v1' & lag == '4')
res_te_LIST[[8]] <- res_te %>% filter(direction == 'v2 -> v1' & lag == '6')

# Calculate sensitivity/specificity (overall):
for (i in 1:length(res_te_LIST)) {
  
  print(names(res_te_LIST)[i])
  
  print('Sensitivity:')
  print((res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true != 'none') %>% nrow()))
  
  print('Specificity:')
  print((res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()))
  
  print('-----------')
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_te_LIST <- vector('list', length = 8)
names(acc_te_LIST) <- names(res_te_LIST)

for (i in 1:length(acc_te_LIST)) {
  acc_te_LIST[[i]] <- calculate_accuracy_matrix(res_te_LIST[[i]])
}

# Keep only best-performing lag for each direction:
best_v1xv2 <- acc_te_LIST[1:4] %>% lapply(., function (ix) {mean(ix$perc_correct)}) %>% bind_rows() %>% which.max() %>% names()
best_v2xv1 <- acc_te_LIST[5:8] %>% lapply(., function (ix) {mean(ix$perc_correct)}) %>% bind_rows() %>% which.max() %>% names()

res_te_LIST <- res_te_LIST[c(best_v1xv2, best_v2xv1)]
acc_te_LIST <- acc_te_LIST[c(best_v1xv2, best_v2xv1)]

rm(best_v1xv2, best_v2xv1)

# Are higher values of te associated with higher true interaction strength?:
assoc_te_LIST <- vector('list', length = 2)
names(assoc_te_LIST) <- names(res_te_LIST)

for (i in 1:length(assoc_te_LIST)) {
  assoc_te_LIST[[i]] <- calculate_assoc_true_strength(res_te_LIST[[i]] %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')),
                                                      method = 'te', met = 'te')
}
rm(i)

# Plot:
p.te.1.1 <- res_te_LIST[[1]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = te, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Transfer Entropy (', names(res_te_LIST)[1], ')'))
p.te.1.2 <- res_te_LIST[[2]] %>%
  mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')) %>%
  ggplot(aes(x = as.character(theta_lambda), y = te, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = sig)) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('Transfer Entropy (', names(res_te_LIST)[2], ')'))
grid.arrange(p.te.1.1, p.te.1.2,  ncol = 1)

p.te.2.1 <- ggplot(data = acc_te_LIST[[1]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Transfer Entropy (', names(acc_te_LIST)[1], ')'))
p.te.2.2 <- ggplot(data = acc_te_LIST[[2]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('Transfer Entropy (', names(acc_te_LIST)[2], ')'))
# grid.arrange(p.te.2.1, p.te.2.2, ncol = 2)

rm(p.te.1.1, p.te.1.2, res_te_LIST)#, res_te)
# rm(p.te.1.1, p.te.1.2, p.te.2.1, p.te.2.2, res_te, res_te_LIST, acc_te_LIST, assoc_te_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (CCM)

# Remove unneeded columns:
res_ccm <- res_ccm %>%
  group_by(run, .id, direction, theta_lambda, delta) %>%
  select(run:direction, rho, MannK:p_surr_alt, theta_lambda:delta) %>%
  summarise(rho_mean = mean(rho), rho_max = rho[LibSize == max(LibSize)], MannK = unique(MannK), tp_opt = unique(tp_opt), p_surr = unique(p_surr_alt)) %>%
  ungroup()

# # Or instead use median rhos:
# res_ccm <- res_ccm %>%
#   group_by(run, .id, direction, theta_lambda, delta) %>%
#   select(run:direction, rho_median, MannK:p_surr_alt, theta_lambda:delta) %>%
#   summarise(rho_mean = mean(rho_median), rho_max = rho_median[LibSize == max(LibSize)], MannK = unique(MannK), tp_opt = unique(tp_opt), p_surr = unique(p_surr_alt)) %>%
#   ungroup()
# # very little difference in results

# Determine significance/direction of true interaction:
res_ccm <- res_ccm %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'))

# Determine significance of detected correlation:
res_ccm <- res_ccm %>%
  mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
         int_est_2 = if_else(MannK < 0.05, 'interaction', 'none'), # method 2: check convergence
         int_est_3 = if_else(MannK < 0.05 & tp_opt < 0, 'interaction', 'none')) # method 3: check convergence + ideal tp negative
# res_ccm <- res_ccm %>%
#   mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-value
#          int_est_2 = if_else(p_surr < 0.05 & MannK < 0.05, 'interaction', 'none'), # method 2: method 1 + convergence
#          int_est_3 = if_else(p_surr < 0.05 & MannK < 0.05 & tp_opt < 0, 'interaction', 'none')) # method 3: method 2 + ideal tp negative

# Get results by direction and method of significance calculation:
res_ccm_LIST <- vector('list', length = 6)
names(res_ccm_LIST) <- c('v1 -> v2 (Method 1)', 'v1 -> v2 (Method 2)', 'v1 -> v2 (Method 3)', 'v2 -> v1 (Method 1)', 'v2 -> v1 (Method 2)', 'v2 -> v1 (Method 3)')

res_ccm_LIST[[1]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[2]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[3]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')
res_ccm_LIST[[4]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[5]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[6]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')

# Calculate sensitivity/specificity (overall):
for (i in 1:length(res_ccm_LIST)) {
  
  print(names(res_ccm_LIST)[i])
  
  print('Sensitivity:')
  print((res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true != 'none') %>% nrow()))
  
  print('Specificity:')
  print((res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()))
  
  print('-----------')
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_ccm_LIST <- vector('list', length = 6)
names(acc_ccm_LIST) <- names(res_ccm_LIST)

for (i in 1:length(acc_ccm_LIST)) {
  acc_ccm_LIST[[i]] <- calculate_accuracy_matrix(res_ccm_LIST[[i]])
}

# Are higher values of rho associated with higher true interaction strength?:
assoc_ccm_LIST_mean = assoc_ccm_LIST_max = vector('list', length = 6)
names(assoc_ccm_LIST_mean) <- names(res_ccm_LIST)
names(assoc_ccm_LIST_max) <- names(res_ccm_LIST)

res_ccm_LIST <- lapply(res_ccm_LIST, function(ix) {
  ix %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no'))
})

for (i in 1:length(assoc_ccm_LIST_mean)) {
  assoc_ccm_LIST_mean[[i]] <- calculate_assoc_true_strength(res_ccm_LIST[[i]], method = 'ccm', met = 'rho_mean')
}
for (i in 1:length(assoc_ccm_LIST_max)) {
  assoc_ccm_LIST_max[[i]] <- calculate_assoc_true_strength(res_ccm_LIST[[i]], method = 'ccm', met = 'rho_max')
}
rm(i)

# Check whether there are any places were surrogate rhos are very different from observed rho
# (Could suggest that the surrogates are not representative of the null distribution accounting
# for underlying seasonality)
p.ccm.surr <- res_ccm_surr %>%
  select(run:.id, direction:rho_ci_upper, theta_lambda:delta) %>%
  inner_join(res_ccm %>% select(run:direction, rho_mean:rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_mean, group = direction, color = direction)) +
  theme_classic() +
  facet_grid(theta_lambda ~ delta, scales = 'free')
# print(p.ccm.surr)
rm(res_ccm_surr, p.ccm.surr)

# Plot:
p.ccm.1.1 <- res_ccm_LIST[[1]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[1], ')'), col = '')
p.ccm.1.2 <- res_ccm_LIST[[2]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[2], ')'), col = '')
p.ccm.1.3 <- res_ccm_LIST[[3]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[3], ')'), col = '')

p.ccm.1.4 <- res_ccm_LIST[[4]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[4], ')'), col = '')
p.ccm.1.5 <- res_ccm_LIST[[5]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[5], ')'), col = '')
p.ccm.1.6 <- res_ccm_LIST[[6]] %>%
  ggplot(aes(x = as.character(theta_lambda), y = rho_mean, group = theta_lambda)) +
  geom_violin() +
  geom_jitter(aes(col = factor(int_est, levels = c('none', 'interaction')))) +
  theme_classic() +
  facet_wrap(~ 1 / delta, nrow = 1) +
  labs(title = paste0('CCM (', names(res_ccm_LIST)[6], ')'), col = '')

grid.arrange(p.ccm.1.1, p.ccm.1.2, p.ccm.1.3, ncol = 1)
grid.arrange(p.ccm.1.4, p.ccm.1.5, p.ccm.1.6, ncol = 1)

p.ccm.3.1 <- ggplot(data = acc_ccm_LIST[[1]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[1], ')'))
p.ccm.3.2 <- ggplot(data = acc_ccm_LIST[[2]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[2], ')'))
p.ccm.3.3 <- ggplot(data = acc_ccm_LIST[[3]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[3], ')'))
p.ccm.3.4 <- ggplot(data = acc_ccm_LIST[[4]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[4], ')'))
p.ccm.3.5 <- ggplot(data = acc_ccm_LIST[[5]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[5], ')'))
p.ccm.3.6 <- ggplot(data = acc_ccm_LIST[[6]], aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1)) +
  labs(title = paste0('CCM (', names(acc_ccm_LIST)[6], ')'))
# grid.arrange(p.ccm.3.1, p.ccm.3.4, p.ccm.3.2, p.ccm.3.5, p.ccm.3.3, p.ccm.3.6, ncol = 2)

rm(p.ccm.1.1, p.ccm.1.2, p.ccm.1.3, p.ccm.1.4, p.ccm.1.5, p.ccm.1.6, res_ccm, res_ccm_LIST)
# rm(p.ccm.1.1, p.ccm.1.2, p.ccm.1.3, p.ccm.1.4, p.ccm.1.5, p.ccm.1.6,
#    p.ccm.3.1, p.ccm.3.2, p.ccm.3.3, p.ccm.3.4, p.ccm.3.5, p.ccm.3.6,
#    res_ccm, res_ccm_LIST, acc_ccm_LIST, assoc_ccm_LIST_mean, assoc_ccm_LIST_max)

# ------------------------------------------------------------------------------

# Plot results of all methods

# Heatmaps:
p.comb.1 <- arrangeGrob(p.corr.3, p.gam.2.1, p.gam.2.2, p.granger.3.1, p.granger.3.2, p.granger.3.3, p.granger.3.4, p.te.2.1, p.te.2.2, p.ccm.3.1, p.ccm.3.4, p.ccm.3.3, p.ccm.3.6,
                        layout_matrix = rbind(c(NA, NA, 1, 1, NA, NA),
                                              c(NA, 2, 2, 3, 3, NA),
                                              c(NA, 4, 4, 5, 5, NA),
                                              c(NA, 6, 6, 7, 7, NA),
                                              c(NA, 8, 8, 9, 9, NA),
                                              c(NA, 10, 10, 11, 11, NA),
                                              c(NA, 12, 12, 13, 13, NA)))
plot(p.comb.1)
rm(p.corr.3, p.gam.2.1, p.gam.2.2, p.granger.3.1, p.granger.3.2, p.granger.3.3, p.granger.3.4, p.te.2.1, p.te.2.2, p.ccm.3.1, p.ccm.3.4, p.ccm.3.2, p.ccm.3.5, p.ccm.3.3, p.ccm.3.6)

# Plot percent accurate for all methods:
res_acc <- acc_corr %>% mutate(method = 'Corr. Coef.', direction = 'v1 -> v2') %>%
  bind_rows(acc_corr %>% mutate(method = 'Corr. Coef.', direction = 'v2 -> v1')) %>%
  bind_rows(acc_gam %>% mutate(method = 'GAMs', direction = 'v1 -> v2')) %>%
  bind_rows(acc_gam %>% mutate(method = 'GAMs', direction = 'v2 -> v1')) %>%
  bind_rows(acc_gam_confound %>% mutate(method = 'GAMs (w/ Seas.)', direction = 'v1 -> v2')) %>%
  bind_rows(acc_gam_confound %>% mutate(method = 'GAMs (w/ Seas.)', direction = 'v2 -> v1')) %>%
  bind_rows(acc_granger_LIST[[1]] %>% mutate(method = 'GC', direction = 'v1 -> v2')) %>%
  bind_rows(acc_granger_LIST[[2]] %>% mutate(method = 'GC', direction = 'v2 -> v1')) %>%
  bind_rows(acc_granger_LIST[[3]] %>% mutate(method = 'GC (w/ Seas.)', direction = 'v1 -> v2')) %>%
  bind_rows(acc_granger_LIST[[4]] %>% mutate(method = 'GC (w/ Seas.)', direction = 'v2 -> v1')) %>%
  bind_rows(acc_te_LIST[[1]] %>% mutate(method = 'TE', direction = 'v1 -> v2')) %>%
  bind_rows(acc_te_LIST[[2]] %>% mutate(method = 'TE', direction = 'v2 -> v1')) %>%
  bind_rows(acc_ccm_LIST[[1]] %>% mutate(method = 'CCM (Method 1)', direction = 'v1 -> v2')) %>%
  bind_rows(acc_ccm_LIST[[2]] %>% mutate(method = 'CCM (Method 1)', direction = 'v2 -> v1')) %>%
  bind_rows(acc_ccm_LIST[[5]] %>% mutate(method = 'CCM (Method 3)', direction = 'v1 -> v2')) %>%
  bind_rows(acc_ccm_LIST[[6]] %>% mutate(method = 'CCM (Method 3)', direction = 'v2 -> v1')) %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GAMs (w/ Seas.)', 'GC', 'GC (w/ Seas.)', 'TE', 'CCM (Method 1)', 'CCM (Method 3)'))) %>%
  mutate(s_proxy = case_match(strength, '0' ~ -2, '0.25' ~ -1, '0.5' ~ 0, '4' ~ 3, .default = as.numeric(strength)),
         method_num = as.numeric(method),
         x_use = s_proxy + 0.1 * (method_num - 5))

p.comb.2 <- ggplot(res_acc, aes(x = x_use, y = perc_correct, group = method, shape = method, col = method)) +
  geom_rect(aes(xmin = -Inf, xmax = -1.5, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_rect(aes(xmin = -1.5, xmax = -0.5, ymin = -Inf, ymax = Inf), fill = 'gray95', col = 'gray95') +
  geom_rect(aes(xmin = -0.5, xmax = 0.5, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_rect(aes(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = 'gray95', col = 'gray95') +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'gray95', col = 'gray95') +
  geom_point(size = 2.5) + #geom_line() +
  facet_grid(duration ~ direction) +
  theme_bw() +
  scale_x_continuous(limits = c(-2.4, 3.3), breaks = -2:3, labels = c('0', '0.25', '0.5', '1', '2', '4')) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_shape_manual(values = c(16, 17, 17, 15, 15, 3, 8, 8)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = 'Interaction Strength', y = '% Correct', shape = 'Method', col = 'Method')
print(p.comb.2)
rm(acc_corr, acc_gam, acc_gam_confound, acc_granger_LIST, acc_te_LIST, acc_ccm_LIST)

# Plot Spearman's rho between metrics and true interaction strength for all methods:
# Note: Values are calculated such that we expect POSITIVE rho for all analyses
res_assoc <- assoc_corr %>% mutate(method = 'Corr. Coef.', direction = 'v1 -> v2') %>%
  bind_rows(assoc_corr %>% mutate(method = 'Corr. Coef.', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_gam %>% mutate(method = 'GAMs', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_gam %>% mutate(method = 'GAMs', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_gam_confound %>% mutate(method = 'GAMs (w/ Seas.)', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_gam_confound %>% mutate(method = 'GAMs (w/ Seas.)', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_granger_LIST[[1]] %>% mutate(method = 'GC', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_granger_LIST[[2]] %>% mutate(method = 'GC', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_granger_LIST[[3]] %>% mutate(method = 'GC (w/ Seas.)', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_granger_LIST[[4]] %>% mutate(method = 'GC (w/ Seas.)', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_te_LIST[[1]] %>% mutate(method = 'TE', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_te_LIST[[2]] %>% mutate(method = 'TE', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_ccm_LIST_max[[1]] %>% mutate(method = 'CCM (Method 1)', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_ccm_LIST_max[[2]] %>% mutate(method = 'CCM (Method 1)', direction = 'v2 -> v1')) %>%
  bind_rows(assoc_ccm_LIST_max[[5]] %>% mutate(method = 'CCM (Method 3)', direction = 'v1 -> v2')) %>%
  bind_rows(assoc_ccm_LIST_max[[6]] %>% mutate(method = 'CCM (Method 3)', direction = 'v2 -> v1')) %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GAMs (w/ Seas.)', 'GC', 'GC (w/ Seas.)', 'TE', 'CCM (Method 1)', 'CCM (Method 3)'))) %>%
  mutate(duration = 1 / delta) %>%
  mutate(d_proxy = case_match(duration, 1 ~ 1, 4 ~ 2, 13 ~ 3),
         method_num = as.numeric(method),
         x_use = d_proxy + 0.1 * (method_num - 5))

p.comb.3 <- ggplot(res_assoc, aes(x = x_use, y = rho, group = method, shape = method, col = method)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = 'gray95', col = 'gray95') +
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_point(size = 2.5) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'white', col = NA, alpha = 0.075) +
  geom_point(data = res_assoc %>% filter(rho > 0 & p_value < 0.05), size = 2.5) +
  geom_hline(yintercept = 0) +
  facet_grid(true_int ~ direction) +
  theme_bw() +
  scale_x_continuous(limits = c(0.6, 3.3), breaks = 1:3, labels = c('1', '4', '13')) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_shape_manual(values = c(16, 17, 17, 15, 15, 3, 18, 18)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = 'Interaction Duration', y = "Spearman's Rho", shape = 'Method', col = 'Method')
print(p.comb.3)
rm(assoc_corr, assoc_gam, assoc_gam_confound, assoc_granger_LIST, assoc_te_LIST, assoc_ccm_LIST_mean, assoc_ccm_LIST_max)

# Close pdf:
dev.off()

# Clean up:
rm(list = ls())
