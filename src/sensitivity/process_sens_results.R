# ------------------------------------------------------------------------------
# Code to extract and process results of sensitivity analyses
# ------------------------------------------------------------------------------

# Setup

# Load packages
library(tidyverse)
library(gridExtra)
library(testthat)
library(viridis)

# Load functions:
source('src/02_dependencies/fxns_process_results.R')

# ------------------------------------------------------------------------------ 

# Read in all results
# To process: forcing low/high, process noise low, obs noise high/low,
# reporting high, 20y, asymmetric (GAM only)

# Get file names:
res_filenames <- list.files(path = 'results/sens/', full.names = TRUE)
expect_true(length(res_filenames) == 18 * 2 * 8)

# Read in results:
data_LIST = corr_LIST = gam_LIST = granger_LIST = te_LIST = ccm_LIST = list()
to_remove_LIST = list()
for (i in 1:length(res_filenames)) {
  print(i)
  
  # Which results are being read?:
  components <- str_split(res_filenames[i], '_')[[1]]
  
  which_run <- as.numeric(components[2])
  run_local <- components[3]
  which_sens <- str_remove(str_flatten(components[4:length(components)], collapse = '_'), '.rds')
  
  # Load results:
  res_temp <- read_rds(res_filenames[i])
  
  # Get true interaction parameters:
  true_params <- res_temp$true_param[c('theta_lambda1', 'delta1'), 1]
  
  # Extract results:
  if (run_local == 'TRUE') {
    
    # Get correlation results:
    if (which_sens != 'asymmetric') {
      
      corr_LIST[[length(corr_LIST) + 1]] <- res_temp$cor %>%
        as_tibble() %>%
        mutate(sens = which_sens,
               run = which_run,
               theta_lambda = true_params[1],
               delta = true_params[2]) %>%
        mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
               int_true = if_else(theta_lambda == 1, 'none', int_true)) %>%
        mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
               int_est = if_else(p_value < 0.05, int_est, 'none')) %>%
        select(sens:run, .id:p_value, theta_lambda:int_est)
      
    }
    
    # Get Granger causality results:
    gc_temp <- res_temp$granger %>%
      as_tibble() %>%
      mutate(sens = which_sens,
             run = which_run,
             theta_lambda = true_params[1],
             delta = true_params[2]) %>%
      mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
             int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
             int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
      mutate(int_est = if_else(ftest_p < 0.05, 'interaction', 'none')) %>%
      select(sens:run, .id, direction:confounding, logRSS, ftest_p, theta_lambda:int_est, adf_p:kpss_p)
    
    to_remove_LIST[[length(to_remove_LIST) + 1]] <- gc_temp %>%
      filter(adf_p >= 0.05 | kpss_p < 0.05) %>%
      select(sens:.id) %>%
      distinct() %>%
      mutate(delete = TRUE)
    
    if (which_sens != 'asymmetric') {
      
      if (which_sens == 'forcing_low') {
        granger_LIST[[length(granger_LIST) + 1]] <- gc_temp %>%
          select(sens:int_est) %>%
          filter(confounding == 'none') %>%
          select(-confounding)
      } else {
        granger_LIST[[length(granger_LIST) + 1]] <- gc_temp %>%
          select(sens:int_est) %>%
          filter(confounding == 'seasonal') %>%
          select(-confounding)
      }
      
    }
    rm(gc_temp)
    
    # Get transfer entropy results:
    if (which_sens != 'asymmetric') {
      
      te_LIST[[length(te_LIST) + 1]] <- res_temp$transfer_entropy %>%
        as_tibble() %>%
        mutate(sens = which_sens,
               run = which_run,
               theta_lambda = true_params[1],
               delta = true_params[2]) %>%
        rename_with(~ str_remove(.x, '_confound')) %>%
        mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
               int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
               int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
        mutate(int_est = if_else(p_value < 0.05 & te > 0, 'interaction', 'none')) %>%
        group_by(direction, .id) %>%
        filter(te == max(te)) %>%
        ungroup() %>%
        select(sens:run, .id, direction, lag, te, p_value, theta_lambda:int_est)
      
    }
    
  } else if (run_local == 'FALSE') {
    
    # Get data:
    data_LIST[[length(data_LIST) + 1]] <- res_temp$data %>%
      as_tibble() %>%
      select(time, date, .id, V1_obs, V2_obs) %>%
      mutate(sens = which_sens,
             run = which_run,
             theta_lambda = true_params[1],
             delta = true_params[2])
    
    # Get GAM results:
    gam_LIST[[length(gam_LIST) + 1]] <- res_temp$gam_cor %>%
      as_tibble() %>%
      mutate(sens = which_sens,
             run = which_run,
             theta_lambda = true_params[1],
             delta = true_params[2]) %>%
      mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
             int_true = if_else(theta_lambda == 1, 'none', int_true)) %>%
      mutate(int_est = if_else(cor_median > 0, 'pos', 'neg'),
             int_est = if_else(CI_lower95 > 0 | CI_upper95 < 0, int_est, 'none')) %>%
      filter(!if_any(b_V1obsln_Intercept:lp__, ~ . >1.05)) %>%
      filter(n_div == 0) %>%
      select(sens:run, .id, cor_median:CI_upper95, theta_lambda:int_est)
    
    # Get CCM results:
    if (which_sens != 'asymmetric') {
      
      ccm_LIST[[length(ccm_LIST) + 1]] <- res_temp$CCM %>%
        as_tibble() %>%
        mutate(sens = which_sens,
               run = which_run,
               theta_lambda = true_params[1],
               delta = true_params[2]) %>%
        filter(data == 'obs') %>%
        group_by(sens, run, .id, direction, theta_lambda, delta) %>%
        summarise(rho_max = rho_median[LibSize == max(LibSize)], tp_opt = unique(tp_opt), max_cmc = unique(max_cmc), p_conv = unique(p_conv), p_surr = unique(p_surr)) %>%
        ungroup() %>%
        mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
               int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
               int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
        mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
               int_est_2 = if_else(p_conv < 0.05 & max_cmc > 0, 'interaction', 'none')) %>% # method 2: check convergence
        select(sens:direction, rho_max, p_surr, p_conv, tp_opt, theta_lambda:delta, int_true:int_est_2)
      
    }
    
  }
  rm(res_temp, true_params, components, which_run, run_local, which_sens)
  
}
rm(i, res_filenames)

# Compile results:
res_corr <- bind_rows(corr_LIST)
res_gam <- bind_rows(gam_LIST)
res_granger <- bind_rows(granger_LIST)
res_te <- bind_rows(te_LIST)
res_ccm <- bind_rows(ccm_LIST)

rm(corr_LIST, gam_LIST, granger_LIST, te_LIST, ccm_LIST)

# Check where GAMs had trouble:
res_gam %>% split(~ .$sens) %>% lapply(function(ix) table(ix$run))

# Remove results where data were not stationary:
df_to_remove <- bind_rows(to_remove_LIST)

res_corr <- res_corr %>% left_join(df_to_remove %>% mutate(.id = as.integer(.id)),
                                   by = c('sens', 'run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(!delete)
res_gam <- res_gam %>% left_join(df_to_remove %>% mutate(.id = as.integer(.id)),
                                 by = c('sens', 'run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(!delete)
res_granger <- res_granger %>% left_join(df_to_remove,
                                         by = c('sens', 'run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(!delete)
res_te <- res_te %>% left_join(df_to_remove %>% mutate(.id = as.integer(.id)),
                               by = c('sens', 'run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(!delete)
res_ccm <- res_ccm %>% left_join(df_to_remove,
                                 by = c('sens', 'run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(!delete)

rm(to_remove_LIST, df_to_remove)

# Get "main" results too:
res_MAIN <- read_rds('results/res_compiled.rds')

res_corr <- res_corr %>% bind_rows(res_MAIN[[1]] %>% select(!delete) %>% mutate(sens = 'Main'))
res_gam <- res_gam %>% bind_rows(res_MAIN[[2]] %>% select(!delete) %>% mutate(sens = 'Main'))
res_granger <- res_granger %>% bind_rows(res_MAIN[[3]] %>% select(!confounding) %>% mutate(sens = 'Main'))
res_te <- res_te %>% bind_rows(res_MAIN[[4]] %>% mutate(sens = 'Main'))
res_ccm <- res_ccm %>% bind_rows(res_MAIN[[5]] %>% select(!int_est_3) %>% mutate(sens = 'Main'))

rm(res_MAIN)

# ------------------------------------------------------------------------------

# Calculate accuracy of results by sensitivity analysis

# Correlation coefficients:
acc_corr <- res_corr %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter((int_true == 'pos' & int_est == 'pos') | (int_true == 'neg' & int_est == 'neg')) %>% nrow()) / (ix %>% filter(int_true == 'pos' | int_true == 'neg') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_byparam_corr <- res_corr %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')
assoc_corr <- res_corr %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'corr', met = 'cor')) %>% bind_rows(.id = 'which_sens')

# GAMs:
acc_gam <- res_gam %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter((int_true == 'pos' & int_est == 'pos') | (int_true == 'neg' & int_est == 'neg')) %>% nrow()) / (ix %>% filter(int_true == 'pos' | int_true == 'neg') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_byparam_gam <- res_gam %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')
assoc_gam <- res_gam %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'gam', met = 'cor_median')) %>% bind_rows(.id = 'which_sens')

# Granger causality:
res_granger_v1v2 <- res_granger %>% filter(direction == 'v1 -> v2')
res_granger_v2v1 <- res_granger %>% filter(direction == 'v2 -> v1')

acc_granger_v1v2 <- res_granger_v1v2 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')
acc_granger_v2v1 <- res_granger_v2v1 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_byparam_granger_v1v2 <- res_granger_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')
acc_byparam_granger_v2v1 <- res_granger_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')

assoc_granger_v1v2 <- res_granger_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'granger', met = 'logRSS')) %>% bind_rows(.id = 'which_sens')
assoc_granger_v2v1 <- res_granger_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'granger', met = 'logRSS')) %>% bind_rows(.id = 'which_sens')

# Transfer entropy:
res_te_v1v2 <- res_te %>% filter(direction == 'v1 -> v2')
res_te_v2v1 <- res_te %>% filter(direction == 'v2 -> v1')

acc_te_v1v2 <- res_te_v1v2 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')
acc_te_v2v1 <- res_te_v2v1 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_byparam_te_v1v2 <- res_te_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')
acc_byparam_te_v2v1 <- res_te_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix)) %>% bind_rows(.id = 'which_sens')

assoc_te_v1v2 <- res_te_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'te', met = 'te')) %>% bind_rows(.id = 'which_sens')
assoc_te_v2v1 <- res_te_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'te', met = 'te')) %>% bind_rows(.id = 'which_sens')

# CCM:
res_ccm_v1v2 <- res_ccm %>% filter(direction == 'v1 -> v2')
res_ccm_v2v1 <- res_ccm %>% filter(direction == 'v2 -> v1')

acc_ccm1_v1v2 <- res_ccm_v1v2 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est_1 == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')
acc_ccm1_v2v1 <- res_ccm_v2v1 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est_1 == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_ccm2_v1v2 <- res_ccm_v1v2 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est_2 == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')
acc_ccm2_v2v1 <- res_ccm_v2v1 %>%
  split(~ .$sens) %>%
  lapply(function(ix) {
    as_tibble(data.frame(sens = (ix %>% filter(int_true != 'none' & int_est_2 == 'interaction') %>% nrow()) / (ix %>% filter(int_true != 'none') %>% nrow()),
                         spec = (ix %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()) / (ix %>% filter(int_true == 'none') %>% nrow())))
  }) %>%
  bind_rows(.id = 'which_sens')

acc_byparam_ccm1_v1v2 <- res_ccm_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix %>% mutate(int_est = int_est_1))) %>% bind_rows(.id = 'which_sens')
acc_byparam_ccm1_v2v1 <- res_ccm_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix %>% mutate(int_est = int_est_1))) %>% bind_rows(.id = 'which_sens')
acc_byparam_ccm2_v1v2 <- res_ccm_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix %>% mutate(int_est = int_est_2))) %>% bind_rows(.id = 'which_sens')
acc_byparam_ccm2_v2v1 <- res_ccm_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_accuracy_matrix(ix %>% mutate(int_est = int_est_2))) %>% bind_rows(.id = 'which_sens')

assoc_ccm_v1v2 <- res_ccm_v1v2 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'ccm', met = 'rho_max')) %>% bind_rows(.id = 'which_sens')
assoc_ccm_v2v1 <- res_ccm_v2v1 %>% split(~ .$sens) %>% lapply(function(ix) calculate_assoc_true_strength(ix, method = 'ccm', met = 'rho_max')) %>% bind_rows(.id = 'which_sens')

# Clean up:
rm(res_corr, res_gam, res_granger, res_granger_v1v2, res_granger_v2v1, res_te, res_te_v1v2, res_te_v2v1, res_ccm, res_ccm_v1v2, res_ccm_v2v1)

# ------------------------------------------------------------------------------

# Plot results

# Plot overall accuracy:
df_acc <- acc_corr %>% mutate(method = 'Corr. Coef.') %>%
  bind_rows(acc_gam %>% mutate(method = 'GAMs')) %>%
  bind_rows(acc_granger_v1v2 %>% mutate(method = 'GC (V1 -> V2)')) %>%
  bind_rows(acc_granger_v2v1 %>% mutate(method = 'GC (V2 -> V1)')) %>%
  bind_rows(acc_te_v1v2 %>% mutate(method = 'TE (V1 -> V2)')) %>%
  bind_rows(acc_te_v2v1 %>% mutate(method = 'TE (V2 -> V1)')) %>%
  bind_rows(acc_ccm1_v1v2 %>% mutate(method = 'CCM (Method 1) (V1 -> V2)')) %>%
  bind_rows(acc_ccm1_v2v1 %>% mutate(method = 'CCM (Method 1) (V2 -> V1)')) %>%
  bind_rows(acc_ccm2_v1v2 %>% mutate(method = 'CCM (Method 2) (V1 -> V2)')) %>%
  bind_rows(acc_ccm2_v2v1 %>% mutate(method = 'CCM (Method 2) (V2 -> V1)'))

df_acc <- df_acc %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC (V1 -> V2)', 'GC (V2 -> V1)', 'TE (V1 -> V2)', 'TE (V2 -> V1)',
                                            'CCM (Method 1) (V1 -> V2)', 'CCM (Method 1) (V2 -> V1)', 'CCM (Method 2) (V1 -> V2)', 'CCM (Method 2) (V2 -> V1)')))

df_acc <- df_acc %>%
  mutate(which_sens = case_when(which_sens == '20y' ~ '20 years',
                                which_sens == 'forcing_high' ~ 'High forcing',
                                which_sens == 'forcing_low' ~ 'No forcing',
                                which_sens == 'obs_noise_high' ~ 'High obs noise',
                                which_sens == 'obs_noise_low' ~ 'Low obs noise',
                                which_sens == 'process_noise_low' ~ 'Low process noise',
                                .default = which_sens)) %>%
  mutate(which_sens = factor(which_sens, levels = c('Main', 'High forcing', 'No forcing', 'Low process noise', 'High obs noise', 'Low obs noise', '20 years'))) %>%
  filter(!is.na(which_sens))

p_acc <- ggplot(data = df_acc %>% filter(which_sens != 'Main')) +
  # geom_point(data = df_acc %>% filter(which_sens == 'Main'),
  #            aes(x = sens, y = spec, shape = method), size = 3) +
  geom_vline(data = df_acc %>% filter(which_sens == 'Main'),
             aes(xintercept = sens), lty = 2) +
  geom_hline(data = df_acc %>% filter(which_sens == 'Main'),
             aes(yintercept = spec), lty = 2) +
  geom_point(aes(x = sens, y = spec, shape = method, col = which_sens, alpha = method), size = 3) +
  facet_wrap(~ method, ncol = 2) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'right') +
  scale_x_continuous(limits = c(0, 1), n.breaks = 10) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 10) +
  scale_shape_manual(values = c(18, 17, 15, 15, 3, 3, 8, 8, 8, 8), guide = 'none') +
  # scale_color_manual(values = c('#fb6a4a', '#fa9fb5', '#8dd3c7', '#bc80bd', '#bebada', '#fdb462')) +
  scale_color_manual(values = c('#d95f02', '#e6ab02', '#7570b3', '#1b9e77', '#66a61e', '#e7298a')) +
  scale_alpha_manual(values = c(0.7, 0.7, 0.7, 0.7, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0), guide = 'none') +
  labs(x = 'Sensitivity', y = 'Specificity', col = 'Analysis')

# Plot accuracy by true interaction parameters:
p_acc_corr <- ggplot(data = acc_byparam_corr %>%
                       filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                       mutate(duration = paste0(duration, ' week'),
                              duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                       mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                       group_by(which_sens) %>%
                       mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                       ungroup(),
                     aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, nrow = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Correlation Coefficients', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

p_acc_gam <- ggplot(data = acc_byparam_gam %>%
                      filter(!(which_sens %in% c('Main', 'reporting_high', 'asymmetric'))) %>%
                      mutate(duration = paste0(duration, ' week'),
                             duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                      mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                      group_by(which_sens) %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                      ungroup(),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, nrow = 2) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'GAMs', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

p_acc_granger12 <- ggplot(data = acc_byparam_granger_v1v2 %>%
                            filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                            mutate(duration = paste0(duration, ' week'),
                                   duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                            mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                            group_by(which_sens) %>%
                            mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                            ungroup(),
                          aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Granger Causality (V1 -> V2)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')
p_acc_granger21 <- ggplot(data = acc_byparam_granger_v2v1 %>%
                            filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                            mutate(duration = paste0(duration, ' week'),
                                   duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                            mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                            group_by(which_sens) %>%
                            mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                            ungroup(),
                          aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Granger Causality (V2 -> V1)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

p_acc_te12 <- ggplot(data = acc_byparam_te_v1v2 %>%
                       filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                       mutate(duration = paste0(duration, ' week'),
                              duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                       mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                       group_by(which_sens) %>%
                       mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                       ungroup(),
                     aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Transfer Entropy (V1 -> V2)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')
p_acc_te21 <- ggplot(data = acc_byparam_te_v2v1 %>%
                       filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                       mutate(duration = paste0(duration, ' week'),
                              duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                       mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                       group_by(which_sens) %>%
                       mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                       ungroup(),
                     aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'Transfer Entropy (V2 -> V1)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

p_acc_ccm1_12 <- ggplot(data = acc_byparam_ccm1_v1v2 %>%
                          filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                          mutate(duration = paste0(duration, ' week'),
                                 duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                          mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                          group_by(which_sens) %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                          ungroup(),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'CCM (Method 1) (V1 -> V2)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')
p_acc_ccm1_21 <- ggplot(data = acc_byparam_ccm1_v2v1 %>%
                          filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                          mutate(duration = paste0(duration, ' week'),
                                 duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                          mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                          group_by(which_sens) %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                          ungroup(),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'CCM (Method 1) (V2 -> V1)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

p_acc_ccm2_12 <- ggplot(data = acc_byparam_ccm2_v1v2 %>%
                          filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                          mutate(duration = paste0(duration, ' week'),
                                 duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                          mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                          group_by(which_sens) %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                          ungroup(),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'CCM (Method 2) (V1 -> V2)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')
p_acc_ccm2_21 <- ggplot(data = acc_byparam_ccm2_v2v1 %>%
                          filter(!(which_sens %in% c('Main', 'reporting_high'))) %>%
                          mutate(duration = paste0(duration, ' week'),
                                 duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                          mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                          group_by(which_sens) %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                          ungroup(),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ which_sens, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(title = 'CCM (Method 2) (V2 -> V1)', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration')

# Plot association between strength and point estimates:
df_assoc <- assoc_corr %>% mutate(method = 'Corr. Coef.') %>%
  bind_rows(assoc_gam %>% mutate(method = 'GAMs')) %>%
  bind_rows(assoc_granger_v1v2 %>% mutate(method = 'GC (V1 -> V2)')) %>%
  bind_rows(assoc_granger_v2v1 %>% mutate(method = 'GC (V2 -> V1)')) %>%
  bind_rows(assoc_te_v1v2 %>% mutate(method = 'TE (V1 -> V2)')) %>%
  bind_rows(assoc_te_v2v1 %>% mutate(method = 'TE (V2 -> V1)')) %>%
  bind_rows(assoc_ccm_v1v2 %>% mutate(method = 'CCM (V1 -> V2)')) %>%
  bind_rows(assoc_ccm_v2v1 %>% mutate(method = 'CCM (V2 -> V1)'))

df_assoc <- df_assoc %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC (V1 -> V2)', 'GC (V2 -> V1)',
                                            'TE (V1 -> V2)', 'TE (V2 -> V1)', 'CCM (V1 -> V2)', 'CCM (V2 -> V1)')))

df_assoc <- df_assoc %>%
  mutate(which_sens = case_when(which_sens == '20y' ~ '20 years',
                                which_sens == 'forcing_high' ~ 'High forcing',
                                which_sens == 'forcing_low' ~ 'No forcing',
                                which_sens == 'obs_noise_high' ~ 'High obs noise',
                                which_sens == 'obs_noise_low' ~ 'Low obs noise',
                                which_sens == 'process_noise_low' ~ 'Low process noise',
                                .default = which_sens)) %>%
  mutate(which_sens = factor(which_sens, levels = c('Main', 'High forcing', 'No forcing', 'Low process noise', 'High obs noise', 'Low obs noise', '20 years'))) %>%
  filter(!is.na(which_sens))

p_assoc <- df_assoc %>%
  mutate(delta = factor(delta)) %>%
  ggplot(aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 1.0, lty = 2) +
  geom_pointrange() +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_grid(which_sens ~ method) +
  scale_x_continuous(n.breaks = 6) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Relative Change in Point Estimate Due to 2x Change in True Strength', y = 'Duration')

# Plot results of asymmetric analysis:
p_asym <- ggplot(data = acc_byparam_gam %>%
                   filter(which_sens %in% c('Main', 'asymmetric')) %>%
                   mutate(duration = paste0(duration, ' week'),
                          duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                   mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))) %>%
                   group_by(which_sens) %>%
                   mutate(strength_proxy = rank(strength, ties.method = 'min')) %>%
                   ungroup(),
                 aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, linetype = which_sens, alpha = which_sens)) +
  geom_line() +
  geom_point(size = 3.5) +
  facet_wrap(~ duration, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_linetype_manual(values = c(2, 1)) +
  scale_alpha_manual(values = c(0.35, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct')

# Print all plots:
pdf(file = 'results/plots/plot_sensitivity.pdf', width = 16, height = 12)
print(p_acc)
print(p_acc_corr)
print(p_acc_gam)
grid.arrange(p_acc_granger12, p_acc_granger21, nrow = 1)
grid.arrange(p_acc_te12, p_acc_te21, nrow = 1)
grid.arrange(p_acc_ccm1_12, p_acc_ccm1_21, nrow = 1)
grid.arrange(p_acc_ccm2_12, p_acc_ccm2_21, nrow = 1)
print(p_assoc)
print(p_asym)
dev.off()

# ------------------------------------------------------------------------------

# Clean up:
rm(list = ls())
