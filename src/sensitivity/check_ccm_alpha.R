# ------------------------------------------------------------------------------
# Code to test the effect of parameter alpha when generated surrogates for CCM
# ------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(rEDM)
library(gridExtra)
library(viridis)

# Read in functions:
source('src/functions_etc/fxns_process_results.R')

# Save plots:
pdf('results/plots/sens_ccm_alpha.pdf', width = 15, height = 8)

# ------------------------------------------------------------------------------

# Read in all results

# Read in "main" results (alpha=20):
source('src/functions_etc/load_main_results.R')

# Clean up unneeded results:
rm(res_corr, res_gam, res_granger, res_granger_LIST, res_te, res_te_LIST,
   acc_weighted_te, best_v1xv2, best_v2xv1, res_ccm_LIST)

# Read in results using various values for alpha (0, 10, 30, 50):
res_filenames_alpha0 <- list.files(path = 'results/ccm_alpha_sens/', pattern = '_0.rds', full.names = TRUE)
res_filenames_alpha10 <- list.files(path = 'results/ccm_alpha_sens/', pattern = '_10.rds', full.names = TRUE)
res_filenames_alpha30 <- list.files(path = 'results/ccm_alpha_sens/', pattern = '_30.rds', full.names = TRUE)
res_filenames_alpha50 <- list.files(path = 'results/ccm_alpha_sens/', pattern = '_50.rds', full.names = TRUE)

res_0 = res_10 = res_30 = res_50 = vector('list', length = length(res_filenames_alpha0))

for (i in 1:length(res_filenames_alpha0)) {
  
  res_0[[i]] <- read_rds(res_filenames_alpha0[i])$CCM
  res_10[[i]] <- read_rds(res_filenames_alpha10[i])$CCM
  res_30[[i]] <- read_rds(res_filenames_alpha30[i])$CCM
  res_50[[i]] <- read_rds(res_filenames_alpha50[i])$CCM
  
}

# Rename lists with correct run numbers:
where_run <- which(!is.na(as.numeric(str_split(res_filenames_alpha0, '_')[[1]])))
names(res_0) = names(res_10) = names(res_30) = names(res_50) <- unlist(map(str_split(res_filenames_alpha0, '_'), where_run))

# Clean up:
rm(i, res_filenames_alpha0, res_filenames_alpha10, res_filenames_alpha30, res_filenames_alpha50, where_run)

# Get true interaction parameter values:
sens_res_LIST <- vector('list', length = 4)
names(sens_res_LIST) <- c('alpha_0', 'alpha_10', 'alpha_30', 'alpha_50')

sens_res_LIST[[1]] <- res_0
sens_res_LIST[[2]] <- res_10
sens_res_LIST[[3]] <- res_30
sens_res_LIST[[4]] <- res_50
rm(res_0, res_10, res_30, res_50)

sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  
  ix %>%
    bind_rows(.id = 'run') %>%
    mutate(run = as.numeric(run)) %>%
    ungroup() %>%
    inner_join(dat %>%
                 select(run:delta) %>%
                 distinct(),
               by = 'run')
  
})

# Get surrogate "data":
sens_res_surr_LIST <- lapply(sens_res_LIST, function(ix) {
  ix %>% filter(data == 'surr')
})

# Get only results using "observed" data:
sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  ix %>% filter(data == 'obs')
})

# Format results of sensitivity analysis:
sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  
  ix %>%
    group_by(run, .id, direction, theta_lambda, delta) %>%
    select(run:direction, rho, MannK:p_surr, theta_lambda, delta) %>%
    summarise(rho_mean = mean(rho), rho_max = rho[LibSize == max(LibSize)], MannK = unique(MannK), tp_opt = unique(tp_opt), max_cmc = unique(max_cmc), p_surr = unique(p_surr)) %>%
    ungroup()
  
})

# Get true and inferred interaction information:
sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  
  ix %>%
    mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
           int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
           int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
    mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
           int_est_3 = if_else(MannK < 0.05 & max_cmc > 0 & tp_opt < 0, 'interaction', 'none')) # method 3: check convergence + ideal tp negative
  
})

# Compare plots of surrogate data:
p.ccm.surr_1 <- res_ccm_surr %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper, theta_lambda:delta) %>%
  inner_join(res_ccm %>% select(run:direction, rho_mean:rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 20')

p.ccm.surr_2 <- sens_res_surr_LIST[[1]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper, theta_lambda:delta) %>%
  inner_join(sens_res_LIST[[1]] %>% select(run:direction, rho_mean:rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 0')

p.ccm.surr_3 <- sens_res_surr_LIST[[4]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper, theta_lambda:delta) %>%
  inner_join(sens_res_LIST[[4]] %>% select(run:direction, rho_mean:rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 50')

grid.arrange(p.ccm.surr_1, p.ccm.surr_2, p.ccm.surr_3, nrow = 1)

rm(res_ccm_surr, sens_res_surr_LIST, p.ccm.surr_1, p.ccm.surr_2, p.ccm.surr_3)

# Keep only relevant columns:
sens_res_LIST <- lapply(1:length(sens_res_LIST), function(ix) {
  sens_res_LIST[[ix]] %>%
    select(run:delta, p_surr:int_est_3) %>%
    mutate(alpha = as.numeric(str_remove(names(sens_res_LIST)[ix], 'alpha_')))
})

# Add "main" results to list:
res_LIST <- vector('list', length = length(sens_res_LIST) + 1)

res_LIST[[1]] <- sens_res_LIST[[1]]
res_LIST[[2]] <- sens_res_LIST[[2]]
res_LIST[[3]] <- res_ccm %>%
  select(run:delta, p_surr:int_est_1, int_est_3) %>%
  mutate(alpha = 20)
res_LIST[[4]] <- sens_res_LIST[[3]]
res_LIST[[5]] <- sens_res_LIST[[4]]

names(res_LIST) <- c('alpha_0', 'alpha_10', 'alpha_20', 'alpha_30', 'alpha_50')
rm(sens_res_LIST, res_ccm)

# ------------------------------------------------------------------------------

# Visualize the effect of different alpha values on surrogate data

# Choose a few simulations to use as examples:
set.seed(28403)
ids_to_use <- sample(1:100, size = 3)

# Get list of data for these simulations:
dat_LIST <- vector('list', length = length(ids_to_use))

dat_LIST[[1]] <- dat %>% filter(run == 4, .id == ids_to_use[1])
dat_LIST[[2]] <- dat %>% filter(run == 4, .id == ids_to_use[2])
dat_LIST[[3]] <- dat %>% filter(run == 4, .id == ids_to_use[3])
rm(ids_to_use)

# Generate surrogates with a range of alpha values:
dat_surr_LIST <- lapply(dat_LIST, function(ix) {
  
  surr_v1 = surr_v2 = vector('list', length = 5)
  
  alpha_vec <- c(0, 10, 20, 30, 50)
  
  for (i in 1:length(alpha_vec)) {
    
    alpha_val <- alpha_vec[i]
    
    surr_v1[[i]] <- SurrogateData(ix$V1_obs, method = 'seasonal', num_surr = 1, T_period = 52.25, alpha = alpha_val) %>%
      as_tibble() %>%
      mutate(time = 1:522) %>%
      pivot_longer(-time, names_to = 'sim', values_to = 'V1_obs') %>%
      mutate(V1_obs = if_else(V1_obs < 0, 0, V1_obs),
             alpha = alpha_val)
    surr_v2[[i]] <- SurrogateData(ix$V2_obs, method = 'seasonal', num_surr = 1, T_period = 52.25, alpha = alpha_val) %>%
      as_tibble() %>%
      mutate(time = 1:522) %>%
      pivot_longer(-time, names_to = 'sim', values_to = 'V2_obs') %>%
      mutate(V2_obs = if_else(V2_obs < 0, 0, V2_obs),
             alpha = alpha_val)
    
  }
  
  dat_surr <- bind_rows(surr_v1) %>%
    left_join(bind_rows(surr_v2),
              by = c('time', 'sim', 'alpha')) %>%
    select(time:V1_obs, V2_obs, alpha)
  
  return(dat_surr)
  
})

# Visualize "observed" data and surrogates:
p_obs_V1_1 <- ggplot() +
  geom_line(data = dat_LIST[[1]], aes(x = time - 104, y = V1_obs)) +
  geom_point(data = dat_LIST[[1]], aes(x = time - 104, y = V1_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs')
p_obs_V1_2 <- ggplot() +
  geom_line(data = dat_LIST[[2]], aes(x = time - 104, y = V1_obs)) +
  geom_point(data = dat_LIST[[2]], aes(x = time - 104, y = V1_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs')
p_obs_V1_3 <- ggplot() +
  geom_line(data = dat_LIST[[3]], aes(x = time - 104, y = V1_obs)) +
  geom_point(data = dat_LIST[[3]], aes(x = time - 104, y = V1_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs')

p_obs_V2_1 <- ggplot() +
  geom_line(data = dat_LIST[[1]], aes(x = time - 104, y = V2_obs)) +
  geom_point(data = dat_LIST[[1]], aes(x = time - 104, y = V2_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs')
p_obs_V2_2 <- ggplot() +
  geom_line(data = dat_LIST[[2]], aes(x = time - 104, y = V2_obs)) +
  geom_point(data = dat_LIST[[2]], aes(x = time - 104, y = V2_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs')
p_obs_V2_3 <- ggplot() +
  geom_line(data = dat_LIST[[3]], aes(x = time - 104, y = V2_obs)) +
  geom_point(data = dat_LIST[[3]], aes(x = time - 104, y = V2_obs)) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs')

p_surr_V1_1 <- ggplot(data = dat_surr_LIST[[1]]) +
  geom_line(aes(x = time, y = V1_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs') +
  scale_color_viridis()
p_surr_V1_2 <- ggplot(data = dat_surr_LIST[[2]]) +
  geom_line(aes(x = time, y = V1_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs') +
  scale_color_viridis()
p_surr_V1_3 <- ggplot(data = dat_surr_LIST[[3]]) +
  geom_line(aes(x = time, y = V1_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs') +
  scale_color_viridis()

p_surr_V2_1 <- ggplot(data = dat_surr_LIST[[1]]) +
  geom_line(aes(x = time, y = V2_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs') +
  scale_color_viridis()
p_surr_V2_2 <- ggplot(data = dat_surr_LIST[[2]]) +
  geom_line(aes(x = time, y = V2_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs') +
  scale_color_viridis()
p_surr_V2_3 <- ggplot(data = dat_surr_LIST[[3]]) +
  geom_line(aes(x = time, y = V2_obs, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs') +
  scale_color_viridis()

# p_surr_1 <- arrangeGrob(arrangeGrob(p_obs_V1_1, p_surr_V1_1, ncol = 1, heights = c(1, 5)),
#                         arrangeGrob(p_obs_V1_2, p_surr_V1_2, ncol = 1, heights = c(1, 5)),
#                         arrangeGrob(p_obs_V1_3, p_surr_V1_3, ncol = 1, heights = c(1, 5)),
#                         ncol = 3)
# p_surr_2 <- arrangeGrob(arrangeGrob(p_obs_V2_1, p_surr_V2_1, ncol = 1, heights = c(1, 5)),
#                         arrangeGrob(p_obs_V2_2, p_surr_V2_2, ncol = 1, heights = c(1, 5)),
#                         arrangeGrob(p_obs_V2_3, p_surr_V2_3, ncol = 1, heights = c(1, 5)),
#                         ncol = 3)

p_surr_1 <- arrangeGrob(arrangeGrob(p_obs_V1_1, p_surr_V1_1, ncol = 1, heights = c(1, 5)),
                        arrangeGrob(p_obs_V2_1, p_surr_V2_1, ncol = 1, heights = c(1, 5)),
                        ncol = 2)
p_surr_2 <- arrangeGrob(arrangeGrob(p_obs_V1_2, p_surr_V1_2, ncol = 1, heights = c(1, 5)),
                        arrangeGrob(p_obs_V2_2, p_surr_V2_2, ncol = 1, heights = c(1, 5)),
                        ncol = 2)
p_surr_3 <- arrangeGrob(arrangeGrob(p_obs_V1_3, p_surr_V1_3, ncol = 1, heights = c(1, 5)),
                        arrangeGrob(p_obs_V2_3, p_surr_V2_3, ncol = 1, heights = c(1, 5)),
                        ncol = 2)

plot(p_surr_1)
plot(p_surr_2)
plot(p_surr_3)

rm(dat_LIST, dat_surr_LIST, p_surr_1, p_surr_2, p_surr_3, p_obs_V1_1, p_obs_V1_2, p_obs_V1_3, p_obs_V2_1,
   p_obs_V2_2, p_obs_V2_3, p_surr_V1_1, p_surr_V1_2, p_surr_V1_3, p_surr_V2_1, p_surr_V2_2, p_surr_V2_3)

# for flu,  smaller value (5-10) is probably fine, as this already adds quite a lot of noise; for RSV, there is a much stronger signal, so a higher value of alpha (~30) might be better

# ------------------------------------------------------------------------------

# Compare results using different values of alpha

# Compare range of p-values:
p1_1 <- res_LIST %>%
  bind_rows() %>%
  filter(direction == 'v1 -> v2') %>%
  ggplot(aes(x = alpha, y = p_surr, group = alpha)) +
  geom_boxplot(aes(fill = alpha)) +
  facet_wrap(~ run, scales = 'free_y') +
  theme_classic() +
  scale_fill_viridis(option = 'G')
p1_2 <- res_LIST %>%
  bind_rows() %>%
  filter(direction == 'v2 -> v1') %>%
  ggplot(aes(x = alpha, y = p_surr, group = alpha)) +
  geom_boxplot(aes(fill = alpha)) +
  facet_wrap(~ run, scales = 'free_y') +
  theme_classic() +
  scale_fill_viridis(option = 'G')
print(p1_1)
print(p1_2)
rm(p1_1, p1_2)

# in general, p-values are lower for higher values of alpha - in other words, high alpha means more likely to be significant
# additionally, there are some runs for which all p-values are pretty much the same - this happens more with positive or short-lived interactions

# Calculate overall sensitivity, specificity, weighted accuracy, and MCC:
res_acc <- lapply(1:length(res_LIST), function(ix) {
  
  sens_pos_1 <- (res_LIST[[ix]] %>% filter(int_true_dir == 'pos' & int_est_1 == 'interaction') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_1 <- (res_LIST[[ix]] %>% filter(int_true_dir == 'neg' & int_est_1 == 'interaction') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_1 <- (res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  sens_pos_3 <- (res_LIST[[ix]] %>% filter(int_true_dir == 'pos' & int_est_3 == 'interaction') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_3 <- (res_LIST[[ix]] %>% filter(int_true_dir == 'neg' & int_est_3 == 'interaction') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_3 <- (res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_3 == 'none') %>% nrow()) / (res_LIST[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  acc_weighted_1 <- (sens_pos_1 * weight_pos + sens_neg_1 * weight_neg + spec_1 * weight_null) / (weight_pos + weight_neg + weight_null)
  acc_weighted_3 <- (sens_pos_3 * weight_pos + sens_neg_3 * weight_neg + spec_3 * weight_null) / (weight_pos + weight_neg + weight_null)
  
  tp1 <- res_LIST[[ix]] %>% filter(int_true != 'none' & int_est_1 != 'none') %>% nrow()
  tn1 <- res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()
  fp1 <- res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_1 != 'none') %>% nrow()
  fn1 <- res_LIST[[ix]] %>% filter(int_true != 'none' & int_est_1 == 'none') %>% nrow()
  
  tp3 <- res_LIST[[ix]] %>% filter(int_true != 'none' & int_est_3 != 'none') %>% nrow()
  tn3 <- res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_3 == 'none') %>% nrow()
  fp3 <- res_LIST[[ix]] %>% filter(int_true == 'none' & int_est_3 != 'none') %>% nrow()
  fn3 <- res_LIST[[ix]] %>% filter(int_true != 'none' & int_est_3 == 'none') %>% nrow()
  
  mcc_1 <- mcc(tp1, tn1, fp1, fn1)
  mcc_3 <- mcc(tp3, tn3, fp3, fn3)
  
  bind_rows(setNames(c(acc_weighted_1, mcc_1, sens_pos_1, sens_neg_1, spec_1), nm = c('acc_weight', 'mcc', 'sens_pos', 'sens_neg', 'spec')),
            setNames(c(acc_weighted_3, mcc_3, sens_pos_3, sens_neg_3, spec_3), nm = c('acc_weight', 'mcc', 'sens_pos', 'sens_neg', 'spec'))) %>%
    mutate(method = c('method_1', 'method_3')) %>%
    mutate(alpha = unique(res_LIST[[ix]]$alpha))
  
}) %>%
  bind_rows()

p2 <- ggplot(data = res_acc %>%
               pivot_longer(acc_weight:spec,
                            names_to = 'metric',
                            values_to = 'value') %>%
               mutate(metric = factor(metric, levels = c('acc_weight', 'mcc', 'sens_pos', 'sens_neg', 'spec'))),
             aes(x = alpha, y = value, group = method, col = method)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ metric, nrow = 1) +
  theme_classic() +
  theme(legend.position = 'right') +
  labs(x = 'alpha', y = 'Value', col = '') +
  scale_color_brewer(palette = 'Set1')

# higher alpha leads to better sens_neg, but worse spec and acc_weight; little impact on mcc or sens_pos; method_3 same for all since doesn't rely on surrogates
# method_1 closest to method_3 when alpha = 0

# Calculate accuracy by interaction parameter values:
acc_by_param <- lapply(1:length(res_LIST), function(ix) {
  res_LIST[[ix]] %>%
    mutate(int_est = int_est_1) %>%
    calculate_accuracy_matrix() %>%
    mutate(alpha = as.numeric(str_remove(names(res_LIST)[ix], 'alpha_')))
}) %>%
  bind_rows()

p3 <- ggplot(data = acc_by_param, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  facet_wrap(~ alpha, nrow = 1) +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1))
grid.arrange(p2, p3, ncol = 1)
rm(p2, p3)

# again, increased alpha leads to higher sensitivity but lower specificity

dev.off()

# ------------------------------------------------------------------------------

# Clean up:
rm(list = ls())
