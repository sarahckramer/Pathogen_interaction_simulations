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
source('src/02_dependencies/fxns_process_results.R')

# Save plots:
pdf('results/plots/sens_ccm_alpha_NEW.pdf', width = 15, height = 8)

# ------------------------------------------------------------------------------

# Read in all results

# Get synthetic data and true parameter values:

res_filenames_F <- list.files(path = 'results/', pattern = 'FALSE', full.names = TRUE) # run locally
results_F <- vector('list', length = length(res_filenames_F))

for (i in 1:length(res_filenames_F)) {
  results_F[[i]] <- read_rds(res_filenames_F[i])
}
rm(i)

where_run <- which(!is.na(as.numeric(str_split(res_filenames_F, '_')[[1]])))
names(results_F) <- unlist(map(str_split(res_filenames_F, '_'), where_run))

int_params <- c('theta_lambda1', 'delta1')
res_trueparams <- lapply(1:length(results_F), function(ix) {
  results_F[[ix]]$true_param[int_params, 1]
}) %>%
  bind_rows() %>%
  rename('theta_lambda' = 'theta_lambda1',
         'delta' = 'delta1') %>%
  mutate(run = as.numeric(names(results_F))) %>%
  arrange(run)

data_list <- vector('list', length = length(results_F))
for (i in 1:length(results_F)) {
  data_list[[i]] <- results_F[[i]]$data %>%
    select(time, date, .id, V1_obs, V2_obs) %>%
    mutate(run = as.numeric(names(results_F)[i]))
}

dat <- bind_rows(data_list) %>%
  arrange(run) %>%
  as_tibble() %>%
  inner_join(res_trueparams, by = 'run')

rm(res_filenames_F, where_run, data_list, i, int_params, res_trueparams)

# Also get CCM results for alpha = 0:
res_0 <- lapply(results_F, getElement, 'CCM')
rm(results_F)

# Read in results using various values for alpha (0, 0.1, 0.25, 0.5):
res_filenames_TENTH <- list.files(path = 'results/sens_ccm_alpha/', pattern = 'TENTH', full.names = TRUE)
res_filenames_QUARTER <- list.files(path = 'results/sens_ccm_alpha/', pattern = 'QUARTER', full.names = TRUE)
res_filenames_HALF <- list.files(path = 'results/sens_ccm_alpha/', pattern = 'HALF', full.names = TRUE)

res_10 = res_25 = res_50 = vector('list', length = length(res_filenames_HALF))

for (i in 1:length(res_filenames_HALF)) {
  
  res_10[[i]] <- read_rds(res_filenames_TENTH[i])$CCM
  res_25[[i]] <- read_rds(res_filenames_QUARTER[i])$CCM
  res_50[[i]] <- read_rds(res_filenames_HALF[i])$CCM
  
}

# Rename lists with correct run numbers:
where_run <- which(!is.na(as.numeric(str_split(res_filenames_HALF, '_')[[1]])))
names(res_10) = names(res_25) = names(res_50) <- unlist(map(str_split(res_filenames_HALF, '_'), where_run))

# Clean up:
rm(i, res_filenames_TENTH, res_filenames_QUARTER, res_filenames_HALF, where_run)

# Get true interaction parameter values:
sens_res_LIST <- vector('list', length = 4)
names(sens_res_LIST) <- c('alpha_0', 'alpha_10', 'alpha_25', 'alpha_50')

sens_res_LIST[[1]] <- res_0
sens_res_LIST[[2]] <- res_10
sens_res_LIST[[3]] <- res_25
sens_res_LIST[[4]] <- res_50
rm(res_0, res_10, res_25, res_50)

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
    select(run:direction, rho_median, tp_opt:p_surr, theta_lambda, delta) %>%
    summarise(rho_mean = mean(rho_median), rho_max = rho_median[LibSize == max(LibSize)], tp_opt = unique(tp_opt), max_cmc = unique(max_cmc), p_conv = unique(p_conv), p_surr = unique(p_surr)) %>%
    ungroup()
  
})

# Get true and inferred interaction information:
sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  
  ix %>%
    mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
           int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
           int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
    mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
           int_est_2 = if_else(p_conv < 0.05 & max_cmc > 0, 'interaction', 'none')) # method 2: check convergence
  
})

# Limit to columns of interest:
sens_res_LIST <- lapply(sens_res_LIST, function(ix) {
  
  ix %>%
    select(run:direction, rho_max, p_surr, p_conv, tp_opt, theta_lambda:delta, int_true:int_est_2)
  
})

# Compare plots of surrogate data:
p.ccm.surr_1 <- sens_res_surr_LIST[[1]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper) %>%
  inner_join(sens_res_LIST[[1]] %>% select(run:direction, rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho_median, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 0')

p.ccm.surr_2 <- sens_res_surr_LIST[[2]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper) %>%
  inner_join(sens_res_LIST[[2]] %>% select(run:direction, rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho_median, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 0.1')

p.ccm.surr_3 <- sens_res_surr_LIST[[3]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper) %>%
  inner_join(sens_res_LIST[[3]] %>% select(run:direction, rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho_median, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 0.25')

p.ccm.surr_4 <- sens_res_surr_LIST[[4]] %>%
  filter(theta_lambda == 0.25) %>%
  select(run:.id, direction:rho_ci_upper) %>%
  inner_join(sens_res_LIST[[4]] %>% select(run:direction, rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho_median, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  theme(legend.position = 'none') +
  facet_wrap(~ delta, scales = 'free', ncol = 1) +
  labs(title = 'alpha = 0.5')

# grid.arrange(p.ccm.surr_1, p.ccm.surr_2, p.ccm.surr_3, p.ccm.surr_4, nrow = 1)

rm(sens_res_surr_LIST, p.ccm.surr_1, p.ccm.surr_2, p.ccm.surr_3, p.ccm.surr_4)

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

# Log-transform and center data:
dat_LIST <- lapply(dat_LIST, function(ix) {
  
  ix %>% mutate(V1_obs_ln = scale(log(V1_obs + 1), scale = FALSE),
                V2_obs_ln = scale(log(V2_obs + 1), scale = FALSE))
  
})

# Generate surrogates with a range of alpha values:
dat_surr_LIST <- lapply(dat_LIST, function(ix) {
  
  surr_v1 = surr_v2 = vector('list', length = 4)
  
  alpha_vec <- c(0, 0.1, 0.25, 0.5)
  
  for (i in 1:length(alpha_vec)) {
    
    alpha_val <- alpha_vec[i]
    
    surr_v1[[i]] <- SurrogateData(ix$V1_obs_ln, method = 'seasonal', num_surr = 1, T_period = 52.25, alpha = alpha_val) %>%
      as_tibble() %>%
      mutate(time = 1:522) %>%
      pivot_longer(-time, names_to = 'sim', values_to = 'V1_obs_ln') %>%
      mutate(alpha = alpha_val)
    surr_v2[[i]] <- SurrogateData(ix$V2_obs_ln, method = 'seasonal', num_surr = 1, T_period = 52.25, alpha = alpha_val) %>%
      as_tibble() %>%
      mutate(time = 1:522) %>%
      pivot_longer(-time, names_to = 'sim', values_to = 'V2_obs_ln') %>%
      mutate(alpha = alpha_val)
    
  }
  
  dat_surr <- bind_rows(surr_v1) %>%
    left_join(bind_rows(surr_v2),
              by = c('time', 'sim', 'alpha')) %>%
    select(time:V1_obs_ln, V2_obs_ln, alpha)
  
  return(dat_surr)
  
})

# Visualize "observed" data and surrogates:
p_obs_V1_1 <- ggplot() +
  geom_line(data = dat_LIST[[1]], aes(x = time - 104, y = V1_obs_ln)) +
  geom_point(data = dat_LIST[[1]], aes(x = time - 104, y = V1_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln')
p_obs_V1_2 <- ggplot() +
  geom_line(data = dat_LIST[[2]], aes(x = time - 104, y = V1_obs_ln)) +
  geom_point(data = dat_LIST[[2]], aes(x = time - 104, y = V1_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln')
p_obs_V1_3 <- ggplot() +
  geom_line(data = dat_LIST[[3]], aes(x = time - 104, y = V1_obs_ln)) +
  geom_point(data = dat_LIST[[3]], aes(x = time - 104, y = V1_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln')

p_obs_V2_1 <- ggplot() +
  geom_line(data = dat_LIST[[1]], aes(x = time - 104, y = V2_obs_ln)) +
  geom_point(data = dat_LIST[[1]], aes(x = time - 104, y = V2_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln')
p_obs_V2_2 <- ggplot() +
  geom_line(data = dat_LIST[[2]], aes(x = time - 104, y = V2_obs_ln)) +
  geom_point(data = dat_LIST[[2]], aes(x = time - 104, y = V2_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln')
p_obs_V2_3 <- ggplot() +
  geom_line(data = dat_LIST[[3]], aes(x = time - 104, y = V2_obs_ln)) +
  geom_point(data = dat_LIST[[3]], aes(x = time - 104, y = V2_obs_ln), size = 0.5) +
  theme_classic() +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln')

p_surr_V1_1 <- ggplot(data = dat_surr_LIST[[1]]) +
  geom_line(aes(x = time, y = V1_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln') +
  scale_color_viridis()
p_surr_V1_2 <- ggplot(data = dat_surr_LIST[[2]]) +
  geom_line(aes(x = time, y = V1_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln') +
  scale_color_viridis()
p_surr_V1_3 <- ggplot(data = dat_surr_LIST[[3]]) +
  geom_line(aes(x = time, y = V1_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V1_obs_ln') +
  scale_color_viridis()

p_surr_V2_1 <- ggplot(data = dat_surr_LIST[[1]]) +
  geom_line(aes(x = time, y = V2_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln') +
  scale_color_viridis()
p_surr_V2_2 <- ggplot(data = dat_surr_LIST[[2]]) +
  geom_line(aes(x = time, y = V2_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln') +
  scale_color_viridis()
p_surr_V2_3 <- ggplot(data = dat_surr_LIST[[3]]) +
  geom_line(aes(x = time, y = V2_obs_ln, col = alpha, group = sim)) +
  facet_wrap(~ alpha, ncol = 1) +
  theme_classic() +
  theme(legend.position = 'none') +
  labs(x = 'Time (Weeks)', y = 'V2_obs_ln') +
  scale_color_viridis()

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

# ------------------------------------------------------------------------------

# Compare results using different values of alpha

# Compare range of p-values:
p1_1 <- sens_res_LIST %>%
  bind_rows(.id = 'alpha') %>%
  mutate(alpha = as.numeric(str_remove(alpha, 'alpha_')) / 100) %>%
  filter(direction == 'v1 -> v2') %>%
  ggplot(aes(x = alpha, y = p_surr, group = alpha)) +
  geom_boxplot(aes(fill = alpha)) +
  facet_wrap(~ run, scales = 'free_y') +
  theme_classic() +
  scale_fill_viridis(option = 'G')
p1_2 <- sens_res_LIST %>%
  bind_rows(.id = 'alpha') %>%
  mutate(alpha = as.numeric(str_remove(alpha, 'alpha_')) / 100) %>%
  filter(direction == 'v2 -> v1') %>%
  ggplot(aes(x = alpha, y = p_surr, group = alpha)) +
  geom_boxplot(aes(fill = alpha)) +
  facet_wrap(~ run, scales = 'free_y') +
  theme_classic() +
  scale_fill_viridis(option = 'G')
print(p1_1)
print(p1_2)
rm(p1_1, p1_2)

# Calculate overall sensitivity, specificity, weighted accuracy, and MCC:
res_v1xv2 <- lapply(sens_res_LIST, function(ix) {
  ix %>% filter(direction == 'v1 -> v2')
})
res_v2xv1 <- lapply(sens_res_LIST, function(ix) {
  ix %>% filter(direction == 'v2 -> v1')
})

res_acc_v1xv2 <- lapply(1:length(res_v1xv2), function(ix) {
  
  sens_pos_1 <- (res_v1xv2[[ix]] %>% filter(int_true_dir == 'pos' & int_est_1 == 'interaction') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_1 <- (res_v1xv2[[ix]] %>% filter(int_true_dir == 'neg' & int_est_1 == 'interaction') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_1 <- (res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  sens_pos_3 <- (res_v1xv2[[ix]] %>% filter(int_true_dir == 'pos' & int_est_2 == 'interaction') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_3 <- (res_v1xv2[[ix]] %>% filter(int_true_dir == 'neg' & int_est_2 == 'interaction') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_3 <- (res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()) / (res_v1xv2[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  tp1 <- res_v1xv2[[ix]] %>% filter(int_true != 'none' & int_est_1 != 'none') %>% nrow()
  tn1 <- res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()
  fp1 <- res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_1 != 'none') %>% nrow()
  fn1 <- res_v1xv2[[ix]] %>% filter(int_true != 'none' & int_est_1 == 'none') %>% nrow()
  
  tp3 <- res_v1xv2[[ix]] %>% filter(int_true != 'none' & int_est_2 != 'none') %>% nrow()
  tn3 <- res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()
  fp3 <- res_v1xv2[[ix]] %>% filter(int_true == 'none' & int_est_2 != 'none') %>% nrow()
  fn3 <- res_v1xv2[[ix]] %>% filter(int_true != 'none' & int_est_2 == 'none') %>% nrow()
  
  mcc_1 <- mcc(tp1, tn1, fp1, fn1)
  mcc_3 <- mcc(tp3, tn3, fp3, fn3)
  
  bind_rows(setNames(c(mcc_1, sens_pos_1, sens_neg_1, spec_1), nm = c('mcc', 'sens_pos', 'sens_neg', 'spec')),
            setNames(c(mcc_3, sens_pos_3, sens_neg_3, spec_3), nm = c('mcc', 'sens_pos', 'sens_neg', 'spec'))) %>%
    mutate(method = c('method_1', 'method_2')) %>%
    mutate(alpha = names(res_v1xv2)[ix])
  
}) %>%
  bind_rows() %>%
  mutate(direction = 'v1 -> v2')

res_acc_v2xv1 <- lapply(1:length(res_v2xv1), function(ix) {
  
  sens_pos_1 <- (res_v2xv1[[ix]] %>% filter(int_true_dir == 'pos' & int_est_1 == 'interaction') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_1 <- (res_v2xv1[[ix]] %>% filter(int_true_dir == 'neg' & int_est_1 == 'interaction') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_1 <- (res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  sens_pos_3 <- (res_v2xv1[[ix]] %>% filter(int_true_dir == 'pos' & int_est_2 == 'interaction') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg_3 <- (res_v2xv1[[ix]] %>% filter(int_true_dir == 'neg' & int_est_2 == 'interaction') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec_3 <- (res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()) / (res_v2xv1[[ix]] %>% filter(int_true == 'none') %>% nrow())
  
  tp1 <- res_v2xv1[[ix]] %>% filter(int_true != 'none' & int_est_1 != 'none') %>% nrow()
  tn1 <- res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_1 == 'none') %>% nrow()
  fp1 <- res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_1 != 'none') %>% nrow()
  fn1 <- res_v2xv1[[ix]] %>% filter(int_true != 'none' & int_est_1 == 'none') %>% nrow()
  
  tp3 <- res_v2xv1[[ix]] %>% filter(int_true != 'none' & int_est_2 != 'none') %>% nrow()
  tn3 <- res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_2 == 'none') %>% nrow()
  fp3 <- res_v2xv1[[ix]] %>% filter(int_true == 'none' & int_est_2 != 'none') %>% nrow()
  fn3 <- res_v2xv1[[ix]] %>% filter(int_true != 'none' & int_est_2 == 'none') %>% nrow()
  
  mcc_1 <- mcc(tp1, tn1, fp1, fn1)
  mcc_3 <- mcc(tp3, tn3, fp3, fn3)
  
  bind_rows(setNames(c(mcc_1, sens_pos_1, sens_neg_1, spec_1), nm = c('mcc', 'sens_pos', 'sens_neg', 'spec')),
            setNames(c(mcc_3, sens_pos_3, sens_neg_3, spec_3), nm = c('mcc', 'sens_pos', 'sens_neg', 'spec'))) %>%
    mutate(method = c('method_1', 'method_2')) %>%
    mutate(alpha = names(res_v2xv1)[ix])
  
}) %>%
  bind_rows() %>%
  mutate(direction = 'v2 -> v1')

res_acc <- res_acc_v1xv2 %>%
  bind_rows(res_acc_v2xv1)
rm(res_acc_v1xv2, res_acc_v2xv1)

p2 <- ggplot(data = res_acc %>%
               pivot_longer(mcc:spec,
                            names_to = 'metric',
                            values_to = 'value') %>%
               mutate(metric = factor(metric, levels = c('mcc', 'sens_pos', 'sens_neg', 'spec'))),
             aes(x = alpha, y = value, group = method, col = method)) +
  geom_line() +
  geom_point() +
  facet_grid(direction ~ metric) +
  theme_classic() +
  theme(legend.position = 'right') +
  labs(x = 'alpha', y = 'Value', col = '') +
  scale_color_brewer(palette = 'Set1')
print(p2)
rm(p2)

# Calculate accuracy by interaction parameter values:
acc_by_param_v1xv2 <- lapply(1:length(res_v1xv2), function(ix) {
  res_v1xv2[[ix]] %>%
    mutate(int_est = int_est_1) %>%
    calculate_accuracy_matrix() %>%
    mutate(alpha = as.numeric(str_remove(names(res_v1xv2)[ix], 'alpha_')))
}) %>%
  bind_rows()
acc_by_param_v2xv1 <- lapply(1:length(res_v2xv1), function(ix) {
  res_v2xv1[[ix]] %>%
    mutate(int_est = int_est_1) %>%
    calculate_accuracy_matrix() %>%
    mutate(alpha = as.numeric(str_remove(names(res_v2xv1)[ix], 'alpha_')))
}) %>%
  bind_rows()

p3.1 <- ggplot(data = acc_by_param_v1xv2, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  facet_wrap(~ alpha, nrow = 1) +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1))
p3.2 <- ggplot(data = acc_by_param_v2xv1, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  facet_wrap(~ alpha, nrow = 1) +
  theme_classic() +
  scale_fill_viridis(option = 'G', limits = c(0, 1))
grid.arrange(p3.1, p3.2, ncol = 1)
rm(p3.1, p3.2)

dev.off()

# ------------------------------------------------------------------------------

# Clean up:
rm(list = ls())
