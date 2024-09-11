# ------------------------------------------------------------------------------
# Code to read in and format all "main" results
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

# Calculate weights for weighted accuracy measurement

# Classify runs as positive, negative, or null:
df_true_int <- dat %>%
  select(run, .id, theta_lambda:delta) %>%
  distinct() %>%
  mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
         int_true = if_else(theta_lambda == 1, 'none', int_true))

# Get weights based on relative amount of simulations with pos/neg/no interaction:
weight_pos <- min(table(df_true_int$int_true)) / (df_true_int %>% filter(int_true == 'pos') %>% nrow())
weight_neg <- min(table(df_true_int$int_true)) / (df_true_int %>% filter(int_true == 'neg') %>% nrow())
weight_null <- min(table(df_true_int$int_true)) / (df_true_int %>% filter(int_true == 'none') %>% nrow())

# Clean up:
rm(df_true_int)

# ------------------------------------------------------------------------------

# Determine significance/direction of true and detected interactions

# Correlation coefficients:
res_corr <- res_corr %>%
  mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
         int_true = if_else(theta_lambda == 1, 'none', int_true)) %>%
  mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
         int_est = if_else(p_value < 0.05, int_est, 'none'))

# GAMs:
res_gam <- res_gam %>%
  mutate(int_true = if_else(theta_lambda > 1, 'pos', 'neg'),
         int_true = if_else(theta_lambda == 1, 'none', int_true)) %>%
  mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
         int_est = if_else(CI_lower95 > 0 | CI_upper95 < 0, int_est, 'none')) %>%
  mutate(int_est_confound = if_else(cor_confound > 0, 'pos', 'neg'),
         int_est_confound = if_else(CI_lower95_confound > 0 | CI_upper95_confound < 0, int_est_confound, 'none'))

# Granger causality:
res_granger <- res_granger %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  ungroup() %>%
  mutate(int_est = if_else(ftest_p < 0.05, 'interaction', 'none'))

res_granger_LIST <- vector('list', length = 4)
names(res_granger_LIST) <- c('v1 -> v2 (No confounding)', 'v2 -> v1 (No confounding)', 'v1 -> v2 (Seasonality Controlled)', 'v2 -> v1 (Seasonality Controlled)')

res_granger_LIST[[1]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'none')
res_granger_LIST[[2]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'none')
res_granger_LIST[[3]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'seasonal')
res_granger_LIST[[4]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'seasonal')

# Transfer entropy:
res_te <- res_te %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  ungroup() %>%
  mutate(int_est = if_else(p_value < 0.05, 'interaction', 'none'))

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

acc_weighted_te <- vector('list', length = length(res_te_LIST))
names(acc_weighted_te) <- names(res_te_LIST)

for (i in 1:length(res_te_LIST)) {
  
  sens_pos <- (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg <- (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec <- (res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  acc_weighted_temp <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  acc_weighted_te[[i]] <- acc_weighted_temp
  rm(acc_weighted_temp, sens_pos, sens_neg, spec)
  
}
rm(i)

# Keep only best-performing lag for each direction:
best_v1xv2 <- acc_weighted_te[1:4] %>% bind_rows() %>% which.max() %>% names()
best_v2xv1 <- acc_weighted_te[5:8] %>% bind_rows() %>% which.max() %>% names()

res_te_LIST <- res_te_LIST[c(best_v1xv2, best_v2xv1)]

# CCM:
res_ccm <- res_ccm %>%
  group_by(run, .id, direction, theta_lambda, delta) %>%
  select(run:direction, rho, MannK:p_surr, theta_lambda:delta) %>%
  summarise(rho_mean = mean(rho), rho_max = rho[LibSize == max(LibSize)], MannK = unique(MannK), tp_opt = unique(tp_opt), max_cmc = unique(max_cmc), p_surr = unique(p_surr)) %>%
  ungroup()

res_ccm <- res_ccm %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
         int_est_2 = if_else(MannK < 0.05 & max_cmc > 0, 'interaction', 'none'), # method 2: check convergence
         int_est_3 = if_else(MannK < 0.05 & max_cmc > 0 & tp_opt < 0, 'interaction', 'none')) # method 3: check convergence + ideal tp negative

res_ccm_LIST <- vector('list', length = 6)
names(res_ccm_LIST) <- c('v1 -> v2 (Method 1)', 'v1 -> v2 (Method 2)', 'v1 -> v2 (Method 3)', 'v2 -> v1 (Method 1)', 'v2 -> v1 (Method 2)', 'v2 -> v1 (Method 3)')

res_ccm_LIST[[1]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[2]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[3]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')
res_ccm_LIST[[4]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[5]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[6]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')
