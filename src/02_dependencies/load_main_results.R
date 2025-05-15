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
  mutate(int_est = if_else(cor_median > 0, 'pos', 'neg'),
         int_est = if_else(CI_lower95 > 0 | CI_upper95 < 0, int_est, 'none'))

res_gam <- res_gam %>%
  # filter(!if_any(b_V1obsln_Intercept:rescor__V1obsln__V2obsln, ~ . > 1.01)) %>%
  filter(!if_any(b_V1obsln_Intercept:lp__, ~ . > 1.01)) %>%
  filter(n_div == 0) %>%
  select(run:.id, cor_median:CI_upper95, theta_lambda:int_est)
print(table(res_gam$run))

# Granger causality:
to_remove <- res_granger %>%
  filter(adf_p >= 0.05 | kpss_p < 0.05) %>%
  select(run:.id) %>%
  distinct() %>%
  mutate(delete = TRUE)

res_granger <- res_granger %>% left_join(to_remove,
                                         by = c('run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(-delete)

res_corr<- res_corr %>% left_join(to_remove %>% mutate(.id = as.integer(.id)),
                       by = c('run', '.id')) %>%
  filter(is.na(delete))
res_gam <- res_gam %>% left_join(to_remove %>% mutate(.id = as.integer(.id)),
                       by = c('run', '.id')) %>%
  filter(is.na(delete))

res_granger <- res_granger %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  ungroup() %>%
  mutate(int_est = if_else(ftest_p < 0.05, 'interaction', 'none'))

res_granger <- res_granger %>%
  select(run:.id, direction:confounding, logRSS, ftest_p, theta_lambda:int_est)

res_granger_LIST <- vector('list', length = 4)
names(res_granger_LIST) <- c('v1 -> v2 (No confounding)', 'v2 -> v1 (No confounding)', 'v1 -> v2 (Seasonality Controlled)', 'v2 -> v1 (Seasonality Controlled)')

res_granger_LIST[[1]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'none')
res_granger_LIST[[2]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'none')
res_granger_LIST[[3]] <- res_granger %>% filter(direction == 'v1 -> v2', confounding == 'seasonal')
res_granger_LIST[[4]] <- res_granger %>% filter(direction == 'v2 -> v1', confounding == 'seasonal')

rm(res_granger)

# Transfer entropy:
res_te <- res_te %>% left_join(to_remove %>% mutate(.id = as.numeric(.id)),
                               by = c('run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(-delete)

res_te <- res_te %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  ungroup() %>%
  mutate(int_est = if_else(p_value < 0.05 & te > 0, 'interaction', 'none'),
         int_est_confound = if_else(p_value_confound < 0.05 & te_confound > 0, 'interaction', 'none'))

# Keep only lag with highest TE for each synthetic dataset:
res_te_raw <- res_te %>% group_by(direction, run, .id) %>% filter(te == max(te)) %>% ungroup() %>% select(-contains('confound'))
res_te_confound1 <- res_te %>% group_by(direction, run, .id) %>% filter(te_confound == max(te_confound)) %>% ungroup() %>% select(-c(te, te_confound2, sd_null, sd_null_confound2, p_value, p_value_confound2, int_est))

res_te_LIST <- vector('list', length = 4)
names(res_te_LIST) <- c('v1 -> v2 (No confounding)', 'v2 -> v1 (No confounding)', 'v1 -> v2 (Seasonality Controlled)', 'v2 -> v1 (Seasonality Controlled)')

res_te_LIST[[1]] <- res_te_raw %>% filter(direction == 'v1 -> v2')
res_te_LIST[[2]] <- res_te_raw %>% filter(direction == 'v2 -> v1')

res_te_confound1 <- res_te_confound1 %>% rename_with(~ str_remove(.x, pattern = '_confound'))
res_te_LIST[[3]] <- res_te_confound1 %>% filter(direction == 'v1 -> v2')
res_te_LIST[[4]] <- res_te_confound1 %>% filter(direction == 'v2 -> v1')

# Clean up:
res_te_LIST <- lapply(res_te_LIST, function(ix) {
  ix %>% select(run, .id, lag, direction, te, p_value, theta_lambda:int_est)
})

rm(res_te, res_te_raw, res_te_confound1)

# CCM:
res_ccm <- res_ccm %>%
  group_by(run, .id, direction, theta_lambda, delta) %>%
  select(run:direction, rho_median, tp_opt:p_surr, theta_lambda:delta) %>%
  summarise(rho_mean = mean(rho_median), rho_max = rho_median[LibSize == max(LibSize)], tp_opt = unique(tp_opt), max_cmc = unique(max_cmc), p_conv = unique(p_conv), p_surr = unique(p_surr)) %>%
  ungroup()

res_ccm <- res_ccm %>%
  mutate(int_true = if_else(theta_lambda == 1, 'none', 'interaction'),
         int_true_dir = if_else(theta_lambda > 1, 'pos', int_true),
         int_true_dir = if_else(theta_lambda < 1, 'neg', int_true_dir)) %>%
  mutate(int_est_1 = if_else(p_surr < 0.05, 'interaction', 'none'), # method 1: check p-values based on surrogates
         int_est_2 = if_else(p_conv < 0.05 & max_cmc > 0, 'interaction', 'none'), # method 2: check convergence
         int_est_3 = if_else(p_conv < 0.05 & max_cmc > 0 & tp_opt < 0, 'interaction', 'none')) # method 3: check convergence + ideal tp negative

res_ccm <- res_ccm %>%
  select(run:direction, rho_max, p_surr, p_conv, tp_opt, theta_lambda:delta, int_true:int_est_3)

res_ccm <- res_ccm %>% left_join(to_remove, by = c('run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(-delete)

res_ccm_LIST <- vector('list', length = 6)
names(res_ccm_LIST) <- c('v1 -> v2 (Method 1)', 'v1 -> v2 (Method 2)', 'v1 -> v2 (Method 3)', 'v2 -> v1 (Method 1)', 'v2 -> v1 (Method 2)', 'v2 -> v1 (Method 3)')

res_ccm_LIST[[1]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[2]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[3]] <- res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')
res_ccm_LIST[[4]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_2, int_est_3)) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[5]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_3)) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[6]] <- res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-c(int_est_1, int_est_2)) %>% rename('int_est' = 'int_est_3')
