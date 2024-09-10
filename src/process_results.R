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

# Load functions:
source('src/functions_etc/fxns_process_results.R')

# Open pdf to save plots:
date <- format(Sys.Date(), '%d%m%y')
# pdf(file = paste0('results/plots/plot_accuracy_by_method_', date, '_UPDATED.pdf'),
#     width = 16, height = 12)

# ------------------------------------------------------------------------------

# Read in all results
source('src/functions_etc/load_main_results.R')

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

# Process accuracy of results (correlation coefficients)

# Calculate sensitivity/specificity (overall):
sens_pos <- (res_corr %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_corr %>% filter(int_true == 'pos') %>% nrow())
sens_neg <- (res_corr %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_corr %>% filter(int_true == 'neg') %>% nrow())
spec <- (res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_corr %>% filter(int_true == 'none') %>% nrow())

print('Sensitivity (Any Interaction) (Overall):')
print((res_corr %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()) / (res_corr %>% filter(int_true != 'none') %>% nrow()))

print('Sensitivity (Correct Direction) (Overall):')
print(sens_pos)
print(sens_neg)

print('Specificity (Overall):')
print(spec)

# Calculate overall accuracy (weighted):
acc_weighted_corr <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)

print('Overall accuracy (weighted):')
print(acc_weighted_corr)

# Calculate Matthews correlation coefficient (MCC):
tp <- res_corr %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()
tn <- res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()
fp <- res_corr %>% filter(int_true == 'none' & int_est != 'none') %>% nrow()
fn <- res_corr %>% filter(int_true != 'none' & int_est == 'none') %>% nrow()

print('MCC:')
print(mcc(tp, tn, fp, fn))

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
rm(p.corr.1, p.corr.2, sens_pos, sens_neg, spec, res_corr)
# rm(p.corr.1, p.corr.2, p.corr.3, res_corr, acc_corr, assoc_corr)

# ------------------------------------------------------------------------------

# Process accuracy of results (GAMs)

# CHECK NEW CONFIDENCE INTERVALS:
res_gam %>%
  select(run, .id, cor, CI_lower95, CI_upper95, cor_confound, CI_lower95_confound, CI_upper95_confound) %>%
  mutate(cor_within = (cor > CI_lower95 & cor < CI_upper95),
         cor_confound_within = (cor_confound > CI_lower95_confound & cor_confound < CI_upper95_confound)) %>%
  summarise(cor_within = sum(cor_within),
            cor_confound_within = sum(cor_confound_within))

# CHECK RESULTS ON LOG-TRANSFORMED DATA:
filelist_temp <- list.files(path = 'results/gam_log_transform/', full.names = TRUE)

res_gam_log <- vector('list', length = length(filelist_temp))
for (i in 1:length(filelist_temp)) {
  res_gam_log[[i]] <- read_rds(filelist_temp[i])
}
rm(filelist_temp, i)

res_gam_log <- bind_rows(res_gam_log) %>%
  mutate(int_est = if_else(cor > 0, 'pos', 'neg'),
         int_est = if_else(CI_lower95 > 0 | CI_upper95 < 0, int_est, 'none')) %>%
  mutate(int_est_confound = if_else(cor_confound > 0, 'pos', 'neg'),
         int_est_confound = if_else(CI_lower95_confound > 0 | CI_upper95_confound < 0, int_est_confound, 'none')) %>%
  mutate(cor_within = (cor > CI_lower95 & cor < CI_upper95),
         cor_confound_within = (cor_confound > CI_lower95_confound & cor_confound < CI_upper95_confound))

summary(res_gam_log$cor_within)
summary(res_gam_log$cor_confound_within)

res_gam_check <- res_gam %>%
  left_join(res_gam_log %>% select(run, .id, int_est:int_est_confound),
            by = c('run', '.id'))

table(res_gam_check$int_true, res_gam_check$int_est.x)
table(res_gam_check$int_true, res_gam_check$int_est.y)

table(res_gam_check$int_true, res_gam_check$int_est_confound.x)
table(res_gam_check$int_true, res_gam_check$int_est_confound.y)

sens_pos <- (res_gam_check %>% filter(int_true == 'pos' & int_est.y == 'pos') %>% nrow()) / (res_gam_check %>% filter(int_true == 'pos') %>% nrow())
sens_neg <- (res_gam_check %>% filter(int_true == 'neg' & int_est.y == 'neg') %>% nrow()) / (res_gam_check %>% filter(int_true == 'neg') %>% nrow())
sens_pos_confound <- (res_gam_check %>% filter(int_true == 'pos' & int_est_confound.y == 'pos') %>% nrow()) / (res_gam_check %>% filter(int_true == 'pos') %>% nrow())
sens_neg_confound <- (res_gam_check %>% filter(int_true == 'neg' & int_est_confound.y == 'neg') %>% nrow()) / (res_gam_check %>% filter(int_true == 'neg') %>% nrow())
spec <- (res_gam_check %>% filter(int_true == 'none' & int_est.y == 'none') %>% nrow()) / (res_gam_check %>% filter(int_true == 'none') %>% nrow())
spec_confound <- (res_gam_check %>% filter(int_true == 'none' & int_est_confound.y == 'none') %>% nrow()) / (res_gam_check %>% filter(int_true == 'none') %>% nrow())

print((sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null))
print((sens_pos_confound * weight_pos + sens_neg_confound * weight_neg + spec_confound * weight_null) / (weight_pos + weight_neg + weight_null))

rm(res_gam_log, res_gam_check)

# Calculate sensitivity/specificity (overall):
sens_pos <- (res_gam %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_gam %>% filter(int_true == 'pos') %>% nrow())
sens_neg <- (res_gam %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_gam %>% filter(int_true == 'neg') %>% nrow())
sens_pos_confound <- (res_gam %>% filter(int_true == 'pos' & int_est_confound == 'pos') %>% nrow()) / (res_gam %>% filter(int_true == 'pos') %>% nrow())
sens_neg_confound <- (res_gam %>% filter(int_true == 'neg' & int_est_confound == 'neg') %>% nrow()) / (res_gam %>% filter(int_true == 'neg') %>% nrow())
spec <- (res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_gam %>% filter(int_true == 'none') %>% nrow())
spec_confound <- (res_gam %>% filter(int_true == 'none' & int_est_confound == 'none') %>% nrow()) / (res_gam %>% filter(int_true == 'none') %>% nrow())

print('Sensitivity (Any Interaction) (Overall):')
print((res_gam %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()) / (res_gam %>% filter(int_true != 'none') %>% nrow()))
print((res_gam %>% filter(int_true != 'none' & int_est_confound != 'none') %>% nrow()) / (res_gam %>% filter(int_true != 'none') %>% nrow()))

print('Sensitivity (Correct Direction) (Overall):')
print(sens_pos)
print(sens_neg)
print(sens_pos_confound)
print(sens_neg_confound)

print('Specificity (Overall):')
print(spec)
print(spec_confound)

# Calculate overall accuracy (weighted):
acc_weighted_gam <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
acc_weighted_gam_confound <- (sens_pos_confound * weight_pos + sens_neg_confound * weight_neg + spec_confound * weight_null) / (weight_pos + weight_neg + weight_null)

print('Overall accuracy (weighted):')
print(acc_weighted_gam)
print(acc_weighted_gam_confound)

# Calculate Matthews correlation coefficient (MCC):
tp <- res_gam %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()
tn <- res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()
fp <- res_gam %>% filter(int_true == 'none' & int_est != 'none') %>% nrow()
fn <- res_gam %>% filter(int_true != 'none' & int_est == 'none') %>% nrow()

tp_confound <- res_gam %>% filter(int_true != 'none' & int_est_confound != 'none') %>% nrow()
tn_confound <- res_gam %>% filter(int_true == 'none' & int_est_confound == 'none') %>% nrow()
fp_confound <- res_gam %>% filter(int_true == 'none' & int_est_confound != 'none') %>% nrow()
fn_confound <- res_gam %>% filter(int_true != 'none' & int_est_confound == 'none') %>% nrow()

print('MCC:')
print(mcc(tp, tn, fp, fn))
print(mcc(tp_confound, tn_confound, fp_confound, fn_confound))

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

rm(p.gam.1.1, p.gam.1.2, sens_pos, sens_pos_confound,
   sens_neg, sens_neg_confound, spec, spec_confound, res_gam)
# rm(p.gam.1.1, p.gam.1.2, p.gam.2.1, p.gam.2.2, res_gam, acc_gam, acc_gam_confound, assoc_gam, assoc_gam_confound)

# ------------------------------------------------------------------------------

# Process accuracy of results (Granger causality)

# Check whether any datasets are not stationary:
res_granger %>% filter(adf_p >= 0.05) %>% nrow() %>% print()
res_granger %>% filter(kpss_p < 0.05) %>% nrow() %>% print()

# Calculate sensitivity/specificity (overall):
acc_weighted_granger <- vector('list', length = length(res_granger_LIST))
names(acc_weighted_granger) <- names(res_granger_LIST)
for (i in 1:length(res_granger_LIST)) {
  
  print(names(res_granger_LIST)[i])
  
  sens_pos <- (res_granger_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg <- (res_granger_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec <- (res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  print('Sensitivity:')
  print(sens_pos)
  print(sens_neg)
  
  print('Specificity:')
  print(spec)
  
  acc_weighted_temp <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  
  print('Overall accuracy (weighted):')
  print(acc_weighted_temp)
  
  acc_weighted_granger[[i]] <- acc_weighted_temp
  rm(acc_weighted_temp, sens_pos, sens_neg, spec)
  
  # Calculate Matthews correlation coefficient (MCC):
  tp <- res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()
  tn <- res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()
  fp <- res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow()
  fn <- res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow()
  
  print('MCC:')
  print(mcc(tp, tn, fp, fn))
  
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

# Calculate sensitivity/specificity (overall):
acc_weighted_te <- vector('list', length = length(res_te_LIST))
names(acc_weighted_te) <- names(res_te_LIST)

for (i in 1:length(res_te_LIST)) {
  
  print(names(res_te_LIST)[i])
  
  sens_pos <- (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg <- (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec <- (res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  print('Sensitivity:')
  print(sens_pos)
  print(sens_neg)
  
  print('Specificity:')
  print(spec)
  
  acc_weighted_temp <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  
  print('Overall accuracy (weighted):')
  print(acc_weighted_temp)
  
  acc_weighted_te[[i]] <- acc_weighted_temp
  rm(acc_weighted_temp, sens_pos, sens_neg, spec)
  
  # Calculate Matthews correlation coefficient (MCC):
  tp <- res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()
  tn <- res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()
  fp <- res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow()
  fn <- res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow()
  
  print('MCC:')
  print(mcc(tp, tn, fp, fn))
  
  print('-----------')
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_te_LIST <- vector('list', length = 2)
names(acc_te_LIST) <- names(res_te_LIST)

for (i in 1:length(acc_te_LIST)) {
  acc_te_LIST[[i]] <- calculate_accuracy_matrix(res_te_LIST[[i]])
}

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

rm(p.te.1.1, p.te.1.2, res_te_LIST, res_te, best_v1xv2, best_v2xv1)
# rm(p.te.1.1, p.te.1.2, p.te.2.1, p.te.2.2, res_te, res_te_LIST, acc_te_LIST, assoc_te_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (CCM)

# Calculate sensitivity/specificity (overall):
acc_weighted_ccm <- vector('list', length = length(res_ccm_LIST))
names(acc_weighted_ccm) <- names(res_ccm_LIST)

for (i in 1:length(res_ccm_LIST)) {
  
  print(names(res_ccm_LIST)[i])
  
  sens_pos <- (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow())
  sens_neg <- (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow())
  spec <- (res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  print('Sensitivity:')
  print(sens_pos)
  print(sens_neg)
  
  print('Specificity:')
  print(spec)
  
  acc_weighted_temp <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  
  print('Overall accuracy (weighted):')
  print(acc_weighted_temp)
  
  acc_weighted_ccm[[i]] <- acc_weighted_temp
  rm(acc_weighted_temp, sens_pos, sens_neg, spec)
  
  # Calculate Matthews correlation coefficient (MCC):
  tp <- res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow()
  tn <- res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()
  fp <- res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow()
  fn <- res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow()
  
  print('MCC:')
  print(mcc(tp, tn, fp, fn))
  
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
assoc_ccm_LIST_max <- vector('list', length = 6)
names(assoc_ccm_LIST_max) <- names(res_ccm_LIST)

res_ccm_LIST <- lapply(res_ccm_LIST, function(ix) {
  ix %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no'))
})

for (i in 1:length(assoc_ccm_LIST_max)) {
  assoc_ccm_LIST_max[[i]] <- calculate_assoc_true_strength(res_ccm_LIST[[i]], method = 'ccm', met = 'rho_max')
}
rm(i)

# CHECK USING MEDIANS:
res_ccm_LIST_check = assoc_ccm_LIST_max_check = vector('list', length = 4)

res_ccm_LIST_check[[1]] <- res_ccm_check %>% filter(direction == 'v1 -> v2') %>% select(-int_est_3) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST_check[[2]] <- res_ccm_check %>% filter(direction == 'v1 -> v2') %>% select(-int_est_1) %>% rename('int_est' = 'int_est_3')
res_ccm_LIST_check[[3]] <- res_ccm_check %>% filter(direction == 'v2 -> v1') %>% select(-int_est_3) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST_check[[4]] <- res_ccm_check %>% filter(direction == 'v2 -> v1') %>% select(-int_est_1) %>% rename('int_est' = 'int_est_3')

res_ccm_LIST_check <- lapply(res_ccm_LIST_check, function(ix) {
  ix %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no'))
})

for (i in 1:length(assoc_ccm_LIST_max_check)) {
  assoc_ccm_LIST_max_check[[i]] <- calculate_assoc_true_strength(res_ccm_LIST_check[[i]], method = 'ccm', met = 'rho_max')
}
rm(i)

print(assoc_ccm_LIST_max[c(1, 3, 4, 6)] %>% bind_rows() %>% filter(rho > 0 & p_value < 0.05))
print(assoc_ccm_LIST_max_check %>% bind_rows() %>% filter(rho > 0 & p_value < 0.05))

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

rm(p.ccm.1.1, p.ccm.1.2, p.ccm.1.3, p.ccm.1.4, p.ccm.1.5, p.ccm.1.6, res_ccm_LIST,
   weight_pos, weight_neg, weight_null, res_ccm)
# rm(p.ccm.1.1, p.ccm.1.2, p.ccm.1.3, p.ccm.1.4, p.ccm.1.5, p.ccm.1.6,
#    p.ccm.3.1, p.ccm.3.2, p.ccm.3.3, p.ccm.3.4, p.ccm.3.5, p.ccm.3.6,
#    res_ccm, res_ccm_LIST, acc_ccm_LIST, assoc_ccm_LIST_mean, assoc_ccm_LIST_max)

# ------------------------------------------------------------------------------

# Compile/plot results of all methods

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
ggsave(filename = 'results/plots/heatmaps.svg', p.comb.1, width = 13, height = 18)
rm(p.corr.3,p.gam.2.1, p.gam.2.2, p.granger.3.1, p.granger.3.2, p.granger.3.3, p.granger.3.4,
   p.te.2.1, p.te.2.2, p.ccm.3.1, p.ccm.3.4, p.ccm.3.2, p.ccm.3.5, p.ccm.3.3, p.ccm.3.6)

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
  scale_shape_manual(values = c(18, 17, 17, 15, 15, 3, 8, 8)) +
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
  scale_shape_manual(values = c(18, 17, 17, 15, 15, 3, 18, 18)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = 'Interaction Duration', y = "Spearman's Rho", shape = 'Method', col = 'Method')
print(p.comb.3)
rm(assoc_corr, assoc_gam, assoc_gam_confound, assoc_granger_LIST, assoc_te_LIST, assoc_ccm_LIST_mean, assoc_ccm_LIST_max)

# Calculate measure of overall accuracy, weighted by number of simulations per true parameter set:
# Note that correlation, GAMs, and gradient boosting will have some disadvantage here, since it matters
# for them whether the predicted interaction is positive or negative

res_acc_weighted <- bind_cols(method = 'Corr. Coef.', weighted_acc = acc_weighted_corr, direction = 'v1 -> v2') %>%
  bind_rows(bind_cols(method = 'Corr. Coef.', weighted_acc = acc_weighted_corr, direction = 'v2 -> v1')) %>%
  bind_rows(bind_cols(method = 'GAMs', weighted_acc = acc_weighted_gam, direction = 'v1 -> v2')) %>%
  bind_rows(bind_cols(method = 'GAMs', weighted_acc = acc_weighted_gam, direction = 'v2 -> v1')) %>%
  bind_rows(bind_cols(method = 'GAMs (w/ Seas.)', weighted_acc = acc_weighted_gam_confound, direction = 'v1 -> v2')) %>%
  bind_rows(bind_cols(method = 'GAMs (w/ Seas.)', weighted_acc = acc_weighted_gam_confound, direction = 'v2 -> v1')) %>%
  bind_rows(acc_weighted_granger %>%
              bind_rows() %>%
              pivot_longer(everything(), values_to = 'weighted_acc') %>%
              mutate(direction = str_sub(name, 1, 8),
                     method = if_else(str_detect(name, 'Seasonality'), 'GC (w/ Seas.)', 'GC')) %>%
              select(method, weighted_acc, direction)) %>%
  bind_rows(acc_weighted_te %>%
              bind_rows() %>%
              pivot_longer(everything(), names_to = 'direction', values_to = 'weighted_acc') %>%
              mutate(direction = str_sub(direction, 1, 8),
                     method = 'TE') %>%
              select(method, weighted_acc, direction)) %>%
  bind_rows(acc_weighted_ccm %>%
              bind_rows() %>%
              pivot_longer(everything(), values_to = 'weighted_acc') %>%
              mutate(direction = str_sub(name, 1, 8),
                     method = paste0('CCM ', str_sub(name, 10))) %>%
              select(method, weighted_acc, direction) %>%
              filter(!str_detect(method, '2')))

res_acc_weighted <- res_acc_weighted %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GAMs (w/ Seas.)', 'GC', 'GC (w/ Seas.)', 'TE', 'CCM (Method 1)', 'CCM (Method 3)')))

p.comb.4 <- ggplot(res_acc_weighted, aes(x = direction, y = weighted_acc, shape = method, col = method)) +
  geom_point(size = 3) +
  # facet_wrap(~ direction, nrow = 1) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 0.75), n.breaks = 20) +
  scale_shape_manual(values = c(18, 17, 17, 15, 15, 3, 8, 8, 16)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3', '#6a3d9a', '#33a02c', '#b2df8a', '#000000')) +
  labs(x = 'Direction', y = 'Accuracy (Weighted)', shape = 'Method', col = 'Method')
print(p.comb.4)
rm(acc_weighted_corr, acc_weighted_gam, acc_weighted_gam_confound,
   acc_weighted_granger, acc_weighted_te, acc_weighted_ccm)

# Close pdf:
dev.off()

# Clean up:
rm(list = ls())
