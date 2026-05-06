# ------------------------------------------------------------------------------
# Code to compare main analysis results with results using higher values for
# lags, history lengths, and embedding dimensions
# ------------------------------------------------------------------------------

# Setup

# Load packages
library(tidyverse)
library(gridExtra)
library(grid)
library(testthat)
library(viridis)

# Load functions:
source('src/02_dependencies/fxns_process_results.R')

# ------------------------------------------------------------------------------ 

# Read in all results

# Sensitivity analysis results:
sensLagEmbedding <- TRUE
source('src/02_dependencies/load_main_results.R')

# Main results:
res_main <- read_rds('results/res_compiled.rds')

# Combine and format data:
res_granger_LIST[[1]] <- res_main[[3]] %>% filter(direction == 'v1 -> v2')
res_granger_LIST[[2]] <- res_main[[3]] %>% filter(direction == 'v2 -> v1')
names(res_granger_LIST) <- c('V1 -> V2 (Main)', 'V2 -> V1 (Main)', 'V1 -> V2 (Sens)', 'V2 -> V1 (Sens)')

res_te_LIST[[1]] <- res_main[[4]] %>% filter(direction == 'v1 -> v2')
res_te_LIST[[2]] <- res_main[[4]] %>% filter(direction == 'v2 -> v1')
names(res_te_LIST) <- c('V1 -> V2 (Main)', 'V2 -> V1 (Main)', 'V1 -> V2 (Sens)', 'V2 -> V1 (Sens)')

res_ccm_LIST_TEMP <- Filter(Negate(is.null), res_ccm_LIST)
res_ccm_LIST <- vector('list', length = 8)

res_ccm_LIST[[1]] <- res_main[[5]] %>% filter(direction == 'v1 -> v2') %>% select(run:p_conv, theta_lambda:int_true_dir, int_est_1) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[2]] <- res_main[[5]] %>% filter(direction == 'v1 -> v2') %>% select(run:p_conv, theta_lambda:int_true_dir, int_est_2) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[3]] <- res_main[[5]] %>% filter(direction == 'v2 -> v1') %>% select(run:p_conv, theta_lambda:int_true_dir, int_est_1) %>% rename('int_est' = 'int_est_1')
res_ccm_LIST[[4]] <- res_main[[5]] %>% filter(direction == 'v2 -> v1') %>% select(run:p_conv, theta_lambda:int_true_dir, int_est_2) %>% rename('int_est' = 'int_est_2')
res_ccm_LIST[[5]] <- res_ccm_LIST_TEMP[[1]]
res_ccm_LIST[[6]] <- res_ccm_LIST_TEMP[[2]]
res_ccm_LIST[[7]] <- res_ccm_LIST_TEMP[[3]]
res_ccm_LIST[[8]] <- res_ccm_LIST_TEMP[[4]]
names(res_ccm_LIST) <- c('V1 -> V2 (Method 1) (Main)', 'V1 -> V2 (Method 2) (Main)', 'V2 -> V1 (Method 1) (Main)', 'V2 -> V1 (Method 2) (Main)',
                         'V1 -> V2 (Method 1) (Sens)', 'V1 -> V2 (Method 2) (Sens)', 'V2 -> V1 (Method 1) (Sens)', 'V2 -> V1 (Method 2) (Sens)')

rm(dat, res_corr, res_ccm, res_ccm_surr, res_ccm_LIST_TEMP, to_remove, res_main)

# ------------------------------------------------------------------------------

# Calculate sensitivity/specificity for all methods

# Initiate tibble to hold results:
df_acc <- as_tibble(NULL)

# Granger causality:
for (i in 1:length(res_granger_LIST)) {
  
  sens_test <- binom.test(res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_granger_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_granger_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('GC ', names(res_granger_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}

acc_granger_LIST <- vector('list', length = length(res_granger_LIST))
names(acc_granger_LIST) <- names(res_granger_LIST)

for (i in 1:length(acc_granger_LIST)) {
  acc_granger_LIST[[i]] <- calculate_accuracy_matrix(res_granger_LIST[[i]])
}
rm(i)

# Transfer entropy:
for (i in 1:length(res_te_LIST)) {
  
  sens_test <- binom.test(res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_te_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('TE ', names(res_te_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}

acc_te_LIST <- vector('list', length = length(res_te_LIST))
names(acc_te_LIST) <- c(names(res_te_LIST))

for (i in 1:length(res_te_LIST)) {
  acc_te_LIST[[i]] <- calculate_accuracy_matrix(res_te_LIST[[i]])
}
rm(i)

# CCM:
for (i in 1:length(res_ccm_LIST)) {
  
  sens_test <- binom.test(res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_ccm_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_ccm_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('CCM ', names(res_ccm_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}

acc_ccm_LIST <- vector('list', length = length(res_ccm_LIST))
names(acc_ccm_LIST) <- names(res_ccm_LIST)

for (i in 1:length(acc_ccm_LIST)) {
  acc_ccm_LIST[[i]] <- calculate_accuracy_matrix(res_ccm_LIST[[i]])
}
rm(i, sens_test, spec_test)

# ------------------------------------------------------------------------------

# Plot results

# Calculate and plot sens/spec differences:
df_diff <- df_acc %>%
  mutate(analysis = if_else(str_detect(method, 'Main'), 'Main', 'Sens'),
         method = str_sub(method, 1, -8),
         .before = method) %>%
  pivot_wider(id_cols = method, names_from = analysis, values_from = sens:upper_spec) %>%
  mutate(sens_diff = sens_Sens - sens_Main,
         sens_diff_lower = lower_sens_Sens - upper_sens_Main,
         sens_diff_upper = upper_sens_Sens - lower_sens_Main,
         spec_diff = spec_Sens - spec_Main,
         spec_diff_lower = lower_spec_Sens - upper_spec_Main,
         spec_diff_upper = upper_spec_Sens - lower_spec_Main) %>%
  select(method, contains('diff')) %>%
  mutate(direction = str_sub(method, 4, 12),
         method = str_remove(method, direction),
         .after = method) %>%
  mutate(method = str_trim(method),
         direction = str_trim(direction)) %>%
  mutate(method = if_else(method == 'CCM (Method 1)', 'CCM2', method),
         method = if_else(method == 'CCM (Method 2)', 'CCM1', method)) %>%
  mutate(method = factor(method, levels = c('GC', 'TE', 'CCM1', 'CCM2')))

p.comp <- ggplot(df_diff %>% mutate(direction = if_else(direction == 'V1 -> V2', 'A %->% B', 'B %->% A')),
                 aes(x = sens_diff, y = spec_diff, shape = method, col = method)) +
  geom_hline(yintercept = 0, lty = 2, col = 'gray20') +
  geom_vline(xintercept = 0, lty = 2, col = 'gray20') +
  geom_segment(aes(x = sens_diff_lower, xend = sens_diff_upper, y = spec_diff)) +
  geom_segment(aes(y = spec_diff_lower, yend = spec_diff_upper, x = sens_diff)) +
  geom_point(size = 2.5) +
  facet_wrap(~ direction, nrow = 1, labeller = label_parsed) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray90'),
        panel.grid.minor = element_blank(),
        legend.position = 'top') +
  scale_x_continuous(limits = c(-0.1, 0.5), n.breaks = 6) +
  scale_y_continuous(limits = c(-0.55, 0.15), n.breaks = 8) +
  scale_shape_manual(values = c(15, 16, 4, 4)) +
  scale_color_manual(values = c('#1f78b4', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = expression(Delta~'Sensitivity'), y = expression(Delta~'Specificity'), shape = 'Method', color = 'Method')
plot(p.comp)
# ggsave(filename = 'results/plots/figures/FigureS4.svg', p.comp, height = 4.6, width = 7.4)
