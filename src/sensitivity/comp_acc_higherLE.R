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

# Plot sensitivity by true interaction strength/duration:
acc_granger <- bind_rows(acc_granger_LIST, .id = 'name') %>%
  mutate(direction = str_sub(name, 1, 8),
         analysis = if_else(str_detect(name, 'Main'), 'Main', 'Sens')) %>%
  select(-name)
acc_te <- bind_rows(acc_te_LIST, .id = 'name') %>%
  mutate(direction = str_sub(name, 1, 8),
         analysis = if_else(str_detect(name, 'Main'), 'Main', 'Sens')) %>%
  select(-name)
acc_ccm <- bind_rows(acc_ccm_LIST, .id = 'name') %>%
  mutate(direction = str_sub(name, 1, 8),
         method = if_else(str_detect(name, 'Method 1'), 'Method 1', 'Method 2'),
         analysis = if_else(str_detect(name, 'Main'), 'Main', 'Sens')) %>%
  select(-name)

x_lab <- textGrob('Strength', gp = gpar(fontsize = 14, hjust = 1))
y_lab <- textGrob('Sensitivity', rot = 90, hjust = -0, gp = gpar(fontsize = 14))

title_1 <- textGrob('Granger Causality (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.85)
title_2 <- textGrob('Granger Causality (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.85)
title_3 <- textGrob('Transfer Entropy (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_4 <- textGrob('Transfer Entropy (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)
title_5 <- textGrob('CCM (Method 1) (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_6 <- textGrob('CCM (Method 1) (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)
title_7 <- textGrob('CCM (Method 2) (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_8 <- textGrob('CCM (Method 2) (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)

p.granger.1a <- ggplot(data = acc_granger %>%
                         filter(direction == 'V1 -> V2') %>%
                         mutate(strength = as.numeric(strength)) %>% 
                         filter(strength > '1'),
                       aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.granger.1b <- ggplot(data = acc_granger %>%
                         filter(direction == 'V1 -> V2') %>%
                         mutate(strength = as.numeric(strength)) %>% 
                         filter(strength < 1) %>%
                         mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                       aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.granger.2a <- ggplot(data = acc_granger %>%
                         filter(direction == 'V2 -> V1') %>%
                         mutate(strength = as.numeric(strength)) %>% 
                         filter(strength > '1'),
                       aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.granger.2b <- ggplot(data = acc_granger %>%
                         filter(direction == 'V2 -> V1') %>%
                         mutate(strength = as.numeric(strength)) %>% 
                         filter(strength < 1) %>%
                         mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                       aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.te.1a <- ggplot(data = acc_te %>%
                    filter(direction == 'V1 -> V2') %>%
                    mutate(strength = as.numeric(strength)) %>% 
                    filter(strength > '1'),
                  aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.te.1b <- ggplot(data = acc_te %>%
                    filter(direction == 'V1 -> V2') %>%
                    mutate(strength = as.numeric(strength)) %>% 
                    filter(strength < 1) %>%
                    mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                  aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.te.2a <- ggplot(data = acc_te %>%
                    filter(direction == 'V2 -> V1') %>%
                    mutate(strength = as.numeric(strength)) %>% 
                    filter(strength > '1'),
                  aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.te.2b <- ggplot(data = acc_te %>%
                    filter(direction == 'V2 -> V1') %>%
                    mutate(strength = as.numeric(strength)) %>% 
                    filter(strength < 1) %>%
                    mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                  aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.ccm.1a <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V1 -> V2', method == 'Method 1') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength > '1'),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.ccm.1b <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V1 -> V2', method == 'Method 1') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength < 1) %>%
                     mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.ccm.2a <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V2 -> V1', method == 'Method 1') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength > '1'),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.ccm.2b <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V2 -> V1', method == 'Method 1') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength < 1) %>%
                     mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.ccm.3a <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V1 -> V2', method == 'Method 2') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength > '1'),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.ccm.3b <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V1 -> V2', method == 'Method 2') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength < 1) %>%
                     mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.ccm.4a <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V2 -> V1', method == 'Method 2') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength > '1'),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(x = NULL, y = '', title = NULL)
p.ccm.4b <- ggplot(data = acc_ccm %>%
                     filter(direction == 'V2 -> V1', method == 'Method 2') %>%
                     mutate(strength = as.numeric(strength)) %>% 
                     filter(strength < 1) %>%
                     mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
  facet_wrap(~ analysis, ncol = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4, 8), labels = c('0.5', '0.25', '0')) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_shape_manual(values = c(1, 2, 0)) +
  labs(x = NULL, y = '', title = NULL)

p.granger <- arrangeGrob(arrangeGrob(p.granger.1a, p.granger.1b, nrow = 1, widths = c(0.8, 1), top = title_1),
                         arrangeGrob(p.granger.2a, p.granger.2b, nrow = 1, widths = c(0.8, 1), top = title_2),
                         ncol = 2, bottom = x_lab, left = y_lab)
p.te <- arrangeGrob(arrangeGrob(p.te.1a, p.te.1b, nrow = 1, widths = c(0.8, 1), top = title_3),
                    arrangeGrob(p.te.2a, p.te.2b, nrow = 1, widths = c(0.8, 1), top = title_4),
                    ncol = 2, bottom = x_lab, left = y_lab)
p.ccm.2 <- arrangeGrob(arrangeGrob(p.ccm.1a, p.ccm.1b, nrow = 1, widths = c(0.8, 1), top = title_7),
                       arrangeGrob(p.ccm.2a, p.ccm.2b, nrow = 1, widths = c(0.8, 1), top = title_8),
                       ncol = 2, bottom = x_lab, left = y_lab)
p.ccm.1 <- arrangeGrob(arrangeGrob(p.ccm.3a, p.ccm.3b, nrow = 1, widths = c(0.8, 1), top = title_5),
                       arrangeGrob(p.ccm.4a, p.ccm.4b, nrow = 1, widths = c(0.8, 1), top = title_6),
                       ncol = 2, bottom = x_lab, left = y_lab)

plot(p.granger)
plot(p.te)
plot(p.ccm.1)
plot(p.ccm.2)
