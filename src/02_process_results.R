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
source('src/02_dependencies/fxns_process_results.R')

# Open pdf to save plots:
date <- format(Sys.Date(), '%d%m%y')
# pdf(file = paste0('results/plots/plot_accuracy_by_method_', date, '.pdf'),
#     width = 16, height = 12)

# ------------------------------------------------------------------------------ 

# Read in all results
source('src/02_dependencies/load_main_results.R')

# ------------------------------------------------------------------------------

# Visualize some datasets

# Choose five random simulations to plot:
set.seed(490275)
ids_to_plot <- sample(1:100, size = 25)

# Format data for plotting:
dat_plot <- dat %>%
  filter(.id %in% ids_to_plot) %>%
  # filter((theta_lambda == 0 & delta == 0.25) | (theta_lambda == 4 & delta == 0.25) | theta_lambda == 1) %>%
  pivot_longer(V1_obs:V2_obs, names_to = 'virus', values_to = 'obs') %>%
  mutate(virus = if_else(virus == 'V1_obs', 'Influenza', 'RSV')) %>%
  mutate(obs = obs / 5000000 * 1000)

# Plot:
year_breaks <- dat_plot %>% filter(str_detect(date, '-07-0[1-7]')) %>% pull(date) %>% unique()
year_breaks <- c(year_breaks, '2022-07-03')

for (i in c(1, 2, 8, 12, 17, 24)) {
  
  p.data <- ggplot(data = dat_plot %>%
                     filter(.id == ids_to_plot[[i]]) %>%
                     mutate(theta_lambda = as.character(theta_lambda),
                            delta = as.character(7 / delta)) %>%
                     filter((theta_lambda %in% c(0, 0.5, 2, 4) & delta %in% c(28, 91)) |
                              theta_lambda == 1) %>%
                     mutate(delta = if_else(theta_lambda == '1', '28', delta),
                            theta_lambda = paste0('Interaction Strength: ', theta_lambda)),
                   aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
    geom_line(linewidth = 0.75) +
    geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
    facet_wrap(~ theta_lambda, ncol = 1) +
    theme_classic() +
    theme(legend.position = 'none',
          axis.title = element_text(size = 13.5),
          axis.text = element_text(size = 11),
          strip.text = element_text(size = 12)) +
    scale_x_continuous(breaks = NULL) +
    scale_color_manual(values = c('#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3')) +
    labs(x = '', y = 'Incidence (per 1000)', col = '', title = i)
  print(p.data)
  
}
# ggsave(filename = 'results/plots/supp_plot1.svg', p.data, width = 11, height = 7)

rm(dat_plot, ids_to_plot)

# ------------------------------------------------------------------------------

# Process accuracy of results (correlation coefficients)

# Calculate sensitivity/specificity:
df_acc <- as_tibble(data.frame(method = 'Corr. Coef.',
                               sens_pos = (res_corr %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_corr %>% filter(int_true == 'pos') %>% nrow()),
                               sens_neg = (res_corr %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_corr %>% filter(int_true == 'neg') %>% nrow()),
                               spec = (res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_corr %>% filter(int_true == 'none') %>% nrow()),
                               mcc = mcc(res_corr %>% filter(int_true != 'none' & int_est != 'none') %>% nrow(),
                                         res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(),
                                         res_corr %>% filter(int_true == 'none' & int_est != 'none') %>% nrow(),
                                         res_corr %>% filter(int_true != 'none' & int_est == 'none') %>% nrow())))

# Calculate sensitivity/specificity (by true params):
acc_corr <- calculate_accuracy_matrix(res_corr)

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_corr <- calculate_assoc_true_strength(res_corr %>% mutate(sig = if_else(int_est == 'none', 'no', 'yes')),
                                            method = 'corr', met = 'cor')

# Plot:
p.corr.1 <- ggplot(data = acc_corr %>%
                     mutate(strength_proxy = rank(strength, ties.method = 'min')),
                   aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'A', x = 'Strength', y = '% Correct   ')

p.legend.1 <- ggplot(data = acc_corr %>%
                       mutate(duration = paste0(duration, ' week'),
                              duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                       mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))),
                     aes(x = as.numeric(strength), y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  theme_classic() +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_brewer(palette = 'Set1') +
  labs(shape = 'Duration', color = 'Duration', lty = 'Duration')
p.legend.1 <- ggplotGrob(p.legend.1)$grobs[[which(sapply(ggplotGrob(p.legend.1)$grobs, function(x) x$name) == 'guide-box')]]

res_corr_sum <- res_corr %>%
  group_by(theta_lambda, delta) %>%
  summarise(median = median(cor),
            lower = quantile(cor, p = 0.1),
            upper = quantile(cor, p = 0.9)) %>%
  ungroup()

p.corr.2 <- ggplot(res_corr_sum %>%
                     mutate(delta = factor(7 / delta)) %>%
                     mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                     mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                            strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                            strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_corr_sum %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  # geom_line(aes(x = strength_proxy, y = median, group = delta, col = delta)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'A', x = 'True Interaction Strength', y = expression(paste("Pearson's  ", rho)))

p.legend.2 <- ggplot(res_corr_sum %>%
                       mutate(delta = factor(1 / delta)) %>%
                       mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                       mutate(delta = paste0(delta, ' week'),
                              delta = if_else(str_detect(delta, '1 '), delta, paste0(delta, 's'))) %>%
                       mutate(delta = factor(delta, levels = c('1 week', '4 weeks', '13 weeks')))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_brewer(palette = 'Set1') +
  labs(col = 'True Interaction Duration')
p.legend.2 <- ggplotGrob(p.legend.2)$grobs[[which(sapply(ggplotGrob(p.legend.2)$grobs, function(x) x$name) == 'guide-box')]]

rm(res_corr, res_corr_sum, acc_corr)

# ------------------------------------------------------------------------------

# Process accuracy of results (GAMs)

# Calculate sensitivity/specificity:
df_acc <- bind_rows(df_acc, data.frame(method = 'GAMs',
                                       sens_pos = (res_gam %>% filter(int_true == 'pos' & int_est == 'pos') %>% nrow()) / (res_gam %>% filter(int_true == 'pos') %>% nrow()),
                                       sens_neg = (res_gam %>% filter(int_true == 'neg' & int_est == 'neg') %>% nrow()) / (res_gam %>% filter(int_true == 'neg') %>% nrow()),
                                       spec = (res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_gam %>% filter(int_true == 'none') %>% nrow()),
                                       mcc = mcc(res_gam %>% filter(int_true != 'none' & int_est != 'none') %>% nrow(),
                                                 res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(),
                                                 res_gam %>% filter(int_true == 'none' & int_est != 'none') %>% nrow(),
                                                 res_gam %>% filter(int_true != 'none' & int_est == 'none') %>% nrow())))

# Calculate sensitivity/specificity (by true params):
acc_gam <- calculate_accuracy_matrix(res_gam)

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_gam <- calculate_assoc_true_strength(res_gam %>%
                                             mutate(sig = if_else(int_est == 'none', 'no', 'yes')),
                                           method = 'gam', met = 'cor_median')

# Plot:
p.gam.1 <- ggplot(data = acc_gam %>%
                    mutate(strength_proxy = rank(strength, ties.method = 'min')),
                  aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'B', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')

res_gam_sum <- res_gam %>%
  group_by(theta_lambda, delta) %>%
  summarise(median = median(cor_median),
            lower = quantile(cor_median, p = 0.1),
            upper = quantile(cor_median, p = 0.9)) %>%
  ungroup()

p.gam.2 <- ggplot(res_gam_sum %>%
                    mutate(delta = factor(7 / delta)) %>%
                    mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                    mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                           strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                           strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_gam_sum %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'B', x = 'True Interaction Strength', y = 'Residual Correlation')

rm(res_gam, res_gam_sum, acc_gam)

# ------------------------------------------------------------------------------

# Process accuracy of results (Granger causality)

# Check whether any datasets are not stationary:
res_granger %>% filter(adf_p >= 0.05) %>% nrow() %>% print()
res_granger %>% filter(kpss_p < 0.05) %>% nrow() %>% print()

# Calculate sensitivity/specificity:
for (i in 1:length(res_granger_LIST)) {
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('GC ', names(res_granger_LIST)[i]),
                                         sens_pos = (res_granger_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow()),
                                         sens_neg = (res_granger_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow()),
                                         spec = (res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_granger_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()),
                                         mcc = mcc(res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow(),
                                                   res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(),
                                                   res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow(),
                                                   res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow())))
  
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
p.granger.1.1 <- ggplot(data = acc_granger_LIST[[1]] %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'C', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.granger.1.2 <- ggplot(data = acc_granger_LIST[[2]] %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.granger.1.3 <- ggplot(data = acc_granger_LIST[[3]] %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'C', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.granger.1.4 <- ggplot(data = acc_granger_LIST[[4]] %>%
                          mutate(strength_proxy = rank(strength, ties.method = 'min')),
                        aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')

res_granger_sum <- lapply(res_granger_LIST, function(ix) {
  ix %>% group_by(theta_lambda, delta) %>%
    summarise(median = median(logRSS),
              lower = quantile(logRSS, p = 0.1),
              upper = quantile(logRSS, p = 0.9)) %>%
    ungroup()
})

p.granger.2.1 <- ggplot(res_granger_sum[[1]] %>%
                          mutate(delta = factor(7 / delta)) %>%
                          mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                          mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                                 strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                                 strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_granger_sum[[1]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'C', x = 'True Interaction Strength', y = expression(G[x %->% y]))
p.granger.2.2 <- ggplot(res_granger_sum[[2]] %>%
                          mutate(delta = factor(7 / delta)) %>%
                          mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                          mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                                 strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                                 strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_granger_sum[[2]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'True Interaction Strength', y = expression(G[y %->% x]))

rm(res_granger, res_granger_LIST, res_granger_sum, acc_granger_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (transfer entropy)

# Calculate sensitivity/specificity:
for (i in 1:length(res_te_LIST)) {
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('TE ', names(res_te_LIST)[i]),
                                         sens_pos = (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow()),
                                         sens_neg = (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow()),
                                         spec = (res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()),
                                         mcc = mcc(res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow(),
                                                   res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(),
                                                   res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow(),
                                                   res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow())))
  
}
for (i in 1:length(res_te_LIST)) {
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('TE ', names(res_te_LIST_confound)[i]),
                                         sens_pos = (res_te_LIST_confound[[i]] %>% filter(int_true_dir == 'pos' & int_est_confound2 == 'interaction') %>% nrow()) / (res_te_LIST_confound[[i]] %>% filter(int_true_dir == 'pos') %>% nrow()),
                                         sens_neg = (res_te_LIST_confound[[i]] %>% filter(int_true_dir == 'neg' & int_est_confound2 == 'interaction') %>% nrow()) / (res_te_LIST_confound[[i]] %>% filter(int_true_dir == 'neg') %>% nrow()),
                                         spec = (res_te_LIST_confound[[i]] %>% filter(int_true == 'none' & int_est_confound2 == 'none') %>% nrow()) / (res_te_LIST_confound[[i]] %>% filter(int_true == 'none') %>% nrow()),
                                         mcc = mcc(res_te_LIST_confound[[i]] %>% filter(int_true != 'none' & int_est_confound2 != 'none') %>% nrow(),
                                                   res_te_LIST_confound[[i]] %>% filter(int_true == 'none' & int_est_confound2 == 'none') %>% nrow(),
                                                   res_te_LIST_confound[[i]] %>% filter(int_true == 'none' & int_est_confound2 != 'none') %>% nrow(),
                                                   res_te_LIST_confound[[i]] %>% filter(int_true != 'none' & int_est_confound2 == 'none') %>% nrow())))
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_te_LIST <- vector('list', length = 4)
names(acc_te_LIST) <- c(names(res_te_LIST), names(res_te_LIST_confound))

for (i in 1:length(res_te_LIST)) {
  acc_te_LIST[[i]] <- calculate_accuracy_matrix(res_te_LIST[[i]])
}
acc_te_LIST[[3]] <- calculate_accuracy_matrix(res_te_LIST_confound[[1]] %>% rename('int_est' = 'int_est_confound2'))
acc_te_LIST[[4]] <- calculate_accuracy_matrix(res_te_LIST_confound[[2]] %>% rename('int_est' = 'int_est_confound2'))

# Are higher values of te associated with higher true interaction strength?:
assoc_te_LIST <- vector('list', length = 4)
names(assoc_te_LIST) <- c(names(res_te_LIST), names(res_te_LIST_confound))

for (i in 1:length(res_te_LIST)) {
  assoc_te_LIST[[i]] <- calculate_assoc_true_strength(res_te_LIST[[i]] %>% mutate(sig = if_else(int_est == 'interaction', 'yes', 'no')),
                                                      method = 'te', met = 'te')
}
assoc_te_LIST[[3]] <- calculate_assoc_true_strength(res_te_LIST_confound[[1]] %>% mutate(sig = if_else(int_est_confound2 == 'interaction', 'yes', 'no')),
                                                    method = 'te', met = 'te_confound2')
assoc_te_LIST[[4]] <- calculate_assoc_true_strength(res_te_LIST_confound[[2]] %>% mutate(sig = if_else(int_est_confound2 == 'interaction', 'yes', 'no')),
                                                    method = 'te', met = 'te_confound2')
rm(i)

# Plot:
p.te.1.1 <- ggplot(data = acc_te_LIST[[1]] %>%
                     mutate(strength_proxy = rank(strength, ties.method = 'min')),
                   aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'D', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.te.1.2 <- ggplot(data = acc_te_LIST[[2]] %>%
                     mutate(strength_proxy = rank(strength, ties.method = 'min')),
                   aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.te.1.3 <- ggplot(data = acc_te_LIST[[3]] %>%
                     mutate(strength_proxy = rank(strength, ties.method = 'min')),
                   aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'D', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.te.1.4 <- ggplot(data = acc_te_LIST[[4]] %>%
                     mutate(strength_proxy = rank(strength, ties.method = 'min')),
                   aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')

res_te_sum <- lapply(res_te_LIST_confound, function(ix) {
  ix %>%
    group_by(theta_lambda, delta) %>%
    summarise(median = median(te_confound2),
              lower = quantile(te_confound2, p = 0.1),
              upper = quantile(te_confound2, p = 0.9)) %>%
    ungroup()
})

p.te.2.1 <- ggplot(res_te_sum[[1]] %>%
                     mutate(delta = factor(7 / delta)) %>%
                     mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                     mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                            strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                            strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_te_sum[[1]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'D', x = 'True Interaction Strength', y = expression(T[x %->% y]))
p.te.2.2 <- ggplot(res_te_sum[[2]] %>%
                     mutate(delta = factor(7 / delta)) %>%
                     mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                     mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                            strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                            strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_te_sum[[2]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'True Interaction Strength', y = expression(T[y %->% x]))

rm(res_te_LIST, res_te_LIST_confound, res_te_sum, acc_te_LIST, best_v1xv2, best_v2xv1, best_v1xv2_confound, best_v2xv1_confound)

# ------------------------------------------------------------------------------

# Process accuracy of results (CCM)

# Calculate sensitivity/specificity:
for (i in 1:length(res_ccm_LIST)) {
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('CCM ', names(res_ccm_LIST)[i]),
                                         sens_pos = (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'pos' & int_est == 'interaction') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'pos') %>% nrow()),
                                         sens_neg = (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'neg' & int_est == 'interaction') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true_dir == 'neg') %>% nrow()),
                                         spec = (res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow()) / (res_ccm_LIST[[i]] %>% filter(int_true == 'none') %>% nrow()),
                                         mcc = mcc(res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est != 'none') %>% nrow(),
                                                   res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(),
                                                   res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est != 'none') %>% nrow(),
                                                   res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'none') %>% nrow())))
  
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
p.ccm.1.1 <- ggplot(data = acc_ccm_LIST[[1]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'E', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.ccm.1.2 <- ggplot(data = acc_ccm_LIST[[2]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'F', x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.ccm.1.3 <- ggplot(data = acc_ccm_LIST[[3]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.995),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.ccm.1.4 <- ggplot(data = acc_ccm_LIST[[4]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.ccm.1.5 <- ggplot(data = acc_ccm_LIST[[5]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')
p.ccm.1.6 <- ggplot(data = acc_ccm_LIST[[6]] %>%
                      mutate(strength_proxy = rank(strength, ties.method = 'min')),
                    aes(x = strength_proxy, y = perc_correct, shape = duration, color = duration, lty = duration)) +
  geom_line() +
  geom_point(size = 3) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 10, 13, 16), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'Strength', y = '% Correct   ', shape = 'Duration', color = 'Duration', lty = 'Duration')

res_ccm_sum <- lapply(res_ccm_LIST, function(ix) {
  ix %>%
    group_by(theta_lambda, delta) %>%
    summarise(median = median(rho_max),
              lower = quantile(rho_max, p = 0.1),
              upper = quantile(rho_max, p = 0.9)) %>%
    ungroup()
})

p.ccm.2.1 <- ggplot(res_ccm_sum[[1]] %>%
                      mutate(delta = factor(7 / delta)) %>%
                      mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                      mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                             strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                             strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_ccm_sum[[1]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(tag = 'E', x = 'True Interaction Strength', y = 'Cross-Map Skill')
p.ccm.2.2 <- ggplot(res_ccm_sum[[4]] %>%
                      mutate(delta = factor(7 / delta)) %>%
                      mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                      mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                             strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                             strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.75, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_ccm_sum[[4]] %>% filter(theta_lambda == 1),
                  size = 0.75, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 16),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 10, limits = c(-0.4, 1.0)) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = 'True Interaction Strength', y = 'Cross-Map Skill')

rm(res_ccm_LIST, res_ccm_sum, res_ccm, acc_ccm_LIST)

# ------------------------------------------------------------------------------

# Compile/plot results of all methods

# Plot overall accuracy:
df_acc <- df_acc %>%
  mutate(direction = if_else(str_detect(method, 'v1 -> v2'), 'v1 -> v2', 'v2 -> v1')) %>%
  mutate(method = str_remove(method, 'v1 -> v2 ')) %>% mutate(method = str_remove(method, 'v2 -> v1 ')) %>%
  mutate(confounding = if_else(str_detect(method, 'Seasonality') | str_detect(method, 'lag 13'), ' (w/ Seas)', '')) %>%
  mutate(method = if_else(!(method %in% c('Corr. Coef.', 'GAMs') | str_detect(method, 'CCM')), str_sub(method, 1, 3), method),
         method = if_else(str_detect(method, 'CCM') | str_detect(method, 'Corr. Coef.'), method, str_remove(method, ' ')),
         method = paste0(method, confounding)) %>%
  select(-confounding)

df_acc <- df_acc %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC', 'TE (w/ Seas)', 'CCM (Method 1)', 'CCM (Method 2)')))

df_acc <- df_acc %>% bind_rows(df_acc[1:2, ] %>%
                                 mutate(direction = 'v1 -> v2'))

df_acc %>%
  filter(!is.na(method)) %>%
  select(method, direction, mcc, sens_pos:spec) %>%
  print()

p.comb.1 <- ggplot(df_acc %>% filter(!is.na(method)), aes(x = direction, y = mcc, shape = method, col = method)) +
  geom_point(size = 3) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'right') +
  scale_shape_manual(values = c(18, 17, 15, 3, 8, 8)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#1f78b4', '#6a3d9a', '#33a02c', '#b2df8a')) +
  scale_y_continuous(limits = c(-0.10, 0.55), n.breaks = 15) +
  labs(x = 'Direction', y = 'Matthews correlation coefficient', shape = 'Method', col = 'Method')
plot(p.comb.1)
# ggsave(filename = 'results/plots/overall_accuracy_by_method.svg', p.comb.1, height = 6, width = 11)

# Plot accuracy by method, strength, and duration:
p.comb.2 <- arrangeGrob(p.corr.1, p.gam.1, p.granger.1.1, p.granger.1.2, p.te.1.3, p.te.1.4, p.ccm.1.1, p.ccm.1.4, p.ccm.1.2, p.ccm.1.5, p.legend.1,
                        layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 10), c(11, 11)),
                        heights = c(1, 1, 1, 1, 1, 0.25))
plot(p.comb.2)
# ggsave(filename = 'results/plots/accuracy_by_method_and_true_int_params.svg', p.comb.2, width = 12, height = 9)
rm(p.corr.1, p.gam.1, p.granger.1.1, p.granger.1.2, p.granger.1.3, p.granger.1.4, p.te.1.1, p.te.1.2,
   p.te.1.3, p.te.1.4, p.ccm.1.1, p.ccm.1.4, p.ccm.1.2, p.ccm.1.5, p.ccm.1.3, p.ccm.1.6, p.legend.1)

# Table of correlation coefficients between inferred and true interaction strength magnitude:
df_assoc <- assoc_corr %>%
  mutate(method = 'Corr. Coef.', .before = delta) %>%
  filter(true_int == 'all') %>%
  bind_rows(assoc_gam %>%
              mutate(method = 'GAMs', .before = delta) %>%
              filter(true_int == 'all')) %>%
  bind_rows(assoc_granger_LIST[1:2] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'GC', .before = delta) %>%
              filter(true_int == 'all') %>%
              mutate(direction = str_sub(direction, 1, 8))) %>%
  bind_rows(assoc_te_LIST[3:4] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'TE (w/ Seas)', .before = delta) %>%
              filter(true_int == 'all') %>%
              mutate(direction = str_sub(direction, 1, 8))) %>%
  bind_rows(assoc_ccm_LIST_max[c(1, 4)] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'CCM', .before = delta) %>%
              filter(true_int == 'all') %>%
              mutate(direction = str_sub(direction, 1, 8)))
print(df_assoc %>% filter(delta == 1))
print(df_assoc %>% filter(delta == 7 / 28))
print(df_assoc %>% filter(delta == 7 / 91))
# rm(assoc_corr, assoc_gam, assoc_granger_LIST, assoc_te_LIST, assoc_ccm_LIST_max)

# Plot metric values vs. true interaction strength by duration:
p.comb.3 <- arrangeGrob(p.corr.2, p.gam.2, p.granger.2.1, p.granger.2.2, p.te.2.1, p.te.2.2, p.ccm.2.1, p.ccm.2.2, p.legend.2,
                        layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9)),
                        heights = c(1, 1, 1, 1, 0.25))
plot(p.comb.3)
# ggsave(filename = 'results/plots/true_strength_vs_point_estimates.svg', p.comb.3, width = 12, height = 9)
rm(p.corr.2, p.gam.2, p.granger.2.1, p.granger.2.2, p.te.2.1, p.te.2.2, p.ccm.2.1, p.ccm.2.2, p.legend.2)

# Plot Spearman's rho between metrics and true interaction strength for all methods:
# Note: Values are calculated such that we expect POSITIVE rho for all analyses
df_assoc <- df_assoc %>%
  mutate(direction = if_else(is.na(direction), 'v1 -> v2', direction)) %>%
  bind_rows(df_assoc %>% filter(method == 'Corr. Coef.') %>% mutate(direction = 'v2 -> v1')) %>%
  bind_rows(df_assoc %>% filter(method == 'GAMs') %>% mutate(direction = 'v2 -> v1')) %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC', 'TE (w/ Seas)', 'CCM'))) %>%
  mutate(duration = 1 / delta) %>%
  mutate(d_proxy = case_match(duration, 1 ~ 1, 4 ~ 2, 13 ~ 3),
         method_num = as.numeric(method),
         x_use = d_proxy + 0.15 * (method_num - 3))

p.comb.4 <- ggplot(df_assoc, aes(x = x_use, y = rho, group = method, shape = method, col = method)) +
  geom_rect(aes(xmin = -Inf, xmax = 1.5, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_rect(aes(xmin = 1.5, xmax = 2.5, ymin = -Inf, ymax = Inf), fill = 'gray95', col = 'gray95') +
  geom_rect(aes(xmin = 2.5, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'white', col = 'white') +
  geom_point(size = 2.5) +
  geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), fill = 'white', col = NA, alpha = 0.075) +
  geom_point(data = df_assoc %>% filter(rho > 0 & p_value < 0.05), size = 2.5) +
  geom_hline(yintercept = 0) +
  facet_wrap(~ direction) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'right') +
  scale_x_continuous(limits = c(0.6, 3.4), breaks = 1:3, labels = c('1', '4', '13')) +
  scale_y_continuous(limits = c(-1, 1)) +
  scale_shape_manual(values = c(18, 17, 15, 3, 8)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#1f78b4', '#6a3d9a', '#33a02c')) +
  labs(x = 'Interaction Duration', y = "Spearman's Rho", shape = 'Method', col = 'Method')
print(p.comb.4)
# ggsave(filename = 'results/plots/correlation_true_strength_vs_point_estimates.svg', height = 5.5, width = 11)
rm(assoc_corr, assoc_gam, assoc_granger_LIST, assoc_te_LIST, assoc_ccm_LIST_max)

# Close pdf:
dev.off()

# Clean up:
rm(list = ls())
