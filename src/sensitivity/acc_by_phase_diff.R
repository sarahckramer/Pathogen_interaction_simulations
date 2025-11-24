# ------------------------------------------------------------------------------
# Check whether method accuracy varies by median peak timing difference
# ------------------------------------------------------------------------------

# Setup

# Load packages
library(tidyverse)
library(testthat)
library(mgcv)
library(gridExtra)
library(grid)

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

# Calculate outbreak metrics by interaction strength/duration

# Get season information:
year_breaks <- dat %>% filter(str_detect(date, '-07-0[1-7]')) %>% pull(date) %>% unique()
year_breaks <- c(year_breaks, '2022-07-03')

dat <- dat %>%
  mutate(season = cut(date, breaks = year_breaks, labels = 1:10, include.lowest = TRUE))

# Calculate seasonal attack rates, peak timings, and phase differences for each simulation:
ph_diff <- dat %>%
  group_by(run, theta_lambda, delta, .id, season) %>%
  summarise(pt1 = which.max(V1_obs), pt2 = which.max(V2_obs),
            pt_diff = pt1 - pt2) %>%
  ungroup() %>%
  filter(theta_lambda == 1) %>%
  group_by(.id) %>%
  summarise(pt_diff = median(pt_diff)) %>%
  ungroup()

# ------------------------------------------------------------------------------

# Get accuracy by phase difference for all methods
res_corr <- res_corr %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff %>% mutate(.id = as.numeric(.id)), by = '.id') %>%
  mutate(.id = factor(.id))
res_gam <- res_gam %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff %>% mutate(.id = as.numeric(.id)), by = '.id') %>%
  mutate(.id = factor(.id))
res_granger <- res_granger_LIST[[3]] %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff, by = '.id') %>%
  mutate(.id = as.numeric(.id)) %>%
  mutate(.id = factor(.id))
res_te <- res_te_LIST[[3]] %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff %>% mutate(.id = as.numeric(.id)), by = '.id') %>%
  mutate(.id = factor(.id))
res_ccm1 <- res_ccm_LIST[[1]] %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff, by = '.id') %>%
  mutate(.id = as.numeric(.id)) %>%
  mutate(.id = factor(.id))
res_ccm2 <- res_ccm_LIST[[2]] %>%
  mutate(correct = if_else(int_true == int_est, 1, 0)) %>%
  select(.id, correct) %>%
  left_join(ph_diff, by = '.id') %>%
  mutate(.id = as.numeric(.id)) %>%
  mutate(.id = factor(.id))
rm(res_granger_LIST, res_te_LIST, res_ccm_LIST, res_ccm_surr, res_ccm)

# ------------------------------------------------------------------------------

# Fit GAMs and check how phase difference affects chance of correct result

# Fit GAMs for all methods:
m_corr <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_corr, family = 'binomial', method = 'ML')
m_gam <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_gam, family = 'binomial', method = 'ML')
m_gc <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_granger, family = 'binomial', method = 'ML')
m_te <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_te, family = 'binomial', method = 'ML')
m_ccm1 <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_ccm1, family = 'binomial', method = 'ML')
m_ccm2 <- gam(correct ~ s(pt_diff, bs = 'cr') + s(.id, bs = 're', k = 100), data = res_ccm2, family = 'binomial', method = 'ML')

# Get predictions from all models:
pred_data <- with(res_corr,
                  expand_grid(pt_diff = seq(-4, 15, by = 0.25))) %>%
  mutate(.id = '1')

pred_data <- pred_data %>%
  bind_cols(as.data.frame(predict(m_corr, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_Corr'), fit:se.fit) %>%
  bind_cols(as.data.frame(predict(m_gam, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_GAM'), fit:se.fit) %>%
  bind_cols(as.data.frame(predict(m_gc, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_GC'), fit:se.fit) %>%
  bind_cols(as.data.frame(predict(m_te, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_TE'), fit:se.fit) %>%
  bind_cols(as.data.frame(predict(m_ccm1, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_CCM1'), fit:se.fit) %>%
  bind_cols(as.data.frame(predict(m_ccm2, pred_data, type = 'link', se.fit = TRUE, exclude = 's(.id)'))) %>%
  rename_with(~ paste0(.x, '_CCM2'), fit:se.fit) %>%
  pivot_longer(fit_Corr:se.fit_CCM2) %>%
  separate_wider_delim(name, delim = '_', names = c('meas', 'method')) %>%
  pivot_wider(names_from = meas, values_from = value)

# Convert to probability correct:
ilink <- family(m_corr)$linkinv

pred_data <- pred_data %>%
  mutate(fitted = ilink(fit),
         lower = ilink(fit - (2 * se.fit)),
         upper = ilink(fit + (2 * se.fit))) %>%
  select(-c(fit:se.fit))

# ------------------------------------------------------------------------------

# Plot impact of peak timing difference on method accuracy

p1 <- ggplot(data = pred_data %>% filter(method == 'Corr'),
             aes(x = pt_diff, y = fitted)) +
  # geom_rug(data = res_corr %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
  #          aes(pt_diff), sides = 'b', inherit.aes = FALSE, alpha = 0.2, linewidth = 1) +
  geom_density(data = res_corr %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = '#ff7f00', alpha = 0.1) +
  geom_line(col = '#ff7f00', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'Corr. Coef.')
p2 <- ggplot(data = pred_data %>% filter(method == 'GAM'),
             aes(x = pt_diff, y = fitted)) +
  geom_density(data = res_gam %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = '#e31a1c', alpha = 0.1) +
  geom_line(col = '#e31a1c', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'GAMs')
p3 <- ggplot(data = pred_data %>% filter(method == 'GC'),
             aes(x = pt_diff, y = fitted - 0.4)) +
  geom_density(data = res_granger %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower - 0.4, ymax = upper - 0.4), fill = '#1f78b4', alpha = 0.1) +
  geom_line(col = '#1f78b4', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), labels = c('0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'), limits = c(0, 0.6), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'GC')
p4 <- ggplot(data = pred_data %>% filter(method == 'TE'),
             aes(x = pt_diff, y = fitted - 0.4)) +
  geom_density(data = res_te %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower - 0.4, ymax = upper - 0.4), fill = '#6a3d9a', alpha = 0.1) +
  geom_line(col = '#6a3d9a', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), labels = c('0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'), limits = c(0, 0.6), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'TE')
p5 <- ggplot(data = pred_data %>% filter(method == 'CCM2'),
             aes(x = pt_diff, y = fitted - 0.4)) +
  geom_density(data = res_ccm1 %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower - 0.4, ymax = upper - 0.4), fill = '#33a02c', alpha = 0.1) +
  geom_line(col = '#33a02c', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), labels = c('0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'), limits = c(0, 0.6), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'CCM (Method 1)')
p6 <- ggplot(data = pred_data %>% filter(method == 'CCM1'),
             aes(x = pt_diff, y = fitted - 0.4)) +
  geom_density(data = res_ccm2 %>% group_by(.id) %>% summarise(pt_diff = unique(pt_diff)),
               aes(x = pt_diff), inherit.aes = FALSE, adjust = 0.5, fill = 'gray95', col = 'gray80') +
  geom_vline(xintercept = 0, lty = 2) +
  geom_ribbon(aes(ymin = lower - 0.4, ymax = upper - 0.4), fill = '#b2df8a', alpha = 0.2) +
  geom_line(col = '#a5da76', linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-4, 15, by = 1), limits = c(-4, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.6, by = 0.1), labels = c('0.4', '0.5', '0.6', '0.7', '0.8', '0.9', '1.0'), limits = c(0, 0.6), expand = c(0, 0)) +
  theme(title = element_text(size = 11),
        axis.text = element_text(size = 12),
        panel.grid.minor = element_blank(),
        plot.margin = unit(c(5.5, 8, 5.5, 5.5), 'pt')) +
  labs(x = NULL, y = NULL, title = 'CCM (Method 2)')

y_lab <- textGrob('Probability Correct', rot = 90, gp = gpar(fontsize = 14))
x_lab <- textGrob('Peak Timing Difference', gp = gpar(fontsize = 14, hjust = 1))

p <- arrangeGrob(p3, p4, p5, p6, ncol = 2, left = y_lab, bottom = x_lab)
plot(p)
ggsave(filename = 'results/plots/figures/FigureS8.svg', p, width = 10, height = 6)
