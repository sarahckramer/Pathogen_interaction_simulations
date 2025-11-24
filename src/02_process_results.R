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
library(grid)
library(testthat)
library(viridis)
library(lme4)

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

# # Store results (for comparison with sensitivity analyses)
# 
# res_LIST <- list(res_corr, res_gam, bind_rows(res_granger_LIST[3:4]),
#                  bind_rows(res_te_LIST[3:4]), res_ccm)
# write_rds(res_LIST, file = 'results/res_compiled.rds')
# rm(res_LIST)

# ------------------------------------------------------------------------------

# Visualize some datasets

# Choose simulations to plot:
ids_to_plot <- c(54, 70, 96, 89, 27)

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

for (i in ids_to_plot) {
  
  p.data <- ggplot(data = dat_plot %>%
                     filter(.id == i) %>%
                     mutate(theta_lambda = as.character(theta_lambda),
                            delta = as.character(7 / delta)) %>%
                     filter((theta_lambda %in% c(0, 4)) |
                              (theta_lambda == 1 & delta == 7)) %>%
                     mutate(theta_lambda = paste0('Interaction Strength: ', theta_lambda)),
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
    # scale_color_manual(values = c('#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3')) +
    scale_color_manual(values = c('#fc9272', '#de2d26', '#fee0d2', '#6baed6', '#2171b5', '#bdd7e7')) +
    labs(x = '', y = 'Incidence (per 1000)', col = '', title = i)
  print(p.data)
  
}

dat_plot_TEMP <- dat_plot %>%
  mutate(theta_lambda = as.character(theta_lambda),
         delta = as.character(7 / delta)) %>%
  filter((theta_lambda %in% c(0, 4)) |
           (theta_lambda == 1 & delta == 7))

p1 <- ggplot(data = dat_plot_TEMP %>% filter(.id == 54),
             aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  # scale_color_manual(values = c('#e31a1c', '#fb9a99', '#1f78b4', '#a6cee3')) +
  # scale_color_manual(values = c('#fc9272', '#de2d26', '#fee0d2', '#6baed6', '#2171b5', '#bdd7e7')) +
  scale_color_manual(values = c('red', '#99000d', '#fc9272', '#6baed6', '#08519c', '#bdd7e7')) +
  labs(x = 'Time', y = 'Incidence (per 1000)', col = '', tag = 'A')

s_p1 <- ggplot(data = dat_plot_TEMP %>% filter(.id == 70),
               aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 32)) +
  scale_color_manual(values = c('red', '#99000d', '#fc9272', '#6baed6', '#08519c', '#bdd7e7')) +
  labs(x = 'Time', y = 'Incidence (per 1000)', col = '', tag = 'A')
s_p2 <- ggplot(data = dat_plot_TEMP %>% filter(.id == 96),
               aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 32)) +
  scale_color_manual(values = c('red', '#99000d', '#fc9272', '#6baed6', '#08519c', '#bdd7e7')) +
  labs(x = 'Time', y = 'Incidence (per 1000)', col = '', tag = 'B')
s_p3 <- ggplot(data = dat_plot_TEMP %>% filter(.id == 89),
               aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 32)) +
  scale_color_manual(values = c('red', '#99000d', '#fc9272', '#6baed6', '#08519c', '#bdd7e7')) +
  labs(x = 'Time', y = 'Incidence (per 1000)', col = '', tag = 'C')
s_p4 <- ggplot(data = dat_plot_TEMP %>% filter(.id == 27),
               aes(x = date, y = obs, group = paste(virus, delta), color = paste(virus, delta))) +
  geom_line(linewidth = 0.75) +
  geom_vline(xintercept = year_breaks, lty = 2, col = 'gray60') +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(legend.position = 'none',
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  scale_y_continuous(limits = c(0, 32)) +
  scale_color_manual(values = c('red', '#99000d', '#fc9272', '#6baed6', '#08519c', '#bdd7e7')) +
  labs(x = 'Time', y = 'Incidence (per 1000)', col = '', tag = 'D')
p.rep.data <- arrangeGrob(s_p1, s_p2, s_p3, s_p4, nrow = 2)
plot(p.rep.data)
# ggsave(filename = 'results/plots/figures/FigureS1.svg', p.rep.data, width = 18, height = 9.75)

rm(dat_plot, ids_to_plot)

# ------------------------------------------------------------------------------

# Calculate outbreak metrics by interaction strength/duration

# Get season information:
dat <- dat %>%
  mutate(season = cut(date, breaks = year_breaks, labels = 1:10, include.lowest = TRUE))

# Set population size:
N <- 5000000

# Calculate seasonal attack rates, peak timings, and phase differences for each simulation:
met <- dat %>%
  group_by(run, theta_lambda, delta, .id, season) %>%
  summarise(ar1 = sum(V1_obs) / N, ar2 = sum(V2_obs) / N,
            pt1 = which.max(V1_obs), pt2 = which.max(V2_obs),
            pt_diff = abs(pt1 - pt2)) %>%
  ungroup()

# Calculate difference in metrics from scenario with no interaction:
met <- met %>%
  group_by(.id, season, delta) %>%
  mutate(change_ar1 = ar1 / ar1[theta_lambda == 1],
         change_ar2 = ar2 / ar2[theta_lambda == 1],
         change_pt1 = pt1 - pt1[theta_lambda == 1],
         change_pt2 = pt2 - pt2[theta_lambda == 1],
         change_ptdiff = pt_diff - pt_diff[theta_lambda == 1]) %>%
  ungroup()

# Get average and st.dev. across all seasons:
met <- met %>%
  group_by(run, theta_lambda, delta, .id) %>%
  summarise(across(ar1:change_ptdiff, mean, .names = '{.col}_mean'),
            ar1_sd = sd(ar1), ar2_sd = sd(ar2),
            pt1_sd = sd(pt1), pt2_sd = sd(pt2),
            pt_diff_sd = sd(pt_diff)) %>%
  ungroup() %>%
  select(run:pt_diff_mean, ar1_sd:pt_diff_sd, change_ar1_mean:change_ptdiff_mean)

# Get change in st.dev. from scenario with no interaction:
met <- met %>%
  group_by(.id, delta) %>%
  mutate(change_ar1_sd = ar1_sd / ar1_sd[theta_lambda == 1],
         change_ar2_sd = ar2_sd / ar2_sd[theta_lambda == 1],
         change_pt1_sd = pt1_sd - pt1_sd[theta_lambda == 1],
         change_pt2_sd = pt2_sd - pt2_sd[theta_lambda == 1],
         change_ptdiff_sd = pt_diff_sd - pt_diff_sd[theta_lambda == 1]) %>%
  ungroup()

# Plot:
met <- met %>%
  mutate(delta = factor(1 / delta, levels = c('1', '4', '13')))

# met %>%
#   select(theta_lambda:delta, contains('change')) %>%
#   filter(theta_lambda != 1) %>%
#   pivot_longer(contains('change')) %>%
#   ggplot(aes(x = theta_lambda, y = value, group = paste(theta_lambda, delta), fill = delta)) +
#   geom_boxplot() + theme_classic() + facet_wrap(~ name, scales = 'free')

p2a <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
              # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
              aes(x = change_ar2_mean, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 1.0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.03, 0.975),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_fill_viridis(discrete = TRUE) + scale_color_viridis(discrete = TRUE) +
  scale_x_continuous(n.breaks = 5) +
  labs(x = 'Relative AR', y = 'Density', tag = 'B')
p2b <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
              # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
              aes(x = change_ar2_sd, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 1.0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Relative sd(AR)', y = '')
p2c <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
              # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
              aes(x = change_pt2_mean, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Change in PT', y = '')
p2d <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
              # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
              aes(x = change_pt2_sd, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Change in sd(PT)', y = '')

p.legend <- ggplot(met %>%
                     filter(theta_lambda %in% c(0, 4)) %>%
                     mutate(delta = if_else(delta == '1', paste0(delta, ' week'), paste0(delta, ' weeks'))) %>%
                     mutate(delta = factor(delta, levels = c('1 week', '4 weeks', '13 weeks'))),
                   aes(x = change_ar1_mean, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  theme_classic() +
  theme(legend.title = element_text(size = 14, margin = margin(b = 8)),
        legend.text = element_text(size = 14),
        legend.key.spacing.y = unit(0.1, 'cm'),
        legend.key.height = unit(1.0, 'cm'),
        legend.key.width = unit(1.0, 'cm'),
        legend.position = 'right') +
  # guides(fill = guide_legend(byrow = TRUE)) +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(fill = expression(7/delta))
p.legend <- ggplotGrob(p.legend)$grobs[[which(sapply(ggplotGrob(p.legend)$grobs, function(x) x$name) == 'guide-box')]]

p.int.impact <- arrangeGrob(p1, arrangeGrob(p2a, p2b, p2c, p2d, p.legend, nrow = 1, widths = c(1, 1, 1, 1, 0.5)), heights = c(1, 0.92), ncol = 1)
plot(p.int.impact)
# ggsave(filename = 'results/plots/figures/Figure2.svg', p.int.impact, width = 12.5, height = 9)

# Also plot impact on pathogen A (for supplement):
p2a_supp <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
                   # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
                   aes(x = change_ar1_mean, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 1.0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        # plot.tag = element_text(size = 24),
        # plot.tag.position = c(0.03, 0.975),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 5) +
  labs(x = 'Relative AR', y = 'Density')
p2b_supp <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
                   # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
                   aes(x = change_ar1_sd, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 1.0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Relative sd(AR)', y = '')
p2c_supp <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
                   # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
                   aes(x = change_pt1_mean, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Change in PT', y = '')
p2d_supp <- ggplot(met %>% filter(theta_lambda %in% c(0, 4)),
                   # met %>% filter(theta_lambda != 1 & theta_lambda != 0.25),
                   aes(x = change_pt1_sd, group = delta, fill = delta)) +
  geom_density(adjust = 0.75, linewidth = 0.5, alpha = 0.3, color = 'white') +
  geom_vline(xintercept = 0, lty = 2) +
  facet_wrap(~ theta_lambda, ncol = 1, labeller = label_bquote(theta[lambda]==.(theta_lambda))) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 11),
        strip.text = element_text(size = 12),
        legend.position = 'none') +
  scale_fill_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  scale_x_continuous(n.breaks = 10) +
  labs(x = 'Change in sd(PT)', y = '')

p.int.impact_supp <- arrangeGrob(p2a_supp, p2b_supp, p2c_supp, p2d_supp, p.legend, nrow = 1, widths = c(1, 1, 1, 1, 0.45))
plot(p.int.impact_supp)
# ggsave(filename = 'results/plots/figures/FigureS2.svg', p.int.impact_supp, width = 12.5, height = 4.3)

# ------------------------------------------------------------------------------

# Process accuracy of results (correlation coefficients)

# Calculate sensitivity/specificity:
sens_test <- binom.test(res_corr %>% filter((int_true == 'pos' & int_est == 'pos') | (int_true == 'neg' & int_est == 'neg')) %>% nrow(), res_corr %>% filter(int_true == 'pos' | int_true == 'neg') %>% nrow())
spec_test <- binom.test(res_corr %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_corr %>% filter(int_true == 'none') %>% nrow())

df_acc <- as_tibble(data.frame(method = 'Corr. Coef.',
                               sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                               spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))

# Sensitivity analysis - remove any runs where GAMs did not fit successfully:
res_CHECK <- res_corr %>%
  left_join(to_remove_GAM, by = c('run', '.id')) %>%
  filter(is.na(delete)) %>%
  select(-delete)

sens_test <- binom.test(res_CHECK %>% filter((int_true == 'pos' & int_est == 'pos') | (int_true == 'neg' & int_est == 'neg')) %>% nrow(), res_CHECK %>% filter(int_true == 'pos' | int_true == 'neg') %>% nrow())
spec_test <- binom.test(res_CHECK %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_CHECK %>% filter(int_true == 'none') %>% nrow())

df_acc_CHECK <- as_tibble(data.frame(method = 'Corr. Coef.',
                                     sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                     spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))

# Calculate sensitivity/specificity (by true params):
acc_corr <- calculate_accuracy_matrix(res_corr)

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_corr <- calculate_assoc_true_strength(res_corr, method = 'corr', met = 'cor')

# Plot:
p.corr.1a <- ggplot(data = acc_corr %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength > '1'),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'A', x = NULL, y = '', title = NULL)
p.corr.1b <- ggplot(data = acc_corr %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength < 1) %>%
                      mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.legend.1 <- ggplot(data = acc_corr %>%
                       mutate(duration = paste0(duration, ' week'),
                              duration = if_else(str_detect(duration, '1 '), duration, paste0(duration, 's'))) %>%
                       mutate(duration = factor(duration, levels = c('1 week', '4 weeks', '13 weeks'))),
                     aes(x = as.numeric(strength), y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5) +
  theme_classic() +
  theme(legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  # scale_color_viridis(discrete = TRUE) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(shape = 'Duration', color = 'Duration')
p.legend.1 <- ggplotGrob(p.legend.1)$grobs[[which(sapply(ggplotGrob(p.legend.1)$grobs, function(x) x$name) == 'guide-box')]]

res_corr_sum <- res_corr %>%
  mutate(delta = if_else(theta_lambda == 1, 1, delta)) %>%
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
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_corr_sum %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  # geom_line(aes(x = strength_proxy, y = median, group = delta, col = delta)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 6) +#, limits = c(-1.0, 1.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'A', x = 'True Interaction Strength', y = expression(r[Pearson]), title = 'Corr. Coef.')

p.legend.2 <- ggplot(res_corr_sum %>%
                       mutate(delta = factor(1 / delta)) %>%
                       mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                       mutate(delta = paste0(delta, ' week'),
                              delta = if_else(str_detect(delta, '1 '), delta, paste0(delta, 's'))) %>%
                       mutate(delta = factor(delta, levels = c('1 week', '4 weeks', '13 weeks')))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  theme_classic() +
  theme(legend.title = element_text(size = 13.5),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(col = 'True Interaction Duration')
p.legend.2 <- ggplotGrob(p.legend.2)$grobs[[which(sapply(ggplotGrob(p.legend.2)$grobs, function(x) x$name) == 'guide-box')]]

rm(res_corr, res_corr_sum, acc_corr)

# ------------------------------------------------------------------------------

# Process accuracy of results (GAMs)

# Calculate sensitivity/specificity:
sens_test <- binom.test(res_gam %>% filter((int_true == 'pos' & int_est == 'pos') | (int_true == 'neg' & int_est == 'neg')) %>% nrow(), res_gam %>% filter(int_true == 'pos' | int_true == 'neg') %>% nrow())
spec_test <- binom.test(res_gam %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_gam %>% filter(int_true == 'none') %>% nrow())

df_acc <- bind_rows(df_acc, data.frame(method = 'GAMs',
                                       sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                       spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))

# Calculate sensitivity/specificity (by true params):
acc_gam <- calculate_accuracy_matrix(res_gam)

# Are correlation coefficient magnitudes associated with true interaction strength?:
assoc_gam <- calculate_assoc_true_strength(res_gam, method = 'gam', met = 'cor_median')

# Plot:
p.gam.1a <- ggplot(data = acc_gam %>%
                     mutate(strength = as.numeric(strength)) %>%
                     filter(strength > '1'),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'B', x = NULL, y = '', title = NULL)
p.gam.1b <- ggplot(data = acc_gam %>%
                     mutate(strength = as.numeric(strength)) %>%
                     filter(strength < 1) %>%
                     mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                   aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

res_gam_sum <- res_gam %>%
  mutate(delta = if_else(theta_lambda == 1, 1, delta)) %>%
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
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_gam_sum %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 6) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'B', x = 'True Interaction Strength', y = expression(r[GAM]), title = 'GAMs')

rm(res_gam, res_gam_sum, acc_gam)

# ------------------------------------------------------------------------------

# Process accuracy of results (Granger causality)

# Calculate sensitivity/specificity:
for (i in 1:length(res_granger_LIST)) {
  
  sens_test <- binom.test(res_granger_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_granger_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_granger_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_granger_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('GC ', names(res_granger_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}
rm(i)

# Sensitivity analysis - remove any runs where GAMs did not fit successfully:
res_granger_LIST_CHECK <- res_granger_LIST[3:4]
res_granger_LIST_CHECK <- lapply(res_granger_LIST_CHECK, function(ix) {
  ix %>%
    mutate(.id = as.numeric(.id)) %>%
    left_join(to_remove_GAM, by = c('run', '.id')) %>%
    filter(is.na(delete)) %>%
    select(-delete)
})

for (i in 1:length(res_granger_LIST_CHECK)) {
  
  sens_test <- binom.test(res_granger_LIST_CHECK[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_granger_LIST_CHECK[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_granger_LIST_CHECK[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_granger_LIST_CHECK[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc_CHECK <- bind_rows(df_acc_CHECK, data.frame(method = paste0('GC ', names(res_granger_LIST_CHECK)[i]),
                                                     sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                                     spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
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
  assoc_granger_LIST[[i]] <- calculate_assoc_true_strength(res_granger_LIST[[i]], method = 'granger', met = 'logRSS')
}
rm(i)

# Plot:
p.granger.1.3a <- ggplot(data = acc_granger_LIST[[3]] %>%
                           mutate(strength = as.numeric(strength)) %>%
                           filter(strength > '1'),
                         aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'C', x = NULL, y = '', title = NULL)
p.granger.1.3b <- ggplot(data = acc_granger_LIST[[3]] %>%
                           mutate(strength = as.numeric(strength)) %>%
                           filter(strength < 1) %>%
                           mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                         aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.granger.1.4a <- ggplot(data = acc_granger_LIST[[4]] %>%
                           mutate(strength = as.numeric(strength)) %>%
                           filter(strength > '1'),
                         aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'D', x = NULL, y = '', title = NULL)
p.granger.1.4b <- ggplot(data = acc_granger_LIST[[4]] %>%
                           mutate(strength = as.numeric(strength)) %>%
                           filter(strength < 1) %>%
                           mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                         aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

res_granger_sum <- lapply(res_granger_LIST, function(ix) {
  ix %>% mutate(delta = if_else(theta_lambda == 1, 1, delta)) %>% group_by(theta_lambda, delta) %>%
    summarise(median = median(logRSS),
              lower = quantile(logRSS, p = 0.1),
              upper = quantile(logRSS, p = 0.9)) %>%
    ungroup()
})

p.granger.2.1 <- ggplot(res_granger_sum[[3]] %>%
                          mutate(delta = factor(7 / delta)) %>%
                          mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                          mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                                 strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                                 strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_granger_sum[[3]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 7) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'C', x = 'True Interaction Strength', y = expression(G[x %->% y]), title = 'Granger Causality (A \u279E B)')
p.granger.2.2 <- ggplot(res_granger_sum[[4]] %>%
                          mutate(delta = factor(7 / delta)) %>%
                          mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                          mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                                 strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                                 strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_granger_sum[[4]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_y_continuous(n.breaks = 4) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'D', x = 'True Interaction Strength', y = expression(G[y %->% x]), title = 'Granger Causality (B \u279E A)')

rm(res_granger_LIST, res_granger_sum, acc_granger_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (transfer entropy)

# Calculate sensitivity/specificity:
for (i in 1:length(res_te_LIST)) {
  
  sens_test <- binom.test(res_te_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_te_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_te_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_te_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('TE ', names(res_te_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}
rm(i)

# Sensitivity analysis - remove any runs where GAMs did not fit successfully:
res_te_LIST_CHECK <- res_te_LIST[3:4]
res_te_LIST_CHECK <- lapply(res_te_LIST_CHECK, function(ix) {
  ix %>%
    mutate(.id = as.numeric(.id)) %>%
    left_join(to_remove_GAM, by = c('run', '.id')) %>%
    filter(is.na(delete)) %>%
    select(-delete)
})

for (i in 1:length(res_te_LIST_CHECK)) {
  
  sens_test <- binom.test(res_te_LIST_CHECK[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_te_LIST_CHECK[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_te_LIST_CHECK[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_te_LIST_CHECK[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc_CHECK <- bind_rows(df_acc_CHECK, data.frame(method = paste0('TE ', names(res_te_LIST_CHECK)[i]),
                                                     sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                                     spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}
rm(i)

# Calculate accuracy by true parameter values:
acc_te_LIST <- vector('list', length = 4)
names(acc_te_LIST) <- c(names(res_te_LIST))

for (i in 1:length(res_te_LIST)) {
  acc_te_LIST[[i]] <- calculate_accuracy_matrix(res_te_LIST[[i]])
}

# Are higher values of te associated with higher true interaction strength?:
assoc_te_LIST <- vector('list', length = 4)
names(assoc_te_LIST) <- c(names(res_te_LIST))

for (i in 1:length(res_te_LIST)) {
  assoc_te_LIST[[i]] <- calculate_assoc_true_strength(res_te_LIST[[i]], method = 'te', met = 'te')
}
rm(i)

# Plot:
p.te.1.3a <- ggplot(data = acc_te_LIST[[3]] %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength > '1'),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'E', x = NULL, y = '', title = NULL)
p.te.1.3b <- ggplot(data = acc_te_LIST[[3]] %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength < 1) %>%
                      mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.te.1.4a <- ggplot(data = acc_te_LIST[[4]] %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength > '1'),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'F', x = NULL, y = '', title = NULL)
p.te.1.4b <- ggplot(data = acc_te_LIST[[4]] %>%
                      mutate(strength = as.numeric(strength)) %>%
                      filter(strength < 1) %>%
                      mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                    aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

res_te_sum <- lapply(res_te_LIST, function(ix) {
  ix %>%
    mutate(delta = if_else(theta_lambda == 1, 1, delta)) %>%
    group_by(theta_lambda, delta) %>%
    summarise(median = median(te),
              lower = quantile(te, p = 0.1),
              upper = quantile(te, p = 0.9)) %>%
    ungroup()
})

p.te.2.1 <- ggplot(res_te_sum[[3]] %>%
                     mutate(delta = factor(7 / delta)) %>%
                     mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                     mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                            strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                            strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_te_sum[[3]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'E', x = 'True Interaction Strength', y = expression(T[x %->% y]), title = 'Transfer Entropy (A \u279E B)')
p.te.2.2 <- ggplot(res_te_sum[[4]] %>%
                     mutate(delta = factor(7 / delta)) %>%
                     mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                     mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                            strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                            strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_te_sum[[4]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'F', x = 'True Interaction Strength', y = expression(T[y %->% x]), title = 'Transfer Entropy (B \u279E A)')

rm(res_te_LIST, res_te_sum, acc_te_LIST)

# ------------------------------------------------------------------------------

# Process accuracy of results (CCM)

# Calculate sensitivity/specificity:
for (i in 1:length(res_ccm_LIST)) {
  
  sens_test <- binom.test(res_ccm_LIST[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_ccm_LIST[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_ccm_LIST[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_ccm_LIST[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc <- bind_rows(df_acc, data.frame(method = paste0('CCM ', names(res_ccm_LIST)[i]),
                                         sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                         spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
}
rm(i)

# Sensitivity analysis - remove any runs where GAMs did not fit successfully:
res_ccm_LIST_CHECK <- res_ccm_LIST[c(1:2, 4:5)]
res_ccm_LIST_CHECK <- lapply(res_ccm_LIST_CHECK, function(ix) {
  ix %>%
    mutate(.id = as.numeric(.id)) %>%
    left_join(to_remove_GAM, by = c('run', '.id')) %>%
    filter(is.na(delete)) %>%
    select(-delete)
})

for (i in 1:length(res_ccm_LIST_CHECK)) {
  
  sens_test <- binom.test(res_ccm_LIST_CHECK[[i]] %>% filter(int_true != 'none' & int_est == 'interaction') %>% nrow(), res_ccm_LIST_CHECK[[i]] %>% filter(int_true != 'none') %>% nrow())
  spec_test <- binom.test(res_ccm_LIST_CHECK[[i]] %>% filter(int_true == 'none' & int_est == 'none') %>% nrow(), res_ccm_LIST_CHECK[[i]] %>% filter(int_true == 'none') %>% nrow())
  
  df_acc_CHECK <- bind_rows(df_acc_CHECK, data.frame(method = paste0('CCM ', names(res_ccm_LIST_CHECK)[i]),
                                                     sens = sens_test$estimate, lower_sens = sens_test$conf.int[1], upper_sens = sens_test$conf.int[2],
                                                     spec = spec_test$estimate, lower_spec = spec_test$conf.int[1], upper_spec = spec_test$conf.int[2]))
  
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

for (i in 1:length(assoc_ccm_LIST_max)) {
  assoc_ccm_LIST_max[[i]] <- calculate_assoc_true_strength(res_ccm_LIST[[i]], method = 'ccm', met = 'rho_max')
}
rm(i)

# Check whether there are any places were surrogate rhos are very different from observed rho
# (Could suggest that the surrogates are not representative of the null distribution accounting
# for underlying seasonality)
p.ccm.surr <- res_ccm_surr %>%
  select(run:.id, direction:rho_ci_upper) %>%
  inner_join(res_ccm %>% select(run:direction, rho_max, p_surr:int_est_1),
             by = c('run', '.id', 'direction')) %>%
  ggplot() +
  geom_violin(aes(x = .id, y = rho, group = paste(.id, direction), fill = direction), alpha = 0.1) +
  geom_point(aes(x = .id, y = rho_max, group = direction, color = direction)) +
  theme_classic() +
  facet_grid(theta_lambda ~ delta, scales = 'free')
# print(p.ccm.surr)
rm(res_ccm_surr, p.ccm.surr)

# Plot:
p.ccm.1.1a <- ggplot(data = acc_ccm_LIST[[1]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength > '1'),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'I', x = NULL, y = '', title = NULL)
p.ccm.1.1b <- ggplot(data = acc_ccm_LIST[[1]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength < 1) %>%
                       mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.ccm.1.2a <- ggplot(data = acc_ccm_LIST[[2]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength > '1'),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'G', x = NULL, y = '', title = NULL)
p.ccm.1.2b <- ggplot(data = acc_ccm_LIST[[2]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength < 1) %>%
                       mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.ccm.1.4a <- ggplot(data = acc_ccm_LIST[[4]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength > '1'),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'J', x = NULL, y = '', title = NULL)
p.ccm.1.4b <- ggplot(data = acc_ccm_LIST[[4]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength < 1) %>%
                       mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

p.ccm.1.5a <- ggplot(data = acc_ccm_LIST[[5]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength > '1'),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line() +
  geom_point(size = 3.5, alpha = 0.7) +
  theme_classic() +
  theme(axis.title = element_text(size = 13),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 1.12),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(2, 4), limits = c(1.75, 4.25)) +
  scale_y_continuous(limits = c(0, 1.01), n.breaks = 5) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_viridis(option = 'C', discrete = TRUE) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(tag = 'H', x = NULL, y = '', title = NULL)
p.ccm.1.5b <- ggplot(data = acc_ccm_LIST[[5]] %>%
                       mutate(strength = as.numeric(strength)) %>%
                       filter(strength < 1) %>%
                       mutate(strength = if_else(strength == 0, 8, 1 / as.numeric(strength))),
                     aes(x = strength, y = perc_correct, shape = duration, color = duration)) +
  geom_line(linetype = 2) +
  geom_point(size = 3.5) +
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

res_ccm_sum <- lapply(res_ccm_LIST, function(ix) {
  ix %>%
    mutate(delta = if_else(theta_lambda == 1, 1, delta)) %>%
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
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_ccm_sum[[1]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'G', x = 'True Interaction Strength', y = expression(rho), title = 'CCM (A \u279E B)')
p.ccm.2.2 <- ggplot(res_ccm_sum[[4]] %>%
                      mutate(delta = factor(7 / delta)) %>%
                      mutate(strength_proxy = rank(theta_lambda, ties.method = 'min')) %>%
                      mutate(strength_proxy = if_else(delta == 7, strength_proxy - 0.25, strength_proxy),
                             strength_proxy = if_else(delta == 91, strength_proxy + 0.25, strength_proxy),
                             strength_proxy = if_else(theta_lambda > 1, strength_proxy + 1.5, strength_proxy))) +
  geom_pointrange(size = 0.5, linewidth = 1.0, aes(x = strength_proxy, y = median, ymin = lower, ymax = upper, col = delta)) +
  geom_pointrange(data = res_ccm_sum[[4]] %>% filter(theta_lambda == 1),
                  size = 0.5, linewidth = 1.0, x = 9.75, aes(y = median, ymin = lower, ymax = upper)) +
  geom_hline(yintercept = 0, lty = 2, linewidth = 1.0, col = 'gray70') +
  theme_classic() +
  theme(axis.title = element_text(size = 13.5),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.0075, 0.975),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1, 4, 7, 9.75, 12.5, 15.5), labels = c(0, 0.25, 0.5, 1.0, 2.0, 4.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'H', x = 'True Interaction Strength', y = expression(rho), title = 'CCM (B \u279E A)')

rm(res_ccm_LIST, res_ccm_sum, res_ccm, acc_ccm_LIST)

# ------------------------------------------------------------------------------

# Compile/plot results of all methods

# Plot overall accuracy:
df_acc <- df_acc %>%
  mutate(direction = if_else(str_detect(method, 'v1 -> v2'), 'v1 -> v2', 'v2 -> v1')) %>%
  mutate(method = str_remove(method, 'v1 -> v2 ')) %>% mutate(method = str_remove(method, 'v2 -> v1 ')) %>%
  mutate(confounding = if_else(str_detect(method, 'Seasonality'), '(w/ Seas)', '')) %>%
  mutate(method = if_else(!(method %in% c('Corr. Coef.', 'GAMs') | str_detect(method, 'CCM')), str_sub(method, 1, 2), method)) %>%
  mutate(sig = str_sub(method, 13, 13)) %>%
  mutate(method = if_else(str_detect(method, 'CCM'), paste0('CCM', sig), method)) %>%
  select(-sig)
df_acc_STORE <- df_acc

df_acc <- df_acc %>%
  filter((method %in% c('GC', 'TE') & confounding == '(w/ Seas)') | !(method %in% c('GC', 'TE')),
         method != 'CCM3') %>%
  select(-confounding) %>%
  mutate(method_new = method,
         method_new = if_else(method == 'CCM1', 'CCM2', method_new),
         method_new = if_else(method == 'CCM2', 'CCM1', method_new)) %>%
  mutate(method = method_new) %>%
  select(-method_new) %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC', 'TE', 'CCM1', 'CCM2')))

df_acc <- df_acc %>% bind_rows(df_acc[1:2, ] %>%
                                 mutate(direction = 'v1 -> v2'))

df_acc %>%
  select(method, direction, sens:upper_spec) %>%
  print()

p.comb.1a <- ggplot(df_acc %>% filter(!is.na(method)) %>% filter(direction == 'v1 -> v2') %>%
                      mutate(lab_x = if_else(str_detect(method, 'TE'), sens + 0.03, sens)) %>%
                      mutate(lab_y = if_else(str_detect(method, 'TE'), spec + 0.04, spec)),
                    aes(x = sens, y = spec, shape = method, col = method)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = 'gray70') +
  geom_segment(aes(x = lower_sens, xend = upper_sens, y = spec)) +
  geom_segment(aes(y = lower_spec, yend = upper_spec, x = sens)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(limits = c(0.4, 1), n.breaks = 4) +
  scale_y_continuous(limits = c(0, 0.8), n.breaks = 6) +
  scale_shape_manual(values = c(18, 17, 15, 16, 4, 4)) +
  # scale_shape_manual(values = c(5, 2, 0, 1, 4, 4)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#1f78b4', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = 'Sensitivity', y = 'Specificity', tag = 'A')
p.comb.1b <- ggplot(df_acc %>% filter(!is.na(method)) %>% filter(direction == 'v2 -> v1') %>%
                      mutate(lab_x = if_else(str_detect(method, 'TE'), sens + 0.03, sens)) %>%
                      mutate(lab_y = if_else(str_detect(method, 'TE'), spec + 0.04, spec)),
                    aes(x = sens, y = spec, shape = method, col = method)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = 'gray70') +
  geom_segment(aes(x = lower_sens, xend = upper_sens, y = spec)) +
  geom_segment(aes(y = lower_spec, yend = upper_spec, x = sens)) +
  geom_point(size = 2.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(limits = c(0.4, 1), n.breaks = 4) +
  scale_y_continuous(limits = c(0, 0.8), n.breaks = 6) +
  scale_shape_manual(values = c(18, 17, 15, 16, 4, 4)) +
  # scale_shape_manual(values = c(5, 2, 0, 1, 4, 4)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#1f78b4', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(x = 'Sensitivity', y = 'Specificity', tag = 'B')

p.comb.1.legend <- ggplot(df_acc %>% filter(!is.na(method)),
                          aes(x = sens, y = spec, shape = method, col = method)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'top') +
  guides(shape = guide_legend(nrow = 1)) +
  scale_shape_manual(values = c(18, 17, 15, 16, 4, 4)) +
  # scale_shape_manual(values = c(5, 2, 0, 1, 4, 4)) +
  scale_color_manual(values = c('#ff7f00', '#e31a1c', '#1f78b4', '#6a3d9a', '#33a02c', '#b2df8a')) +
  labs(shape = 'Method', col = 'Method')
p.comb.1.legend <- ggplotGrob(p.comb.1.legend)$grobs[[which(sapply(ggplotGrob(p.comb.1.legend)$grobs, function(x) x$name) == 'guide-box')]]

p.comb.1 <- arrangeGrob(p.comb.1.legend, arrangeGrob(p.comb.1a, p.comb.1b, nrow = 1), heights = c(0.15, 1))
plot(p.comb.1)
# ggsave(filename = 'results/plots/figures/Figure3.svg', p.comb.1, height = 4.8, width = 7.5)

# Plot accuracy by method, strength, and duration:
x_lab <- textGrob('Strength', gp = gpar(fontsize = 14, hjust = 1))
y_lab <- textGrob('Sensitivity', rot = 90, hjust = -0, gp = gpar(fontsize = 14))

title_1 <- textGrob('Corr. Coef.', gp = gpar(fontsize = 13), hjust = 2.03)
title_2 <- textGrob('GAMs', gp = gpar(fontsize = 13), hjust = 3.55)
title_3 <- textGrob('Granger Causality (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.85)
title_4 <- textGrob('Granger Causality (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.85)
title_5 <- textGrob('Transfer Entropy (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_6 <- textGrob('Transfer Entropy (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)
title_7 <- textGrob('CCM (Method 1) (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_8 <- textGrob('CCM (Method 1) (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)
title_9 <- textGrob('CCM (Method 2) (A \u279E B)', gp = gpar(fontsize = 13), hjust = 0.9)
title_0 <- textGrob('CCM (Method 2) (B \u279E A)', gp = gpar(fontsize = 13), hjust = 0.9)

p.comb.2 <- arrangeGrob(arrangeGrob(arrangeGrob(p.corr.1a, p.corr.1b, nrow = 1, widths = c(0.8, 1), top = title_1),
                                    arrangeGrob(p.gam.1a, p.gam.1b, nrow = 1, widths = c(0.8, 1), top = title_2),
                                    arrangeGrob(p.granger.1.3a, p.granger.1.3b, nrow = 1, widths = c(0.8, 1), top = title_3),
                                    arrangeGrob(p.granger.1.4a, p.granger.1.4b, nrow = 1, widths = c(0.8, 1), top = title_4),
                                    arrangeGrob(p.te.1.3a, p.te.1.3b, nrow = 1, widths = c(0.8, 1), top = title_5),
                                    arrangeGrob(p.te.1.4a, p.te.1.4b, nrow = 1, widths = c(0.8, 1), top = title_6),
                                    arrangeGrob(p.ccm.1.2a, p.ccm.1.2b, nrow = 1, widths = c(0.8, 1), top = title_7),
                                    arrangeGrob(p.ccm.1.5a, p.ccm.1.5b, nrow = 1, widths = c(0.8, 1), top = title_8),
                                    arrangeGrob(p.ccm.1.1a, p.ccm.1.1b, nrow = 1, widths = c(0.8, 1), top = title_9),
                                    arrangeGrob(p.ccm.1.4a, p.ccm.1.4b, nrow = 1, widths = c(0.8, 1), top = title_0),
                                    ncol = 2, bottom = x_lab, left = y_lab),
                        p.legend.1, heights = c(1, 0.05))
plot(p.comb.2)
# ggsave(filename = 'results/plots/figures/Figure4.svg', p.comb.2, width = 10, height = 10)
rm(p.corr.1, p.gam.1, p.granger.1.1, p.granger.1.2, p.granger.1.3, p.granger.1.4, p.te.1.1, p.te.1.2,
   p.te.1.3, p.te.1.4, p.ccm.1.1, p.ccm.1.4, p.ccm.1.2, p.ccm.1.5, p.ccm.1.3, p.ccm.1.6, p.legend.1)

# Plot accuracy for methods accounting for vs. not accounting for seasonality:
p.seas <- ggplot(df_acc_STORE %>% filter(str_detect(method, 'GC') | str_detect(method, 'TE')) %>%
                   mutate(seas = str_detect(confounding, 'Seas'), method = str_sub(method, 1, 2),
                          direction = if_else(direction == 'v1 -> v2', 'A %->% B', 'B %->% A')),
                 aes(x = sens, y = spec, shape = method, col = method, alpha = seas)) +
  geom_abline(slope = 1, intercept = 0, lty = 2, col = 'gray70') +
  geom_segment(aes(x = lower_sens, xend = upper_sens, y = spec)) +
  geom_segment(aes(y = lower_spec, yend = upper_spec, x = sens)) +
  geom_point(size = 2.5) +
  facet_wrap(~ direction, nrow = 1, labeller = label_parsed) +
  theme_bw() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_rect(fill = 'gray90'),
        panel.grid.minor = element_blank(),
        legend.position = 'top') +
  scale_x_continuous(limits = c(0.4, 1), n.breaks = 4) +
  scale_y_continuous(limits = c(0, 1), n.breaks = 6) +
  scale_shape_manual(values = c(15, 16)) +
  scale_color_manual(values = c('#1f78b4', '#6a3d9a')) +
  scale_alpha_manual(values = c(0.4, 1.0), guide = 'none') +
  labs(x = 'Sensitivity', y = 'Specificity', shape = 'Method', col = 'Method')
plot(p.seas)
# ggsave(filename = 'results/plots/figures/FigureS3.svg', p.seas, height = 4.5, width = 5)

# Plot effect of true interaction strength on point estimates:
df_assoc <- assoc_corr %>%
  mutate(method = 'Corr. Coef.', .before = coef) %>%
  bind_rows(assoc_gam %>%
              mutate(method = 'GAMs')) %>%
  bind_rows(assoc_granger_LIST[3:4] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'GC (w/ Seas)') %>%
              mutate(direction = str_sub(direction, 1, 8))) %>%
  bind_rows(assoc_te_LIST[3:4] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'TE (w/ Seas)') %>%
              mutate(direction = str_sub(direction, 1, 8))) %>%
  bind_rows(assoc_ccm_LIST_max[c(1, 4)] %>%
              bind_rows(.id = 'direction') %>%
              mutate(method = 'CCM') %>%
              mutate(direction = str_sub(direction, 1, 8)))

df_assoc <- df_assoc %>%
  mutate(direction = if_else(is.na(direction), 'v1 -> v2', direction)) %>%
  bind_rows(df_assoc %>% filter(method == 'Corr. Coef.') %>% mutate(direction = 'v2 -> v1')) %>%
  bind_rows(df_assoc %>% filter(method == 'GAMs') %>% mutate(direction = 'v2 -> v1'))

df_assoc <- df_assoc %>%
  mutate(method = factor(method, levels = c('Corr. Coef.', 'GAMs', 'GC (w/ Seas)', 'TE (w/ Seas)', 'CCM')))

df_assoc <- df_assoc %>%
  mutate(delta = factor(delta))

p3.a <- ggplot(df_assoc %>% filter(method == 'Corr. Coef.'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  # scale_color_brewer(palette = 'Set1') +
  labs(tag = 'A', x = NULL, y = '', title = 'Corr. Coef.')
p3.b <- ggplot(df_assoc %>% filter(method == 'GAMs'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'B', x = NULL, y = '', title = 'GAMs')
p3.c <- ggplot(df_assoc %>% filter(method == 'GC (w/ Seas)', direction == 'v1 -> v2'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'C', x = NULL, y = '', title = 'Granger Causality (A \u279E B)')
p3.d <- ggplot(df_assoc %>% filter(method == 'GC (w/ Seas)', direction == 'v2 -> v1'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'D', x = NULL, y = '', title = 'Granger Causality (B \u279E A)')
p3.e <- ggplot(df_assoc %>% filter(method == 'TE (w/ Seas)', direction == 'v1 -> v2'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'E', x = NULL, y = '', title = 'Transfer Entropy (A \u279E B)')
p3.f <- ggplot(df_assoc %>% filter(method == 'TE (w/ Seas)', direction == 'v2 -> v1'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'F', x = NULL, y = '', title = 'Transfer Entropy (B \u279E A)')
p3.g <- ggplot(df_assoc %>% filter(method == 'CCM', direction == 'v1 -> v2'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'G', x = NULL, y = '', title = 'CCM (A \u279E B)')
p3.h <- ggplot(df_assoc %>% filter(method == 'CCM', direction == 'v2 -> v1'),
               aes(x = 2**coef, xmin = 2**lower, xmax = 2**upper, y = delta, col = delta)) +
  geom_vline(xintercept = 2.0, lty = 1, col = 'red3') +
  geom_vline(xintercept = 1.0, lty = 2, col = 'gray70') +
  geom_pointrange() +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 20),
        plot.tag.position = c(0.025, 0.975),
        legend.position = 'none',
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  scale_x_continuous(breaks = c(0.75, 1.0, 1.25, 1.5, 1.75, 2.0), limits = c(0.715, 2.0)) +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(tag = 'H', x = NULL, y = '', title = 'CCM (B \u279E A)')

p3.legend <- ggplot(df_assoc %>%
                      mutate(delta = if_else(delta == '1', paste0(delta, ' week'), paste0(delta, ' weeks'))) %>%
                      mutate(delta = factor(delta, levels = c('1 week', '4 weeks', '13 weeks'))),
                    aes(x = 2**coef, y = delta, col = delta)) +
  geom_point() +
  theme_classic() +
  theme(legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'bottom') +
  scale_color_manual(values = viridis(6, option = 'C')[c(1, 3, 5)]) +
  labs(col = 'Duration')
p3.legend <- ggplotGrob(p3.legend)$grobs[[which(sapply(ggplotGrob(p3.legend)$grobs, function(x) x$name) == 'guide-box')]]

x_lab <- textGrob('Change in Point Estimate', gp = gpar(fontsize = 14, hjust = 1))

p.comb.3 <- arrangeGrob(arrangeGrob(p3.a, p3.b, p3.c, p3.d, p3.e, p3.f, p3.g, p3.h, ncol = 2, bottom = x_lab), p3.legend, heights = c(1, 0.085))
plot(p.comb.3)
# ggsave(filename = 'results/plots/figures/Figure5.svg', p.comb.3, width = 8, height = 7.5)

# Plot metric values vs. true interaction strength by duration:
p.comb.4 <- arrangeGrob(p.corr.2, p.gam.2, p.granger.2.1, p.granger.2.2, p.te.2.1, p.te.2.2, p.ccm.2.1, p.ccm.2.2, p.legend.2,
                        layout_matrix = rbind(c(1, 2), c(3, 4), c(5, 6), c(7, 8), c(9, 9)),
                        heights = c(1, 1, 1, 1, 0.25))
plot(p.comb.4)
# ggsave(filename = 'results/plots/figures/FigureS4.svg', p.comb.4, width = 10, height = 8.75)
rm(p.corr.2, p.gam.2, p.granger.2.1, p.granger.2.2, p.te.2.1, p.te.2.2, p.ccm.2.1, p.ccm.2.2, p.legend.2)

# Close pdf:
dev.off()

# Clean up:
rm(list = ls())
