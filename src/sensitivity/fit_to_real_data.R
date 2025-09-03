# ---------------------------------------------------------------------------------------------------------------------
# Code to check various methods using real-world data
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(gridExtra)
library(viridis)
library(mgcv)
library(gratia)
library(brms)

# Read in and format data:
hk_dat <- read_rds('../resp_virus_interactions/data/formatted/dat_hk_byOutbreak.rds')
can_dat <- read_csv('../resp_virus_interactions/data/formatted/dat_canada.csv')

hk_dat <- hk_dat$h1_plus_b_rsv %>%
  select(time, Year, Week, n_T:n_P2) %>%
  rename('h1b' = 'n_P1', 'rsv' = 'n_P2') %>%
  left_join(hk_dat$h3_rsv %>% select(time, Year, Week, n_P1) %>% rename('h3' = 'n_P1'),
            by = c('time', 'Year', 'Week')) %>%
  select(time:n_T, h1b, h3, rsv) %>%
  mutate(flu = h1b + h3, .after = h3) %>%
  drop_na()

can_dat <- can_dat %>%
  select(time, year, week, n_T1, n_T2, n_P1, n_P2) %>%
  rename('flu' = 'n_P1',
         'rsv' = 'n_P2',
         'Year' = 'year',
         'Week' = 'week')

# Visualize data:
p1a <- ggplot(hk_dat, aes(x = Week, y = flu / n_T, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.41)) + scale_color_viridis()
p1b <- ggplot(hk_dat, aes(x = Week, y = h1b / n_T, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.41)) + scale_color_viridis()
p1c <- ggplot(hk_dat, aes(x = Week, y = h3 / n_T, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.41)) + scale_color_viridis()
p1d <- ggplot(hk_dat, aes(x = Week, y = rsv / n_T, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.15)) + scale_color_viridis()

p1 <- arrangeGrob(p1a, p1b, p1c, p1d, layout_matrix = rbind(c(1, 4), c(2, 4), c(3, 4)))
plot(p1)

p2a <- ggplot(can_dat, aes(x = Week, y = flu / n_T1, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.35)) + scale_color_viridis()
p2b <- ggplot(can_dat, aes(x = Week, y = rsv / n_T2, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.35)) + scale_color_viridis()

p2 <- arrangeGrob(p2a, p2b, nrow = 1)
plot(p2)

# Prep data for input into models:
hk_dat <- hk_dat %>%
  mutate(h1b = h1b / n_T,
         h3 = h3 / n_T,
         flu = flu / n_T,
         rsv = rsv / n_T) %>%
  select(-n_T) %>%
  pivot_longer(h1b:flu) %>%
  rename('flu' = 'value',
         'flu_type' = 'name') %>%
  select(Year, Week, flu_type, flu, rsv) %>%
  split(f = .$flu_type)

can_dat <- can_dat %>%
  mutate(flu = flu / n_T1,
         rsv = rsv / n_T2) %>%
  select(Year:Week, flu:rsv)

# ---------------------------------------------------------------------------------------------------------------------

# GAMs

# Fit various GAMs (untransformed data):
res <- NULL

mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(Week, bs = 'cc', k = 53)) + set_rescor(TRUE)
for (nm in names(hk_dat)) {
  print(nm)
  
  tic <- Sys.time()
  mvn_mod <- brm(mvn_mod_form, data = hk_dat[[nm]], chains = 4, cores = 4, warmup = 2000, iter = 3000)
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)
  
  corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__flu__rsv')$rescor__flu__rsv
  n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
    sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
  }) %>% unlist() %>% sum()
  
  res <- bind_rows(res, data.frame(cbind(loc = 'hk', flu_type = nm, trans = 'no', cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                                   row.names = ''))
}

mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(Week, bs = 'cc', k = 52)) + set_rescor(TRUE)

tic <- Sys.time()
mvn_mod <- brm(mvn_mod_form, data = can_dat, chains = 4, cores = 4, warmup = 2000, iter = 3000)
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'mins'
print(etime)

corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__flu__rsv')$rescor__flu__rsv
n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
  sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
}) %>% unlist() %>% sum()

res <- bind_rows(res, data.frame(cbind(loc = 'can', flu_type = 'flu', trans = 'no', cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                                 row.names = ''))

# Fit various GAMs (transformed data):
hk_dat <- lapply(hk_dat, function(ix) {
  ix %>%
    mutate(flu = scale(log(flu), scale = FALSE)[, 1],
           rsv = scale(log(rsv), scale = FALSE)[, 1])
})
can_dat <- can_dat %>%
  mutate(flu = scale(log(flu), scale = FALSE)[, 1],
         rsv = scale(log(rsv), scale = FALSE)[, 1])

mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(Week, bs = 'cc', k = 53)) + set_rescor(TRUE)
for (nm in names(hk_dat)) {
  print(nm)
  
  tic <- Sys.time()
  mvn_mod <- brm(mvn_mod_form, data = hk_dat[[nm]], chains = 4, cores = 4, warmup = 2000, iter = 3000)
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'mins'
  print(etime)
  
  corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__flu__rsv')$rescor__flu__rsv
  n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
    sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
  }) %>% unlist() %>% sum()
  
  res <- bind_rows(res, data.frame(cbind(loc = 'hk', flu_type = nm, trans = 'yes', cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                                   row.names = ''))
}

mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(Week, bs = 'cc', k = 52)) + set_rescor(TRUE)

tic <- Sys.time()
mvn_mod <- brm(mvn_mod_form, data = can_dat, chains = 4, cores = 4, warmup = 2000, iter = 3000)
toc <- Sys.time()
etime <- toc - tic
units(etime) <- 'mins'
print(etime)

corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__flu__rsv')$rescor__flu__rsv
n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
  sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
}) %>% unlist() %>% sum()

res <- bind_rows(res, data.frame(cbind(loc = 'can', flu_type = 'flu', trans = 'yes', cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                                 row.names = ''))

# Format and summarize results:
res %>%
  as_tibble() %>%
  filter(if_any(b_flu_Intercept:lp__, ~ . > 1.05)) %>%
  print()
res %>%
  as_tibble() %>%
  filter(n_div > 0) %>%
  print()

res %>%
  as_tibble() %>%
  select(loc:trans, cor_median:CI_upper95, n_div) %>%
  mutate(across(cor_median:CI_upper95, as.numeric)) %>%
  ggplot() +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(aes(x = trans, y = cor_median, ymin = CI_lower95, ymax = CI_upper95, col = flu_type),
                  size = 1) +
  facet_wrap(~ loc, ncol = 1) +
  theme_classic() +
  scale_color_brewer(palette = 'Set1')

