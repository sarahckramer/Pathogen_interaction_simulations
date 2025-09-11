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
library(vars)
library(tseries)
library(VARtests)
library(rJava)
library(rEDM)
library(Kendall)
library(pracma)
library(parallel)
library(doSNOW)

# Read in and format data:
hk_dat <- read_rds('../resp_virus_interactions/data/formatted/dat_hk_byOutbreak.rds')
can_dat <- read_csv('../resp_virus_interactions/data/formatted/dat_canada.csv')

hk_dat <- hk_dat$h1_plus_b_rsv %>%
  dplyr::select(time, Year, Week, n_T:n_P2) %>%
  rename('h1b' = 'n_P1', 'rsv' = 'n_P2') %>%
  left_join(hk_dat$h3_rsv %>% dplyr::select(time, Year, Week, n_P1) %>% rename('h3' = 'n_P1'),
            by = c('time', 'Year', 'Week')) %>%
  dplyr::select(time:n_T, h1b, h3, rsv) %>%
  mutate(flu = h1b + h3, .after = h3) %>%
  drop_na() %>%
  mutate(h1b = h1b / n_T,
         flu = flu / n_T,
         rsv = rsv / n_T) %>%
  dplyr::select(-c(n_T, h3)) %>%
  pivot_longer(h1b:flu) %>%
  rename('flu' = 'value',
         'flu_type' = 'name') %>%
  dplyr::select(Year, Week, flu_type, flu, rsv) %>%
  split(f = .$flu_type)

can_dat <- can_dat %>%
  dplyr::select(time, year, week, n_T1, n_T2, n_P1, n_P2) %>%
  rename('flu' = 'n_P1',
         'rsv' = 'n_P2',
         'Year' = 'year',
         'Week' = 'week') %>%
  mutate(flu = flu / n_T1,
         rsv = rsv / n_T2) %>%
  dplyr::select(Year:Week, flu:rsv)

# Visualize data:
p1a <- ggplot(hk_dat$flu, aes(x = Week, y = flu, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.41)) + scale_color_viridis()
p1b <- ggplot(hk_dat$h1b, aes(x = Week, y = flu, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.41)) + scale_color_viridis()
p1c <- ggplot(hk_dat$flu, aes(x = Week, y = rsv, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.15)) + scale_color_viridis()

p1 <- arrangeGrob(p1a, p1b, p1c, layout_matrix = rbind(c(1, 3), c(2, 3)))

p2a <- ggplot(can_dat, aes(x = Week, y = flu, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.35)) + scale_color_viridis()
p2b <- ggplot(can_dat, aes(x = Week, y = rsv, group = Year, col = Year)) + geom_line() + geom_point() + theme_classic() + scale_y_continuous(limits = c(0, 0.35)) + scale_color_viridis()

p2 <- arrangeGrob(p2a, p2b, nrow = 1)

grid.arrange(p1, p2, ncol = 1)
rm(p1a, p1b, p1c, p2a, p2b)

# ---------------------------------------------------------------------------------------------------------------------

# Transform data

hk_dat <- lapply(hk_dat, function(ix) {
  ix %>%
    mutate(flu = scale(log(flu), scale = FALSE)[, 1],
           rsv = scale(log(rsv), scale = FALSE)[, 1])
})
can_dat <- can_dat %>%
  mutate(flu = scale(log(flu), scale = FALSE)[, 1],
         rsv = scale(log(rsv), scale = FALSE)[, 1])

# ---------------------------------------------------------------------------------------------------------------------

# Pearson correlation

# Fit correlations:
cor_hk_flu <- cor.test(hk_dat$flu$flu, hk_dat$flu$rsv)
cor_hk_h1b <- cor.test(hk_dat$h1b$flu, hk_dat$h1b$rsv)
cor_can <- cor.test(can_dat$flu, can_dat$rsv)

# Format output:
res_cor_LIST <- list(cor_hk_flu, cor_hk_h1b, cor_can)

res_cor <- bind_cols(loc = c('hk', 'hk', 'can'),
                     flu_type = c('flu', 'h1b', 'flu'),
                     cor = lapply(res_cor_LIST, getElement, 'estimate') %>%
                       unlist(),
                     CI_lower_95 = lapply(lapply(res_cor_LIST, getElement, 'conf.int'), '[[', 1) %>%
                       unlist(),
                     CI_upper_95 = lapply(lapply(res_cor_LIST, getElement, 'conf.int'), '[[', 2) %>%
                       unlist(),
                     p_value = lapply(res_cor_LIST, getElement, 'p.value') %>%
                       unlist())

# Clean up:
rm(cor_hk_flu, cor_hk_h1b, cor_can, res_cor_LIST)

# ---------------------------------------------------------------------------------------------------------------------

# GAMs

# Fit models:
res_gam <- NULL

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
  
  res_gam <- bind_rows(res_gam, data.frame(cbind(loc = 'hk', flu_type = nm, cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
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

res_gam <- bind_rows(res_gam, data.frame(cbind(loc = 'can', flu_type = 'flu', cor = mean(corrs_alt), cor_median = median(corrs_alt), CI_lower95 = quantile(corrs_alt, p = 0.025), CI_upper95 = quantile(corrs_alt, p = 0.975), t(rhat(mvn_mod)), n_div = n_divergent),
                                         row.names = ''))

# Check for convergence:
res_gam %>%
  as_tibble() %>%
  filter(if_any(b_flu_Intercept:lp__, ~ . > 1.05)) %>%
  print()
res_gam %>%
  as_tibble() %>%
  filter(n_div > 0) %>%
  print()

# Format results:
res_gam <- res_gam %>%
  as_tibble() %>%
  dplyr::select(loc:flu_type, cor_median:CI_upper95, n_div) %>%
  mutate(across(cor_median:CI_upper95, as.numeric))

# Clean up:
rm(mvn_mod, mvn_mod_form, corrs_alt, n_divergent, nm, tic, toc, etime)

# ---------------------------------------------------------------------------------------------------------------------

# Granger causality

# Get climate info:
dat_clim_hk <- read_csv('../resp_virus_interactions/data/formatted/clim_dat_hk.csv')
dat_clim_can <- read_csv('../resp_virus_interactions/data/formatted/clim_dat_can.csv')

hk_dat <- lapply(hk_dat, function(ix) {
  ix %>%
    left_join(dat_clim_hk %>%
                dplyr::select(year, week, temp, ah) %>%
                rename('Year' = 'year', 'Week' = 'week'),
              by = c('Year', 'Week'))
})
can_dat <- can_dat %>%
  left_join(dat_clim_can %>%
              dplyr::select(YEAR, WEEK, temp, ah) %>%
              rename('Year' = 'YEAR', 'Week' = 'WEEK'),
            by = c('Year', 'Week'))

rm(dat_clim_hk, dat_clim_can)

# Log-transform and scale climate data:
hk_dat <- lapply(hk_dat, function(ix) {
  ix %>%
    mutate(temp = scale(log(temp), scale = FALSE)[, 1],
           ah = scale(log(ah), scale = FALSE)[, 1])
})
can_dat <- can_dat %>%
  mutate(temp = scale(temp)[, 1],
         ah = scale(log(ah), scale = FALSE)[, 1])

# Run GC analysis:
res_gc <- NULL

for (i in 1:3) {
  
  if (i == 3) {
    dat <- can_dat
    loc_temp <- 'can'
  } else {
    dat <- hk_dat[[i]]
    loc_temp <- 'hk'
  }
  
  if (i == 2) {
    type_temp <- 'h1b'
  } else {
    type_temp <- 'flu'
  }
  
  # determine best lags:
  df <- dat %>% dplyr::select(flu, rsv)
  
  lags <- lapply(df, VARselect, lag.max = 20)
  
  lag_v1 <- as.numeric(lags$flu$selection[3])
  lag_v2 <- as.numeric(lags$rsv$selection[3])
  rm(lags)
  
  # check stationary:
  adf_v1 <- adf.test(dat$flu, k = lag_v1)
  adf_v2 <- adf.test(dat$rsv, k = lag_v2)
  
  kpss_v1 <- kpss.test(dat$flu)
  kpss_v2 <- kpss.test(dat$rsv)
  
  # get lag value for further analysis:
  p <- as.numeric(VARselect(df, lag.max = 20)$selection[3])
  rm(df)
  
  # estimate GC:
  var1 <- VAR(y = dat[, c('flu', 'rsv')], p = p, exogen = dat[, c('temp', 'ah')])
  
  sink(file = nullfile())
  ar_v1 <- VARfit(dat$flu, p = p, exogen = dat[, c('temp', 'ah')])
  ar_v2 <- VARfit(dat$rsv, p = p, exogen = dat[, c('temp', 'ah')])
  sink()
  
  logRSS_v1 <- log(sum(ar_v1$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'flu'] ** 2))
  logRSS_v2 <- log(sum(ar_v2$resid ** 2, na.rm = TRUE) / sum(residuals(var1)[, 'rsv'] ** 2))
  
  # get p-values
  p1_ftest <- causality(var1, cause = 'rsv')$Granger$p.value
  p2_ftest <- causality(var1, cause = 'flu')$Granger$p.value
  
  # output results
  res_gc <- bind_rows(res_gc, data.frame(cbind(logRSS = c(logRSS_v1, logRSS_v2),
                                               direction = c('rsv -> flu', 'flu -> rsv'),
                                               ftest_p = c(p1_ftest, p2_ftest),
                                               adf_p = c(adf_v1$p.value, adf_v2$p.value),
                                               kpss_p = c(kpss_v1$p.value, kpss_v2$p.value))) %>%
                        as_tibble() %>%
                        mutate(loc = loc_temp,
                               flu_type = type_temp,
                               logRSS = as.numeric(logRSS),
                               ftest_p = as.numeric(ftest_p),
                               adf_p = as.numeric(adf_p),
                               kpss_p = as.numeric(kpss_p)))
  
}

# Clean up:
rm(lag_v1, lag_v2, p, var, var1, ar_v1, ar_v2, logRSS_v1, logRSS_v2, p1_ftest, p2_ftest, adf_v1, adf_v2, kpss_v1, kpss_v2,
   i, dat, loc_temp, type_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Transfer entropy

# Start Java:
.jinit()
.jaddClassPath('infodynamics-dist-1.6.1/infodynamics.jar')

# Run TE:
res_te <- NULL

for (i in 1:3) {
  
  if (i == 3) {
    dat <- can_dat
    loc_temp <- 'can'
  } else {
    dat <- hk_dat[[i]]
    loc_temp <- 'hk'
  }
  
  if (i == 2) {
    type_temp <- 'h1b'
  } else {
    type_temp <- 'flu'
  }
  
  # get relevant data:
  sourceArray <- dat$flu
  destArray <- dat$rsv
  condArray <- dat$temp # dat[, c('temp', 'ah')]
  
  # create and set up TE calculator:
  teCalc <- .jnew('infodynamics/measures/continuous/kraskov/ConditionalTransferEntropyCalculatorKraskov')
  
  .jcall(teCalc, 'V', 'setProperty', 'k', '4')
  .jcall(teCalc, 'V', 'setProperty', 'delay', '4')
  
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray, condArray)
  
  # set embedding dimensions:
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '5')
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20')
  
  # run analysis:
  result_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  nullDist_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                             'computeSignificance', 500L)
  p_value_v2_x_v1 <- nullDist_v2_x_v1$pValue
  
  # set up reverse direction:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray, condArray)
  
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '2')
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20')
  
  # run analysis:
  result_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  nullDist_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                             'computeSignificance', 500L)
  p_value_v1_x_v2 <- nullDist_v1_x_v2$pValue
  
  # output results:
  res_te <- bind_rows(res_te, data.frame(cbind(loc = loc_temp,
                                               flu_type = type_temp,
                                               te = c(result_v1_x_v2, result_v2_x_v1),
                                               direction = c('rsv -> flu', 'flu -> rsv'),
                                               p_value = c(p_value_v1_x_v2, p_value_v2_x_v1))) %>%
                        as_tibble() %>%
                        mutate(across(c(te, p_value), as.numeric)))
  
}

# Clean up:
rm(sourceArray, destArray, condArray, teCalc, result_v2_x_v1, result_v1_x_v2, nullDist_v2_x_v1, nullDist_v1_x_v2,
   p_value_v2_x_v1, p_value_v1_x_v2, i, dat, loc_temp, type_temp)

# ---------------------------------------------------------------------------------------------------------------------

# CCM

# Run CCM:
res_ccm <- NULL

for (i in 1:3) {
  
  if (i == 3) {
    dat <- can_dat
    loc_temp <- 'can'
  } else {
    dat <- hk_dat[[i]]
    loc_temp <- 'hk'
  }
  
  if (i == 2) {
    type_temp <- 'h1b'
  } else {
    type_temp <- 'flu'
  }
  
  # Get only relevant columns and add column for time:
  dat <- dat %>%
    mutate(time = 1:nrow(dat)) %>%
    dplyr::select(time, flu:rsv)
  
  # set library and prediction sets:
  lib <- paste('1', nrow(dat), sep = ' ')
  pred <- lib
  
  # get embedding dimensions:
  E_v1 <- EmbedDimension(dataFrame = dat, columns = 'flu', target = 'flu',
                         lib = lib, pred = pred, maxE = 2, showPlot = FALSE)
  E_v1 <- E_v1 %>% slice_max(rho) %>% pull(E) # keep the row with max prediction skill
  
  # get E for v2
  E_v2 <- EmbedDimension(dataFrame = dat, columns = "rsv", target = "rsv",
                         lib = lib, pred = pred, maxE = 5, showPlot = FALSE)
  E_v2 <- E_v2 %>% slice_max(rho) %>% pull(E)
  
  # get best negative tp:
  vars <- c('flu', 'rsv')
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -20:0) %>%
    filter(lib_column != target_column) %>%
    mutate(E = if_else(lib_column == 'flu', E_v1, E_v2))
  lib_size_tp <- nrow(dat) - (20 - 1) - (max(params$E) - 1)
  
  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = dat, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]), 
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))
  
  optimal_tp_v1xv2 <- output %>% filter(target_column == 'rsv') %>% filter(`flu:rsv` == max(`flu:rsv`)) %>% pull(tp)
  optimal_tp_v2xv1 <- output %>% filter(target_column == 'flu') %>% filter(`rsv:flu` == max(`rsv:flu`)) %>% pull(tp)
  
  highest_crossmap_cor_v1xv2 <- output %>% filter(target_column == 'rsv') %>% pull(`flu:rsv`) %>% max()
  highest_crossmap_cor_v2xv1 <- output %>% filter(target_column == 'flu') %>% pull(`rsv:flu`) %>% max()
  
  # run CCM:
  lib_max <- nrow(dat) - (max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1)) - 1) - (max(abs(E_v1), abs(E_v2)) - 1)
  if (optimal_tp_v1xv2 == 0 & optimal_tp_v2xv1 == 0) lib_max <- lib_max - 1
  
  v1_xmap_v2 <- CCM(dataFrame = dat, E = E_v1, columns = 'flu', target = 'rsv', 
                    libSizes = c(E_v1 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v1xv2, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  v2_xmap_v1 <- CCM(dataFrame = dat, E = E_v2, columns = 'rsv', target = 'flu', 
                    libSizes = c(E_v2 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v2xv1, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  # extract and format results:
  res <- v1_xmap_v2$LibMeans %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize)) %>%
    dplyr::select(1:2) %>%
    inner_join(v2_xmap_v1$LibMeans %>%
                 mutate(LibSize = if_else(LibSize < 20, 12, LibSize)) %>%
                 dplyr::select(1:2),
               by = 'LibSize') %>%
    rename('mean_v1xv2' = 'flu:rsv',
           'mean_v2xv1' = 'rsv:flu') %>%
    pivot_longer(-LibSize, names_to = 'direction', values_to = 'rho') %>%
    mutate(direction = if_else(direction == 'mean_v1xv2', 'v2 -> v1', 'v1 -> v2'))
  
  # get p-value (method 2):
  all_predictions_v1 <- v1_xmap_v2$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  all_predictions_v2 <- v2_xmap_v1$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  
  corrs_Lmin_v1 <- all_predictions_v1 %>% filter(LibSize == min(LibSize)) %>% pull(rho)
  corrs_Lmax_v1 <- all_predictions_v1 %>% filter(LibSize == max(LibSize)) %>% pull(rho)
  
  corrs_Lmin_v2 <- all_predictions_v2 %>% filter(LibSize == min(LibSize)) %>% pull(rho)
  corrs_Lmax_v2 <- all_predictions_v2 %>% filter(LibSize == max(LibSize)) %>% pull(rho)
  
  inverse_quantile <- function(x, y) {
    # Based on: https://github.com/cobeylab/pyembedding/blob/master/statutils.py
    
    x <- sort(x)
    
    if (x[1] == x[length(x)]) {
      print('Error: Same initial and final value for Lmin')
    } else {
      quantiles <- seq(0, 1, length.out = length(x))
    }
    
    y[y < x[1]] <- x[1]
    y[y > x[length(x)]] <- x[length(x)]
    
    inv_q_y <- interp1(x, quantiles, xi = y)
    
    return(inv_q_y)
    
  }
  
  p_conv_v1xv2 <- 1 - mean(inverse_quantile(corrs_Lmin_v1, corrs_Lmax_v1))
  p_conv_v2xv1 <- 1 - mean(inverse_quantile(corrs_Lmin_v2, corrs_Lmax_v2))
  
  # get p-value (method 1):
  num_surr <- 100
  
  surr_v1 <- SurrogateData(dat$flu, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  surr_v2 <- SurrogateData(dat$rsv, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  
  surr_dat_list_v1xv2 <- vector('list', length = num_surr)
  surr_dat_list_v2xv1 <- vector('list', length = num_surr)
  
  for (i in 1:num_surr) {
    
    surr_dat_list_v1xv2[[i]] <- dat %>%
      dplyr::select(time, flu) %>%
      mutate(rsv = surr_v2[, i])
    
    surr_dat_list_v2xv1[[i]] <- dat %>%
      dplyr::select(time, rsv) %>%
      mutate(flu = surr_v1[, i])
    
  }
  
  lib_max_use <- max(res$LibSize)
  
  cl <- makeCluster(10, type = 'SOCK')
  registerDoSNOW(cl)
  surr_res <- foreach(i = 1:num_surr, .packages = c('rEDM', 'tidyverse')) %dopar% {
    
    ccm_v1xv2 <- CCM(dataFrame = surr_dat_list_v1xv2[[i]], E = E_v1, columns = 'flu', target = 'rsv', 
                     libSizes = lib_max_use, Tp = optimal_tp_v1xv2, random = TRUE, sample = 100, includeData = TRUE)
    ccm_v2xv1 <- CCM(dataFrame = surr_dat_list_v2xv1[[i]], E = E_v2, columns = 'rsv', target = 'flu', 
                     libSizes = lib_max_use, Tp = optimal_tp_v2xv1, random = TRUE, sample = 100, includeData = TRUE)
    
    temp_means <- ccm_v1xv2$LibMeans %>%
      dplyr::select(1:2) %>%
      rename('rho' = 'flu:rsv') %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$LibMeans %>%
                  dplyr::select(1:2) %>%
                  rename('rho' = 'rsv:flu') %>%
                  mutate(direction = 'v1 -> v2'))
    
    temp_ci <- ccm_v1xv2$CCM1_PredictStat %>%
      summarise(rho_median = quantile(rho, probs = 0.5),
                rho_ci_lower = quantile(rho, probs = 0.025),
                rho_ci_upper = quantile(rho, probs = 0.975)) %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$CCM1_PredictStat %>%
                  summarise(rho_median = quantile(rho, probs = 0.5),
                            rho_ci_lower = quantile(rho, probs = 0.025),
                            rho_ci_upper = quantile(rho, probs = 0.975)) %>%
                  mutate(direction = 'v1 -> v2'))
    
    res_temp <- temp_means %>%
      left_join(temp_ci, by = 'direction')
    
    return(res_temp)
    
  }
  stopCluster(cl)
  
  # combine all results
  surr_res <- bind_rows(surr_res)
  
  # estimate p-value using empirical cumulative distribution
  rho_v1xv2_surr <- surr_res %>% filter(direction == 'v2 -> v1') %>% pull(rho)
  rho_v2xv1_surr <- surr_res %>% filter(direction == 'v1 -> v2') %>% pull(rho)
  
  rho_v1xv2 <- res %>% filter(direction == 'v2 -> v1', LibSize == lib_max_use) %>% pull(rho)
  rho_v2xv1 <- res %>% filter(direction == 'v1 -> v2', LibSize == lib_max_use) %>% pull(rho)
  
  # as in ha0ye tutorial:
  p_surr_v1xv2 <- (sum(rho_v1xv2 < rho_v1xv2_surr) + 1) / (length(rho_v1xv2_surr) + 1)
  p_surr_v2xv1 <- (sum(rho_v2xv1 < rho_v2xv1_surr) + 1) / (length(rho_v2xv1_surr) + 1)
  
  # output results:
  res <- res %>%
    filter(LibSize == lib_max_use) %>%
    mutate(E = if_else(direction == 'v2 -> v1', E_v1, E_v2),
           tp_use = if_else(direction == 'v2 -> v1', optimal_tp_v1xv2, optimal_tp_v2xv1),
           max_cmc = if_else(direction == 'v2 -> v1', highest_crossmap_cor_v1xv2, highest_crossmap_cor_v2xv1),
           p_conv = if_else(direction == 'v2 -> v1', p_conv_v1xv2, p_conv_v2xv1),
           p_surr = if_else(direction == 'v2 -> v1', p_surr_v1xv2, p_surr_v2xv1))
  res_ccm <- bind_rows(res_ccm, res %>%
                         mutate(loc = loc_temp,
                                flu_type = type_temp,
                                .before = LibSize))
  
}

# Clean up:
rm(lib, pred, E_v1, E_v2, vars, params, lib_size_tp, output, optimal_tp_v1xv2, optimal_tp_v2xv1,
   highest_crossmap_cor_v1xv2, highest_crossmap_cor_v2xv1, lib_max, v1_xmap_v2, v2_xmap_v1,
   all_predictions_v1, all_predictions_v2, corrs_Lmin_v1, corrs_Lmax_v1, corrs_Lmin_v2, corrs_Lmax_v2,
   inverse_quantile, p_conv_v1xv2, p_conv_v2xv1, num_surr, surr_v1, surr_v2, surr_dat_list_v1xv2,
   surr_dat_list_v2xv1, lib_max_use, cl, surr_res, rho_v1xv2, rho_v2xv1, rho_v1xv2_surr, rho_v2xv1_surr,
   p_surr_v1xv2, p_surr_v2xv1, i, res, dat, loc_temp, type_temp)

# ---------------------------------------------------------------------------------------------------------------------

# Compile and visualize results

print(res_gc) # only flu -> rsv for Canada sig
print(res_te) # varies by run? all HK sig, none Canada
print(res_ccm) # flu->rsv, h1b->rsv sig for HK (former only with method 2)

res <- res_cor %>%
  dplyr::select(-p_value) %>%
  mutate(met = 'cor') %>%
  rename('CI_lower95' = 'CI_lower_95',
         'CI_upper95' = 'CI_upper_95') %>%
  bind_rows(res_gam %>%
              dplyr::select(-n_div) %>%
              mutate(met = 'gam') %>%
              rename('cor' = 'cor_median'))

p3 <- ggplot(res %>% mutate(x_art = c(0.95, 1.95, 2.95, 1.05, 2.05, 3.05)) %>%
               mutate(met = if_else(met == 'cor', 'Pearson', 'GAMs')) %>%
               mutate(met = factor(met, levels = c('Pearson', 'GAMs'))),
             aes(x = x_art, y = cor, ymin = CI_lower95, ymax = CI_upper95, col = met)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_pointrange(size = 0.7) +
  theme_classic() +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c('Hong Kong (All Flu)', 'Hong Kong (H1+B)', 'Canada (All Flu)')) +
  scale_y_continuous(breaks = seq(-0.3, 0.9, by = 0.1), labels = c(-0.3, -0.2, -0.1, 0, seq(0.1, 0.9, by = 0.1))) +
  scale_color_brewer(palette = 'Set1') +
  labs(x = NULL, y = 'Corr. Coef.', col = 'Method')
print(p3)

p4a <- ggplot(res_gc %>% mutate(x_art = c(1.05, 0.95, 2.05, 1.95, 3.05, 2.95)),
              aes(x = x_art, y = logRSS, col = direction, alpha = ftest_p < 0.05)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c('Hong Kong (All Flu)', 'Hong Kong (H1+B)', 'Canada (All Flu)')) +
  # scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  scale_color_brewer(palette = 'Set1') +
  scale_alpha_discrete(range = c(0.15, 1)) +
  labs(x = NULL, y = 'logRSS (GC)', col = 'Direction', alpha = 'Sig.')
p4b <- ggplot(res_te %>% mutate(x_art = c(1.05, 0.95, 2.05, 1.95, 3.05, 2.95)),
              aes(x = x_art, y = te, col = direction, alpha = p_value < 0.05)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c('Hong Kong (All Flu)', 'Hong Kong (H1+B)', 'Canada (All Flu)')) +
  # scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  scale_color_brewer(palette = 'Set1') +
  scale_alpha_discrete(range = c(0.15, 1)) +
  labs(x = NULL, y = 'TE', col = 'Direction', alpha = 'Sig.')
p4c <- ggplot(res_ccm %>% mutate(x_art = c(1.05, 0.95, 2.05, 1.95, 3.05, 2.95)) %>% pivot_longer(p_conv:p_surr, names_to = 'p_type', values_to = 'p_val') %>% mutate(x_art = if_else(p_type == 'p_conv', x_art - 0.025, x_art + 0.025)),
              aes(x = x_art, y = rho, col = direction, alpha = p_val < 0.05, shape = p_type)) +
  geom_point(size = 3) +
  geom_point(size = 3) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_x_continuous(breaks = c(1, 2, 3),
                     labels = c('Hong Kong (All Flu)', 'Hong Kong (H1+B)', 'Canada (All Flu)')) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.1)) +
  scale_color_brewer(palette = 'Set1') +
  scale_alpha_discrete(range = c(0.15, 1)) +
  scale_shape_manual(values = c(16, 15)) +
  labs(x = NULL, y = 'rho (CCM)', col = 'Direction', alpha = 'Sig.', shape = 'Method')
grid.arrange(p4a, p4b, p4c, ncol = 1)

# Output results:
pdf('results/plots/fit_to_real_data.pdf', width = 15, height = 7)
grid.arrange(p1, p2, ncol = 1)
print(p3)
grid.arrange(p4a, p4b, p4c, ncol = 1)
dev.off()
