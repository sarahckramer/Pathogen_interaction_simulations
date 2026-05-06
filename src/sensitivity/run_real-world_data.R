# ---------------------------------------------------------------------------------------------------------------------
# Code to check various methods using real-world data
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tidyverse)
library(testthat)
library(mgcv)
library(gratia)
library(brms)
library(tseries)
library(vars)
library(VARtests)
library(rJava)
library(rEDM)
library(Kendall)
library(pracma)
library(parallel)
library(doSNOW)
library(gridExtra)

# ---------------------------------------------------------------------------------------------------------------------

# Prep all data

#---- Read in virologic data ----
hk_dat <- read_rds('data/dat_hk_byOutbreak.rds')
can_dat <- read_csv('data/dat_canada.csv')

#---- Read in climate data ----
dat_clim_hk <- read_csv('data/clim_dat_hk.csv')
dat_clim_can <- read_csv('data/clim_dat_can_NEW.csv')

#---- Get proportion positive at each time point ----
hk_dat <- hk_dat$h1_plus_b_rsv %>%
  mutate(flu = n_P1 / n_T,
         rsv = n_P2 / n_T) %>%
  dplyr::select(Year, Week, flu:rsv) %>%
  rename_with(tolower) %>%
  drop_na()

can_dat <- can_dat %>%
  mutate(flu = n_P1 / n_T1,
         rsv = n_P2 / n_T2) %>%
  dplyr::select(year:week, flu:rsv)

#---- Join climate data ----
hk_dat <- hk_dat %>%
  left_join(dat_clim_hk %>%
              dplyr::select(year:week, temp, ah),
            by = c('year', 'week'))
can_dat <- can_dat %>%
  left_join(dat_clim_can %>%
              dplyr::select(year:week, temp, ah),
            by = c('year', 'week'))
rm(dat_clim_hk, dat_clim_can)

#---- Log-transform and scale data ----
hk_dat <- hk_dat %>%
  mutate(across(flu:ah, ~ scale(log(.x), scale = FALSE)[, 1]))
can_dat <- can_dat %>%
  mutate(across(flu:ah, ~ scale(log(.x), scale = FALSE)[, 1]))

#---- Create list of data ----
dat_list <- list(hk_dat, can_dat)
rm(hk_dat, can_dat)

# ---------------------------------------------------------------------------------------------------------------------

# Pearson correlation

#---- Create list of results ----
cor_list <- vector('list', length(dat_list))

#---- Fit correlations ----
for (i in 1:length(dat_list)) {
  
  cor_temp <- cor.test(dat_list[[i]]$flu, dat_list[[i]]$rsv)
  
  cor_list[[i]] <- bind_cols(cor = cor_temp$estimate,
                             CI_lower_95 = cor_temp$conf.int[1],
                             CI_upper_95 = cor_temp$conf.int[2],
                             p_value = cor_temp$p.value)
  
}

#---- Clean up ----
rm(cor_temp, i)

# ---------------------------------------------------------------------------------------------------------------------

# GAMs

#---- Create list of results ----
gam_list <- vector('list', length(dat_list))

#---- Fit GAMs ----
for (i in 1:length(dat_list)) {
  
  if (i == 1) {
    mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(week, bs = 'cc', k = 53)) + set_rescor(TRUE)
  } else if (i == 2) {
    mvn_mod_form <- bf(mvbind(flu, rsv) ~ s(week, bs = 'cc', k = 52)) + set_rescor(TRUE)
  }
  
  n_divergent <- 1
  
  while (n_divergent > 0) {
    mvn_mod <- brm(mvn_mod_form, data = dat_list[[i]], chains = 4, cores = 4, warmup = 2000, iter = 3000)
    
    corrs_alt <- as_draws_df(mvn_mod, variable = 'rescor__flu__rsv')$rescor__flu__rsv
    n_divergent <- lapply(mvn_mod$fit@sim$samples, function(ix) {
      sum(attr(ix, 'sampler_params')[['divergent__']][2001:3000])
    }) %>% unlist() %>% sum()
  }
  
  gam_list[[i]] <- bind_cols(cor_median = median(corrs_alt),
                             CI_lower95 = quantile(corrs_alt, p = 0.025),
                             CI_upper95 = quantile(corrs_alt, p = 0.975),
                             t(rhat(mvn_mod)),
                             n_div = n_divergent)
  
}

#---- Check for convergence ----
for (i in 1:length(gam_list)) {
  expect_equal(gam_list[[i]] %>% filter(if_any(b_flu_Intercept:lp__, ~ . > 1.05)) %>% nrow(), 0)
  expect_equal(gam_list[[i]] %>% filter(n_div > 0) %>% nrow(), 0)
}

#---- Limit to results columns ----
gam_list <- lapply(gam_list, function(ix) {
  ix %>% dplyr::select(cor_median:CI_upper95)
}
)

#---- Clean up ----
rm(mvn_mod, mvn_mod_form, corrs_alt, n_divergent, i)

# ---------------------------------------------------------------------------------------------------------------------

# Granger causality

#---- Create list of results ----
gc_list <- vector('list', length(dat_list))

#---- Run Granger causality analysis ----
for (i in 1:length(dat_list)) {
  
  # Determine best lags:
  df <- dat_list[[i]] %>% dplyr::select(flu, rsv)
  lags <- lapply(df, VARselect, lag.max = 20)
  lag_v1 <- as.numeric(lags$flu$selection[3])
  lag_v2 <- as.numeric(lags$rsv$selection[3])
  rm(lags)
  
  # Check for stationary data:
  adf_v1 <- adf.test(dat_list[[i]]$flu, k = lag_v1)
  adf_v2 <- adf.test(dat_list[[i]]$rsv, k = lag_v2)
  
  kpss_v1 <- kpss.test(dat_list[[i]]$flu)
  kpss_v2 <- kpss.test(dat_list[[i]]$rsv)
  
  # Get lag for analysis:
  p <- as.numeric(VARselect(df, lag.max = 20)$selection[3])
  
  # Estimate GC:
  var_confound <- VAR(y = dat_list[[i]][, c('flu', 'rsv')], p = p, exogen = dat_list[[i]][, 'ah'])
  
  sink(file = nullfile())
  ar_v1_confound <- VARfit(dat_list[[i]]$flu, p = var_confound$p, exogen = dat_list[[i]]$ah)
  ar_v2_confound <- VARfit(dat_list[[i]]$rsv, p = var_confound$p, exogen = dat_list[[i]]$ah)
  sink()
  
  logRSS_v1 <- log(sum(ar_v1_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var_confound)[, 'flu'] ** 2))
  logRSS_v2 <- log(sum(ar_v2_confound$resid ** 2, na.rm = TRUE) / sum(residuals(var_confound)[, 'rsv'] ** 2))
  
  # Get p-values:
  p1_ftest <- causality(var_confound, cause = 'rsv')$Granger$p.value
  p2_ftest <- causality(var_confound, cause = 'flu')$Granger$p.value
  
  # Store results:
  gc_list[[i]] <- bind_cols(logRSS = c(logRSS_v1, logRSS_v2),
                            direction = c('rsv -> flu', 'flu -> rsv'),
                            ftest_p = c(p1_ftest, p2_ftest),
                            adf_p = c(adf_v1$p.value, adf_v2$p.value),
                            kpss_p = c(kpss_v1$p.value, kpss_v2$p.value))
  
}

#---- Check for stationary data ----
for (i in 1:length(gc_list)) {
  gc_list[[i]] %>% filter(adf_p >= 0.05) %>% print()
  gc_list[[i]] %>% filter(kpss_p < 0.05) %>% print()
}
# Note: By ADF, HK data (flu only) not stationary, but by KPSS fine

#---- Limit to results columns ----
gc_list <- lapply(gc_list, function(ix) {
  ix %>% dplyr::select(logRSS:ftest_p)
})

#---- Clean up ----
rm(df, lag_v1, lag_v2, p, var_confound, ar_v1_confound, ar_v2_confound, logRSS_v1, logRSS_v2,
   p1_ftest, p2_ftest, adf_v1, adf_v2, kpss_v1, kpss_v2, i)

# ---------------------------------------------------------------------------------------------------------------------

# Transfer entropy

#---- Create list of results ----
te_list <- vector('list', length(dat_list))

#---- Run TE ----
.jinit()
.jaddClassPath('infodynamics-dist-1.6.1/infodynamics.jar')

for (i in 1:length(te_list)) {
  
  # Get relevant data:
  sourceArray <- dat_list[[i]]$flu
  destArray <- dat_list[[i]]$rsv
  condArray <- dat_list[[i]]$ah
  
  # Create TE calculator:
  teCalc <- .jnew('infodynamics/measures/continuous/kraskov/ConditionalTransferEntropyCalculatorKraskov')
  .jcall(teCalc, 'V', 'setProperty', 'k', '4') # use Kraskov parameter k = 4 for nearest 4 points
  
  #---- Loop through different lags ----
  res_list_temp <- vector('list', length = 4)
  lags <- c('1', '2', '4', '13')
  
  for (j in 1:length(lags)) {
    
    lag <- lags[j]
    .jcall(teCalc, 'V', 'setProperty', 'delay', lag) # lag between source and destination
    
    #---- Calculation for flu -> RSV ----
    
    # Set relevant variables:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray, condArray)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '5') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
    .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '1') # set confounding embedding
    
    # Run analysis:
    result_confound_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
    
    nullDist_confound_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                        'computeSignificance', 500L)
    p_value_confound_v2_x_v1 <- nullDist_confound_v2_x_v1$pValue
    
    #---- Calculation for RSV -> flu ----
    
    # Set relevant variables:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray, condArray)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
    .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '1') # set confounding embedding
    
    # Run analysis:
    result_confound_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
    
    nullDist_confound_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                        'computeSignificance', 500L)
    p_value_confound_v1_x_v2 <- nullDist_confound_v1_x_v2$pValue
    
    # Store results (TEMP):
    res_list_temp[[j]] <- bind_cols(te = c(result_confound_v1_x_v2, result_confound_v2_x_v1),
                                    direction = c('v2 -> v1', 'v1 -> v2'),
                                    p_value = c(p_value_confound_v1_x_v2, p_value_confound_v2_x_v1),
                                    lag = lag)
    
  }
  
  # Store results (keep only lag with highest TE):
  te_list[[i]] <- bind_rows(res_list_temp) %>%
    group_by(direction) %>%
    filter(te == max(te))
  
}

#---- Clean up ----
rm(sourceArray, destArray, condArray, teCalc, result_confound_v2_x_v1, result_confound_v1_x_v2,
   nullDist_confound_v2_x_v1, nullDist_confound_v1_x_v2, p_value_confound_v2_x_v1, p_value_confound_v1_x_v2,
   lags, lag, res_list_temp, i, j)

# ---------------------------------------------------------------------------------------------------------------------

# CCM

#---- Create list of results ----
ccm_list <- vector('list', length(dat_list))

#---- Run CCM ----

for (i in 1:length(ccm_list)) {
  
  #---- Setup ----
  
  # Get data and add 'time' column:
  dat_temp <- dat_list[[i]] %>%
    mutate(time = 1:nrow(dat_list[[i]])) %>%
    dplyr::select(time, flu:rsv)
  
  # Get library and prediction sets:
  lib <- paste('1', nrow(dat_temp), sep = ' ')
  pred <- paste('1', nrow(dat_temp), sep = ' ')
  
  # Determine embedding dimensions:
  E_v1 <- EmbedDimension(dataFrame = dat_temp, columns = 'flu', target = 'flu',
                         lib = lib, pred = pred, maxE = 2, showPlot = FALSE)
  E_v1 <- E_v1 %>% slice_max(rho) %>% pull(E) # keep the row with max prediction skill
  
  E_v2 <- EmbedDimension(dataFrame = dat_temp, columns = "rsv", target = "rsv",
                         lib = lib, pred = pred, maxE = 5, showPlot = FALSE)
  E_v2 <- E_v2 %>% slice_max(rho) %>% pull(E)
  
  # Determine time delay:
  vars <- c('flu', 'rsv')
  params <- expand.grid(lib_column = vars, target_column = vars, tp = -20:0) %>%
    filter(lib_column != target_column) %>%
    mutate(E = if_else(lib_column == 'flu', E_v1, E_v2))
  
  lib_size_tp <- nrow(dat_temp) - (20 - 1) - (max(params$E) - 1)
  
  output <- do.call(rbind, lapply(seq_len(nrow(params)), function(i) {
    CCM(dataFrame = dat_temp, E = params$E[i], libSizes = lib_size_tp, random = FALSE,
        columns = as.character(params$lib_column[i]), target = as.character(params$target_column[i]), 
        Tp = params$tp[i], verbose = FALSE) %>%
      bind_cols(params[i, ])
  }))
  
  optimal_tp_v1xv2 <- output %>% filter(target_column == 'rsv') %>% filter(`flu:rsv` == max(`flu:rsv`)) %>% pull(tp)
  optimal_tp_v2xv1 <- output %>% filter(target_column == 'flu') %>% filter(`rsv:flu` == max(`rsv:flu`)) %>% pull(tp)
  
  highest_crossmap_cor_v1xv2 <- output %>% filter(target_column == 'rsv') %>% pull(`flu:rsv`) %>% max()
  highest_crossmap_cor_v2xv1 <- output %>% filter(target_column == 'flu') %>% pull(`rsv:flu`) %>% max()
  
  expect_true(highest_crossmap_cor_v1xv2 > 0)
  expect_true(highest_crossmap_cor_v2xv1 > 0)
  
  #---- Run CCM ----
  
  # Get max library size:
  lib_max <- nrow(dat_temp) - (max(abs(optimal_tp_v1xv2), abs(optimal_tp_v2xv1)) - 1) - (max(abs(E_v1), abs(E_v2)) - 1)
  if (optimal_tp_v1xv2 == 0 & optimal_tp_v2xv1 == 0) lib_max <- lib_max - 1
  
  # Run:
  v1_xmap_v2 <- CCM(dataFrame = dat_temp, E = E_v1, columns = 'flu', target = 'rsv', 
                    libSizes = c(E_v1 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v1xv2, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  v2_xmap_v1 <- CCM(dataFrame = dat_temp, E = E_v2, columns = 'rsv', target = 'flu', 
                    libSizes = c(E_v2 + 2, seq(30, 100, by = 10), seq(125, lib_max - 1, 25), lib_max), Tp = optimal_tp_v2xv1, random = TRUE,
                    sample = 100, includeData = TRUE, showPlot = FALSE)
  
  mean_rho_v1_xmap_v2 <- v1_xmap_v2$LibMeans %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  mean_rho_v2_xmap_v1 <- v2_xmap_v1$LibMeans %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  
  res <- mean_rho_v1_xmap_v2 %>%
    dplyr::select(1:2) %>%
    inner_join(mean_rho_v2_xmap_v1 %>%
                 dplyr::select(1:2),
               by = 'LibSize') %>%
    rename('mean_v1xv2' = 'flu:rsv',
           'mean_v2xv1' = 'rsv:flu')
  
  # Get all predictions:
  all_predictions_v1 <- v1_xmap_v2$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  all_predictions_v2 <- v2_xmap_v1$CCM1_PredictStat %>%
    mutate(LibSize = if_else(LibSize < 20, 12, LibSize))
  
  # Get medians and CIs for all library sizes:
  intervals_perc_v1 <- all_predictions_v1 %>%
    group_by(LibSize) %>%
    summarise(rho_median = quantile(rho, probs = 0.5),
              rho_ci_lower = quantile(rho, probs = 0.025),
              rho_ci_upper = quantile(rho, probs = 0.975)) %>%
    mutate(direction = 'v2 -> v1')
  intervals_perc_v2 <- all_predictions_v2 %>%
    group_by(LibSize) %>%
    summarise(rho_median = quantile(rho, probs = 0.5),
              rho_ci_lower = quantile(rho, probs = 0.025),
              rho_ci_upper = quantile(rho, probs = 0.975)) %>%
    mutate(direction = 'v1 -> v2')
  
  #---- Evaluate significance ----
  
  # Method 1: Check for convergence
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
  
  # Method 2: Seasonal surrogates
  num_surr <- 100
  
  surr_v1 <- SurrogateData(dat_temp$flu, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  surr_v2 <- SurrogateData(dat_temp$rsv, method = "seasonal", num_surr = num_surr, T_period = 52.25, alpha = 0)
  
  surr_dat_list_v1xv2 <- vector('list', length = num_surr)
  surr_dat_list_v2xv1 <- vector('list', length = num_surr)
  
  for (j in 1:num_surr) {
    
    surr_dat_list_v1xv2[[i]] <- dat_temp %>%
      dplyr::select(time, flu) %>%
      mutate(rsv = surr_v2[, i])
    
    surr_dat_list_v2xv1[[i]] <- dat_temp %>%
      dplyr::select(time, rsv) %>%
      mutate(flu = surr_v1[, i])
    
  }
  
  lib_max_use <- max(intervals_perc_v1$LibSize)
  
  cl <- makeCluster(10, type = 'SOCK')
  registerDoSNOW(cl)
  
  surr_res <- foreach(j = 1:num_surr, .packages = c('rEDM', 'tidyverse')) %dopar% {
    
    ccm_v1xv2 <- CCM(dataFrame = surr_dat_list_v1xv2[[i]], E = E_v1, columns = 'flu', target = 'rsv', 
                     libSizes = lib_max_use, Tp = optimal_tp_v1xv2, random = TRUE, sample = 100, includeData = TRUE)
    ccm_v2xv1 <- CCM(dataFrame = surr_dat_list_v2xv1[[i]], E = E_v2, columns = 'rsv', target = 'flu', 
                     libSizes = lib_max_use, Tp = optimal_tp_v2xv1, random = TRUE, sample = 100, includeData = TRUE)
    
    res_temp <- ccm_v1xv2$LibMeans %>%
      dplyr::select(1:2) %>%
      rename('rho' = 'flu:rsv') %>%
      mutate(direction = 'v2 -> v1') %>%
      bind_rows(ccm_v2xv1$LibMeans %>%
                  dplyr::select(1:2) %>%
                  rename('rho' = 'rsv:flu') %>%
                  mutate(direction = 'v1 -> v2'))
    
    return(res_temp)
    
  }
  
  stopCluster(cl)
  rm(ccm_v1xv2, ccm_v2xv1, res_temp, cl, j)
  
  surr_res <- bind_rows(surr_res)
  
  rho_v1xv2_surr <- surr_res %>% filter(direction == 'v2 -> v1') %>% pull(rho)
  rho_v2xv1_surr <- surr_res %>% filter(direction == 'v1 -> v2') %>% pull(rho)
  
  rho_v1xv2 <- res %>% filter(LibSize == lib_max_use) %>% pull(mean_v1xv2)
  rho_v2xv1 <- res %>% filter(LibSize == lib_max_use) %>% pull(mean_v2xv1)
  
  p_surr_v1xv2 <- (sum(rho_v1xv2 < rho_v1xv2_surr) + 1) / (length(rho_v1xv2_surr) + 1)
  p_surr_v2xv1 <- (sum(rho_v2xv1 < rho_v2xv1_surr) + 1) / (length(rho_v2xv1_surr) + 1)
  
  #---- Store results ----
  ccm_list[[i]] <- bind_rows(intervals_perc_v1, intervals_perc_v2) %>%
    group_by(direction) %>%
    summarise(rho_max = rho_median[LibSize == max(LibSize)]) %>%
    ungroup() %>%
    mutate(E = if_else(direction == 'v2 -> v1', E_v1, E_v2),
           tp_use = if_else(direction == 'v2 -> v1', optimal_tp_v1xv2, optimal_tp_v2xv1),
           p_conv = if_else(direction == 'v2 -> v1', p_conv_v1xv2, p_conv_v2xv1),
           p_surr = if_else(direction == 'v2 -> v1', p_surr_v1xv2, p_surr_v2xv1))
  
}

#---- Clean up ----
rm(lib, pred, E_v1, E_v2, vars, params, lib_size_tp, output, optimal_tp_v1xv2, optimal_tp_v2xv1, highest_crossmap_cor_v1xv2, highest_crossmap_cor_v2xv1,
   lib_max, v1_xmap_v2, v2_xmap_v1, mean_rho_v2_xmap_v1, mean_rho_v1_xmap_v2, all_predictions_v1, all_predictions_v2, intervals_perc_v1, intervals_perc_v2,
   corrs_Lmin_v1, corrs_Lmax_v1, corrs_Lmin_v2, corrs_Lmax_v2, inverse_quantile, p_conv_v1xv2, p_conv_v2xv1, num_surr, surr_v1, surr_v2, surr_dat_list_v1xv2,
   surr_dat_list_v2xv1, lib_max_use, surr_res, rho_v1xv2, rho_v2xv1, rho_v1xv2_surr, rho_v2xv1_surr, p_surr_v1xv2, p_surr_v2xv1, res, dat_temp, i)

# ---------------------------------------------------------------------------------------------------------------------

# Visualize results

p6a <- bind_rows(cor_list) %>%
  mutate(loc = c('Hong Kong', 'Canada')) %>%
  mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
  ggplot(aes(x = loc, y = cor, ymin = CI_lower_95, ymax = CI_upper_95, color = loc)) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_point(size = 4, shape = 17) +
  geom_pointrange(size = 0.8, lwd = 1, shape = 17) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(limits = c(-0.4, 1), breaks = seq(-0.4, 1, by = 0.2)) +
  scale_color_manual(values = c('#762a83', '#1b7837')) +
  labs(x = NULL, y = expression(r[Pearson]), tag = 'A')
p6b <- bind_rows(gam_list) %>%
  mutate(loc = c('Hong Kong', 'Canada')) %>%
  mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
  ggplot(aes(x = loc, y = cor_median, ymin = CI_lower95, ymax = CI_upper95, color = loc)) +
  geom_hline(yintercept = 0, lty = 2) +
  # geom_point(size = 4, shape = 17) +
  geom_pointrange(size = 0.8, lwd = 1, shape = 17) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_y_continuous(limits = c(-0.3, 0.02), breaks = c(-0.3, -0.2, -0.1, 0)) +
  scale_color_manual(values = c('#762a83', '#1b7837')) +
  labs(x = NULL, y = expression(r[GAM]), tag = 'B')
p6c <- bind_rows(gc_list) %>%
  mutate(x_pos = c(1.15 + c(0.075, -0.075), 1.85 + c(0.075, -0.075))) %>%
  # mutate(x_pos = c(1.2, 1.1, 1.9, 1.8)) %>%
  mutate(loc = c('Hong Kong', 'Hong Kong', 'Canada', 'Canada')) %>%
  mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
  mutate(sig = if_else(ftest_p < 0.05, 'yes', 'no')) %>%
  mutate(sig = paste(direction, sig)) %>%
  mutate(sig = factor(sig, levels = c('flu -> rsv yes', 'rsv -> flu yes', 'flu -> rsv no', 'rsv -> flu no'))) %>%
  ggplot(aes(x = x_pos, y = logRSS, col = loc, shape = sig)) +
  geom_hline(yintercept = 0, lty = 2) +
  geom_point(size = 3.5, stroke = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1.15, 1.85), labels = c('Hong Kong', 'Canada'), limits = c(0.8, 2.2)) +
  scale_y_continuous(limits = c(0, 0.042), n.breaks = 5) +
  scale_color_manual(values = c('#762a83', '#1b7837')) +
  scale_shape_manual(values = c(16, 0, 15, 1)) +
  labs(x = NULL, y = expression(G[y %->% x]), tag = 'C')
p6d <- bind_rows(te_list) %>%
  ungroup() %>%
  mutate(x_pos = c(1.15 + c(0.075, -0.075), 1.85 + c(0.075, -0.075))) %>%
  mutate(loc = c('Hong Kong', 'Hong Kong', 'Canada', 'Canada')) %>%
  mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
  mutate(sig = if_else(p_value < 0.05, 'yes', 'no')) %>%
  mutate(sig = paste(direction, sig)) %>%
  mutate(sig = factor(sig, levels = c('v1 -> v2 yes', 'v2 -> v1 yes', 'v1 -> v2 no', 'v2 -> v1 no'))) %>%
  ggplot(aes(x = x_pos, y = te, col = loc, shape = sig)) +
  geom_hline(yintercept = 0, lty = 2, col = 'gray70') +
  geom_point(size = 3.5, stroke = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1.15, 1.85), labels = c('Hong Kong', 'Canada'), limits = c(0.8, 2.2)) +
  scale_y_continuous(limits = c(-0.001, 0.08), n.breaks = 5) +
  scale_color_manual(values = c('#762a83', '#1b7837')) +
  scale_shape_manual(values = c(16, 15, 0, 1)) +
  labs(x = NULL, y = expression(T[y %->% x]), tag = 'D')
p6e <- bind_rows(ccm_list) %>%
  mutate(x_pos1 = c(1.00, 1.22, 1.70, 1.92),
         x_pos2 = c(1.08, 1.30, 1.78, 2.00)) %>%
  mutate(loc = c('Hong Kong', 'Hong Kong', 'Canada', 'Canada')) %>%
  mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
  mutate(sig1 = if_else(p_conv < 0.05, 'yes', 'no'),
         sig2 = if_else(p_surr < 0.05, 'yes', 'no')) %>%
  mutate(sig1 = paste(direction, sig1),
         sig2 = paste(direction, sig2)) %>%
  mutate(sig1 = factor(sig1, levels = c('v1 -> v2 yes', 'v2 -> v1 yes', 'v1 -> v2 no', 'v2 -> v1 no')),
         sig2 = factor(sig2, levels = c('v1 -> v2 yes', 'v2 -> v1 yes', 'v1 -> v2 no', 'v2 -> v1 no'))) %>%
  ggplot(aes(y = rho_max, col = loc)) +
  geom_point(aes(x_pos1, shape = sig1), size = 3.5, stroke = 1) +
  geom_point(aes(x_pos2, shape = sig2), size = 3.5, stroke = 1) +
  theme_classic() +
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        plot.tag = element_text(size = 24),
        plot.tag.position = c(0.008, 0.975),
        panel.grid.minor = element_blank(),
        legend.position = 'none') +
  scale_x_continuous(breaks = c(1.15, 1.85), labels = c('Hong Kong', 'Canada'), limits = c(0.8, 2.2)) +
  scale_y_continuous(limits = c(0.2, 0.9), n.breaks = 8) +
  scale_color_manual(values = c('#762a83', '#1b7837')) +
  scale_shape_manual(values = c(16, 1, 15, 0)) +
  labs(x = NULL, y = expression(rho), tag = 'E')

# bind_rows(ccm_list) %>%
#   mutate(x_pos = c(0.95, 1.05, 1.95, 2.05)) %>%
#   mutate(loc = c('Hong Kong', 'Hong Kong', 'Canada', 'Canada')) %>%
#   mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
#   mutate(sig = if_else(p_conv < 0.05, 'yes', 'no')) %>%
#   mutate(sig = paste(direction, sig)) %>%
#   mutate(sig = factor(sig, levels = c('v1 -> v2 yes', 'v2 -> v1 yes', 'v1 -> v2 no', 'v2 -> v1 no'))) %>%
#   ggplot(aes(x = x_pos, y = rho_max, col = loc, shape = sig)) +
#   geom_point(size = 3.5, stroke = 1) +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.tag = element_text(size = 24),
#         plot.tag.position = c(0.008, 0.975),
#         panel.grid.minor = element_blank(),
#         legend.position = 'none') +
#   scale_x_continuous(breaks = c(1, 2), labels = c('Hong Kong', 'Canada'), limits = c(0.8, 2.2)) +
#   scale_y_continuous(limits = c(0.2, 0.9), n.breaks = 8) +
#   scale_color_manual(values = c('#762a83', '#1b7837')) +
#   scale_shape_manual(values = c(16, 1, 0, 15)) +
#   labs(x = NULL, y = expression(rho), tag = 'E')
# bind_rows(ccm_list) %>%
#   mutate(x_pos = c(0.95, 1.05, 1.95, 2.05)) %>%
#   mutate(loc = c('Hong Kong', 'Hong Kong', 'Canada', 'Canada')) %>%
#   mutate(loc = factor(loc, levels = c('Hong Kong', 'Canada'))) %>%
#   mutate(sig = if_else(p_surr < 0.05, 'yes', 'no')) %>%
#   mutate(sig = paste(direction, sig)) %>%
#   mutate(sig = factor(sig, levels = c('v1 -> v2 yes', 'v2 -> v1 yes', 'v1 -> v2 no', 'v2 -> v1 no'))) %>%
#   ggplot(aes(x = x_pos, y = rho_max, col = loc, shape = sig)) +
#   geom_point(size = 3.5, stroke = 1) +
#   theme_classic() +
#   theme(axis.title = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         plot.tag = element_text(size = 24),
#         plot.tag.position = c(0.008, 0.975),
#         panel.grid.minor = element_blank(),
#         legend.position = 'none') +
#   scale_x_continuous(breaks = c(1, 2), labels = c('Hong Kong', 'Canada'), limits = c(0.8, 2.2)) +
#   scale_y_continuous(limits = c(0.2, 0.9), n.breaks = 8) +
#   scale_color_manual(values = c('#762a83', '#1b7837')) +
#   scale_shape_manual(values = c(16, 15, 1, 0)) +
#   labs(x = NULL, y = expression(rho), tag = 'F')

p6 <- arrangeGrob(arrangeGrob(p6a, p6b, p6c, nrow = 1, widths = c(0.7, 0.7, 1)),
                  arrangeGrob(p6d, p6e, nrow = 1, widths = c(1, 1.4)))
p6.alt <- arrangeGrob(p6a, p6b, p6c, p6d, p6e, nrow = 1, widths = c(0.6, 0.6, 0.8, 0.8, 1.1))

plot(p6)
plot(p6.alt)

ggsave('results/plots/figures/Figure6.svg', p6, width = 10, height = 6)
