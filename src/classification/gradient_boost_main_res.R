# ------------------------------------------------------------------------------
# Gradient boosting to compile results from all statistical methods
#
# Created by: Sarah Kramer
# Creation date: 09 September 2024
# ------------------------------------------------------------------------------

# Setup

# Load packages
library(tidyverse)
library(gridExtra)
library(testthat)
library(viridis)
library(xgboost)

# Load functions:
source('src/functions_etc/fxns_process_results.R')

# ------------------------------------------------------------------------------

# Read in all results
source('src/functions_etc/load_main_results.R')

# ------------------------------------------------------------------------------

# Further formatting of results

# Get combined tibble of results and significance levels for all methods:
res_corr <- res_corr %>%
  select(run:.id, theta_lambda:int_true, cor, int_est) %>%
  rename('cor_val' = cor, 'cor_sig' = int_est) %>%
  mutate(cor_sig = case_match(cor_sig, 'none' ~ 0, 'neg' ~ 1, 'pos' ~ 2))

res_gam <- res_gam %>%
  select(run:.id, cor, int_est, cor_confound, int_est_confound) %>%
  rename('gam_val' = cor, 'gam_sig' = int_est, 'gam_confound_val' = cor_confound, 'gam_confound_sig' = int_est_confound) %>%
  mutate(gam_sig = case_match(gam_sig, 'none' ~ 0, 'neg' ~ 1, 'pos' ~ 2),
         gam_confound_sig = case_match(gam_confound_sig, 'none' ~ 0, 'neg' ~ 1, 'pos' ~ 2))

res_granger <- res_granger %>%
  select(run:.id, direction:confounding, logRSS) %>%
  pivot_wider(names_from = confounding, values_from = logRSS) %>%
  rename('granger_val' = none, 'granger_confound_val' = seasonal) %>%
  inner_join(res_granger %>%
               select(run:.id, direction:confounding, int_est) %>%
               pivot_wider(names_from = confounding, values_from = int_est) %>%
               rename('granger_sig' = none, 'granger_confound_sig' = seasonal),
             by = c('run', '.id', 'direction')) %>%
  select(run:granger_val, granger_sig, granger_confound_val, granger_confound_sig) %>%
  mutate(.id = as.integer(.id),
         granger_sig = case_match(granger_sig, 'none' ~ 0, 'interaction' ~ 1),
         granger_confound_sig = case_match(granger_confound_sig, 'none' ~ 0, 'interaction' ~ 1))

res_te <- res_te %>%
  select(run:.id, direction, lag, te, int_est) %>%
  filter((str_detect(best_v1xv2, direction) & str_detect(best_v1xv2, paste0('lag ', lag))) |
           (str_detect(best_v2xv1, direction) & str_detect(best_v2xv1, paste0('lag ', lag)))) %>%
  select(-lag) %>%
  rename('te_val' = te, 'te_sig' = int_est) %>%
  mutate(.id = as.integer(.id),
         te_sig = case_match(te_sig, 'none' ~ 0, 'interaction' ~ 1))
rm(best_v1xv2, best_v2xv1)

res_ccm <- res_ccm %>%
  select(run:direction, rho_mean, int_est_1) %>%
  rename('ccm_val' = rho_mean, 'ccm_sig' = int_est_1) %>%
  mutate(.id = as.integer(.id),
         ccm_sig = case_match(ccm_sig, 'none' ~ 0, 'interaction' ~ 1))

# Separate out information for determining effect of V1->V2 vs. V2->V1:
res_v1tov2 <- res_corr %>%
  inner_join(res_gam, by = c('run', '.id')) %>%
  inner_join(res_granger %>% filter(direction == 'v1 -> v2') %>% select(-direction), by = c('run', '.id')) %>%
  inner_join(res_te %>% filter(direction == 'v1 -> v2') %>% select(-direction), by = c('run', '.id')) %>%
  inner_join(res_ccm %>% filter(direction == 'v1 -> v2') %>% select(-direction), by = c('run', '.id'))

res_v2tov1 <- res_corr %>%
  inner_join(res_gam, by = c('run', '.id')) %>%
  inner_join(res_granger %>% filter(direction == 'v2 -> v1') %>% select(-direction), by = c('run', '.id')) %>%
  inner_join(res_te %>% filter(direction == 'v2 -> v1') %>% select(-direction), by = c('run', '.id')) %>%
  inner_join(res_ccm %>% filter(direction == 'v2 -> v1') %>% select(-direction), by = c('run', '.id'))

rm(res_corr, res_gam, res_granger, res_te, res_ccm)

# Get also tibble including all results (so, for interactions in both directions):
res_all <- res_v1tov2 %>%
  rename_with(~ paste0(.x, '_v1xv2'), granger_val:ccm_sig) %>%
  inner_join(res_v2tov1 %>%
               select(run:.id, granger_val:ccm_sig) %>%
               rename_with(~ paste0(.x, '_v2xv1'), granger_val:ccm_sig),
             by = c('run', '.id'))

# Split all tibbles into training and testing sets:
set.seed(3857)
train_ind <- sample(nrow(res_all), size = round(0.7 * nrow(res_all)))

res_all_train <- res_all[train_ind, ]
res_v1tov2_train <- res_v1tov2[train_ind, ]
res_v2tov1_train <- res_v2tov1[train_ind, ]

res_all_test <- res_all[-train_ind, ]
res_v1tov2_test <- res_v1tov2[-train_ind, ]
res_v2tov1_test <- res_v2tov1[-train_ind, ]

rm(train_ind)

# Get features and targets for prediction:
x_all_train <- res_all_train %>% select(cor_val:ccm_sig_v2xv1)
x_all_test <- res_all_test %>% select(cor_val:ccm_sig_v2xv1)

x_v1tov2_train <- res_v1tov2_train %>% select(cor_val:ccm_sig)
x_v1tov2_test <- res_v1tov2_test %>% select(cor_val:ccm_sig)

x_v2tov1_train <- res_v2tov1_train %>% select(cor_val:ccm_sig)
x_v2tov1_test <- res_v2tov1_test %>% select(cor_val:ccm_sig)

y_train <- res_all_train %>% select(int_true) %>% mutate(int_true = case_match(int_true, 'none' ~ 0, 'neg' ~ 1, 'pos' ~ 2)) %>% pull(int_true)
y_test <- res_all_test %>% select(int_true) %>% mutate(int_true = case_match(int_true, 'none' ~ 0, 'neg' ~ 1, 'pos' ~ 2)) %>% pull(int_true)

# Weight training data to deal with unbalanced data:
training_weights <- c(1.0, 0.1056662, 0.1733668)[y_train + 1]

# Perform gradient boosting (first, on all metrics in both directions):
# Source: https://www.r-bloggers.com/2021/02/machine-learning-with-r-a-complete-guide-to-gradient-boosting-and-xgboost/

xgb_train <- xgb.DMatrix(data = as.matrix(x_all_train), label = y_train, weight = training_weights)
xgb_test <- xgb.DMatrix(data = as.matrix(x_all_test), label = y_test)

xgb_params <- list(
  booster = 'gbtree',
  eta = 0.01, # learning rate
  max_depth = 4, # maximum tree depth
  gamma = 4, # minimum loss reduction to make further partition on a leaf node
  subsample = 0.75, # proportion of data randomly selected and used to grow trees
  colsample_bytree = 1, # proportion of columns used when constructing each tree
  objective = 'multi:softprob',
  eval_metric = 'mlogloss',
  num_class = length(unique(y_train))
)

set.seed(5829)
model_all <- xgb.train(params = xgb_params, data = xgb_train, nrounds = 2500, verbose = 1)

# # Check accuracy on training data itself:
# # If this improves with a greater number of rounds, but testing accuracy
# # does not, there may be evidence of overfitting
# res_all_pred_train <- res_all_train %>%
#   select(run:int_true) %>%
#   mutate(int_pred = predict(model_all, as.matrix(x_all_train), reshape = TRUE) %>%
#            apply(1, which.max),
#          int_pred = int_pred - 1) %>%
#   mutate(int_pred = case_match(int_pred, 0 ~ 'none', 1 ~ 'neg', 2 ~ 'pos'))
# 
# print('Overall accuracy:')
# print((res_all_pred_train %>% filter(int_true == int_pred) %>% nrow()) / nrow(res_all_pred_train))
# 
# print('Sensitivity (Correct Direction) (Overall):')
# print((res_all_pred_train %>% filter(int_true == 'neg' & int_pred == 'neg') %>% nrow()) / (res_all_pred_train %>% filter(int_true == 'neg') %>% nrow()))
# print((res_all_pred_train %>% filter(int_true == 'pos' & int_pred == 'pos') %>% nrow()) / (res_all_pred_train %>% filter(int_true == 'pos') %>% nrow()))
# 
# print('Specificity (Overall):')
# print((res_all_pred_train %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()) / (res_all_pred_train %>% filter(int_true == 'none') %>% nrow()))
# rm(res_all_pred_train)

# Now check accuracy on test data:
res_all_pred <- res_all_test %>%
  select(run:int_true) %>%
  mutate(int_pred = predict(model_all, as.matrix(x_all_test), reshape = TRUE) %>%
           apply(1, which.max),
         int_pred = int_pred - 1) %>%
  mutate(int_pred = case_match(int_pred, 0 ~ 'none', 1 ~ 'neg', 2 ~ 'pos'))

# Calculate overall accuracy:
print('Overall accuracy:')
print((res_all_pred %>% filter(int_true == int_pred) %>% nrow()) / nrow(res_all_pred))

# Calculate sensitivity/specificity overall:
print('Sensitivity (Any Interaction) (Overall):')
print((res_all_pred %>% filter(int_true != 'none' & int_pred != 'none') %>% nrow()) / (res_all_pred %>% filter(int_true != 'none') %>% nrow()))

print('Sensitivity (Correct Direction) (Overall):')
sens_pos <- (res_all_pred %>% filter(int_true == 'neg' & int_pred == 'neg') %>% nrow()) / (res_all_pred %>% filter(int_true == 'neg') %>% nrow())
sens_neg <- (res_all_pred %>% filter(int_true == 'pos' & int_pred == 'pos') %>% nrow()) / (res_all_pred %>% filter(int_true == 'pos') %>% nrow())
print(sens_pos)
print(sens_neg)

print('Specificity (Overall):')
spec <- (res_all_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()) / (res_all_pred %>% filter(int_true == 'none') %>% nrow())
print(spec)

# Calculate accuracy weighted by prevalence of true interaction values:
weight_pos <- min(table(res_all_test$int_true)) / (res_all_test %>% filter(int_true == 'pos') %>% nrow())
weight_neg <- min(table(res_all_test$int_true)) / (res_all_test %>% filter(int_true == 'neg') %>% nrow())
weight_null <- min(table(res_all_test$int_true)) / (res_all_test %>% filter(int_true == 'none') %>% nrow())

acc_weighted_gb <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)

print('Overall accuracy (weighted):')
print(acc_weighted_gb)

# Calculate Matthews correlation coefficient (MCC):
tp <- res_all_pred %>% filter(int_true != 'none' & int_pred != 'none') %>% nrow()
tn <- res_all_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()
fp <- res_all_pred %>% filter(int_true == 'none' & int_pred != 'none') %>% nrow()
fn <- res_all_pred %>% filter(int_true != 'none' & int_pred == 'none') %>% nrow()

print('MCC:')
print(mcc(tp, tn, fp, fn))

# Calculate accuracy by true param values:
calculate_accuracy_matrix <- function(df) {
  
  all_tl <- sort(unique(df$theta_lambda))
  all_delta <- sort(unique(1 / df$delta))
  
  mat_correct <- matrix(nrow = length(all_tl), ncol = length(all_delta))
  rownames(mat_correct) <- all_tl
  colnames(mat_correct) <- all_delta
  
  for (tl in all_tl) {
    
    if (tl != 1) {
      for (d in all_delta) {
        
        df_temp <- df %>% filter(theta_lambda == tl & delta == 1 / d)
        mat_correct[rownames(mat_correct) == tl, colnames(mat_correct) == d] <- (df_temp %>% filter(int_true == int_pred) %>% nrow()) / nrow(df_temp)
        
      }
    } else {
      
      df_temp <- df %>% filter(theta_lambda == tl)
      mat_correct[rownames(mat_correct) == tl, ] <- (df_temp %>% filter(int_true == int_pred) %>% nrow()) / nrow(df_temp)
      
    }
    
  }
  
  mat_correct <- mat_correct %>%
    as_tibble(rownames = 'strength') %>%
    pivot_longer(-strength, names_to = 'duration', values_to = 'perc_correct') %>%
    mutate(duration = factor(duration, levels = c(1, 4, 13)))
  return(mat_correct)
  
}

acc_gb <- calculate_accuracy_matrix(res_all_pred)

p.combine.metrics.all <- ggplot(data = acc_gb, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Combine Methods')

# Perform gradient boosting (identify interactions in one direction):
xgb_train <- xgb.DMatrix(data = as.matrix(x_v1tov2_train), label = y_train, weight = training_weights)
xgb_test <- xgb.DMatrix(data = as.matrix(x_v1tov2_test), label = y_test)

set.seed(5829)
model_v1tov2 <- xgb.train(params = xgb_params, data = xgb_train, nrounds = 2500, verbose = 1)

xgb_train <- xgb.DMatrix(data = as.matrix(x_v2tov1_train), label = y_train, weight = training_weights)
xgb_test <- xgb.DMatrix(data = as.matrix(x_v2tov1_test), label = y_test)

set.seed(5829)
model_v2tov1 <- xgb.train(params = xgb_params, data = xgb_train, nrounds = 2500, verbose = 1)

res_v1tov2_pred <- res_v1tov2_test %>%
  select(run:int_true) %>%
  mutate(int_pred = predict(model_v1tov2, as.matrix(x_v1tov2_test), reshape = TRUE) %>%
           apply(1, which.max),
         int_pred = int_pred - 1) %>%
  mutate(int_pred = case_match(int_pred, 0 ~ 'none', 1 ~ 'neg', 2 ~ 'pos'))

res_v2tov1_pred <- res_v2tov1_test %>%
  select(run:int_true) %>%
  mutate(int_pred = predict(model_v2tov1, as.matrix(x_v2tov1_test), reshape = TRUE) %>%
           apply(1, which.max),
         int_pred = int_pred - 1) %>%
  mutate(int_pred = case_match(int_pred, 0 ~ 'none', 1 ~ 'neg', 2 ~ 'pos'))

print('V1 -> V2')
print('Sensitivity (Correct Direction) (Overall):')
sens_pos <- (res_v1tov2_pred %>% filter(int_true == 'neg' & int_pred == 'neg') %>% nrow()) / (res_v1tov2_pred %>% filter(int_true == 'neg') %>% nrow())
sens_neg <- (res_v1tov2_pred %>% filter(int_true == 'pos' & int_pred == 'pos') %>% nrow()) / (res_v1tov2_pred %>% filter(int_true == 'pos') %>% nrow())
print(sens_pos)
print(sens_neg)

print('Specificity (Overall):')
spec <- (res_v1tov2_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()) / (res_v1tov2_pred %>% filter(int_true == 'none') %>% nrow())
print(spec)

acc_weighted_gb_v1tov2 <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)

print('Overall accuracy (weighted):')
print(acc_weighted_gb_v1tov2)

tp <- res_v1tov2_pred %>% filter(int_true != 'none' & int_pred != 'none') %>% nrow()
tn <- res_v1tov2_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()
fp <- res_v1tov2_pred %>% filter(int_true == 'none' & int_pred != 'none') %>% nrow()
fn <- res_v1tov2_pred %>% filter(int_true != 'none' & int_pred == 'none') %>% nrow()

print('MCC:')
print(mcc(tp, tn, fp, fn))

print('V2 -> V1')
print('Sensitivity (Correct Direction) (Overall):')
sens_pos <- (res_v2tov1_pred %>% filter(int_true == 'neg' & int_pred == 'neg') %>% nrow()) / (res_v2tov1_pred %>% filter(int_true == 'neg') %>% nrow())
sens_neg <- (res_v2tov1_pred %>% filter(int_true == 'pos' & int_pred == 'pos') %>% nrow()) / (res_v2tov1_pred %>% filter(int_true == 'pos') %>% nrow())
print(sens_pos)
print(sens_neg)

print('Specificity (Overall):')
spec <- (res_v2tov1_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()) / (res_v2tov1_pred %>% filter(int_true == 'none') %>% nrow())
print(spec)

acc_weighted_gb_v2tov1 <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)

print('Overall accuracy (weighted):')
print(acc_weighted_gb_v2tov1)

tp <- res_v2tov1_pred %>% filter(int_true != 'none' & int_pred != 'none') %>% nrow()
tn <- res_v2tov1_pred %>% filter(int_true == 'none' & int_pred == 'none') %>% nrow()
fp <- res_v2tov1_pred %>% filter(int_true == 'none' & int_pred != 'none') %>% nrow()
fn <- res_v2tov1_pred %>% filter(int_true != 'none' & int_pred == 'none') %>% nrow()

print('MCC:')
print(mcc(tp, tn, fp, fn))

acc_gb_v1tov2 <- calculate_accuracy_matrix(res_v1tov2_pred)
acc_gb_v2tov1 <- calculate_accuracy_matrix(res_v2tov1_pred)

p.combine.metrics.v1tov2 <- ggplot(data = acc_gb_v1tov2, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Combine Methods (V1 -> V2)')
p.combine.metrics.v2tov1 <- ggplot(data = acc_gb_v2tov1, aes(x = strength, y = duration, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Combine Methods (V2 -> V1)')
grid.arrange(p.combine.metrics.all, p.combine.metrics.v1tov2, p.combine.metrics.v2tov1, nrow = 1)

# Explore which features most contribute to correct classification:
# https://xgboost.readthedocs.io/en/stable/R-package/xgboostPresentation.html
import_mat_v1tov2 <- xgb.importance(model = model_v1tov2)
import_mat_v2tov1 <- xgb.importance(model = model_v2tov1)

xgb.plot.importance(importance_matrix = import_mat_v1tov2)
xgb.plot.importance(importance_matrix = import_mat_v2tov1)

# Clean up:
rm(model_all, model_v1tov2, model_v2tov1, res_all, res_all_pred, res_all_train,
   res_all_test, res_v1tov2, res_v1tov2_train, res_v1tov2_test, res_v2tov1,
   res_v2tov1_train, res_v2tov1_test, res_v1tov2_pred, res_v2tov1_pred, x_all_train,
   x_all_test, x_v1tov2_train, x_v1tov2_test, x_v2tov1_train, x_v2tov1_test,
   xgb_params, training_weights, xgb_train, xgb_test, y_train, y_test,
   weight_pos, weight_neg, weight_null, sens_pos, sens_neg, spec, fn, fp, tn, tp,
   import_mat_v1tov2, import_mat_v2tov1)
