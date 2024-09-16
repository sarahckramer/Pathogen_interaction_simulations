# ------------------------------------------------------------------------------
# Code to use a convolutional neural network approach to classify pairs of time series
# Source: https://tensorflow.rstudio.com/examples/timeseries_classification_from_scratch
# 
# Created by: Sarah Kramer
# Creation date: 20.08.2024
# ------------------------------------------------------------------------------

# Setup

# Load libraries:
library(tensorflow)
library(keras3)
library(tidyverse)
library(gridExtra)
library(viridis)

# Set seed:
set.seed(18943)

# Set number of times to fit model:
num_runs <- 5

# ------------------------------------------------------------------------------

# Read in all results

# Get file names:
res_filenames_T <- list.files(path = 'results/', pattern = 'TRUE', full.names = TRUE)

# Read in results:
results_T <- vector('list', length = length(res_filenames_T))

for (i in 1:length(res_filenames_T)) {
  results_T[[i]] <- read_rds(res_filenames_T[i])
}
rm(i)

# Label each results list with run number:
where_run <- which(!is.na(as.numeric(str_split(res_filenames_T, '_')[[1]])))

names(results_T) <- unlist(map(str_split(res_filenames_T, '_'), where_run))

rm(res_filenames_T, where_run)

# ------------------------------------------------------------------------------

# Extract true parameter values and data

# Get true parameter values:
int_params <- c('theta_lambda1', 'delta1')

res_trueparams <- lapply(1:length(results_T), function(ix) {
  results_T[[ix]]$true_param[int_params, 1]
}) %>%
  bind_rows() %>%
  rename('theta_lambda' = 'theta_lambda1',
         'delta' = 'delta1') %>%
  mutate(run = as.numeric(names(results_T))) %>%
  arrange(run)

# Get data:
data_list <- vector('list', length = length(results_T))

for (i in 1:length(results_T)) {
  
  data_list[[i]] <- results_T[[i]]$data %>%
    select(time, date, .id, V1_obs, V2_obs) %>%
    mutate(run = as.numeric(names(results_T)[i]))
  
}

dat <- bind_rows(data_list) %>%
  arrange(run)
rm(data_list, i, int_params)

# Join true parameter values:
dat <- dat %>%
  as_tibble() %>%
  inner_join(res_trueparams, by = 'run')
rm(res_trueparams, results_T)

# Get true classification:
dat <- dat %>%
  mutate(class1 = if_else(theta_lambda == 1, 'null', 'int'),
         class2 = if_else(theta_lambda < 1, 'neg', class1),
         class2 = if_else(theta_lambda > 1, 'pos', class2)) %>%
  select(time, .id, run, V1_obs:V2_obs, theta_lambda:class2)

# ------------------------------------------------------------------------------

# Format data for input into model

# Z-normalize data?:
# For now, leave raw

# Set classes as integers:
dat <- dat %>%
  mutate(class1_int = if_else(class1 == 'int', 1, 0),
         class2_int = if_else(class2 == 'null', 0, 1),
         class2_int = if_else(class2 == 'pos', 2, class2_int))

# Reformat tibble:
# For now, use class2
datV1 <- dat %>%
  select(time:V1_obs, class2_int) %>%
  mutate(time = paste0('t_', time)) %>%
  pivot_wider(names_from = time, values_from = V1_obs)
datV2 <- dat %>%
  select(time:run, V2_obs, class2_int) %>%
  mutate(time = paste0('t_', time)) %>%
  pivot_wider(names_from = time, values_from = V2_obs)

# Split into training and testing sets (70/30% split):
set.seed(4820)
train_ind <- sample(nrow(datV1), size = round(0.7 * nrow(datV1)))

datV1_train <- datV1[train_ind, ]
datV2_train <- datV2[train_ind, ]

datV1_test <- datV1[-train_ind, ]
datV2_test <- datV2[-train_ind, ]

# Shuffle the training set:
# This is already done in the previous step

# Reformat data into matrices (training data):
datV1_train_mat <- as.matrix(datV1_train[, 4:ncol(datV1_train)])
datV2_train_mat <- as.matrix(datV2_train[, 4:ncol(datV2_train)])

x_train <- array(NA, dim = c(dim(datV1_train_mat), 2))
x_train[,, 1] <- datV1_train_mat
x_train[,, 2] <- datV2_train_mat
rm(datV1_train_mat, datV2_train_mat)

y_train <- as.matrix(datV1_train$class2_int)

# Reformat data into matrices (testing data):
datV1_test_mat <- as.matrix(datV1_test[, 4:ncol(datV1_test)])
datV2_test_mat <- as.matrix(datV2_test[, 4:ncol(datV2_test)])

x_test <- array(NA, dim = c(dim(datV1_test_mat), 2))
x_test[,, 1] <- datV1_test_mat
x_test[,, 2] <- datV2_test_mat
rm(datV1_test_mat, datV2_test_mat)

y_test <- as.matrix(datV1_test$class2_int)

# ------------------------------------------------------------------------------

# Fit model

# Get number of possible classes:
num_classes <- length(unique(y_train))

# Fit several times, to average over effect of random initializations
if (!file.exists('results/cnn/best_model_1.keras')) {
  for (i in 1:num_runs) {
    
    print(i)
    
    # Build a model
    
    # Create input layer:
    input_shape <- dim(x_train)[-1] # drop batch dim
    input_layer <- layer_input(input_shape)
    
    # Create output layer:
    output_layer <- input_layer %>%
      
      # First convolutional layer
      layer_conv_1d(16, 8, padding = 'same', kernel_regularizer = regularizer_l2(0.001)) %>%
      layer_batch_normalization() %>%
      layer_activation_relu() %>%
      layer_dropout(rate = 0.2) %>%
      
      # Second convolutional layer
      layer_conv_1d(16, 5, padding = 'same', kernel_regularizer = regularizer_l2(0.001)) %>%
      layer_batch_normalization() %>%
      layer_activation_relu() %>%
      layer_dropout(rate = 0.2) %>%
      
      # Third convolutional layer
      layer_conv_1d(16, 3, padding = 'same', kernel_regularizer = regularizer_l2(0.001)) %>%
      layer_batch_normalization() %>%
      layer_activation_relu() %>%
      layer_dropout(rate = 0.2) %>%
      
      layer_global_average_pooling_1d() %>%
      layer_dense(num_classes, activation = "softmax")
    
    model <- keras_model(input_layer, output_layer)
    
    # ------------------------------------------------------------------------------
    
    # Train
    
    # Set parameters for model fitting:
    epochs <- 500
    batch_size <- 32
    
    callbacks <- list(
      callback_model_checkpoint(paste0('results/cnn/best_model_', i, '.keras'), monitor = 'val_loss',
                                save_best_only = TRUE),
      callback_reduce_lr_on_plateau(monitor = 'val_loss', factor = 0.75,
                                    patience = 20, min_lr = 0.0001),
      callback_early_stopping(monitor = 'val_loss', patience = 50,
                              verbose = 1)
    )
    
    # Compile model:
    model %>% compile(
      optimizer = optimizer_adam(learning_rate = 0.0005),#'adam',
      loss = 'sparse_categorical_crossentropy',
      metrics = list('sparse_categorical_accuracy')
    )
    
    # Fit model to training data:
    history <- model %>%
      fit(x_train, y_train,
          batch_size = batch_size,
          epochs = epochs,
          callbacks = callbacks,
          validation_split = 0.2,
          verbose = 1,
          class_weight = list('0' = 1.0, '1' = 0.103, '2' = 0.159))
    
    # Plot training and validation loss:
    print(plot(history))
    
    # Check for signs of overfitting:
    p.overfit <- history %>%
      as_tibble() %>%
      filter(metric != 'learning_rate') %>%
      ggplot(aes(x = epoch, y = value, col = data)) +
      # geom_line() +
      geom_smooth(se = FALSE) +
      facet_wrap(~ metric, scales = 'free_y') +
      theme_classic() +
      scale_x_log10() +
      labs(title = paste0('Run ', i))
    print(p.overfit)
    
  }
  
  # Clean up:
  rm(i, input_shape, input_layer, output_layer, model, epochs, batch_size,
     callbacks, history, p.overfit, num_classes)
  
} else {
  
  # Clean up:
  rm(num_classes)
  
}

# ------------------------------------------------------------------------------

# Check model performance on test data

# Function to calculate accuracy by true interaction parameter values:
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
        mat_correct[rownames(mat_correct) == tl, colnames(mat_correct) == d] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
        
      }
    } else {
      
      df_temp <- df %>% filter(theta_lambda == tl)
      mat_correct[rownames(mat_correct) == tl, ] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
      
    }
    
  }
  
  mat_correct <- mat_correct %>%
    as_tibble(rownames = 'strength') %>%
    pivot_longer(-strength, names_to = 'duration', values_to = 'perc_correct') %>%
    mutate(duration = factor(duration, levels = c(1, 4, 13)))
  return(mat_correct)
  
}

# Loop through all models and calculate accuracy:
dat_pred_LIST <- vector('list', length = num_runs)
sens_pos_vec = sens_neg_vec = spec_vec = acc_weighted_vec = c()
for (i in 1:num_runs) {
  
  # Load fitted model:
  loaded_model <- load_model(paste0('results/cnn/best_model_', i, '.keras'))
  
  # Get loss/accuracy for test data:
  result <- loaded_model %>% evaluate(x_test, y_test)
  
  sprintf("Test loss: %s", result[["loss"]])
  sprintf("Test accuracy: %s", result[["sparse_categorical_accuracy"]])
  
  # Get predictions for test data:
  dat_pred <- datV1_test %>%
    select(.id:class2_int) %>%
    mutate(class2_pred = loaded_model %>%
             predict(x_test) %>%
             apply(1, which.max),
           class2_pred = class2_pred - 1)
  
  dat_pred <- dat %>%
    select(.id:run, theta_lambda:delta, class2) %>%
    unique() %>%
    inner_join(dat_pred,
               by = c('.id', 'run'))
  
  dat_pred_LIST[[i]] <- dat_pred
  
  # Calculate sensitivity/specificity overall:
  sens_pos <- (dat_pred %>% filter(class2 == 'pos' & class2_pred == 2) %>% nrow()) / (dat_pred %>% filter(class2 == 'pos') %>% nrow())
  sens_neg <- (dat_pred %>% filter(class2 == 'neg' & class2_pred == 1) %>% nrow()) / (dat_pred %>% filter(class2 == 'neg') %>% nrow())
  spec <- (dat_pred %>% filter(class2 == 'null' & class2_pred == 0) %>% nrow()) / (dat_pred %>% filter(class2 == 'null') %>% nrow())
  
  print('Sensitivity (Any Interaction) (Overall):')
  print((dat_pred %>% filter(class2 != 'null' & class2_pred != 0) %>% nrow()) / (dat_pred %>% filter(class2 != 'null') %>% nrow()))
  
  print('Sensitivity (Correct Direction) (Overall):')
  print(sens_pos)
  print(sens_neg)
  
  print('Specificity (Overall):')
  print(spec)
  
  # Calculate overall (weighted) accuracy:
  weight_pos <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'pos') %>% nrow())
  weight_neg <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'neg') %>% nrow())
  weight_null <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'null') %>% nrow())
  
  acc_weighted <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  
  print('Overall accuracy (weighted):')
  print(acc_weighted)
  
  sens_pos_vec <- c(sens_pos_vec, sens_pos)
  sens_neg_vec <- c(sens_neg_vec, sens_neg)
  spec_vec <- c(spec_vec, spec)
  acc_weighted_vec <- c(acc_weighted_vec, acc_weighted)
  
}
rm(i, loaded_model, result, dat_pred, sens_pos, sens_neg, spec,
   weight_pos, weight_neg, weight_null, acc_weighted)

# Calculate accuracy by true param values:
acc_by_param_LIST <- lapply(dat_pred_LIST, function(ix) {
  calculate_accuracy_matrix(ix)
})
rm(dat_pred_LIST)

# Get average accuracy over all runs:
summary(sens_pos_vec)
summary(sens_neg_vec)
summary(spec_vec)
summary(acc_weighted_vec)
rm(sens_pos_vec, sens_neg_vec, spec_vec, acc_weighted_vec)

acc_by_param <- acc_by_param_LIST %>%
  bind_rows() %>%
  group_by(strength, duration) %>%
  summarise(perc_correct_mean = mean(perc_correct),
            perc_correct_median = median(perc_correct))
rm(acc_by_param_LIST)

p1 <- ggplot(data = acc_by_param, aes(x = strength, y = duration, fill = perc_correct_median)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Convolutional Neural Network', x = 'True Strength',
       y = 'True Duration (Weeks)', fill = '% Correct (Median)')
# plot(p1)
rm(acc_by_param, datV1, datV1_train, datV1_test, datV2, datV2_train, datV2_test,
   x_train, x_test, y_test, y_train, train_ind)

# ------------------------------------------------------------------------------

# Also explore performance using wider parameter ranges

# Load new test data:
source('src/classification/generate_new_test_data.R')

# Get only data of interest:
dat <- dat %>%
  select(.id, time, V1_obs:V2_obs)

# Join true interaction parameter values:
res_trueparams <- true_params %>%
  t() %>%
  as_tibble() %>%
  select('theta_lambda1', 'delta1', Ri1:Ri2, w2, A1:phi1, R01:R02) %>%
  rename('theta_lambda' = 'theta_lambda1', 'delta' = 'delta1', 'A' = 'A1', 'phi' = 'phi1') %>%
  mutate(.id = 1:ncol(true_params))

dat <- dat %>%
  inner_join(res_trueparams, by = '.id')
rm(true_params, res_trueparams)

# Get true classification:
dat <- dat %>%
  mutate(class1 = if_else(theta_lambda == 1, 'null', 'int'),
         class2 = if_else(theta_lambda < 1, 'neg', class1),
         class2 = if_else(theta_lambda > 1, 'pos', class2)) %>%
  select(.id, time, V1_obs:V2_obs, theta_lambda:class2)

# Format data for use in CNN models:
dat <- dat %>%
  mutate(class1_int = if_else(class1 == 'int', 1, 0),
         class2_int = if_else(class2 == 'null', 0, 1),
         class2_int = if_else(class2 == 'pos', 2, class2_int))

datV1 <- dat %>%
  select(.id:V1_obs, class2_int) %>%
  mutate(time = paste0('t_', time)) %>%
  pivot_wider(names_from = time, values_from = V1_obs)
datV2 <- dat %>%
  select(.id:time, V2_obs, class2_int) %>%
  mutate(time = paste0('t_', time)) %>%
  pivot_wider(names_from = time, values_from = V2_obs)

datV1_mat <- as.matrix(datV1[, 3:ncol(datV1)])
datV2_mat <- as.matrix(datV2[, 3:ncol(datV2)])

x_test <- array(NA, dim = c(dim(datV1_mat), 2))
x_test[,, 1] <- datV1_mat
x_test[,, 2] <- datV2_mat
rm(datV1_mat, datV2_mat)

y_test <- as.matrix(datV1$class2_int)

# Loop through all models and calculate accuracy:
dat_pred_LIST <- vector('list', length = num_runs)
sens_pos_vec = sens_neg_vec = spec_vec = acc_weighted_vec = c()
for (i in 1:num_runs) {
  
  # Load fitted model:
  loaded_model <- load_model(paste0('results/cnn/best_model_', i, '.keras'))
  
  # Get loss/accuracy for test data:
  result <- loaded_model %>% evaluate(x_test, y_test)
  
  sprintf("Test loss: %s", result[["loss"]])
  sprintf("Test accuracy: %s", result[["sparse_categorical_accuracy"]])
  
  # Get predictions for test data:
  dat_pred <- datV1 %>%
    select(.id:class2_int) %>%
    mutate(class2_pred = loaded_model %>%
             predict(x_test) %>%
             apply(1, which.max),
           class2_pred = class2_pred - 1)
  
  dat_pred <- dat %>%
    select(.id, theta_lambda:R02, class2) %>%
    unique() %>%
    inner_join(dat_pred,
               by = '.id')
  
  dat_pred_LIST[[i]] <- dat_pred
  
  # Calculate sensitivity/specificity overall:
  sens_pos <- (dat_pred %>% filter(class2 == 'pos' & class2_pred == 2) %>% nrow()) / (dat_pred %>% filter(class2 == 'pos') %>% nrow())
  sens_neg <- (dat_pred %>% filter(class2 == 'neg' & class2_pred == 1) %>% nrow()) / (dat_pred %>% filter(class2 == 'neg') %>% nrow())
  spec <- (dat_pred %>% filter(class2 == 'null' & class2_pred == 0) %>% nrow()) / (dat_pred %>% filter(class2 == 'null') %>% nrow())
  
  print('Sensitivity (Any Interaction) (Overall):')
  print((dat_pred %>% filter(class2 != 'null' & class2_pred != 0) %>% nrow()) / (dat_pred %>% filter(class2 != 'null') %>% nrow()))
  
  print('Sensitivity (Correct Direction) (Overall):')
  print(sens_pos)
  print(sens_neg)
  
  print('Specificity (Overall):')
  print(spec)
  
  # Calculate overall (weighted) accuracy:
  weight_pos <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'pos') %>% nrow())
  weight_neg <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'neg') %>% nrow())
  weight_null <- min(table(dat_pred$class2)) / (dat_pred %>% filter(class2 == 'null') %>% nrow())
  
  acc_weighted <- (sens_pos * weight_pos + sens_neg * weight_neg + spec * weight_null) / (weight_pos + weight_neg + weight_null)
  
  print('Overall accuracy (weighted):')
  print(acc_weighted)
  
  sens_pos_vec <- c(sens_pos_vec, sens_pos)
  sens_neg_vec <- c(sens_neg_vec, sens_neg)
  spec_vec <- c(spec_vec, spec)
  acc_weighted_vec <- c(acc_weighted_vec, acc_weighted)
  
}
rm(i, loaded_model, result, dat_pred, sens_pos, sens_neg, spec,
   weight_pos, weight_neg, weight_null, acc_weighted,
   datV1, datV2, x_test, y_test)

# Calculate accuracy by true param values:
dat_pred_CAT <- lapply(dat_pred_LIST, function(ix) {
  ix %>%
    mutate(theta_lambda_cat = as.character(cut(theta_lambda, breaks = c(0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0))),
           theta_lambda_cat = if_else(theta_lambda == 1.0, '[1.0]', theta_lambda_cat),
           theta_lambda_cat = factor(theta_lambda_cat, levels = c('(0,0.25]', '(0.25,0.5]', '(0.5,1]', '[1.0]', '(1,2]', '(2,5]', '(5,10]')),
           delta_cat = (cut(1 / delta, breaks = c(1, 4, 8, 13, 17, 21, 26))))
})

acc_by_param_LIST <- lapply(dat_pred_CAT, function(ix) {
  
  all_tl <- levels(ix$theta_lambda_cat)
  all_delta <- levels(ix$delta_cat)
  
  mat_correct <- matrix(nrow = length(all_tl), ncol = length(all_delta))
  rownames(mat_correct) <- all_tl
  colnames(mat_correct) <- all_delta
  
  for (tl in all_tl) {
    for (d in all_delta) {
      
      df_temp <- ix %>% filter(theta_lambda_cat == tl & delta_cat == d)
      mat_correct[rownames(mat_correct) == tl, colnames(mat_correct) == d] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
      
    }
  }
  
  mat_correct <- mat_correct %>%
    as_tibble(rownames = 'strength') %>%
    pivot_longer(-strength, names_to = 'duration', values_to = 'perc_correct') %>%
    mutate(strength = factor(strength, levels = all_tl),
           duration = factor(duration, levels = all_delta))
  return(mat_correct)
  
})
rm(dat_pred_CAT)

# Get average accuracy over all runs:
summary(sens_pos_vec)
summary(sens_neg_vec)
summary(spec_vec)
summary(acc_weighted_vec)
rm(sens_pos_vec, sens_neg_vec, spec_vec, acc_weighted_vec)

acc_by_param <- acc_by_param_LIST %>%
  bind_rows() %>%
  group_by(strength, duration) %>%
  summarise(perc_correct_mean = mean(perc_correct),
            perc_correct_median = median(perc_correct))
rm(acc_by_param_LIST)

p2 <- ggplot(data = acc_by_param, aes(x = strength, y = duration, fill = perc_correct_median)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(title = 'Convolutional Neural Network (WIDE)', x = 'True Strength',
       y = 'True Duration (Weeks)', fill = '% Correct (Median)')
grid.arrange(p1, p2, nrow = 1)
rm(acc_by_param)

# Explore which other parameters might influence accuracy:
acc_by_param_LIST <- lapply(dat_pred_LIST, function(ix) {
  
  ix_cat <- ix %>%
    mutate(Ri1_cat = cut(Ri1, breaks = c(1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0)),
           Ri2_cat = cut(Ri2, breaks = c(1.0, 1.25, 1.5, 2.0, 3.0, 4.0, 5.0)),
           w2_cat = cut(w2, breaks = c(1/(52.25 * 1.0), 1/(52.25 * 0.9), 1/(52.25 * 0.8), 1/(52.25 * 0.7), 1/(52.25 * 0.6))),
           A_cat = cut(A, breaks = seq(0, 0.5, by = 0.1)),
           phi_cat = cut(phi, breaks = c(11, 16, 21, 26, 31, 36, 41)),
           R01_cat = cut(R01, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.45)),
           R02_cat = cut(R02, breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.45)))
  
  mat_correct_Ri1 <- matrix(nrow = 1, ncol = length(levels(ix_cat$Ri1_cat)))
  mat_correct_Ri2 <- matrix(nrow = 1, ncol = length(levels(ix_cat$Ri2_cat)))
  mat_correct_w2 <- matrix(nrow = 1, ncol = length(levels(ix_cat$w2_cat)))
  mat_correct_A <- matrix(nrow = 1, ncol = length(levels(ix_cat$A_cat)))
  mat_correct_phi <- matrix(nrow = 1, ncol = length(levels(ix_cat$phi_cat)))
  mat_correct_R01 <- matrix(nrow = 1, ncol = length(levels(ix_cat$R01_cat)))
  mat_correct_R02 <- matrix(nrow = 1, ncol = length(levels(ix_cat$R02_cat)))
  
  colnames(mat_correct_Ri1) <- levels(ix_cat$Ri1_cat)
  colnames(mat_correct_Ri2) <- levels(ix_cat$Ri2_cat)
  colnames(mat_correct_w2) <- levels(ix_cat$w2_cat)
  colnames(mat_correct_A) <- levels(ix_cat$A_cat)
  colnames(mat_correct_phi) <- levels(ix_cat$phi_cat)
  colnames(mat_correct_R01) <- levels(ix_cat$R01_cat)
  colnames(mat_correct_R02) <- levels(ix_cat$R02_cat)
  
  for (Ri1_val in levels(ix_cat$Ri1_cat)) {
    df_temp <- ix_cat %>% filter(Ri1_cat == Ri1_val)
    mat_correct_Ri1[1, colnames(mat_correct_Ri1) == Ri1_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (Ri2_val in levels(ix_cat$Ri2_cat)) {
    df_temp <- ix_cat %>% filter(Ri2_cat == Ri2_val)
    mat_correct_Ri2[1, colnames(mat_correct_Ri2) == Ri2_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (w2_val in levels(ix_cat$w2_cat)) {
    df_temp <- ix_cat %>% filter(w2_cat == w2_val)
    mat_correct_w2[1, colnames(mat_correct_w2) == w2_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (A_val in levels(ix_cat$A_cat)) {
    df_temp <- ix_cat %>% filter(A_cat == A_val)
    mat_correct_A[1, colnames(mat_correct_A) == A_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (phi_val in levels(ix_cat$phi_cat)) {
    df_temp <- ix_cat %>% filter(phi_cat == phi_val)
    mat_correct_phi[1, colnames(mat_correct_phi) == phi_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (R01_val in levels(ix_cat$R01_cat)) {
    df_temp <- ix_cat %>% filter(R01_cat == R01_val)
    mat_correct_R01[1, colnames(mat_correct_R01) == R01_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  for (R02_val in levels(ix_cat$R02_cat)) {
    df_temp <- ix_cat %>% filter(R02_cat == R02_val)
    mat_correct_R02[1, colnames(mat_correct_R02) == R02_val] <- (df_temp %>% filter(class2_int == class2_pred) %>% nrow()) / nrow(df_temp)
  }
  
  mat_correct_LIST <- vector('list', length = 7)
  
  mat_correct_LIST[[1]] <- mat_correct_Ri1 %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'Ri1',
           value = factor(value, levels = levels(ix_cat$Ri1_cat)))
  mat_correct_LIST[[2]] <- mat_correct_Ri2 %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'Ri2',
           value = factor(value, levels = levels(ix_cat$Ri2_cat)))
  mat_correct_LIST[[3]] <- mat_correct_w2 %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'w2',
           value = factor(value, levels = levels(ix_cat$w2_cat)))
  mat_correct_LIST[[4]] <- mat_correct_A %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'A',
           value = factor(value, levels = levels(ix_cat$A_cat)))
  mat_correct_LIST[[5]] <- mat_correct_phi %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'phi',
           value = factor(value, levels = levels(ix_cat$phi_cat)))
  mat_correct_LIST[[6]] <- mat_correct_R01 %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'R01',
           value = factor(value, levels = levels(ix_cat$R01_cat)))
  mat_correct_LIST[[7]] <- mat_correct_R02 %>%
    as_tibble() %>%
    pivot_longer(everything(), names_to = 'value', values_to = 'perc_correct') %>%
    mutate(param = 'R02',
           value = factor(value, levels = levels(ix_cat$R02_cat)))
  
  return(mat_correct_LIST)
  
})

acc_by_Ri1 <- lapply(acc_by_param_LIST, function(ix) {
  ix[[1]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_Ri2 <- lapply(acc_by_param_LIST, function(ix) {
  ix[[2]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_w2 <- lapply(acc_by_param_LIST, function(ix) {
  ix[[3]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_A <- lapply(acc_by_param_LIST, function(ix) {
  ix[[4]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_phi <- lapply(acc_by_param_LIST, function(ix) {
  ix[[5]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_R01 <- lapply(acc_by_param_LIST, function(ix) {
  ix[[6]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
acc_by_R02 <- lapply(acc_by_param_LIST, function(ix) {
  ix[[7]]
}) %>%
  bind_rows() %>%
  group_by(value) %>%
  summarise(perc_correct = median(perc_correct))
rm(acc_by_param_LIST)

p_ri1 <- ggplot(data = acc_by_Ri1, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'Ri1', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_ri2 <- ggplot(data = acc_by_Ri2, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'Ri2', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_w2 <- ggplot(data = acc_by_w2, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'w2', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_A <- ggplot(data = acc_by_A, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'A', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_phi <- ggplot(data = acc_by_phi, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'phi', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_r01 <- ggplot(data = acc_by_R01, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'none') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'R01', fill = 'perc_correct') +
  scale_y_continuous(breaks = c())
p_r02 <- ggplot(data = acc_by_R02, aes(x = value, y = 1, fill = perc_correct)) +
  geom_tile() +
  theme_classic() +
  theme(legend.position = 'bottom') +
  scale_fill_viridis(limits = c(0, 1), option = 'G') +
  labs(x = '', y = 'R02', fill = '% Correct') +
  scale_y_continuous(breaks = c())

p3 <- arrangeGrob(p_ri1, p_ri2, p_w2, p_A, p_phi, p_r01, p_r02,
                  ncol = 1,
                  heights = c(1, 1, 1, 1, 1, 1, 1.8))
plot(p3)

# Save plots:
pdf(file = 'results/plots/resuts_CNN.pdf', width = 12, height = 7)
grid.arrange(p1, p2, nrow = 1)
plot(p3)
dev.off()

# ------------------------------------------------------------------------------

# Clean up
rm(list = ls())
