# ------------------------------------------------------------------------------
# Functions used for processing "main" results
# ------------------------------------------------------------------------------

# Function to calculate accuracy for all values of theta_lambda/delta:
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
        mat_correct[rownames(mat_correct) == tl, colnames(mat_correct) == d] <- (df_temp %>% filter(int_est == int_true) %>% nrow()) / nrow(df_temp)
        
      }
    } else {
      
      df_temp <- df %>% filter(theta_lambda == tl)
      mat_correct[rownames(mat_correct) == tl, ] <- (df_temp %>% filter(int_est == int_true) %>% nrow()) / nrow(df_temp)
      
    }
    
  }
  
  mat_correct <- mat_correct %>%
    as_tibble(rownames = 'strength') %>%
    pivot_longer(-strength, names_to = 'duration', values_to = 'perc_correct') %>%
    mutate(duration = factor(duration, levels = c(1, 4, 13)))
  return(mat_correct)
  
}

# Function to assess whether higher values of a method's metric are associated with stronger true interaction strengths:
calculate_assoc_true_strength <- function(df, method, met) {
  
  # Format tibble for use in regression model:
  df <- df %>%
    rename('metric' = all_of(met))
  
  df <- df %>%
    mutate(theta_lambda = if_else(theta_lambda == 0, 1/16, theta_lambda)) %>%
    mutate(theta_lambda = if_else(theta_lambda < 1, 1 / theta_lambda, theta_lambda)) %>%
    mutate(delta = factor(1 / delta))
  df <- df %>%
    mutate(ln_x = log(theta_lambda, base = 2)) %>%
    mutate(ln_y = log(abs(metric), base = 2))
  
  # # Plot:
  # p1 <- ggplot(data = df, aes(x = ln_x, y = metric, group = ln_x)) + geom_violin() + facet_wrap(~ delta, ncol = 1) + theme_classic()
  # print(p1)
  
  # Fit model and get confidence intervals:
  # m1 <- lm(ln_y ~ ln_x * delta, data = df)
  m1 <- lmer(ln_y ~ ln_x * delta + (1 | .id), data = df %>% mutate(.id = factor(.id)))
  
  coefs <- summary(m1)$coefficients[, 'Estimate']
  coefs_var <- vcov(m1)
  
  assoc_lmer <- c(coefs['ln_x'], confint(m1)['ln_x', ])
  
  plusminus_4 <- qt(0.975, df = nrow(df) - 7) * sqrt(coefs_var['ln_x', 'ln_x'] + coefs_var['ln_x:delta4', 'ln_x:delta4'] + 2 * coefs_var['ln_x', 'ln_x:delta4'])
  plusminus_13 <- qt(0.975, df = nrow(df) - 7) * sqrt(coefs_var['ln_x', 'ln_x'] + coefs_var['ln_x:delta13', 'ln_x:delta13'] + 2 * coefs_var['ln_x', 'ln_x:delta13'])
  
  # print(c(qt(0.975, df = nrow(df) - 6) * sqrt(coefs_var['ln_x', 'ln_x']), plusminus_4, plusminus_13))
  
  assoc_lmer <- assoc_lmer %>%
    rbind(c(coefs['ln_x'] + coefs['ln_x:delta4'], coefs['ln_x'] + coefs['ln_x:delta4'] - plusminus_4, coefs['ln_x'] + coefs['ln_x:delta4'] + plusminus_4)) %>%
    rbind(c(coefs['ln_x'] + coefs['ln_x:delta13'], coefs['ln_x'] + coefs['ln_x:delta13'] - plusminus_13, coefs['ln_x'] + coefs['ln_x:delta13'] + plusminus_13)) %>%
    as_tibble()
  names(assoc_lmer) <- c('coef', 'lower', 'upper')
  
  assoc_lmer <- assoc_lmer %>% mutate(delta = c(1, 4, 13))
  
  return(assoc_lmer)
  
}

# Function to calculate Matthews correlation coefficient (MCC):
mcc <- function(tp, tn, fp, fn) {
  return((tp * tn - fp * fn) / sqrt(exp(log(tp + fp) + log(tp + fn) + log(tn + fp) + log(tn + fn))))
}
