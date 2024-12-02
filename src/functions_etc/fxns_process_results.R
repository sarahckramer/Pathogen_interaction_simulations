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
  
  df <- df %>%
    rename('metric' = all_of(met))
  
  res_temp <- NULL
  
  for (d in unique(df$delta)) {
    
    if (method %in% c('granger', 'te', 'ccm')) {
      
      cor_temp_overall <- df %>%
        mutate(theta_lambda = if_else(theta_lambda < 1, 1 / theta_lambda, theta_lambda)) %>%
        filter(delta == d | theta_lambda == 1) %>%
        # filter(p_value < 0.05) %>%
        cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
      
    } else {
      
      cor_temp_overall <- df %>%
        filter(delta == d | theta_lambda == 1) %>%
        # filter(p_value < 0.05) %>%
        cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
      
    }
    
    cor_temp_neg <- df %>%
      filter(theta_lambda < 1,
             delta == d) %>%
      filter(sig == 'yes') %>%
      cor.test(~ theta_lambda + metric, data = ., method = 'spearman') # lm(data = ., metric ~ theta_lambda)
    cor_temp_pos <- df %>%
      filter(theta_lambda > 1,
             delta == d) %>%
      filter(sig == 'yes') %>%
      cor.test(~ theta_lambda + metric, data = ., method = 'spearman')
    
    res_temp <- res_temp %>% bind_rows(rbind(c(d, cor_temp_overall$estimate, cor_temp_overall$p.value),
                                             c(d, cor_temp_neg$estimate, cor_temp_neg$p.value),
                                             c(d, cor_temp_pos$estimate, cor_temp_pos$p.value)) %>%
                                         as_tibble() %>%
                                         mutate(true_int = c('all', 'neg', 'pos')))
    
  }
  
  names(res_temp) <- c('delta', 'rho', 'p_value', 'true_int')
  
  if (method %in% c('granger', 'te', 'ccm')) {
    
    res_temp <- res_temp %>%
      mutate(rho = if_else(true_int == 'neg', -1 * rho, rho))
    
  }
  
  return(res_temp)
  
}

# Function to calculate Matthews correlation coefficient (MCC):
mcc <- function(tp, tn, fp, fn) {
  return((tp * tn - fp * fn) / sqrt(exp(log(tp + fp) + log(tp + fn) + log(tn + fp) + log(tn + fn))))
}
