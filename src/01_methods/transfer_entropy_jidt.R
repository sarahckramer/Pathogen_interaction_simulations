# ---------------------------------------------------------------------------------------------------------------------
# Code to run transfer entropy

# Course and short tutorial can be found here:
# https://github.com/jlizier/jidt/wiki/Tutorial
# ---------------------------------------------------------------------------------------------------------------------

# Load the rJava library and start the JVM
library(rJava)
.jinit()

# pointing to where the package folder is sitting
.jaddClassPath('infodynamics-dist-1.6.1/infodynamics.jar')

te_jidt <- function(data, lag){
  
  if (nrow(data) > 0) {
    
    # Get relevant data:
    sourceArray <- data$V1_obs_ln
    destArray <- data$V2_obs_ln
    id <- unique(data$.id)
    
    #---- Analysis NOT accounting for seasonal confounding ----#
    
    # Create a TE calculator:
    teCalc <- .jnew('infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov')
    .jcall(teCalc, 'V', 'setProperty', 'k', '4') # use Kraskov parameter k = 4 for nearest 4 points
    .jcall(teCalc, 'V', 'setProperty', 'delay', lag) # lag between source and destination

    # TE calculation for V1 -> V2:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray)

    # Set embedding dimensions (source: 20, dest: 5)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '5') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding

    result_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')

    nullDist_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                               'computeSignificance', 500L)
    mean_null_v2_x_v1 <- .jcall(nullDist_v2_x_v1, 'D', 'getMeanOfDistribution')
    sd_null_v2_x_v1 <- .jcall(nullDist_v2_x_v1, 'D', 'getStdOfDistribution')
    p_value_v2_x_v1 <- nullDist_v2_x_v1$pValue

    # TE calculation for V2 -> V1:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray)

    # Set embedding dimensions (source: 20, dest: 2)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding

    result_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')

    nullDist_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                               'computeSignificance', 500L)
    mean_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, 'D', 'getMeanOfDistribution')
    sd_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, 'D', 'getStdOfDistribution')
    p_value_v1_x_v2 <- nullDist_v1_x_v2$pValue
    
    #---- Analysis accounting for seasonal forcing (conditional transfer entropy) ----#
    
    # Calculate seasonal component:
    data <- data %>%
      mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52.25 * (time - 26))) %>%
      mutate(seasonal_component = scale(log(seasonal_component), scale = FALSE))
    condArray <- data$seasonal_component
    
    # Create a new TE calculator:
    teCalc <- .jnew('infodynamics/measures/continuous/kraskov/ConditionalTransferEntropyCalculatorKraskov')
    
    .jcall(teCalc, 'V', 'setProperty', 'k', '4') # use Kraskov parameter k = 4 for nearest 4 points
    .jcall(teCalc, 'V', 'setProperty', 'delay', lag) # lag between source and destination
    
    # TE calculation for V1 -> V2:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray, condArray)
    
    # Set embedding dimensions (source: 20, dest: 5)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '5') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
    .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '1') # set confounding embedding
    
    result_confound_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
    
    nullDist_confound_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                        'computeSignificance', 500L)
    mean_null_confound_v2_x_v1 <- .jcall(nullDist_confound_v2_x_v1, 'D', 'getMeanOfDistribution')
    sd_null_confound_v2_x_v1 <- .jcall(nullDist_confound_v2_x_v1, 'D', 'getStdOfDistribution')
    p_value_confound_v2_x_v1 <- nullDist_confound_v2_x_v1$pValue
    
    # TE calculation for V2 -> V1:
    .jcall(teCalc, 'V', 'initialise')
    .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray, condArray)
    
    # Set embedding dimensions (source: 20, dest: 2)
    .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
    .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
    .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '1') # set confounding embedding
    
    result_confound_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
    
    nullDist_confound_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                        'computeSignificance', 500L)
    mean_null_confound_v1_x_v2 <- .jcall(nullDist_confound_v1_x_v2, 'D', 'getMeanOfDistribution')
    sd_null_confound_v1_x_v2 <- .jcall(nullDist_confound_v1_x_v2, 'D', 'getStdOfDistribution')
    p_value_confound_v1_x_v2 <- nullDist_confound_v1_x_v2$pValue
    
    #---- Compile and return results ----#
    res <- data.frame(cbind(te = c(result_v1_x_v2, result_v2_x_v1),
                            te_confound = c(result_confound_v1_x_v2, result_confound_v2_x_v1),
                            direction = c('v2 -> v1', 'v1 -> v2'),
                            sd_null = c(sd_null_v1_x_v2, sd_null_v2_x_v1),
                            sd_null_confound = c(sd_null_confound_v1_x_v2, sd_null_confound_v2_x_v1),
                            p_value = c(p_value_v1_x_v2, p_value_v2_x_v1),
                            p_value_confound = c(p_value_confound_v1_x_v2, p_value_confound_v2_x_v1),
                            lag = rep(lag))) %>%
      as_tibble() %>%
      mutate(across(-direction, as.numeric))
    
  } else {
    res <- data.frame(cbind(lag = as.numeric(lag))) %>% as_tibble()
  }
  
  return(res)
  
}
