###############################################################
#                 Transfer Entropy using jidt package
# 
# Java needs to be installed before this code will work. 
# Installation procedure to just jidt can be found here: 
# https://github.com/jlizier/jidt/tree/master
# Course and short tutorial can be found here:
# https://github.com/jlizier/jidt/wiki/Tutorial
#
# Created by: Sarah Pirikahu
# Creation date: 3 April 2024
###############################################################

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# pointing to where the package folder is sitting
.jaddClassPath('infodynamics-dist-1.6.1/infodynamics.jar')

te_jidt <- function(data, lag){
  
  # Log-transform and center data:
  data <- data %>%
    mutate(V1_obs_ln = scale(log(V1_obs + 1), scale = FALSE),
           V2_obs_ln = scale(log(V2_obs + 1), scale = FALSE))
  
  # Get relevant data:
  sourceArray <- data$V1_obs_ln
  destArray <- data$V2_obs_ln
  id <- unique(data$.id)
  
  #---- Analysis w/o confounding ----#
  
  # Create a TE calculator:
  teCalc <- .jnew('infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov')
  .jcall(teCalc, 'V', 'setProperty', 'k', '4') # use Kraskov parameter k = 4 for nearest 4 points
  .jcall(teCalc, 'V', 'setProperty', 'delay', lag) # lag between source and destination
  
  # TE calculation for V1 -> V2:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray)
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 5)
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
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 2)
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
  
  result_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  
  nullDist_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                             'computeSignificance', 500L)
  mean_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, 'D', 'getMeanOfDistribution')
  sd_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, 'D', 'getStdOfDistribution')
  p_value_v1_x_v2 <- nullDist_v1_x_v2$pValue
  
  #---- Analysis w/ confounding ----#
  
  # Calculate seasonal component:
  data <- data %>%
    mutate(seasonal_component = 1 + 0.2 * cos((2 * pi) / 52.25 * (time - 26)))
  condArray <- data$seasonal_component
  
  # Create a new TE calculator:
  teCalc <- .jnew('infodynamics/measures/continuous/kraskov/ConditionalTransferEntropyCalculatorKraskov')
  
  .jcall(teCalc, 'V', 'setProperty', 'k', '4') # use Kraskov parameter k = 4 for nearest 4 points
  .jcall(teCalc, 'V', 'setProperty', 'delay', lag) # lag between source and destination
  
  # First, use embedding dimension 1 for seasonal component:
  
  # TE calculation for V1 -> V2:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray, condArray)
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 5)
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '5') # set destination embedding
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
  
  result_confound_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  
  nullDist_confound_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                      'computeSignificance', 500L)
  mean_null_confound_v2_x_v1 <- .jcall(nullDist_confound_v2_x_v1, 'D', 'getMeanOfDistribution')
  sd_null_confound_v2_x_v1 <- .jcall(nullDist_confound_v2_x_v1, 'D', 'getStdOfDistribution')
  p_value_confound_v2_x_v1 <- nullDist_confound_v2_x_v1$pValue
  
  # TE calculation for V2 -> V1:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray, condArray)
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 2)
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
  
  result_confound_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  
  nullDist_confound_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                      'computeSignificance', 500L)
  mean_null_confound_v1_x_v2 <- .jcall(nullDist_confound_v1_x_v2, 'D', 'getMeanOfDistribution')
  sd_null_confound_v1_x_v2 <- .jcall(nullDist_confound_v1_x_v2, 'D', 'getStdOfDistribution')
  p_value_confound_v1_x_v2 <- nullDist_confound_v1_x_v2$pValue
  
  # Also try embedding dimension equal to embedding dimension for destination:
  
  # TE calculation for V1 -> V2:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', sourceArray, destArray, condArray)
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 5)
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '5') # set destination embedding
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
  .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '5') # set conditional variable embedding
  
  result_confound2_v2_x_v1 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  
  nullDist_confound2_v2_x_v1 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                      'computeSignificance', 500L)
  mean_null_confound2_v2_x_v1 <- .jcall(nullDist_confound2_v2_x_v1, 'D', 'getMeanOfDistribution')
  sd_null_confound2_v2_x_v1 <- .jcall(nullDist_confound2_v2_x_v1, 'D', 'getStdOfDistribution')
  p_value_confound2_v2_x_v1 <- nullDist_confound2_v2_x_v1$pValue
  
  # TE calculation for V2 -> V1:
  .jcall(teCalc, 'V', 'initialise')
  .jcall(teCalc, 'V', 'setObservations', destArray, sourceArray, condArray)
  
  # SET EMBEDDING DIMENSIONS (source: 20, dest: 2)
  .jcall(teCalc, 'V', 'setProperty', 'k_history', '2') # set destination embedding
  .jcall(teCalc, 'V', 'setProperty', 'l_history', '20') # set source embedding
  .jcall(teCalc, 'V', 'setProperty', 'cond_embed_lengths', '2') # set conditional variable embedding
  
  result_confound2_v1_x_v2 <- .jcall(teCalc, 'D', 'computeAverageLocalOfObservations')
  
  nullDist_confound2_v1_x_v2 <- .jcall(teCalc, 'Linfodynamics/utils/EmpiricalMeasurementDistribution;',
                                      'computeSignificance', 500L)
  mean_null_confound2_v1_x_v2 <- .jcall(nullDist_confound2_v1_x_v2, 'D', 'getMeanOfDistribution')
  sd_null_confound2_v1_x_v2 <- .jcall(nullDist_confound2_v1_x_v2, 'D', 'getStdOfDistribution')
  p_value_confound2_v1_x_v2 <- nullDist_confound2_v1_x_v2$pValue
  
  #---- Compile and return results ----#
  res <- data.frame(cbind(te = c(result_v1_x_v2, result_v2_x_v1),
                          te_confound = c(result_confound_v1_x_v2, result_confound_v2_x_v1),
                          te_confound2 = c(result_confound2_v1_x_v2, result_confound2_v2_x_v1),
                          direction = c('v2 -> v1', 'v1 -> v2'),
                          sd_null = c(sd_null_v1_x_v2, sd_null_v2_x_v1),
                          sd_null_confound = c(sd_null_confound_v1_x_v2, sd_null_confound_v2_x_v1),
                          sd_null_confound2 = c(sd_null_confound2_v1_x_v2, sd_null_confound2_v2_x_v1),
                          p_value = c(p_value_v1_x_v2, p_value_v2_x_v1),
                          p_value_confound = c(p_value_confound_v1_x_v2, p_value_confound_v2_x_v1),
                          p_value_confound2 = c(p_value_confound2_v1_x_v2, p_value_confound2_v2_x_v1),
                          lag = rep(lag))) %>%
    as_tibble() %>%
    mutate(across(-direction, as.numeric))
  
  return(res)
  
}
