###############################################################
#                 Transfer Entropy using jidt package
# 
# Java needs to be installed before this code will work. 
# Installation procedure to just jidt can be found here: 
# https://github.com/jlizier/jidt/tree/master
#
# Created by: Sarah Pirikahu
# Creation date: 3 April 2024
###############################################################

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# pointing to where the package folder is sitting 
.jaddClassPath("/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/infodynamics/infodynamics.jar")

te_jidt <- function(data, lag){
  
  sourceArray <- data$v1_obs
  destArray <- data$v2_obs
  id <- unique(data$.id)
  k_tau <- lag
  
  # Create a TE calculator for V1 -> V2
  teCalc<-.jnew("infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov")
  .jcall(teCalc,"V","setProperty", "k", "4") # Use Kraskov parameter K=4 for 4 nearest points
  .jcall(teCalc,"V","setProperty", "k_tau", k_tau) # lag for destination 
  
  # Perform calculation with correlated source:
  .jcall(teCalc,"V","initialise", 1L) # Use history length 1 (Schreiber k=1)
  .jcall(teCalc,"V","setObservations", sourceArray, destArray)
  result_v2_x_v1 <- .jcall(teCalc,"D","computeAverageLocalOfObservations")
  
  nullDist_v2_x_v1 <- .jcall(teCalc,"Linfodynamics/utils/EmpiricalMeasurementDistribution;",
                             "computeSignificance", 100L)
  mean_null_v2_x_v1 <- .jcall(nullDist_v2_x_v1, "D", "getMeanOfDistribution")
  sd_null_v2_x_v1 <- .jcall(nullDist_v2_x_v1, "D", "getStdOfDistribution")
  p_value_v2_x_v1 <- nullDist_v2_x_v1$pValue
  
  CI_2.5_v2_x_v1 <- result_v2_x_v1 - 1.96*sd_null_v2_x_v1
  CI_97.5_v2_x_v1 <- result_v2_x_v1 + 1.96*sd_null_v2_x_v1
  
  # TE calculation for V2 -> V1 
  .jcall(teCalc,"V","setObservations", destArray, sourceArray)
  result_v1_x_v2 <- .jcall(teCalc,"D","computeAverageLocalOfObservations")
  
  nullDist_v1_x_v2 <- .jcall(teCalc,"Linfodynamics/utils/EmpiricalMeasurementDistribution;",
                             "computeSignificance", 100L)
  mean_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, "D", "getMeanOfDistribution")
  sd_null_v1_x_v2 <- .jcall(nullDist_v1_x_v2, "D", "getStdOfDistribution")
  p_value_v1_x_v2 <- nullDist_v1_x_v2$pValue
  
  CI_2.5_v1_x_v2 <- result_v1_x_v2 - 1.96*sd_null_v1_x_v2
  CI_97.5_v1_x_v2 <- result_v1_x_v2 + 1.96*sd_null_v1_x_v2
  
  res <- data.frame(cbind(id=rep(id,2),te = c(result_v1_x_v2, result_v2_x_v1), 
                          direction=c("v1 -> v2", "v2 -> v1"), 
                          sd_null = c(sd_null_v1_x_v2,sd_null_v2_x_v1), 
                          p_value = c(p_value_v1_x_v2,p_value_v2_x_v1),
                          CI_2.5 = c(CI_2.5_v1_x_v2,CI_2.5_v2_x_v1),
                          CI_97.5 = c(CI_97.5_v1_x_v2,CI_97.5_v2_x_v1), 
                          lag=rep(lag2)))
  res$te <- as.numeric(res$te)
  res$sd_null <- as.numeric(res$sd_null)
  res$p_value <- as.numeric(res$p_value)
  res$CI_2.5 <- as.numeric(res$CI_2.5)
  res$CI_97.5 <- as.numeric(res$CI_97.5)
  
  return(res)
}
