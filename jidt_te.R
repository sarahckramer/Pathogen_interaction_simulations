

# Load the rJava library and start the JVM
library("rJava")
.jinit()

# Change location of jar to match yours:
#  IMPORTANT -- If using the default below, make sure you have set the working directory
#   in R (e.g. with setwd()) to the location of this file (i.e. demos/r) !!
.jaddClassPath("/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/infodynamics/infodynamics.jar")

te_jidt <- function(data){
  
  sourceArray <- data$v1_obs
  destArray <- data$v2_obs
  id <- unique(data$.id)
  
  # Create a TE calculator for V1 -> V2
  teCalc<-.jnew("infodynamics/measures/continuous/kraskov/TransferEntropyCalculatorKraskov")
  .jcall(teCalc,"V","setProperty", "k", "4") # Use Kraskov parameter K=4 for 4 nearest points
  .jcall(teCalc,"V","setProperty", "k_tau", "1") # lag for destination of 4 
  
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
                          CI_97.5 = c(CI_97.5_v1_x_v2,CI_97.5_v2_x_v1)))
  res$te <- as.numeric(res$te)
  res$sd_null <- as.numeric(res$sd_null)
  res$p_value <- as.numeric(res$p_value)
  res$CI_2.5 <- as.numeric(res$CI_2.5)
  res$CI_97.5 <- as.numeric(res$CI_97.5)
  
  return(res)
}


results$transfer_entropy <- results$data %>% group_by(.id) %>% do(te_jidt(.))


hist(results$transfer_entropy$te)
results$transfer_entropy$significant <- ifelse(results$transfer_entropy$p_value<=0.05, "Y","N")
  
ggplot(aes(x=te, y =.id), data=results$transfer_entropy) + geom_point() + 
  geom_errorbar(aes(xmin = CI_2.5, xmax = CI_97.5, colour=significant)) +
  geom_vline(xintercept = 0, lty=2) + facet_grid(.~direction) + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("theta_lambda=0, delta=1,lag=-8") 


table(results$transfer_entropy$direction,results$transfer_entropy$significant)

# plot for RTransferEntropy package
results$transfer_entropy$id <- rep(1:100, each=2)
results$transfer_entropy$CI_2.5 <- results$transfer_entropy$ete - 1.96*results$transfer_entropy$se
results$transfer_entropy$CI_97.5 <- results$transfer_entropy$ete + 1.96*results$transfer_entropy$se
results$transfer_entropy$significant <- ifelse(results$transfer_entropy$p.value<=0.05, "Y","N")
results$transfer_entropy$lag <- as.factor(results$transfer_entropy$lag)
  
ggplot(aes(x=ete, y =id), data=results$transfer_entropy) + geom_point() + 
  geom_errorbar(aes(xmin = CI_2.5, xmax = CI_97.5, colour=lag)) +
  geom_vline(xintercept = 0, lty=2) + facet_grid(.~direction) + 
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("theta_lambda=1, delta=1") 

table(results$transfer_entropy$direction,results$transfer_entropy$significant)
