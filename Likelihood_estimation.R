################################################################################
#                       Maximum likelihood Estimation       
#
# We are going to do a deterministic approximation of the MLE to make the 
# computation quicker and also so we don't end up doing estimation on the model 
# from which we have generated the data from giving it an unfair advantage
#
# Note: remember that the pomp model to generate each dataset is stored in
# results$pomp_model
#
# Created by: Sarah Pirikahu
# Creation date: 22 Aug 2023
################################################################################

lik <- function(pomp_model, data){
  # adding the new simulated data to the pomp model and telling it which parameters we want to estimate
  # this will then determine the likelihood function for our model 
  f_x <- pomp_model %>% 
            traj_objfun(times="time", t0=105,
              data = data,
              est=c("Ri1","Ri2","E01","E02","R01","R02","R12","rho1","rho2","A1","phi1","A2","phi2","delta1","delta2","theta_lambda1","theta_lambda2",
                    "w1","w2"))
  
  # 
  
  
  
}