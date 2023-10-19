################################################################################
#                       Maximum likelihood Estimation       
#
# We are going to do a deterministic approximation of the MLE to make the 
# computation quicker and also so we don't end up doing estimation on the model 
# from which we have generated the data from giving it an unfair advantage
#
# inputs: data = data with time, v1_obs, v2_obs
#         true_params = the true parameter value inputs (simply to specify param names for pomp model)
#         components_l = the csnippets that specify the components of the pomp model 
#         sobol_size = the total number of different initial conditions we want to do for the numerical optimizer
#         jobid  = 
#         maxtime = the maximum amount of time we want to allow the numerical optimizer to run 
#
# Created by: Sarah Pirikahu
# Creation date: 22 Aug 2023
################################################################################

# load packages 
library(nloptr)

lik <- function(data, true_params, components_l = components_l, sobol_size, jobid, maxtime){
  
  # creating new pomp model with the simulated data 
  po <- pomp(data = data,
             times = "time",
             t0 = data$time[1],
             obsnames = c('v1_obs', 'v2_obs'),
             accumvars = c('v1_T', 'v2_T'),
             statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                            'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                            'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                            'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                            'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                            'v1_T', 'v2_T'),
             paramnames = names(true_params),
             params = true_params,
             partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
             globals = components_l[['globs']],
             dmeasure = components_l[['dmeas']],
             rmeasure = components_l[['rmeas']],
             rprocess = euler(step.fun = components_l[['rsim']], delta.t = 1),
             skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
             rinit = components_l[['rinit']]
  )
  
  # parameters to estimate
  est_pars <- c("Ri1","Ri2","E01","E02","R01","R02","R12","rho1","rho2","A1","phi1","A2","phi2","delta1","delta2","theta_lambda1","theta_lambda2",
                "w1","w2")
  
  # Working out the objective function (i.e. working out the likelihood) and
  # specifying the parameters we want to estimate 
  fx <- traj_objfun(data=po,
                    est = est_pars,
                    partrans = po@partrans,
                    verbose=TRUE)
  
  # set up starting values for the numerical optimiser 
  # starting range for each parameter 
  start_range <- data.frame(Ri1 = c(1.0, 2.1),
                            Ri2 = c(1.0, 3),
                            E01 = c(0, 0.0001),
                            E02 = c(0, 0.0001),
                            R01 = c(0.1, 0.5),
                            R02 = c(0.15, 0.3),
                            R12 = c(0, 0.0001),
                            rho1 = c(0.002, 0.004),
                            rho2 =c(0.001, 0.004),
                            A1 = c(0.1, 0.4),
                            phi1 = c(24, 28),
                            A2 = c(0.2, 0.4),
                            phi2 = c(20, 25),
                            delta1 = c(1, 1/3),
                            delta2 = c(1, 1/3),
                            theta_lambda1 = c(0, 2),
                            theta_lambda2 = c(0, 2),
                            w1 = c(1/26,1/78),
                            w2 = c(1/26, 1/29))
  # chosing the ranges of the parameters we are estimating 
  start_range <- start_range[, est_pars]
  
  # coming up with a number of starting values based on the starting ranges specified above
  # the starting values are based on the Latin hypercube sampling (note: all parameters vary - none are fixed)
  start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                               upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                               nseq = sobol_size)
  
  # Get unique identifiers:
  sub_start <- 1:sobol_size
  
  # Loop through start values and perform trajectory matching:
  out <- vector(mode = "list", length = length(sub_start)) # initialise output vector
  for (i in seq_along(sub_start)) {
    # Get param start values:
    x0 <- as.numeric(start_values[sub_start[i], ])
    coef(po, est_pars) <- x0
    x0_trans <- coef(po, est_pars, transform = TRUE)
    
    # check if transform has been done correctly 
    po2 <- po
    p0_trans <- coef(po, transform = T)
    coef(po2, transform = T) <- p0_trans
    stopifnot(all.equal(coef(po), coef(po2)))
    
    # Run trajectory matching using subplex algorithm:
    # http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
    tic <- Sys.time()
    m <- try(
      nloptr(x0 = unname(x0_trans),
             eval_f = fx, # evaluating the pomp models dmeasure component? 
             opts = list(algorithm = 'NLOPT_LN_SBPLX',
                         maxtime = maxtime,
                         maxeval = -1, # disabled
                         xtol_rel = -1, # disabled
                         print_level = 1))
    )
    toc <- Sys.time()
    # estimate run time and print 
    etime <- toc - tic
    units(etime) <- 'mins'
    print(etime)
    
    # If estimation is successful, save results:
    if (!inherits(m, 'try-error')) {
      coef(po, est_pars, transform = TRUE) <- m$solution
      
      # Collect all results:
      out[[i]] <- list(allpars = coef(po),
                       estpars = coef(po, est_pars),
                       ll = -m$objective,
                       conv = m$status,
                       message = m$message,
                       niter = m$iterations,
                       etime = as.numeric(etime))
      # Write to file:
      # saveRDS(out,
      #         file = sprintf('results/res_%s_%s_%d.rds',
      #                        vir1, vir2,
      #                        sub_start[i])
      #)
      
  
    }
    
  }
}

# res <- NULL
# for(i in 1:sobol_size){
#   res <- cbind(res, out[[i]]$estpars) 
# }
# 
# res <- data.frame(t(res), row.names=NULL)
# 
# res_long <- gather(res, param, value, Ri1:w2, factor_key=T)
# 
# ggplot(aes(x=value), data=res_long) + geom_histogram() + 
#   facet_wrap(.~param, scales="free") 
# 
# res_long %>% group_by(param) %>% 
#   summarise(mean=mean(value),median = median(value), sd=sd(value),
#             Q0_025 = quantile(value, 0.025), Q50=quantile(value, 0.5), Q975=quantile(value, 0.975))
# 
