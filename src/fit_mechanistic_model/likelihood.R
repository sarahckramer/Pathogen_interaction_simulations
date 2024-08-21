# ---------------------------------------------------------------------------------------------------------------------
#         Fit model using maximum likelihood approach
# 
# Code to fit a mechanistic model of virus-virus cocirculation
# to the synthetic data; by explicitly accounting for the
# mechanisms underlying both transmission and interaction, this
# approach may be more successful than purely statistical
# methods.
#
# Created by: Sarah Kramer
# Creation date: 24 July 2024
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load packages:
library(nloptr)
library(foreach)
library(doParallel)
library(doSNOW)
library(doMC)

print(detectCores())

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid) # 1:160
# jobid <- 31
run_parallel <- as.logical(Sys.getenv("RUNPARALLEL")); print(run_parallel)

data_id <- (jobid - 1) %% 10 + 1; print(data_id)
jobid <- ceiling(jobid / 10); print(jobid)

# Set relevant parameters:
use_traj_matching <- TRUE
sobol_size <- 500
search_type <- 'broad'

n_cores <- 50 # Number of cores to use
time_max <- 23.75 # Maximal execution time (in hours)

# ---------------------------------------------------------------------------------------------------------------------

# Get model and data

# Create model and generate synthetic data:
source('src/generate_data.R')

# Select a random subset of generated datasets:
set.seed(93859)
ids_to_fit <- sample(1:100, size = 10)
data_id <- ids_to_fit[data_id]

# dat %>%
#   select(time:.id, V1_obs:V2_obs) %>%
#   pivot_longer(V1_obs:V2_obs, names_to = 'vir') %>%
#   filter(.id %in% ids_to_fit) %>%
#   ggplot(aes(x = time, y = value, col = vir)) +
#   geom_line() +
#   facet_wrap(~ .id) +
#   theme_classic()

# Get single dataset, with only columns of interest:
dat <- dat %>%
  filter(.id == data_id) %>%
  select(time, V1_obs:V2_obs) %>%
  as_tibble()
expect_true(nrow(dat) == 522)

# Add data to model object:
resp_mod@data[, 106:627] <- dat %>%
  select(-time) %>%
  t()

# Check that all parameters not being estimated are set to their true values:
coef(resp_mod) <- true_params[, data_id]
expect_true(all(coef(resp_mod) == true_params[, data_id]))

# # Check particle filter:
# pf_res <- c()
# for (i in 1:10) {
#   pf_res <- c(pf_res, resp_mod %>% pfilter(Np = 100))
# }
# logmeanexp(map_vec(pf_res, logLik), se = TRUE)
# 
# 
# p_to_try <- c(1000, 2000, 5000)
# pf_list <- list()
# for (i in 1:length(p_to_try)) {
#   
#   pf_res <- c()
#   for (ix in 1:5) {
#     pf_res <- c(pf_res, resp_mod %>% pfilter(Np = p_to_try[i]))
#   }
#   pf_list[[i]] <- pf_res
#   
# }


# ---------------------------------------------------------------------------------------------------------------------

# Setup fitting exercise

# Choose parameters to be estimated:
if (use_traj_matching) {
  estpars <- c('Ri1', 'Ri2', 'rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'delta2',
               'A1', 'phi1', 'A2', 'phi2', 'k1', 'k2',
               't_si_1', 't_si_2', 't_si_3', 't_si_4', 't_si_5', 't_si_6',
               't_si_7', 't_si_8', 't_si_9', 't_si_10',
               'w_delta_i_1', 'w_delta_i_2', 'w_delta_i_3', 'w_delta_i_4', 'w_delta_i_5', 'w_delta_i_6',
               'w_delta_i_7', 'w_delta_i_8', 'w_delta_i_9', 'w_delta_i_10')
} else {
  estpars <- c('Ri1', 'Ri2', 'rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'delta2',
               'A1', 'phi1', 'A2', 'phi2', 'k1', 'k2', 'beta_sd1', 'beta_sd2',
               't_si_1', 't_si_2', 't_si_3', 't_si_4', 't_si_5', 't_si_6',
               't_si_7', 't_si_8', 't_si_9', 't_si_10',
               'w_delta_i_1', 'w_delta_i_2', 'w_delta_i_3', 'w_delta_i_4', 'w_delta_i_5', 'w_delta_i_6',
               'w_delta_i_7', 'w_delta_i_8', 'w_delta_i_9', 'w_delta_i_10')
}

# Set start ranges for parameters of interest:
start_range <- data.frame(Ri1 = c(1.0, 2.0),
                          Ri2 = c(1.0, 2.0),
                          w1 = c(1 / (52.25 * 5), 1),
                          w2 = c(1 / (52.25 * 5), 1),
                          rho1 = c(0, 1.0),
                          rho2 = c(0, 1.0),
                          theta_lambda1 = c(0, 5.0),
                          theta_lambda2 = c(0, 5.0),
                          delta1 = c(7 / 60, 7),
                          delta2 = c(7 / 60, 7),
                          A1 = c(0.05, 0.5),
                          phi1 = c(0, 52.25),
                          # phi1 = c(20, 42),
                          A2 = c(0.05, 0.5),
                          phi2 = c(0, 52.25),
                          # phi2 = c(20, 42),
                          k1 = c(0, 0.1),
                          k2 = c(0, 0.1),
                          beta_sd1 = c(0, 0.2),
                          beta_sd2 = c(0, 0.2),
                          E10 = c(0, 1e-3),
                          E20 = c(0, 1e-3),
                          R10 = c(0, 0.3),
                          R20 = c(0, 0.3),
                          R120 = c(0, 0.3),
                          t_si_1 = c(105, 156),
                          t_si_2 = c(157, 208),
                          t_si_3 = c(209, 260),
                          t_si_4 = c(261, 312),
                          t_si_5 = c(313, 365),
                          t_si_6 = c(366, 417),
                          t_si_7 = c(418, 469),
                          t_si_8 = c(470, 521),
                          t_si_9 = c(522, 573),
                          t_si_10 = c(574, 626),
                          w_delta_i_1 = c(0, 0.5),
                          w_delta_i_2 = c(0, 0.5),
                          w_delta_i_3 = c(0, 0.5),
                          w_delta_i_4 = c(0, 0.5),
                          w_delta_i_5 = c(0, 0.5),
                          w_delta_i_6 = c(0, 0.5),
                          w_delta_i_7 = c(0, 0.5),
                          w_delta_i_8 = c(0, 0.5),
                          w_delta_i_9 = c(0, 0.5),
                          w_delta_i_10 = c(0, 0.5))
start_range <- start_range[, estpars]

# Draw from start ranges:
if (search_type == 'broad') {
  set.seed(38564478)
  start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                               upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                               nseq = sobol_size)
}

print(estpars)
print(start_range)
print(summary(start_values))

# Create objective function for call to nloptr:
obj_fun <- traj_objfun(data = resp_mod,
                       est = estpars,
                       partrans = resp_mod@partrans,
                       verbose = TRUE)

# Run in parallel?:
if (run_parallel) {
  
  # Set maximal execution time for each estimation:
  nmins_exec <- time_max * 60 / (sobol_size / n_cores)
  print(sprintf("Max estimation time=%.1f min", nmins_exec))
  
  # Get unique identifiers:
  sub_start <- 1:n_cores
  
  # Set up parallelization:
  registerDoMC(n_cores)
  
  print(getDoParRegistered())
  print(getDoParWorkers())
  
  # Loop through start value sets:
  for (i in 1:(sobol_size / n_cores)) {
    
    # Get subset:
    start_values_temp <- start_values[(n_cores * (i - 1) + 1):(n_cores * i), ]
    
    # Transform start values:
    start_values_tran <- t(
      apply(start_values_temp, 1, function(ix) {
        coef(resp_mod, estpars) <- as.numeric(ix)
        coef(resp_mod, estpars, transform = TRUE)
      }, simplify = TRUE)
    )
    x0_trans_names <- colnames(start_values_tran)
    print(x0_trans_names)
    
    if (use_traj_matching) {
      
      # Fit:
      tic <- Sys.time()
      m <- foreach(ix = sub_start, .packages = c('tidyverse', 'testthat', 'pomp', 'nloptr')) %dopar% {
        
        x0_trans <- start_values_tran[ix, ]
        
        return(
          try(
            nloptr(x0 = x0_trans,
                   eval_f = obj_fun,
                   opts = list(algorithm = "NLOPT_LN_SBPLX",
                               maxtime = 60 * nmins_exec,
                               maxeval = -1, # Negative value: criterion is disabled
                               xtol_rel = -1, # Default value: 1e-4
                               print_level = 0,
                               ranseed = 12345))
          )
        )
        
      }
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- 'hours'
      print(etime)
      
      # Process results:
      m <- lapply(m, function(ix) {
        
        if (!inherits(ix, 'try-error')) {
          
          coef(resp_mod, estpars, transform = TRUE) <- ix$solution
          x0_fit_untrans <- coef(resp_mod, estpars)
          
          out <- list(estpars = x0_fit_untrans,
                      ll = -ix$objective,
                      conv = ix$status,
                      message = ix$message,
                      niter = ix$iterations)
          
        } else {
          out <- 'error'
        }
        
        return(out)
        
      })
      
    } else {
      
      # Run MIF to fit to stochastic model:
      tic <- Sys.time()
      
      mf <- foreach(i = sub_start, .packages = c('tidyverse', 'testthat', 'pomp')) %dopar% {
        
        x0_trans <- start_values_tran[i, ]
        coef(resp_mod, estpars, transform = TRUE) <- x0_trans
        
        return(
          try(
            mif2(resp_mod,
                 params = coef(resp_mod),
                 Np = 1000,
                 Nmif = 50,
                 cooling.fraction.50 = 0.5,
                 rw.sd = rw_sd(Ri1 = ifelse(time >= 105, 0.02, 0),
                               Ri2 = ifelse(time >= 105, 0.02, 0),
                               rho1 = ifelse(time >= 105, 0.02, 0),
                               rho2 = ifelse(time >= 105, 0.02, 0),
                               theta_lambda1 = ifelse(time >= 105, 0.02, 0),
                               theta_lambda2 = ifelse(time >= 105, 0.02, 0),
                               delta1 = ifelse(time >= 105, 0.02, 0),
                               delta2 = ifelse(time >= 105, 0.02, 0),
                               A1 = ifelse(time >= 105, 0.02, 0),
                               phi1 = ifelse(time >= 105, 0.02, 0),
                               A2 = ifelse(time >= 105, 0.02, 0),
                               phi2 = ifelse(time >= 105, 0.02, 0),
                               k1 = ifelse(time >= 105, 0.02, 0),
                               k2 = ifelse(time >= 105, 0.02, 0),
                               # beta_sd1 = ifelse(time >= 105, 0.02, 0),
                               # beta_sd2 = ifelse(time >= 105, 0.02, 0),
                               t_si_1 = ifelse(time >= 105 & time <= 156, 0.1, 0),
                               t_si_2 = ifelse(time >= 157 & time <= 208, 0.1, 0),
                               t_si_3 = ifelse(time >= 209 & time <= 260, 0.1, 0),
                               t_si_4 = ifelse(time >= 261 & time <= 312, 0.1, 0),
                               t_si_5 = ifelse(time >= 313 & time <= 365, 0.1, 0),
                               t_si_6 = ifelse(time >= 366 & time <= 417, 0.1, 0),
                               t_si_7 = ifelse(time >= 418 & time <= 469, 0.1, 0),
                               t_si_8 = ifelse(time >= 470 & time <= 521, 0.1, 0),
                               t_si_9 = ifelse(time >= 522 & time <= 573, 0.1, 0),
                               t_si_10 = ifelse(time >= 574 & time <= 626, 0.1, 0),
                               w_delta_i_1 = ifelse(time >= 105 & time <= 156, 0.02, 0),
                               w_delta_i_2 = ifelse(time >= 157 & time <= 208, 0.02, 0),
                               w_delta_i_3 = ifelse(time >= 209 & time <= 260, 0.02, 0),
                               w_delta_i_4 = ifelse(time >= 261 & time <= 312, 0.02, 0),
                               w_delta_i_5 = ifelse(time >= 313 & time <= 365, 0.02, 0),
                               w_delta_i_6 = ifelse(time >= 366 & time <= 417, 0.02, 0),
                               w_delta_i_7 = ifelse(time >= 418 & time <= 469, 0.02, 0),
                               w_delta_i_8 = ifelse(time >= 470 & time <= 521, 0.02, 0),
                               w_delta_i_9 = ifelse(time >= 522 & time <= 573, 0.02, 0),
                               w_delta_i_10 = ifelse(time >= 574 & time <= 626, 0.02, 0))) %>%
              mif2(Nmif = 50,
                   cooling.fraction.50 = 0.25)# %>%
            # mif2(Nmif = 100,
            #      cooling.fraction.50 = 0.1)
          )
        )
        
      }
      
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- 'mins'
      print(etime)

      # Write to file (temporary):
      saveRDS(mf,
            file = sprintf('results/res_%d_%d_%d_PARALLEL_mif_TEMP.rds',
                           jobid,
                           data_id,
                           i)
      )
      
      # Process results:
      m <- lapply(mf, function(ix) {
        
        if (!inherits(ix, 'try-error')) {
          
          ll <- replicate(10, ix %>% pfilter(Np = 2500) %>% logLik()) %>%
            logmeanexp(se = TRUE)
          
          mf.ll <- ix %>%
            coef() %>%
            bind_rows() %>%
            bind_cols(loglik = ll[1],
                      loglik.se = ll[2])
          
          out <- list(mf = ix,
                      res = mf.ll)
          
        } else {
          out <- 'error'
        }
        
        return(out)
        
      })
      
    }
    
    # Write to file:
    saveRDS(m,
            file = sprintf('results/res_%d_%d_%d_PARALLEL.rds',
                           jobid,
                           data_id,
                           i)
    )
    
  }
  
} else {
  
  # Set maximal execution time for each estimation:
  nmins_exec <- time_max * 60 / sobol_size
  print(sprintf("Max estimation time=%.1f min", nmins_exec))
  
  # Loop through start values and perform trajectory matching:
  for (i in seq_along(1:sobol_size)) {
    
    print(paste0('Estimation: ', i))
    
    # Get param start values:
    x0 <- as.numeric(start_values[i, ])
    coef(resp_mod, estpars) <- x0
    x0_trans <- coef(resp_mod, estpars, transform = TRUE)
    print(-1 * obj_fun(x0_trans))
    
    if (use_traj_matching) {
      
      # Run trajectory matching using subplex algorithm:
      # http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
      tic <- Sys.time()
      m <- try(
        nloptr(x0 = unname(x0_trans),
               eval_f = obj_fun,
               opts = list(algorithm = 'NLOPT_LN_SBPLX',
                           maxtime = 60.0 * nmins_exec,
                           maxeval = -1, # disabled
                           xtol_rel = -1, # disabled; default: 1e-4
                           print_level = 0))
      )
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- 'mins'
      print(etime)
      
      # If estimation is successful, compile results:
      if (!inherits(m, 'try-error')) {
        coef(resp_mod, estpars, transform = TRUE) <- m$solution
        
        # Collect all results:
        out <- list(allpars = coef(resp_mod),
                    estpars = coef(resp_mod, estpars),
                    ll = -m$objective,
                    conv = m$status,
                    message = m$message,
                    niter = m$iterations,
                    etime = as.numeric(etime))
        
        # Print results:
        print(out$ll)
        print(out$estpars, digits = 2)
        print(out$conv)
        print(out$message)
      }
      
    } else {
      
      # Run MIF to fit to stochastic model:
      tic <- Sys.time()
      
      mf <- try (
        
        mif2(resp_mod,
             params = coef(resp_mod),
             Np = 50,
             Nmif = 5,
             cooling.fraction.50 = 0.5,
             rw.sd = rw_sd(Ri1 = ifelse(time >= 105, 0.02, 0),
                           Ri2 = ifelse(time >= 105, 0.02, 0),
                           rho1 = ifelse(time >= 105, 0.02, 0),
                           rho2 = ifelse(time >= 105, 0.02, 0),
                           theta_lambda1 = ifelse(time >= 105, 0.02, 0),
                           theta_lambda2 = ifelse(time >= 105, 0.02, 0),
                           delta1 = ifelse(time >= 105, 0.02, 0),
                           delta2 = ifelse(time >= 105, 0.02, 0),
                           A1 = ifelse(time >= 105, 0.02, 0),
                           phi1 = ifelse(time >= 105, 0.02, 0),
                           A2 = ifelse(time >= 105, 0.02, 0),
                           phi2 = ifelse(time >= 105, 0.02, 0),
                           k1 = ifelse(time >= 105, 0.02, 0),
                           k2 = ifelse(time >= 105, 0.02, 0),
                           # beta_sd1 = ifelse(time >= 105, 0.02, 0),
                           # beta_sd2 = ifelse(time >= 105, 0.02, 0),
                           t_si_1 = ifelse(time >= 105 & time <= 156, 0.1, 0),
                           t_si_2 = ifelse(time >= 157 & time <= 208, 0.1, 0),
                           t_si_3 = ifelse(time >= 209 & time <= 260, 0.1, 0),
                           t_si_4 = ifelse(time >= 261 & time <= 312, 0.1, 0),
                           t_si_5 = ifelse(time >= 313 & time <= 365, 0.1, 0),
                           t_si_6 = ifelse(time >= 366 & time <= 417, 0.1, 0),
                           t_si_7 = ifelse(time >= 418 & time <= 469, 0.1, 0),
                           t_si_8 = ifelse(time >= 470 & time <= 521, 0.1, 0),
                           t_si_9 = ifelse(time >= 522 & time <= 573, 0.1, 0),
                           t_si_10 = ifelse(time >= 574 & time <= 626, 0.1, 0),
                           w_delta_i_1 = ifelse(time >= 105 & time <= 156, 0.02, 0),
                           w_delta_i_2 = ifelse(time >= 157 & time <= 208, 0.02, 0),
                           w_delta_i_3 = ifelse(time >= 209 & time <= 260, 0.02, 0),
                           w_delta_i_4 = ifelse(time >= 261 & time <= 312, 0.02, 0),
                           w_delta_i_5 = ifelse(time >= 313 & time <= 365, 0.02, 0),
                           w_delta_i_6 = ifelse(time >= 366 & time <= 417, 0.02, 0),
                           w_delta_i_7 = ifelse(time >= 418 & time <= 469, 0.02, 0),
                           w_delta_i_8 = ifelse(time >= 470 & time <= 521, 0.02, 0),
                           w_delta_i_9 = ifelse(time >= 522 & time <= 573, 0.02, 0),
                           w_delta_i_10 = ifelse(time >= 574 & time <= 626, 0.02, 0))) %>%
          mif2(Nmif = 5,
               cooling.fraction.50 = 0.25)# %>%
        # mif2(Nmif = filtering_number,
        #      cooling.fraction.50 = 0.1)
        
      )
      
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- 'mins'
      print(etime)
      
      # If estimation is successful, compile results:
      if (!inherits(mf, 'try-error')) {
        tic <- Sys.time()
        
        ll <- replicate(10, mf %>% pfilter(Np = 5000) %>% logLik()) %>%
          logmeanexp(se = TRUE)
        
        toc <- Sys.time()
        etime <- toc - tic
        units(etime) <- 'mins'
        print(etime)
        
        mf.ll <- mf %>%
          coef() %>%
          bind_rows() %>%
          bind_cols(loglik = ll[1],
                    loglik.se = ll[2])
        
        out <- list(mf = mf,
                    res = mf.ll)
        
        # Print results:
        print(out$res)
        
      }
      
    }
    
  }
  
}

print('Done!')
