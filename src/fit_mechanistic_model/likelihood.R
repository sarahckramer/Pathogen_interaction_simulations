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
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)

data_id <- (jobid - 1) %% 10 + 1; print(data_id)
jobid <- ceiling(jobid / 10); print(jobid)

# Set relevant parameters:
sobol_size <- 500
# search_type <- 'round1_CIs'

n_cores <- 50 # Number of cores to use
time_max <- 23.75 # Maximal execution time (in hours)

# ---------------------------------------------------------------------------------------------------------------------

# Get model and data

# Create model and generate synthetic data:
source('src/functions_etc/generate_data.R')

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

# ---------------------------------------------------------------------------------------------------------------------

# Setup fitting exercise

# Choose parameters to be estimated:
estpars <- c('Ri1', 'Ri2', 'rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'delta2',
             'A1', 'phi1', 'A2', 'phi2', 'k1', 'k2', 'E01', 'E02', 'R01', 'R02', 'R012',
             't_si_1', 't_si_2', 't_si_3', 't_si_4', 't_si_5', 't_si_6',
             't_si_7', 't_si_8', 't_si_9', 't_si_10',
             'w_delta_i_1', 'w_delta_i_2', 'w_delta_i_3', 'w_delta_i_4', 'w_delta_i_5', 'w_delta_i_6',
             'w_delta_i_7', 'w_delta_i_8', 'w_delta_i_9', 'w_delta_i_10')

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
                          E01 = c(0, 1e-3),
                          E02 = c(0, 1e-3),
                          R01 = c(0, 0.3),
                          R02 = c(0, 0.3),
                          R012 = c(0, 0.3),
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
if (search_type == 'round1_CIs') {
  
  start_range <- read_rds('results/trajectory_matching/round2CI_startvals.rds') %>%
    filter(int_set == jobid, .id == data_id) %>%
    ungroup() %>%
    select(-c(int_set, .id, minmax)) %>%
    as.data.frame()
  
  expect_true(all(estpars %in% names(start_range)))
  expect_true(all(names(start_range) %in% estpars))
  
} else if (search_type == 'round2_CIs') {
  
  start_range <- read_rds('results/trajectory_matching/round3CI_startvals.rds') %>%
    filter(int_set == jobid, .id == data_id) %>%
    ungroup() %>%
    select(-c(int_set, .id, minmax)) %>%
    as.data.frame()
  
  expect_true(all(estpars %in% names(start_range)))
  expect_true(all(names(start_range) %in% estpars))
  
} else if (search_type == 'round3_CIs') {
  
  start_range <- read_rds('results/trajectory_matching/round4CI_startvals.rds') %>%
    filter(int_set == jobid, .id == data_id) %>%
    ungroup() %>%
    select(-c(int_set, .id, minmax)) %>%
    as.data.frame()
  
  expect_true(all(estpars %in% names(start_range)))
  expect_true(all(names(start_range) %in% estpars))
  
}

set.seed(38564478)
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

if (search_type %in% c('round1_CIs', 'round2_CIs', 'round3_CIs')) {
  
  start_values <- start_values %>%
    mutate(phi1 = if_else(phi1 > 52.25, phi1 - 52.25, phi1),
           phi2 = if_else(phi2 > 52.25, phi2 - 52.25, phi2))
  
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
    
  }
  
}

print('Done!')
