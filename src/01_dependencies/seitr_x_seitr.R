# ---------------------------------------------------------------------------------------------------------------------
# Functions to assist with interaction model
# ---------------------------------------------------------------------------------------------------------------------

create_SEITRxSEITR_mod <- function(n_weeks, start_time, parms, debug_bool = FALSE) {
  
  # Function to create pomp object
  # param n_weeks: Number of weeks to run model
  # param start_time: How many weeks of burn-in to run?
  # param parms: True values of all model parameters
  # param debug_bool: Should information for debugging be printed?
  # returns: pomp object with all parameters set to "true" values
  
  # Read model C code:
  mod_code <- readLines('src/01_dependencies/seitr_x_seitr.c')
  
  components_nm <- c('globs', 'toest', 'fromest', 'dmeas', 'rmeas', 'rinit', 'skel', 'rsim')
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for (nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
      str_flatten(collapse = '\n')
    
    if(nm == 'globs') {
      components_l[[nm]] <- paste(components_l[[nm]],
                                  sprintf('static int debug = %d;',
                                          as.integer(debug_bool)),
                                  sep = '\n')
    }
    
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  # Create pomp object:
  po <- pomp(data = data.frame(time = seq(from = start_time, to = n_weeks, by = 1), V1_obs = NA, V2_obs = NA),
             times = "time",
             t0 = start_time,
             obsnames = c('V1_obs', 'V2_obs'),
             accumvars = c('V1', 'V2'),
             statenames = c('X_SS', 'X_ES' , 'X_IS', 'X_TS', 'X_RS', 
                            'X_SE', 'X_EE', 'X_IE', 'X_TE', 'X_RE',
                            'X_SI', 'X_EI' ,'X_II', 'X_TI', 'X_RI', 
                            'X_ST', 'X_ET' ,'X_IT', 'X_TT', 'X_RT',
                            'X_SR', 'X_ER' ,'X_IR', 'X_TR', 'X_RR', 
                            'V1', 'V2'),
             paramnames = names(parms),
             params = parms,
             partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
             globals = components_l[['globs']],
             dmeasure = components_l[['dmeas']],
             rmeasure = components_l[['rmeas']],
             rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
             skeleton = vectorfield(components_l[['skel']]), # putting in deterministic for testing
             rinit = components_l[['rinit']]
  )
  
  # Return pomp object:
  return(po)
  
}


check_transformations <- function(pomp_object) {
  # Function to check that parameters are correctly being transformed
  # params pomp_object: The pomp model object to be checked
  
  obj <- pomp_object
  x <- coef(obj, transform = TRUE)
  obj1 <- obj
  coef(obj1, transform = TRUE) <- x
  
  expect_true(all.equal(coef(obj), coef(obj1)),
              info = 'Parameters not correctly transformed')
}


check_correct_N_CONST <- function(pomp_object, true_n, n_sim = 10) {
  # Function to check that deterministic and stochastic simulations maintain correct N (when pop size constant)
  # params pomp_object: The pomp model object to be checked
  # params true_n: The expected population size
  # params n_sim: The number of stochastic simulations to run; if 0, don't run stochastic model
  
  sim_determ <- trajectory(object = pomp_object, format = 'data.frame') %>%
    rowwise() %>%
    mutate(Ncheck = sum(c_across(contains('X'))),
           Ntrue = true_n)
  expect_true(all.equal(sim_determ$Ncheck, sim_determ$Ntrue))
  
  if (n_sim > 0) {
    sim_stoch <- simulate(object = pomp_object, nsim = n_sim, format = 'data.frame') %>%
      rowwise() %>%
      mutate(Ncheck = sum(c_across(X_SS:X_RR)),
             Ntrue = true_n)
    expect_true(all.equal(sim_stoch$Ncheck, sim_stoch$Ntrue))
  }
  
}


check_independent_dynamics <- function(pomp_object) {
  # Function to check that, when no interaction is specified, virus dynamics are independent
  # params pomp_object: The pomp model object to be checked
  # returns: Plot of epidemic dynamics with and without interaction
  
  p_mat <- parmat(params = coef(pomp_object), nrep = 3)
  p_mat['E01', ] <- c(1e-3, 1e-3, 0)
  p_mat['E02', ] <- c(1e-3, 0, 1e-3)
  p_mat['theta_lambda1', ] <- c(1.0, 1.0, 1.0)
  p_mat['theta_lambda2', ] <- c(1.0, 1.0, 1.0)
  
  sim_determ <- trajectory(object = pomp_object,
                           t0 = 0, times = 0:522,
                           params = p_mat,
                           format = 'data.frame') %>%
    pivot_longer(cols = -c('time', '.id'), names_to = 'var_nm', values_to = 'val') %>%
    filter(var_nm %in% c('V1', 'V2'))
  
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'V1'],
                        sim_determ$val[sim_determ$.id == 2 & sim_determ$var_nm == 'V1']))
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'V2'],
                        sim_determ$val[sim_determ$.id == 3 & sim_determ$var_nm == 'V2']))
  
  p_temp <- ggplot(data = sim_determ %>% filter(.id == 1), aes(x = time, y = 100 * val)) +
    geom_point(aes(color = var_nm)) +
    geom_line(data = sim_determ %>% filter(.id == 2, var_nm == 'V1'), color = 'pink') +
    geom_line(data = sim_determ %>% filter(.id == 3, var_nm == 'V2'), color = 'purple') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)') +
    theme_classic()
  
  return(p_temp)
}
