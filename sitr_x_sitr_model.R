##################################################################################################################
# R code to run pomp model 
# 
# The C code for the pomp model is here: 
# /Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/sitr_x_sitr.cpp
#
# Created by: Sarah Pirikahu 
# Based on code by: Sarah Kramer 
# Creation date: 19 December 
##################################################################################################################

# load libraries
library(tidyverse)
library(pomp)

# read in the c code
mod_code <- readLines('/Volumes/Abt.Domenech/Sarah P/Project 1 - Simulation interaction between influenza and RSV/Analysis/Simulation/sitr_x_sitr.cpp')

# pull out the various components of the C code ready to feed into pomp
# components looking for
components_nm <- c('globs', 'dmeas', 'rmeas', 'rinit', 'skel', 'rsim')
# initialise list
components_l <- vector(mode = 'list', length = length(components_nm))
names(components_l) <- components_nm

for (nm in components_nm) {
  components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
    str_flatten(collapse = '\n')
  components_l[[nm]] <- Csnippet(text = components_l[[nm]])
}


po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
           times = 'time',
           t0 = 0,
           covar = covariate_table(dat[, c('time', 'i_ILI', 'n_T', 'temp', 'ah')], times = 'time'),
           accumvars = c('H1_tot', 'H2_tot', 'H1', 'H2'),
           obsnames = c('n_P1', 'n_P2'),
           statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                          'X_SI', 'X_II', 'X_TI', 'X_RI', 
                          'X_ST', 'X_IT', 'X_TT', 'X_RT',
                          'X_SR', 'X_IR', 'X_TR', 'X_RR', 
                          'H1_tot', 'H2_tot', 
                          'H1', 'H2'),
           paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                          'gamma1', 'gamma2', # 1 / average infectious periods
                          # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                          'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                          'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                          'rho1', 'rho2', # probs. infection leads to ILI consultation
                          'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                          'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                          'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                          'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                          'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                          'N', # population size
                          'I10', 'I20', # props. infectious at outbreak start
                          'R10', 'R20', 'R120'), # props. recovered at outbreak start
           params = c(Ri1 = 1.5, Ri2 = 2,
                      gamma1 = 7 / 5, gamma2 = 7 / 10, # or 4 for flu?
                      # delta = 7 / 5,
                      delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                      theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                      rho1 = 0.5, rho2 = 0.15,
                      alpha = 0, phi = 0,
                      theta_rho1 = 1.0, theta_rho2 = 1.0,
                      eta_temp1 = 0, eta_temp2 = 0,
                      eta_ah1 = 0, eta_ah2 = 0,
                      beta_sd1 = 0, beta_sd2 = 0,
                      N = unique(dat$pop),
                      I10 = 1e-5, I20 = 1e-5,
                      R10 = 0, R20 = 0, R120 = 0),
           globals = components_l[['globs']],
           dmeasure = components_l[['dmeas']],
           rmeasure = components_l[['rmeas']],
           skeleton = vectorfield(components_l[['skel']]),
           rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
           partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
           rinit = components_l[['rinit']]
)
