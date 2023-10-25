# load packages 
library(tidyverse)
library(rstan)

# Stan set up code
rstan_options (auto_write = TRUE) # to avoid recompulation of unchanged stan programs 
options (mc.cores = parallel::detectCores ()) # run this when implementing locally 

# time series of cases
cases <- results$data %>% dplyr::select(time, v1_obs, v2_obs)

# fixed inputs 
mu  <- 0.0002 
nu  <- 0.0007
N <- 3800000;

gamma1 <- 4.8
gamma2 <- 7.8
sigma1 <- 1
sigma2 <- 4.5

# times
n_weeks <- dim(cases)[1] 
t <- seq(0, n_weeks, by = 1)
t0 = 0 
t <- t[-1]



#initial conditions
E01 <- 0.0001
E02 <- 0.0001
R01 <- 0.4
R02 <- 0.2
R12 <- 0.001

# initialising each compartment
X_SS = round((1.0 - E01 - E02 - R01 - R02 - R12) * N)-2
X_ES = round(E01 * N)
X_IS = 1
X_TS = 0
X_RS = round(R01 * N)
X_SE = round(E02 * N)
X_EE = 0
X_IE = 0
X_TE = 0
X_RE = 0
X_SI = 1
X_EI = 0
X_II = 0
X_TI = 0
X_RI = 0
X_ST = 0
X_ET = 0
X_IT = 0
X_TT = 0
X_RT = 0
X_SR = round(R02 * N)
X_ER = 0
X_IR = 0
X_TR = 0
X_RR = round(R12 * N)

y0 = c(X_SS = X_SS, X_ES = X_ES, X_IS = X_IS, X_TS = X_TS, X_RS = X_RS,
       X_SE = X_SE, X_EE = X_EE, X_IE = X_IE, X_TE = X_TE, X_RE = X_RE,
       X_SI = X_SI, X_EI = X_EI, X_II = X_II, X_TI = X_TI, X_RI = X_RI, 
       X_ST = X_ST, X_ET = X_ET, X_IT = X_IT, X_TT = X_TT, X_RT = X_RT, 
       X_SR = X_SR, X_ER = X_ER, X_IR = X_IR, X_TR = X_TR, X_RR = X_RR)


# data for Stan
data_sir <- list(n_weeks = n_weeks, mu = mu, nu = nu,y0 = y0, t0 = t0, ts = t, N = N, 
                 v1_obs = cases$v1_obs, v2_obs = cases$v2_obs,
                 gamma1 = gamma1, gamma2 = gamma2, 
                 sigma1 = sigma1, sigma2 = sigma2)

# number of MCMC steps
niter <- 2000

# reading in the stan model 
model <- stan_model("seitr_x_seitr.stan")

# draw samples from a stan model 
tic <-  Sys.time()
fit_int_model  <- sampling(model,
                         data = data_sir,
                         iter = niter,
                         chains = 1, 
                         seed = 2908)
tock <-  Sys.time()
time_taken <- tock - tic

  
# parameters to give statistics/diagnostics on  
pars <- c("delta1", "delta2", "theta_lambda1","theta_lambda2", "w1",
          "w2", "Ri1", "Ri2", "R01", "R02", "A1", "A2", "phi1", "phi2",
          "rho1", "rho2")

# stats
print(fit_int_model, pars = pars)

# trace plots
traceplot(fit_int_model, pars = pars) 

# plot marginal posteriors for each chain 
stan_dens(fit_int_model, pars = pars, separate_chains = TRUE)

# number of X_SS each day 
params <- lapply(t, function(i){sprintf("y[%s,1]", i)})



