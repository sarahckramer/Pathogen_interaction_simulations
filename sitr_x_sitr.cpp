#include <Rcpp.h>
using namespace Rcpp;

// Implementation of the SITR x SITR model for ciculation of two respiratory 
// viruses

// Created by: Sarah Pirikahu
// Based on code by: Sarah Kramer (https://github.com/sarahckramer/resp_virus_interactions/blob/main/src/resp_interaction_model.c)
// Creation date: 19 December 2022

// for the pomp model we want to have inputs for: 
// rinit - inital conditions
// rmeasure - simulation of f_{Y_n|X_n} for measurement model
// dmeasure - evaluation of f_{Y_n|X_n} for measurement model 
// rprocess - simulation of f_{X_n|X_{n-1}} for the process model (specified in rsim)
// dprocess - evaluation of f_{X_n|X_{n-1}} for the process model (specified in skeleton - skel) 

//start_rinit
// initalising each compartment 
X_SS = nearbyint((1.0 - I10 - I20 - R10 - R20 - R120) * N);
X_IS = nearbyint(I10 * N);
X_TS = 0;
X_RS = nearbyint(R10 * N);
X_SI = nearbyint(I20 * N);
X_II = 0;
X_TI = 0;
X_RI = 0;
X_ST = 0;
X_IT = 0;
X_TT = 0;
X_RT = 0;
X_SR = nearbyint(R20 * N);
X_IR = 0;
X_TR = 0;
X_RR = nearbyint(R120 * N);

// initalising accumulator variables 
H1_tot = 0; // takes into consideration interaction  
H2_tot = 0; 
H1 = 0;     // does not take into consideration the interaction term 
H2 = 0; 

// if the compartments that are assigned some inital value don't sum to
// N (likely due to rounding in the nearbyint function) print the values
// then re calculate the inital value of X_SS by doing N - all other
// inital values for the other compartments
if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
  Rprintf("SS=%f, IS=%f, RS=%f, SI=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_IS, X_RS, X_SI, X_SR, X_RR, X_SS + X_IS + X_RS + X_SI + X_SR + X_RR, N);
  X_SS = nearbyint(N - X_IS - X_RS - X_SI - X_SR - X_RR);
}

// for debugging - printing all variables to figure out what is going on 
if(debug) {
  Rprintf("%f, %f, %f, %f, %f, %f, %f\n", Ri1, Ri2, I10, I20, R10, R20, R120);
}
//end_rinit

//start_dmeas
// initalisation
double fP1, fP2, ll;
// defining the period 
double omega = (2 * M_PI) / 52.25;

// probability of detecting each virus; eqn (2) p24. Kramer (2023)
double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // virus 2

// if probability of detection is less than 0 make it 0.
if (rho1_w < 0) {  // for virus 1 
  rho1_w = 0.0;
}
if (rho2_w < 0) {  // for virus 2 
  rho2_w = 0.0;
}

// calculating components for the likelihood 
fP1 = dbinom(n_P1, n_T, rho1_w, 1); // First likelihood component, natural scale
fP2 = dbinom(n_P2, n_T, rho2_w, 1); // Second likelihood component, natural scale

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(fP1 + fP2, -1e3);

// If data are NA, ll will be NA; in this case, set to zero
ll = ISNA(ll) ? 0.0 : ll;

// for debugging - printing all variables to figure out what is going on 
if(debug) {
  Rprintf("t=%.1f, rho1_w=%.1f, rho2_w=%.1f, n_T=%.1f, fP1=%.1f, fP2=%.1f, sum=%.1f, ll=%.f\n", t, rho1_w, rho2_w, n_T, fP1, fP2, fP1 + fP2, ll);
} 

// calculate likelihood by back transforming the log likelihood 
lik = (give_log) ? ll : exp(ll);
//end_dmeas


//start_rmeas
// defining the period 
double omega = (2 * M_PI) / 52.25;

// probability of detecting each virus
double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); //  virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); //  virus 2

// if probability of detection is less than 0 make it 0
if (rho1_w < 0) {  // for virus 1 
  rho1_w = 0.0;
}
if (rho2_w < 0) {  // for virus 2 
  rho2_w = 0.0;
}

// generate the total number of tests positive to each virus 
n_P1 = rbinom(n_T, rho1_w); // virus 1
n_P2 = rbinom(n_T, rho2_w); // virus 2
//end_rmeas

//start_skel
// calculating the prevelence of each infection 
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) / N; // virus 2

// estimating the transmission rate for each virus taking into consideration climate forcing; eq (1) Kramer (2023)
// note: beta_i = Reff*gamma where Reff is the effective reproductive number at time i in a partially susceptible population  
double beta1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp) * gamma1; // virus 1 
double beta2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp) * gamma2; // virus 2 

// estimating force of infection for each virus 
double lambda1 = beta1 * p1; // virus 1
double lambda2 = beta2 * p2; // virus 2

// relative parameter to describe the rate of cross protection acheived after infection 
double delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// ODEs
DX_SS = -(lambda1 + lambda2) * X_SS; 
DX_IS = lambda1 * X_SS - (gamma1 + theta_lambda1 * lambda2) * X_IS; 
DX_TS = gamma1 * X_IS - (delta1 + theta_lambda1 * lambda2) * X_TS;
DX_RS = delta1 * X_TS - lambda2 * X_RS; 
DX_SI = lambda2 * X_SS - (theta_lambda2 * lambda1 + gamma2) * X_SI;
DX_II = theta_lambda1 * lambda2 * X_IS + theta_lambda2 * lambda1 * X_SI - (gamma1 + gamma2) * X_II; 
DX_TI = theta_lambda1 * lambda2 * X_TS + gamma1 * X_II - (delta1 + gamma2) * X_TI;
DX_RI = lambda2 * X_RS + delta1 * X_TI - gamma2 * X_RI; 
DX_ST = gamma2 * X_SI - (theta_lambda2 * lambda1 + delta2) * X_ST; 
DX_IT = gamma2 * X_II + theta_lambda2 * lambda1 * X_ST - (gamma1 + delta2) * X_IT; 
DX_TT = gamma2 * X_TI + gamma1 * X_IT - (delta1 + delta2) * X_TT;
DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2 * X_RT;
DX_SR = delta2 * X_ST - lambda1 * X_SR; 
DX_IR = delta2 * X_IT + lambda1 * X_SR - gamma1 * X_IR; 
DX_TR = delta2 * X_TT + gamma1 * X_IR - delta1 * X_TR; 
DX_RR = delta2 * X_RT + delta1 * X_TR;

// defining the accumulator variables
// incidence rates of infection overall for each virus not taking into consideration interaction 
DH1_tot = gamma1 * p1; // virus 1 
DH2_tot = gamma2 * p2; // virus 2 
// incidence rates of each virus taking into consideration interaction
DH1 = gamma1 * (X_IS + theta_rho2 * X_II + X_IT + X_IR) / N; // virus 1
DH2 = gamma2 * (X_SI + theta_rho1 * X_II + X_TI + X_RI) / N; // virus 2
//end_skel

//start_rsim


//end_rsim
