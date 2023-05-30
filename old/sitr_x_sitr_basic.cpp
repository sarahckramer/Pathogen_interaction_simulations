#include <Rcpp.h>
using namespace Rcpp;

// Implementation of the SITR x SITR model for circulation of two respiratory 
// viruses
// simplifying first to do without seasonality in transmission and without incorporation of 
// under reporting in measurement model 

// Created by: Sarah Pirikahu
// Creation date: 22 December 2022


//start_globs
// specifying global functions required to make the code run 

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t)
static double pTrans(double r, double delta_t) {
  
  // r: event (instantaneous rate)
  // delta_t: time step
  // Returns: probability of transition
  double out = 1.0 - exp(-r * delta_t);
  return out;
}
//end_globs

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
//if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
//  Rprintf("SS=%f, IS=%f, RS=%f, SI=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_IS, X_RS, X_SI, X_SR, X_RR, X_SS + X_IS + X_RS + X_SI + X_SR + X_RR, N);
//  X_SS = nearbyint(N - X_IS - X_RS - X_SI - X_SR - X_RR);
//}

//end_rinit

//start_dmeas
// initalisation
double fP1, fP2, ll;

// calculating components for the likelihood 
fP1 = dbinom(H1_obs, H1, rho1, 1); // First likelihood component, natural scale
fP2 = dbinom(H2_obs, H2, rho2, 1); // Second likelihood component, natural scale

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(fP1 + fP2, -1e3);

// calculate likelihood by back transforming the log likelihood 
lik = (give_log) ? ll : exp(ll);
//end_dmeas


//start_rmeas
//generate total number of specimens tested for each virus where r1 probability
// of someone been sick and turning up for testing 
//n_T = rbinom(N, rho2);

// generate the total number of tests positive to each virus 
H1_obs = rbinom(H1, rho1); // virus 1
H2_obs = rbinom(H2, rho2); // virus 2
//end_rmeas

//start_skel
// calculating the prevelence of each infection 
double p1 = (X_IS + X_II + X_IT + X_IR) ; // virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) ; // virus 2

// calculate the transmission rate for each virus taking into consideration climate forcing; eq (1) Kramer (2023)
// note: beta_i = Reff*gamma where Reff is the effective reproductive number at time i in a partially susceptible population  
double beta1 = Ri1 / (1.0 - (R10 + R120)) * gamma1; // virus 1 
double beta2 = Ri2 / (1.0 - (R20 + R120)) * gamma2; // virus 2 

// calculate force of infection for each virus 
double lambda1 = beta1 * (p1/N); // virus 1
double lambda2 = beta2 * (p2/N); // virus 2

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
DH1 = gamma1 * (X_IS + theta_rho2 * X_II + X_IT + X_IR) ; // virus 1
DH2 = gamma2 * (X_SI + theta_rho1 * X_II + X_TI + X_RI) ; // virus 2
//end_skel

//start_rsim
// calculate prevalence of each virus 
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) / N; // virus 2

// calculate basic reproductive number R0 for each virus taking into consideration 
double R0_1 = Ri1 / (1.0 - (R10 + R120)); // virus 1
double R0_2 = Ri2 / (1.0 - (R20 + R120)); // virus 2

// initialisation of the transmission terms 
double beta1, beta2;
// incorporate extra demographic stochasticity with the gamma distributed white noise process
// dt is the time step hence it is not defined but is a variable created within pomp
if (p1 > 0.0 && beta_sd1 > 0.0) { 
  beta1 = rgammawn(sqrt(R0_1 / (p1 * N * beta_sd1 * dt)), R0_1 * gamma1);
} else {
  beta1 = R0_1 * gamma1;
}
if (p2 > 0.0 && beta_sd2 > 0.0) {
  beta2 = rgammawn(sqrt(R0_2 / (p2 * N * beta_sd2 * dt)), R0_2 * gamma2);
} else {
  beta2 = R0_2 * gamma2;
}

// calculate force of infection for each virus 
double lambda1 = beta1 * p1; // virus 1
double lambda2 = beta2 * p2; // virus 2

// Calculate duration of refractory period of virus 2:
double delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// initalising transitions 
double rates[18];// vector of length 18
double fromSS[2], fromIS[2], fromTS[2], fromSI[2], fromII[2], fromTI[2], fromST[2], fromIT[2], fromTT[2]; // vectors of length 2
double fromRS, fromRI, fromRT, fromSR, fromIR, fromTR; 

// specifying rate for transition - note: vector indexing starts at 0 for C++ rather than 1 like R 
rates[0] = lambda1; // stochastic force of infection virus 1 (X_SS -> X_IS)
rates[1] = lambda2; // stochastic force of infection virus 2 (X_SS -> X_SI)
rates[2] = gamma1; // rate of ending latent stage for virus 1 whilst susceptible to virus 2 (X_IS -> X_TS)
rates[3] = theta_lambda1 * lambda2; // rate of infection for virus 2 when already infected with virus 1 (X_IS -> X_II)
rates[4] = delta1; // rate of recovery virus 1 whilst susceptible to virus 2 (X_TS -> X_RS) 
rates[5] = theta_lambda1 * lambda2; // rate of infection with virus 2 whilst in the latent stage of virus 1 infection (X_TS -> X_TI)
rates[6] = theta_lambda2 * lambda1; // rate of infection with virus 1 whilst infected with virus 2 (X_SI -> X_II)
rates[7] = gamma2; // rate of ending latent stage for virus 2 when susceptible to virus 1 (X_SI -> X_ST)
rates[8] = gamma1; // rate of ending latent stage for virus 1 whilst co-infected with virus 2 (X_II -> X_TI)
rates[9] = gamma2; // rate of ending latent stage for virus 2 whilst co-infected with virus 1 (X_II -> X_IT)
rates[10] = delta1; // rate of recovery virus 2 whilst infected with virus 1 (X_TI -> X_RI)
rates[11] = gamma2; // rate of ending latent stage for virus 2 whilst co-infected with virus 1 (X_II -> X_IT)
rates[12] = theta_lambda2 * lambda1; // rate of infection with virus 1 whilst in latent infection with virus 2 (X_ST -> X_IT)
rates[13] = delta2; // rate of recovery for virus 2 whilst susceptible to virus 1 (X_ST -> X_SR)
rates[14] = gamma1; // rate of ending latent stage for virus 1 whilst in latent infection with virus 2 (X_IT -> X_TT)
rates[15] = delta2; // recovery rate for virus 2 which co-infected with virus 1  (X_IT -> X_IR)
rates[16] = delta1; // recovery rate of virus 1 whilst in latent infection with virus 2 (X_TT -> X_RT) 
rates[17] = delta2; // recovery rate of virus 2 whilst in latent infection with virus 1 (X_TT -> X_TR)

// drawing sample for each of the compartments to determine who is transitioning.
// This is drawn from the Euler-multinomial distribution and the function
// returns a length(rate[i]) by n matrix where in our case we have 2 columns c1 which we let represent 
// transitions due to virus 1 (i.e. horizontally across compartments) and c2 the transitions
// due to virus 2 (i.e. vertically down compartments)
reulermultinom(2, X_SS, &rates[0], dt, &fromSS[0]);
reulermultinom(2, X_IS, &rates[2], dt, &fromIS[0]);
reulermultinom(2, X_TS, &rates[4], dt, &fromTS[0]);
reulermultinom(2, X_SI, &rates[6], dt, &fromSI[0]);
reulermultinom(2, X_II, &rates[8], dt, &fromII[0]);
reulermultinom(2, X_TI, &rates[10], dt, &fromTI[0]);
reulermultinom(2, X_ST, &rates[12], dt, &fromST[0]);
reulermultinom(2, X_IT, &rates[14], dt, &fromIT[0]);
reulermultinom(2, X_TT, &rates[16], dt, &fromTT[0]);

// drawing samples for each of the recovered compartments from binomial distributions  
fromRS = rbinom(X_RS, pTrans(lambda2, dt));
fromRI = rbinom(X_RI, pTrans(gamma2, dt));
fromRT = rbinom(X_RT, pTrans(delta2, dt));
fromSR = rbinom(X_SR, pTrans(lambda1, dt));
fromIR = rbinom(X_IR, pTrans(gamma1, dt));
fromTR = rbinom(X_TR, pTrans(delta1, dt));

Rprintf("first=%.1f, second=%.f, first=%.1f, second=%.f\n", fromSS[0], fromSS[1], fromIS[0], fromIS[1]);

// balance equations
X_SS += -fromSS[0] - fromSS[1];
X_IS += fromSS[0] - fromIS[0] - fromIS[1];
X_TS += fromIS[0] - fromTS[0] - fromTS[1];
X_RS += fromTS[0] - fromRS;

X_SI += fromSS[1] - fromSI[0] - fromSI[1];
X_II += fromIS[1] + fromSI[0] - fromII[0] - fromII[1];
X_TI += fromTS[1] + fromII[0] - fromTI[0] - fromTI[1];
X_RI += fromRS + fromTI[0] - fromRI;

X_ST += fromSI[1] - fromST[0] - fromST[1];
X_IT += fromII[1] + fromST[0] - fromIT[0] - fromIT[1];
X_TT += fromTI[1] + fromIT[0] - fromTT[0] - fromTT[1];
X_RT += fromRI + fromTT[0] - fromRT;

X_SR += fromST[1] - fromSR;
X_IR += fromIT[1] + fromSR - fromIR;
X_TR += fromTT[1] + fromIR - fromTR;
X_RR += fromRT + fromTR;
  
H1_tot += (fromIS[0] + fromII[0] + fromIT[0] + fromIR);
H2_tot += (fromSI[1] + fromII[1] + fromTI[1] + fromRI);
H1 += (fromIS[0] + theta_rho2 * fromII[0] + fromIT[0] + fromIR);
H2 += (fromSI[1] + theta_rho1 * fromII[1] + fromTI[1] + fromRI);
//end_rsim
