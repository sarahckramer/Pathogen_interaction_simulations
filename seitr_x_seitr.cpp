#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------------------------------//
// Implementation of the SEITR x SEITR model for circulation of two respiratory 
// viruses

// For a pomp model we want to have inputs for: 
// rinit - inital conditions
// rmeasure - simulation of f_{Y_n|X_n} for measurement model
// dmeasure - evaluation of f_{Y_n|X_n} for measurement model 
// rprocess - simulation of f_{X_n|X_{n-1}} for the process model (specified in rsim)

// Created by: Sarah Pirikahu
// Creation date: 22 December 2022

// -----------------------------------------------------------------------------------------------------//

//------- GLOBAL FUNCTIONS -----//
//start_globs

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t). This function transforms our transmission, recovery, etc rates 
// into probabilies for a particular time step delta_t
static double pTrans(double r, double delta_t) {
  
  // r: event (instantaneous rate)
  // delta_t: time step
  // Returns: probability of transition
  double out = 1.0 - exp(-r * delta_t);
  return out;
}

// Scaled logit transform for parameters constrained in interval [a,b]
static double logitCons(double x, double a, double b) {
  x = (x <= a) ? a : (x >= b ? b : x);
  double out = log((x - a) / (b - x));
  return out;
}

// Expit transform for parameters constrained in interval [a,b]
// this is the back transform of the logit function 
static double expitCons(double x, double a, double b) {
  double out = (a + b * exp(x)) / (1.0 + exp(x));
  if(ISNAN(out) | isinf(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
  return out;
}
//end_globs


//----- TRANSFORMATIONS ---//
// here we specify the transformations applied to each parameter to
// constrain the integrator search space. We also specify the corresponding
// back transformation

//start_toest
T_Ri1 = logitCons(Ri1, 1.0, 10); 
T_Ri2 = logitCons(Ri2, 1.0, 10); 

T_sigma1 = sigma1;
T_sigma2 = sigma2;
T_gamma1 = gamma1;
T_gamma2 = gamma2;

T_delta1 = log(delta1);
T_delta2 = log(delta2);

// T_mu = mu;
// T_nu = nu;

T_w1 = log(w1);
T_w2 = log(w2);

T_rho1 = logit(rho1);
T_rho2 = logit(rho2);
T_theta_lambda1 = logitCons(theta_lambda1,0,10);
T_theta_lambda2 = logitCons(theta_lambda2,0,10);

T_A1 = logit(A1);
T_phi1 = logitCons(phi1,20,42); // 4 month period ofter Oct
T_A2 = logit(A2);
T_phi2 = logitCons(phi2,20,42);

T_beta_sd1 = beta_sd1;
T_beta_sd2 = beta_sd2;
T_N = N;
T_nsurges = nsurges;

// int *t_vec = (int *) &t_si_1;
// int *T_t_vec = (int *) &t_si_1;
// for (int i = 0; i < (int) n_surge; i++) {
//  T_t_vec[i] = t_vec[i];
// }

// double T_delta_i_1[n_surge]; 
// double delta_i_1[n_surge];   
// for (int i = 0; i < n_surge; i++) {
//   T_delta_i_1[i] = log(delta_i_1[i]);
// }

T_t_si_1 = t_si_1;
T_t_si_2 = t_si_2;
T_t_si_3 = t_si_3;
T_t_si_4 = t_si_4;
T_t_si_5 = t_si_5;

T_t_si_6 = t_si_6;
T_t_si_7 = t_si_7;
T_t_si_8 = t_si_8;
T_t_si_9 = t_si_9;
T_t_si_10 = t_si_10;

T_t_si_11 = t_si_11;
T_t_si_12 = t_si_12;
T_t_si_13 = t_si_13;
T_t_si_14 = t_si_14;
T_t_si_15 = t_si_15;

T_t_si_16 = t_si_16;
T_t_si_17 = t_si_17;
T_t_si_18 = t_si_18;
T_t_si_19 = t_si_19;
T_t_si_20 = t_si_20;
// 
// T_t_si_21 = t_si_21;
// T_t_si_22 = t_si_22;
// T_t_si_23 = t_si_23;
// T_t_si_24 = t_si_24;
// T_t_si_25 = t_si_25;
// 
// T_t_si_26 = t_si_26;
// T_t_si_27 = t_si_27;
// T_t_si_28 = t_si_28;
// T_t_si_29 = t_si_29;
// T_t_si_30 = t_si_30;
// 
// T_t_si_31 = t_si_31;
// T_t_si_32 = t_si_32;
// T_t_si_33 = t_si_33;
// T_t_si_34 = t_si_34;
// T_t_si_35 = t_si_35;
// 
// T_t_si_36 = t_si_36;
// T_t_si_37 = t_si_37;
// T_t_si_38 = t_si_38;
// T_t_si_39 = t_si_39;
// T_t_si_40 = t_si_40;

T_delta_i_1 = log(delta_i_1);
T_delta_i_2 = log(delta_i_2);
T_delta_i_3 = log(delta_i_3);
T_delta_i_4 = log(delta_i_4);
T_delta_i_5 = log(delta_i_5);

T_delta_i_6 = log(delta_i_6);
T_delta_i_7 = log(delta_i_7);
T_delta_i_8 = log(delta_i_8);
T_delta_i_9 = log(delta_i_9);
T_delta_i_10 = log(delta_i_10);

T_delta_i_11 = log(delta_i_11);
T_delta_i_12 = log(delta_i_12);
T_delta_i_13 = log(delta_i_13);
T_delta_i_14 = log(delta_i_14);
T_delta_i_15 = log(delta_i_15);

T_delta_i_16 = log(delta_i_16);
T_delta_i_17 = log(delta_i_17);
T_delta_i_18 = log(delta_i_18);
T_delta_i_19 = log(delta_i_19);
T_delta_i_20 = log(delta_i_20);
// 
// T_delta_i_21 = log(delta_i_21);
// T_delta_i_22 = log(delta_i_22);
// T_delta_i_23 = log(delta_i_23);
// T_delta_i_24 = log(delta_i_24);
// T_delta_i_25 = log(delta_i_25);
// 
// T_delta_i_6 = log(delta_i_26);
// T_delta_i_7 = log(delta_i_27);
// T_delta_i_8 = log(delta_i_28);
// T_delta_i_9 = log(delta_i_29);
// T_delta_i_10 = log(delta_i_30);
// 
// T_delta_i_31 = log(delta_i_31);
// T_delta_i_32 = log(delta_i_32);
// T_delta_i_33 = log(delta_i_33);
// T_delta_i_34 = log(delta_i_34);
// T_delta_i_35 = log(delta_i_35);
// 
// T_delta_i_36 = log(delta_i_36);
// T_delta_i_37 = log(delta_i_37);
// T_delta_i_38 = log(delta_i_38);
// T_delta_i_39 = log(delta_i_39);
// T_delta_i_40 = log(delta_i_40);

// we need to specify a transform on the sum of E and R
// such that the compartments are not allowed to turn
// negative
double sum_init = 0.0; // initialise
sum_init = E01 + E02 + R01 + R02 + R12;
T_E01 = log(E01 / (1.0 - sum_init));
T_E02 = log(E02 / (1.0 - sum_init));
T_R01 = log(R01 / (1.0 - sum_init));
T_R02 = log(R02 / (1.0 - sum_init));
T_R12 = log(R12 / (1.0 - sum_init));
//end_toest

//start_fromest
Ri1 = expitCons(T_Ri1, 1.0, 10);
Ri2 = expitCons(T_Ri2, 1.0, 10);
sigma1 = T_sigma1;
sigma2 = T_sigma2;
gamma1 = T_gamma1;
gamma2 = T_gamma2;

delta1 = exp(T_delta1);
delta2 = exp(T_delta2);

// mu = T_mu;
// nu = T_nu;

w1 = exp(T_w1);
w2 = exp(T_w2);

rho1 = expit(T_rho1);
rho2 = expit(T_rho2);
theta_lambda1 = expitCons(T_theta_lambda1,0,10);
theta_lambda2 = expitCons(T_theta_lambda2,0,10);

A1 = expit(T_A1);
phi1 = expitCons(T_phi1,20,42);
A2 = expit(T_A2);
phi2 = expitCons(T_phi2,20,42);
 
beta_sd1 = T_beta_sd1;
beta_sd2 = T_beta_sd2;
N = T_N;
nsurges = T_nsurges;

//int *t_vec = (int *) &t_si_1;
//int *T_t_vec = (int *) &t_si_1;
// for (int i = 0; i < (int) n_surge; i++) {
//   t_vec[i] = T_t_vec[i];
// }

// double T_delta_i_1[n_surge]; 
// double delta_i_1[n_surge];   
// for (int i = 0; i < n_surge; i++) {
//   delta_i_1[i] = exp(T_delta_i_1[i]);
// }

t_si_1 = T_t_si_1;
t_si_2 = T_t_si_2;
t_si_3 = T_t_si_3;
t_si_4 = T_t_si_4;
t_si_5 = T_t_si_5;

t_si_6 = T_t_si_6;
t_si_7 = T_t_si_7;
t_si_8 = T_t_si_8;
t_si_9 = T_t_si_9;
t_si_10 = T_t_si_10;

t_si_11 = T_t_si_11;
t_si_12 = T_t_si_12;
t_si_13 = T_t_si_13;
t_si_14 = T_t_si_14;
t_si_15 = T_t_si_15;

t_si_16 = T_t_si_16;
t_si_17 = T_t_si_17;
t_si_18 = T_t_si_18;
t_si_19 = T_t_si_19;
t_si_20 = T_t_si_20;
// 
// t_si_21 = T_t_si_21;
// t_si_22 = T_t_si_22;
// t_si_23 = T_t_si_23;
// t_si_24 = T_t_si_24;
// t_si_25 = T_t_si_25;
// 
// t_si_26 = T_t_si_26;
// t_si_27 = T_t_si_27;
// t_si_28 = T_t_si_28;
// t_si_29 = T_t_si_29;
// t_si_30 = T_t_si_30;
// 
// t_si_31 = T_t_si_31;
// t_si_32 = T_t_si_32;
// t_si_33 = T_t_si_33;
// t_si_34 = T_t_si_34;
// t_si_35 = T_t_si_35;
// 
// t_si_36 = T_t_si_36;
// t_si_37 = T_t_si_37;
// t_si_38 = T_t_si_38;
// t_si_39 = T_t_si_39;
// t_si_40 = T_t_si_40;

delta_i_1 = exp(T_delta_i_1);
delta_i_2 = exp(T_delta_i_2);
delta_i_3 = exp(T_delta_i_3);
delta_i_4 = exp(T_delta_i_4);
delta_i_5 = exp(T_delta_i_5);

delta_i_6 = exp(T_delta_i_6);
delta_i_7 = exp(T_delta_i_7);
delta_i_8 = exp(T_delta_i_8);
delta_i_9 = exp(T_delta_i_9);
delta_i_10 = exp(T_delta_i_10);

delta_i_11 = exp(T_delta_i_11);
delta_i_12 = exp(T_delta_i_12);
delta_i_13 = exp(T_delta_i_13);
delta_i_14 = exp(T_delta_i_14);
delta_i_15 = exp(T_delta_i_15);

delta_i_16 = exp(T_delta_i_16);
delta_i_17 = exp(T_delta_i_17);
delta_i_18 = exp(T_delta_i_18);
delta_i_19 = exp(T_delta_i_19);
delta_i_20 = exp(T_delta_i_20);
// 
// delta_i_21 = exp(T_delta_i_21);
// delta_i_22 = exp(T_delta_i_22);
// delta_i_23 = exp(T_delta_i_23);
// delta_i_24 = exp(T_delta_i_24);
// delta_i_25 = exp(T_delta_i_25);
// 
// delta_i_26 = exp(T_delta_i_26);
// delta_i_27 = exp(T_delta_i_27);
// delta_i_28 = exp(T_delta_i_28);
// delta_i_29 = exp(T_delta_i_29);
// delta_i_30 = exp(T_delta_i_30);
// 
// delta_i_31 = exp(T_delta_i_31);
// delta_i_32 = exp(T_delta_i_32);
// delta_i_33 = exp(T_delta_i_33);
// delta_i_34 = exp(T_delta_i_34);
// delta_i_35 = exp(T_delta_i_35);
// 
// delta_i_36 = exp(T_delta_i_36);
// delta_i_37 = exp(T_delta_i_37);
// delta_i_38 = exp(T_delta_i_38);
// delta_i_39 = exp(T_delta_i_39);
// delta_i_40 = exp(T_delta_i_40);

 
double sum_init = 0.0;
sum_init = exp(E01) + exp(E02) + exp(R01) + exp(R02) + exp(R12);
E01 = exp(T_E01) / (1.0 + sum_init);
E02 = exp(T_E02) / (1.0 + sum_init);
R01 = exp(T_R01) / (1.0 + sum_init);
R02 = exp(T_R02) / (1.0 + sum_init);
R12 = exp(T_R12) / (1.0 + sum_init);
//end_fromest

// -----INITIALISATION-----//
//start_rinit

// E01 = initial proportion of the population infected with virus 1
// E02 = inital proportion of the population infected with virus 2
// R01 = inital proportion of the population immune to virus 1
// R02 = initial proportion of the population immune to virus 2
// R12 = initial proportion of the population immune to both virus 1 and 2 
// N = overall population size
// note: E01 +  E02 + R01 + R02 +  R12 must be <=1

// Each compartment
X_SS = nearbyint((1.0 - E01 - E02 - R01 - R02 - R12) * N);
X_ES = nearbyint(E01 * N);
X_IS = 0;
X_TS = 0;
X_RS = nearbyint(R01 * N);
X_SE = nearbyint(E02 * N);
X_EE = 0;
X_IE = 0;
X_TE = 0;
X_RE = 0;
X_SI = 0;
X_EI = 0;
X_II = 0;
X_TI = 0;
X_RI = 0;
X_ST = 0;
X_ET = 0;
X_IT = 0;
X_TT = 0;
X_RT = 0;
X_SR = nearbyint(R02 * N);
X_ER = 0;
X_IR = 0;
X_TR = 0;
X_RR = nearbyint(R12 * N);

// if the compartments that are assigned some initial value don't sum to
// N (likely due to rounding in the nearbyint function) print the values
// then re calculate the initial value of X_SS by doing N - all other
// inital values for the other compartments
if ((X_SS + X_ES + X_RS + X_SE + X_SR + X_RR) != N) {
  Rprintf("SS=%f, ES=%f, RS=%f, SE=%f, SR=%f, RR=%f, sum=%f, N=%f, E01=%f, E02=%f, R01=%f, R02=%f, R12=%f\n", X_SS, X_ES, X_RS, X_SE, X_SR, X_RR, X_SS + X_ES + X_RS + X_SE + X_SR + X_RR, N, E01, E02,R01,R02,R12);
  X_SS = nearbyint(N - X_ES - X_RS - X_SE - X_SR - X_RR);
}


//end_rinit


//------ MEASUREMENNT MODEL ------//

// SIMULATION 
//start_rmeas
// generate the total number of tests positive to each virus 
v1_obs = rbinom(v1_T, rho1); // virus 1
v2_obs = rbinom(v2_T, rho2); // virus 2
//end_rmeas

// EVALUATION
//start_dmeas
double ll_1, ll_2, ll;

// observation model: v1_obs ~ binomial(size=v1_T, probability=rho1)
// where: v1_obs is the observed number of v1 cases 
//        v1_T the total number of v1 cases in the population  
//        rho1 the probability of testing positive to v1 
// similarly for v2

// calculating components for the likelihood 
ll_1 = dbinom(v1_obs, nearbyint(v1_T), rho1, 1); 
ll_2 = dbinom(v2_obs, nearbyint(v2_T), rho2, 1); 

//Rprintf("v1_obs=%.4f, v2_obs=%.4f, v1_T=%.4f, v2_T=%.4f, ll_1=%.4f, ll_2=%.4f\n",v1_obs, v2_obs, v1_T, v2_T, ll_1, ll_2);

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(ll_1 + ll_2, -1e3);

// calculate likelihood by back transforming the log likelihood 
lik = (give_log) ? ll : exp(ll);
//end_dmeas


//---------- PROCESS MODEL ----------//


// DETERMINISTIC SKELETON
// note: you don't technically need a skeleton when you are doing only a stochastic model
// but it is useful to use the deterministic model for debugging purposes

//start_skel
// calculating the incidence of each infection 
double p1 = (X_IS + X_IE + X_II + X_IT + X_IR); // virus 1
double p2 = (X_SI + X_EI + X_II + X_TI + X_RI); // virus 2

// calculate the transmission rate for each virus 
// where beta_i = Reff*gamma and Reff is the effective reproductive number at time i in a partially susceptible 
// population.   

double beta1 = Ri1 / (1.0 - (R01 + R12)) * gamma1; // virus 1 
double beta2 = Ri2 / (1.0 - (R02 + R12)) * gamma2; // virus 2 

// incorporate seasonality parameter for each virus 
// where A = amplitude, omega = annual angular frequency, t = time and phi = phase
double omega = (2 * M_PI)/52;
double s1 = 1 + A1 * cos(omega * (t - phi1));
double s2 = 1 + A2 * cos(omega * (t - phi2));

// calculate force of infection for each virus - note A = 0 means no seasonality component  
double lambda1 = beta1 * (p1/N) * s1; // virus 1
double lambda2 = beta2 * (p2/N) * s2; // virus 2

// addition of surges for v1
// note: we are using these surges for virus 1 only as it represents influenza 
// which has a number of different strains and new mutations each season      

// the new rate of immunity w1_loss = w1 + delta(t)
// where delta(t) = delta_i for i in [1,n] when t=t_si days since start of season i;
//                  0 otherwise

// note: to see large noticeable changes in the simulated data the
// loss in immunity delta_i needs to be quite large - otherwise the 
// changes are relatively subtle

// initialising vectors for t_si and delta_i
double *t_vec = (double *) &t_si_1;
double *delta_vec = (double *) &delta_i_1;
double w1_s;

// assigning the loss in immunity depending on the number of surges we have
for(int i = 0; i < nsurges + 1; i++){
  if(floor(t) == t_vec[i]) { // if t is a surge time point the add the surge in loss of immunity
    w1_s = w1 + delta_vec[i];
    break; // exit if we find a surge point
  } else{
    w1_s = w1; // if we don't find a surge point then just set the constant immunity loss
  }
}

// ODEs
//column 1 of schematic
DX_SS = -(lambda1 + lambda2) * X_SS + w2 * X_SR + w1_s * X_RS;
DX_SE = lambda2 * X_SS - (lambda1 + sigma2) * X_SE + w1_s * X_RE;
DX_SI = sigma2 * X_SE - (lambda1 * theta_lambda2 + gamma2) * X_SI + w1_s * X_RI;
DX_ST = gamma2 * X_SI - (lambda1 * theta_lambda2 + delta2) * X_ST + w1_s * X_RT;
DX_SR = delta2 * X_ST - lambda1 * X_SR + w1_s * X_RR - w2 * X_SR;

// column 2  of schematic
DX_ES = lambda1 * X_SS - (sigma1 + lambda2) * X_ES + w2 * X_ER;
DX_EE = lambda2 * X_ES + lambda1 * X_SE - (sigma1 + sigma2) * X_EE;
DX_EI = lambda1 * theta_lambda2 * X_SI + sigma2 * X_EE - (sigma1 + gamma2) * X_EI;
DX_ET = gamma2 * X_EI + lambda1 * theta_lambda2 * X_ST - (sigma1 + delta2) * X_ET;
DX_ER = lambda1 * X_SR + delta2 * X_ET - sigma1 * X_ER - w2 * X_ER ;

// column 3  of schematic
DX_IS = sigma1 * X_ES - (gamma1 + lambda2 * theta_lambda1) * X_IS + w2 * X_IR;
DX_IE = lambda2 * theta_lambda1 * X_IS + sigma1 * X_EE - (gamma1 + sigma2) * X_IE;
DX_II = sigma1 * X_EI + sigma2 * X_IE - (gamma1 + gamma2) * X_II;
DX_IT = sigma1 * X_ET + gamma2 * X_II - (gamma1 + delta2) * X_IT;
DX_IR = delta2 * X_IT + sigma1 * X_ER - gamma1 * X_IR -w2 * X_IR;

//column 4  of schematic
DX_TS = gamma1 * X_IS - (delta1 + lambda2 * theta_lambda1) * X_TS + w2 * X_TR;
DX_TE = lambda2 * theta_lambda1 * X_TS + gamma1 * X_IE - (delta1 + sigma2) * X_TE;
DX_TI = sigma2 * X_TE + gamma1 * X_II - (delta1 + gamma2) * X_TI;
DX_TT = gamma1 * X_IT + gamma2 * X_TI - (delta1 + delta2)* X_TT;
DX_TR = gamma1 * X_IR + delta2 * X_TT - delta1 * X_TR - w2 * X_TR;

//column 5  of schematic
DX_RS = delta1 * X_TS - lambda2 * X_RS + w2 * X_RR - w1_s * X_RS;
DX_RE = lambda2 * X_RS + delta1 * X_TE - sigma2 * X_RE - w1_s * X_RE;
DX_RI = sigma2 * X_RE + delta1 * X_TI - gamma2 * X_RI - w1_s * X_RI;
DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2* X_RT - w1_s * X_RT;
DX_RR = delta2 * X_RT + delta1  * X_TR - w1_s * X_RR - w2 * X_RR;

Dv1_T = p1 * gamma1;
Dv2_T = p2 * gamma2;
//end_skel

// Stochastic model SIMULATION

//start_rsim
// calculate prevalence of each virus 
double p1 = (X_IS + X_IE +  X_II + X_IT + X_IR); // virus 1
double p2 = (X_SI + X_EI +  X_II + X_TI + X_RI); // virus 2

// calculate basic reproductive number R0 for each virus 
double R0_1 = Ri1 / (1.0 - (R01 + R12)); // virus 1
double R0_2 = Ri2 / (1.0 - (R02 + R12)); // virus 2

//Rprintf("R0_1=%.4f, Ri1=%.4f, R01=%.4f, R12_1=%.4f\n", R0_1, Ri1, R01, R12);

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

// incorporate seasonality parameter for each virus 
// where A = amplitude, omega = annual angular frequency, t = time and phi = phase
double omega = (2 * M_PI)/52;
double s1 = 1 + A1 * cos(omega * (t - phi1));
double s2 = 1 + A2 * cos(omega * (t - phi2));

// calculate force of infection for each virus 
double lambda1 = beta1 * (p1/N) * s1; // virus 1
double lambda2 = beta2 * (p2/N) * s2; // virus 2

// addition of surges for vir1
// initialising vectors for t_si and delta_i
double *t_vec = (double *) &t_si_1;
double *delta_vec = (double *) &delta_i_1;
double w1_s;

// assigning the loss in immunity depending on the number of surges we have
for(int i = 0; i < nsurges + 1; i++){
    if(t == t_vec[i]) { // if t is a surge time point the add the surge in loss of immunity
      w1_s = w1 + delta_vec[i];
      break; // exit if we find a surge point
   } else{
      w1_s = w1; // if we don't find a surge point then just set the constant immunity loss
   }
}

// specifying the transitions 
double rates[50];// vector of length 50
double fromSS[2], fromES[2], fromIS[2], fromTS[2], fromRS[2];
double fromSE[2], fromEE[2], fromIE[2], fromTE[2], fromRE[2];
double fromSI[2], fromEI[2], fromII[2], fromTI[2], fromRI[2];
double fromST[2], fromET[2], fromIT[2], fromTT[2], fromRT[2]; 
double fromSR[2], fromER[2], fromIR[2], fromTR[2], fromRR[2]; // vectors of length 2

// row 1 of schematic
rates[0] = lambda1; // force of infection virus 1 (X_SS -> X_ES)
rates[1] = lambda2; // force of infection virus 2 (X_SS -> X_SE)
rates[2] = sigma1;  // (X_ES -> X_IS)
rates[3] = lambda2; //(X_ES -> X_EE)
rates[4] = gamma1; // (X_IS -> X_TS)
rates[5] = lambda2 * theta_lambda1; // (X_IS -> X_IE)
rates[6] = delta1; // (X_TS -> X_RS)
rates[7] = lambda2 * theta_lambda1; // (X_TS -> X_TE)
rates[8] = w1_s; // (X_RS -> X_SS)
rates[9] = lambda2; // (X_RS -> X_RE)

// row 2 of schematic
rates[10] = lambda1; // (X_SE -> X_EE)
rates[11] = sigma2;  //  (X_SE -> X_SI)
rates[12] = sigma1; // (X_EE -> X_IE)
rates[13] = sigma2; // (X_EE -> X_EI)
rates[14] = gamma1; // (X_IE -> X_TE)
rates[15] = sigma2; // (X_IE -> X_II)
rates[16] = delta1; // (X_TE -> X_RE)
rates[17] = sigma2; // (X_TE -> X_TI)
rates[18] = w1_s; // (X_RE -> X_SE)
rates[19] = sigma2; // (X_RE -> X_RI)

// row 3 of schematic
rates[20] = lambda1 * theta_lambda2; // (X_SI -> X_EI)
rates[21] = gamma2; // (X_SI -> X_ST)
rates[22] = sigma1; // (X_EI -> X_II)
rates[23] = gamma2; // (X_EI -> X_ET)
rates[24] = gamma1; // (X_II -> X_TI)
rates[25] = gamma2; // (X_II -> X_IT)
rates[26] = delta1; // (X_TI -> X_RI)
rates[27] = gamma2; // (X_TI -> X_TT)
rates[28] = w1_s; // (X_RI -> X_SI)
rates[29] = gamma2; // (X_RI -> X_RT)

// row 4 of schematic
rates[30] = lambda1 * theta_lambda2; // (X_ST -> X_ET)
rates[31] = delta2; // (X_ST -> X_SR)
rates[32] = sigma1; // (X_ET -> X_IT)
rates[33] = delta2; // (X_ET -> X_ER)
rates[34] = gamma1; // (X_IT -> X_TT)
rates[35] = delta2; // (X_IT -> X_IR)
rates[36] = delta1; // (X_TT -> X_RT)
rates[37] = delta2; // (X_TT -> X_TR)
rates[38] = w1_s; // (X_RT -> X_ST)
rates[39] = delta2; // (X_RT -> X_RR)

// row 5 of schematic
rates[40] = lambda1; // (X_SR -> X_ER)
rates[41] = w2; // (X_SR -> X_SS)
rates[42] = sigma1; // (X_ER -> X_IR)
rates[43] = w2; // (X_ER -> X_ES)
rates[44] = gamma1; // (X_IR -> X_TR)
rates[45] = w2; // (X_IR -> X_IS)
rates[46] = delta1; // (X_TR -> X_RR)
rates[47] = w2; // (X_TR -> X_TS)
rates[48] = w1_s; // (X_RR -> X_SR)
rates[49] = w2;// (X_RR -> X_RS)


// drawing sample for each of the compartments from the Euler-multinomial distribution
// returns a length(rate[i]) by n matrix where in our case we have 2 columns c1 (i.e. vec[0]) which we let represent
// transitions due to virus 1 (i.e. horizontally across compartments) and c2 the transitions (i.e. vec[1])
// due to virus 2 (i.e. vertically down compartments)

// row 1 of schematic
reulermultinom(2, X_SS, &rates[0], dt, &fromSS[0]);
reulermultinom(2, X_ES, &rates[2], dt, &fromES[0]);
reulermultinom(2, X_IS, &rates[4], dt, &fromIS[0]);
reulermultinom(2, X_TS, &rates[6], dt, &fromTS[0]);
reulermultinom(2, X_RS, &rates[8], dt, &fromRS[0]);

// row 2
reulermultinom(2, X_SE, &rates[10], dt, &fromSE[0]);
reulermultinom(2, X_EE, &rates[12], dt, &fromEE[0]);
reulermultinom(2, X_IE, &rates[14], dt, &fromIE[0]);
reulermultinom(2, X_TE, &rates[16], dt, &fromTE[0]);
reulermultinom(2, X_RE, &rates[18], dt, &fromRE[0]);

// row 3
reulermultinom(2, X_SI, &rates[20], dt, &fromSI[0]);
reulermultinom(2, X_EI, &rates[22], dt, &fromEI[0]);
reulermultinom(2, X_II, &rates[24], dt, &fromII[0]);
reulermultinom(2, X_TI, &rates[26], dt, &fromTI[0]);
reulermultinom(2, X_RI, &rates[28], dt, &fromRI[0]);

// row 4
reulermultinom(2, X_ST, &rates[30], dt, &fromST[0]);
reulermultinom(2, X_ET, &rates[32], dt, &fromET[0]);
reulermultinom(2, X_IT, &rates[34], dt, &fromIT[0]);
reulermultinom(2, X_TT, &rates[36], dt, &fromTT[0]);
reulermultinom(2, X_RT, &rates[38], dt, &fromRT[0]);

// row 5
reulermultinom(2, X_SR, &rates[40], dt, &fromSR[0]);
reulermultinom(2, X_ER, &rates[42], dt, &fromER[0]);
reulermultinom(2, X_IR, &rates[44], dt, &fromIR[0]);
reulermultinom(2, X_TR, &rates[46], dt, &fromTR[0]);
reulermultinom(2, X_RR, &rates[48], dt, &fromRR[0]);

// balance equations
// row 1 of schematic
X_SS += - fromSS[0] - fromSS[1] + fromRS[0] + fromSR[1];
X_ES += fromSS[0] - fromES[0] - fromES[1] + fromER[1];
X_IS += fromES[0] - fromIS[0] - fromIS[1] + fromIR[1];
X_TS += fromIS[0] - fromTS[0] - fromTS[1] + fromTR[1];
X_RS += fromTS[0] - fromRS[0] - fromRS[1] + fromRR[1];

// row 2
X_SE += fromSS[1] - fromSE[0] - fromSE[1] + fromRE[0];
X_EE += fromES[1] + fromSE[0] - fromEE[0] - fromEE[1];
X_IE += fromIS[1] + fromEE[0] - fromIE[0] - fromIE[1];
X_TE += fromTS[1] + fromIE[0] - fromTE[0] - fromTE[1];
X_RE += fromRS[1] + fromTE[0] - fromRE[0] - fromRE[1];

// row 3
X_SI += fromSE[1] - fromSI[0] - fromSI[1] + fromRI[0];
X_EI += fromEE[1] + fromSI[0] - fromEI[0] - fromEI[1];
X_II += fromIE[1] + fromEI[0] - fromII[0] - fromII[1];
X_TI += fromTE[1] + fromII[0] - fromTI[0] - fromTI[1];
X_RI += fromRE[1] + fromTI[0] - fromRI[0] - fromRI[1];

// row 4
X_ST += fromSI[1] - fromST[0] - fromST[1] + fromRT[0];
X_ET += fromEI[1] + fromST[0] - fromET[0] - fromET[1];
X_IT += fromII[1] + fromET[0] - fromIT[0] - fromIT[1];
X_TT += fromTI[1] + fromIT[0] - fromTT[0] - fromTT[1];
X_RT += fromRI[1] + fromTT[0] - fromRT[0] - fromRT[1];

// row 5
X_SR += fromST[1] - fromSR[0] - fromSR[1] + fromRR[0];
X_ER += fromET[1] + fromSR[0] - fromER[0] - fromER[1];
X_IR += fromIT[1] + fromER[0] - fromIR[0] - fromIR[1];
X_TR += fromTT[1] + fromIR[0] - fromTR[0] - fromTR[1];
X_RR += fromRT[1] + fromTR[0] - fromRR[0] - fromRR[1];

// Total incidence of each virus in the population
v1_T += fromIS[0] + fromIE[0] + fromII[0] + fromIT[0] + fromIR[0];
v2_T += fromSI[1] + fromEI[1] + fromII[1] + fromTI[1] + fromRI[1];

//end_rsim
