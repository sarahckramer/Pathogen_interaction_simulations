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
//end_globs

// -----INITIALISATION-----//
//start_rinit

// E01 = initial proportion of the population infected with virus 1
// E02 = inital proportion of the population infected with virus 2
// R01 = inital proportion of the population immune to virus 1
// R02 = initial proportion of the population immune to virus 2
// R12 = initial proportion of the population immune to both virus 1 and 2 
// N = overall population size

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

// Accumulator variables 
v1_T = 0; // takes into consideration interaction  
v2_T = 0; 

// if the compartments that are assigned some initial value don't sum to
// N (likely due to rounding in the nearbyint function) print the values
// then re calculate the initial value of X_SS by doing N - all other
// inital values for the other compartments
if ((X_SS + X_ES + X_RS + X_SE + X_SR + X_RR) != N) {
  Rprintf("SS=%f, ES=%f, RS=%f, SE=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_ES, X_RS, X_SE, X_SR, X_RR, X_SS + X_ES + X_RS + X_SE + X_SR + X_RR, N);
  X_SS = nearbyint(N - X_ES - X_RS - X_SE - X_SR - X_RR);
}

Rprintf("E01=%f, E02=%f, R01=%f, R02=%f, R12=%f\n", E01, E02, R01, R02, R12);

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
ll_1 = dbinom(v1_obs, v1_T, rho1, 1); 
ll_2 = dbinom(v2_obs, v2_T, rho2, 1); 

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(ll_1 + ll_2, -1e3);

// calculate likelihood by back transforming the log likelihood 
lik = (give_log) ? ll : exp(ll);
//end_dmeas


//---------- PROCESS MODEL ----------//

// Stochastic model SIMULATION
// note: you don't need a skeleton when you are doing only a stochastic model

//start_rsim
// calculate prevalence of each virus 
double p1 = (X_IS + X_IE +  X_II + X_IT + X_IR); // virus 1
double p2 = (X_SI + X_EI +  X_II + X_TI + X_RI); // virus 2

// calculate basic reproductive number R0 for each virus 
double R0_1 = Ri1 / (1.0 - (R01 + R12)); // virus 1
double R0_2 = Ri2 / (1.0 - (R02 + R12)); // virus 2

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
double omega = (2 * M_PI)/365;
double s = 1 + A * cos(omega * (t - phi));

// calculate force of infection for each virus 
double lambda1 = beta1 * (p1/N) * s; // virus 1
double lambda2 = beta2 * (p2/N) * s; // virus 2

Rprintf("p1=%.2f, p2=%.2f, beta1=%.1f, beta2=%.1f lambda1=%.3f, lambda2=%.3f, dt=%.3f\n", 
        p1, p2, beta1, beta2, lambda1, lambda2,dt);
        
// specifying the transitions 
double rates[49];// vector of length 49
double fromSS[2], fromES[2], fromIS[2], fromTS[2], fromRS[2];
double fromSE[2], fromEE[2], fromIE[2], fromTE[2], fromRE[2];
double fromSI[2], fromEI[2], fromII[2], fromTI[2], fromRI[2];
double fromST[2], fromET[2], fromIT[2], fromTT[2], fromRT[2]; 
double fromSR[2], fromER[2], fromIR[2], fromTR[2], fromRR[2]; // vectors of length 2

// specifying rate for transition - note: vector indexing starts at 0 for C++ rather than 1 like R

// row 1 of schematic
rates[0] = lambda1; // force of infection virus 1 (X_SS -> X_ES)
rates[1] = lambda2; // force of infection virus 2 (X_SS -> X_SE)
rates[2] = sigma1;  // (X_ES -> X_IS)
rates[3] = lambda2; //(X_ES -> X_EE)
rates[4] = gamma1; // (X_IS -> X_TS)
rates[5] = lambda2 * theta_lambda1; // (X_IS -> X_IE)
rates[6] = delta1; // (X_TS -> X_RS)
rates[7] = lambda2 * theta_lambda1; // (X_TS -> X_TE)
rates[8] = lambda2; // (X_RS -> X_RE)
rates[9] = w1; // (X_RS -> X_SS)

// row 2 of schematic
rates[10] = lambda1; // (X_SE -> X_EE)
rates[11] = sigma2;  //  (X_SE -> X_SI)
rates[12] = sigma1; // (X_EE -> X_IE)
rates[13] = sigma2; // (X_EE -> X_EI)
rates[14] = gamma1; // (X_IE -> X_TE)
rates[15] = sigma2; // (X_IE -> X_II)
rates[16] = delta1; // (X_TE -> X_RE)
rates[17] = sigma2; // (X_TE -> X_TI)
rates[18] = sigma2; // (X_RE -> X_RI)
rates[19] = w1; // (X_RE -> X_SE)

// row 3 of schematic
rates[20] = lambda1 * theta_lambda2; // (X_SI -> X_EI)
rates[21] = gamma2; // (X_SI -> X_ST)
rates[22] = sigma1; // (X_EI -> X_II)
rates[23] = gamma2; // (X_EI -> X_ET)
rates[24] = gamma1; // (X_II -> X_TI)
rates[25] = gamma2; // (X_II -> X_IT)
rates[26] = delta1; // (X_TI -> X_RI)
rates[27] = gamma2; // (X_TI -> X_TT)
rates[28] = gamma2; // (X_RI -> X_RT)
rates[29] = w1; // (X_RI -> X_SI)

// row 4 of schematic
rates[30] = lambda1 * theta_lambda2; // (X_ST -> X_ET)
rates[31] = delta2; // (X_ST -> X_SR)
rates[32] = sigma1; // (X_ET -> X_IT)
rates[33] = delta2; // (X_ET -> X_ER)
rates[34] = gamma1; // (X_IT -> X_TT)
rates[35] = delta2; // (X_IT -> X_IR)
rates[36] = delta1; // (X_TT -> X_RT)
rates[37] = delta2; // (X_TT -> X_TR)
rates[38] = delta2; // (X_RT -> X_RR)
rates[39] = w1; // (X_RT -> X_ST)

// row 5 of schematic
rates[40] = lambda1; // (X_SR -> X_ER) 
rates[41] = w2; // (X_SR -> X_SS)
rates[42] = sigma1; // (X_ER -> X_IR)
rates[43] = w2; // (X_ER -> X_ES)
rates[44] = gamma1; // (X_IR -> X_TR)
rates[45] = w2; // (X_IR -> X_IS)
rates[46] = delta1; // (X_TR -> X_RR)
rates[47] = w2; // (X_TR -> X_TS)
rates[48] = w1; // (X_RR -> X_SR)
rates[49] = w2;// (X_RR -> X_RS)

// drawing sample for each of the compartments from the Euler-multinomial distribution
// returns a length(rate[i]) by n matrix where in our case we have 2 columns c1 which we let represent
// transitions due to virus 1 (i.e. vertically down compartments) and c2 the transitions
// due to virus 2 (i.e. horizontally across compartments)

// row 1
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

//Rprintf("fromRS=%.1f, fromRE=%.1f, fromRI=%.1f, fromRT=%.1f\n",
//        fromRS, fromRE, fromRI, fromRT);
          
// balance equations

// row 1 of schmematic
X_SS += -fromSS[0] - fromSS[1] + fromRS[1] + fromSR[1];
X_ES += fromSS[0] - fromES[0] - fromES[1] + fromER[1];
X_IS += fromES[0] - fromIS[0] - fromIS[1] + fromIR[1];
X_TS += fromIS[0] - fromTS[0] - fromTS[1] + fromTR[1];
X_RS += fromTS[0] - fromRS[0] - fromRS[1] + fromRR[1];

// row 2
X_SE += fromSS[1] - fromSE[0] - fromSE[1] + fromRE[1];
X_EE += fromES[1] + fromSE[0] - fromEE[0] - fromEE[1];
X_IE += fromIS[1] + fromEE[0] - fromIE[0] - fromIE[1];
X_TE += fromTS[1] + fromIE[0] - fromTE[0] - fromTE[1];
X_RE += fromRS + fromTE[0] - fromRE[0] -fromRE[1];

// row 3
X_SI += fromSE[1] - fromSI[0] - fromSI[1] + fromRI[1];
X_EI += fromEE[1] + fromSI[0] - fromEI[0] - fromEI[1];
X_II += fromIE[1] + fromEI[0] - fromII[0] - fromII[1];
X_TI += fromTE[1] + fromII[0] - fromTI[0] - fromTI[1];
X_RI += fromRE + fromTI[0] - fromRI[0] - fromRI[1];

// row 4
X_ST += fromSI[1] - fromST[0] - fromST[1] + fromRT[1];
X_ET += fromEI[1] + fromST[0] - fromET[0] - fromET[1];
X_IT += fromII[1] + fromET[0] - fromIT[0] - fromIT[1];
X_TT += fromTI[1] + fromIT[0] - fromTT[0] - fromTT[1];
X_RT += fromRI + fromTT[0] - fromRT[0] - fromRT[1];

// row 5
X_SR += fromST[1] - fromSR[0] - fromSR[1] + fromRR[0];
X_ER += fromET[1] + fromSR[0] - fromER[0] - fromER[1];
X_IR += fromIT[1] + fromER[0] - fromIR[0] - fromER[1];
X_TR += fromTT[1] + fromIR[0] - fromTR[0] - fromTR[1];
X_RR += fromRT[1] + fromTR[0] - fromRR[0] - fromRR[1];

// Total number of cases of each virus in the population
v1_T += (fromIS[1] + fromIE[1] + fromII[1] + fromIT[1] + fromIR[1]);
v2_T += (fromSI[0] + fromEI[0] + fromII[0] + fromTI[0] + fromRI[0]);
//end_rsim
