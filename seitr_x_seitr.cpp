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

// addition of other variables to be able to monitor output
// w = 0;
// delta = 0;
// lambda_1 = 0;
// lambda_2 = 0;
// gamma_1 = 0;
// beta_1 = 0;
// s_1 =0;

// if the compartments that are assigned some initial value don't sum to
// N (likely due to rounding in the nearbyint function) print the values
// then re calculate the initial value of X_SS by doing N - all other
// inital values for the other compartments
if ((X_SS + X_ES + X_RS + X_SE + X_SR + X_RR) != N) {
  Rprintf("SS=%f, ES=%f, RS=%f, SE=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_ES, X_RS, X_SE, X_SR, X_RR, X_SS + X_ES + X_RS + X_SE + X_SR + X_RR, N);
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
ll_1 = dbinom(v1_obs, v1_T, rho1, 1); 
ll_2 = dbinom(v2_obs, v2_T, rho2, 1); 

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
// calculating the prevalence of each infection 
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

// ODEs
//column 1 of schematic
DX_SS = -(lambda1 + lambda2) * X_SS + w2 * X_SR + w1 * X_RS + mu * N - nu * X_SS;
DX_SE = lambda2 * X_SS - (lambda1 + sigma2) * X_SE + w1 * X_RE - nu * X_SE;
DX_SI = sigma2 * X_SE - (lambda1 * theta_lambda2 + gamma2) * X_SI + w1 * X_RI - nu * X_SI;
DX_ST = gamma2 * X_SI - (lambda1 * theta_lambda2 + delta2) * X_ST + w1 * X_RT - nu * X_ST;
DX_SR = delta2 * X_ST - lambda1 * X_SR + w1 * X_RR - w2 * X_SR - nu * X_SR;

// column 2  of schematic
DX_ES = lambda1 * X_SS - (sigma1 + lambda2) * X_ES + w2 * X_ER - nu * X_ES;
DX_EE = lambda2 * X_ES + lambda1 * X_SE - (sigma1 + sigma2) * X_EE - nu * X_EE;
DX_EI = lambda1 * theta_lambda2 * X_SI + sigma2 * X_EE - (sigma1 + gamma2) * X_EI - nu * X_EI;
DX_ET = gamma2 * X_EI + lambda1 * theta_lambda2 * X_ST - (sigma1 + delta2) * X_ET - nu * X_ET;
DX_ER = lambda1 * X_SR + delta2 * X_ET - sigma1 * X_ER - w2 * X_ER - nu * X_ER;

// column 3  of schematic
DX_IS = sigma1 * X_ES - (gamma1 + lambda2 * theta_lambda1) * X_IS + w2 * X_IR - nu * X_IS;
DX_IE = lambda2 * theta_lambda1 * X_IS + sigma1 * X_EE - (gamma1 + sigma2) * X_IE - nu * X_IE;
DX_II = sigma1 * X_EI + sigma2 * X_IE - (gamma1 + gamma2) * X_II - nu * X_II;
DX_IT = sigma1 * X_ET + gamma2 * X_II - (gamma1 + delta2) * X_IT - nu * X_IT;
DX_IR = delta2 * X_IT + sigma1 * X_ER - gamma1 * X_IR -w2 * X_IR - nu * X_IR;

//column 4  of schematic
DX_TS = gamma1 * X_IS - (delta1 + lambda2 * theta_lambda1) * X_TS + w2 * X_TR - nu * X_TS;
DX_TE = lambda2 * theta_lambda1 * X_TS + gamma1 * X_IE - (delta1 + sigma2) * X_TE - nu * X_TE;
DX_TI = sigma2 * X_TE + gamma1 * X_II - (delta1 + gamma2) * X_TI - nu * X_TI;
DX_TT = gamma1 * X_IT + gamma2 * X_TI - (delta1 + delta2)* X_TT - nu * X_TT;
DX_TR = gamma1 * X_IR + delta2 * X_TT - delta1 * X_TR - w2 * X_TR - nu * X_TR;

//column 5  of schematic
DX_RS = delta1 * X_TS - lambda2 * X_RS + w2 * X_RR - w1 * X_RS - nu * X_RS;
DX_RE = lambda2 * X_RS + delta1 * X_TE - sigma2 * X_RE - w1 * X_RE - nu * X_RE;
DX_RI = sigma2 * X_RE + delta1 * X_TI - gamma2 * X_RI - w1 * X_RI - nu * X_RI;
DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2* X_RT - w1 * X_RT - nu * X_RT;
DX_RR = delta2 * X_RT + delta1  * X_TR - w1 * X_RR - w2 * X_RR - nu * X_RR;

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
// note: we are using these surges for virus 1 only as it represents influenza 
// which has a number of different strains and new mutations each season      

// the new rate of immunity w1_loss = w1 + delta(t)
// where delta(t) = delta_i for i in [1,n] when t=t_si days since start of season i;
//                  0 otherwise

// it is difficult to make this part of the code dynamic therefore
// will need to change the code if we want more than 5 surges
// (note: to see large noticeable changes in the simulated data the
// loss in immunity delta_i needs to be quite large - otherwise the 
// changes are relatively subtle)
double w1_s;
if(t==t_si_1){
   w1_s = w1 + delta_i_1;
   } else if (t==t_si_2){
     w1_s = w1 + delta_i_2;
   } else if (t==t_si_3){
     w1_s = w1 + delta_i_3;
   } else if(t==t_si_4){
     w1_s = w1 + delta_i_4;
   } else if(t==t_si_5){
     w1_s = w1 + delta_i_5;
   } else{
   w1_s = w1;
}

// specifying the transitions 
double rates[75];// vector of length 75
double fromSS[3], fromES[3], fromIS[3], fromTS[3], fromRS[3];
double fromSE[3], fromEE[3], fromIE[3], fromTE[3], fromRE[3];
double fromSI[3], fromEI[3], fromII[3], fromTI[3], fromRI[3];
double fromST[3], fromET[3], fromIT[3], fromTT[3], fromRT[3]; 
double fromSR[3], fromER[3], fromIR[3], fromTR[3], fromRR[3]; // vectors of length 3

// specifying rate for transition - note: vector indexing starts at 0 for C++ rather than 1 like R

// row 1 of schematic
rates[0] = lambda1; // force of infection virus 1 (X_SS -> X_ES)
rates[1] = lambda2; // force of infection virus 2 (X_SS -> X_SE)
rates[2] = nu; // natural X_SS death rate
rates[3] = sigma1;  // (X_ES -> X_IS)
rates[4] = lambda2; //(X_ES -> X_EE)
rates[5] = nu; // natural X_ES death rate  
rates[6] = gamma1; // (X_IS -> X_TS)
rates[7] = lambda2 * theta_lambda1; // (X_IS -> X_IE)
rates[8] = nu; // natural X_IS death rate
rates[9] = delta1; // (X_TS -> X_RS)
rates[10] = lambda2 * theta_lambda1; // (X_TS -> X_TE)
rates[11] = nu; // natural X_TS death rate
rates[12] = w1_s; // (X_RS -> X_SS)
rates[13] = lambda2; // (X_RS -> X_RE)
rates[14] = nu; // natural X_RS death rate

// row 2 of schematic
rates[15] = lambda1; // (X_SE -> X_EE)
rates[16] = sigma2;  //  (X_SE -> X_SI)
rates[17] = nu; // natural X_SE death rate
rates[18] = sigma1; // (X_EE -> X_IE)
rates[19] = sigma2; // (X_EE -> X_EI)
rates[20] = nu; // natural X_EE death rate
rates[21] = gamma1; // (X_IE -> X_TE)
rates[22] = sigma2; // (X_IE -> X_II)
rates[23] = nu; // natural X_IE death rate
rates[24] = delta1; // (X_TE -> X_RE)
rates[25] = sigma2; // (X_TE -> X_TI)
rates[26] = nu; // natural X_TE death rate
rates[27] = w1_s; // (X_RE -> X_SE)
rates[28] = sigma2; // (X_RE -> X_RI)
rates[29] = nu; // natural X_RE death rate

// row 3 of schematic
rates[30] = lambda1 * theta_lambda2; // (X_SI -> X_EI)
rates[31] = gamma2; // (X_SI -> X_ST)
rates[32] = nu; // natural X_SI death rate 
rates[33] = sigma1; // (X_EI -> X_II)
rates[34] = gamma2; // (X_EI -> X_ET)
rates[35] = nu; // natural X_EI death rate
rates[36] = gamma1; // (X_II -> X_TI)
rates[37] = gamma2; // (X_II -> X_IT)
rates[38] = nu; // natural X_II death rate
rates[39] = delta1; // (X_TI -> X_RI)
rates[40] = gamma2; // (X_TI -> X_TT)
rates[41] = nu; // natural X_TI death rate
rates[42] = w1_s; // (X_RI -> X_SI)
rates[43] = gamma2; // (X_RI -> X_RT)
rates[44] = nu; // natural X_RI death rate

// row 4 of schematic
rates[45] = lambda1 * theta_lambda2; // (X_ST -> X_ET)
rates[46] = delta2; // (X_ST -> X_SR)
rates[47] = nu; // natural X_ST death rate
rates[48] = sigma1; // (X_ET -> X_IT)
rates[49] = delta2; // (X_ET -> X_ER)
rates[50] = nu; // natural X_ET death rate
rates[51] = gamma1; // (X_IT -> X_TT)
rates[52] = delta2; // (X_IT -> X_IR)
rates[53] = nu; // natural X_IT death rate
rates[54] = delta1; // (X_TT -> X_RT)
rates[55] = delta2; // (X_TT -> X_TR)
rates[56] = nu; // natural X_TT death rate
rates[57] = w1_s; // (X_RT -> X_ST)
rates[58] = delta2; // (X_RT -> X_RR)
rates[59] = nu; // natural X_RT death rate

// row 5 of schematic
rates[60] = lambda1; // (X_SR -> X_ER) 
rates[61] = w2; // (X_SR -> X_SS)
rates[62] = nu; // natural X_SR death rate
rates[63] = sigma1; // (X_ER -> X_IR)
rates[64] = w2; // (X_ER -> X_ES)
rates[65] = nu; // natural X_ER death rate
rates[66] = gamma1; // (X_IR -> X_TR)
rates[67] = w2; // (X_IR -> X_IS)
rates[68] = nu; // natural X_IR death rate
rates[69] = delta1; // (X_TR -> X_RR)
rates[70] = w2; // (X_TR -> X_TS)
rates[71] = nu; // natural X_TR death rate
rates[72] = w1_s; // (X_RR -> X_SR)
rates[73] = w2;// (X_RR -> X_RS)
rates[74] = nu; // natural X_RR death rate

// drawing sample for each of the compartments from the Euler-multinomial distribution
// returns a length(rate[i]) by n matrix where in our case we have 2 columns c1 (i.e. vec[0]) which we let represent
// transitions due to virus 1 (i.e. horizontally across compartments) and c2 the transitions (i.e. vec[1])
// due to virus 2 (i.e. vertically down compartments)

// row 1 of schematic 
reulermultinom(3, X_SS, &rates[0], dt, &fromSS[0]);
reulermultinom(3, X_ES, &rates[3], dt, &fromES[0]);
reulermultinom(3, X_IS, &rates[6], dt, &fromIS[0]);
reulermultinom(3, X_TS, &rates[9], dt, &fromTS[0]);
reulermultinom(3, X_RS, &rates[12], dt, &fromRS[0]);

// row 2
reulermultinom(3, X_SE, &rates[15], dt, &fromSE[0]);
reulermultinom(3, X_EE, &rates[18], dt, &fromEE[0]);
reulermultinom(3, X_IE, &rates[21], dt, &fromIE[0]);
reulermultinom(3, X_TE, &rates[24], dt, &fromTE[0]);
reulermultinom(3, X_RE, &rates[27], dt, &fromRE[0]);

// row 3
reulermultinom(3, X_SI, &rates[30], dt, &fromSI[0]);
reulermultinom(3, X_EI, &rates[33], dt, &fromEI[0]);
reulermultinom(3, X_II, &rates[36], dt, &fromII[0]);
reulermultinom(3, X_TI, &rates[39], dt, &fromTI[0]);
reulermultinom(3, X_RI, &rates[42], dt, &fromRI[0]);

// row 4
reulermultinom(3, X_ST, &rates[45], dt, &fromST[0]);
reulermultinom(3, X_ET, &rates[48], dt, &fromET[0]);
reulermultinom(3, X_IT, &rates[51], dt, &fromIT[0]);
reulermultinom(3, X_TT, &rates[54], dt, &fromTT[0]);
reulermultinom(3, X_RT, &rates[57], dt, &fromRT[0]);

// row 5
reulermultinom(3, X_SR, &rates[60], dt, &fromSR[0]);
reulermultinom(3, X_ER, &rates[63], dt, &fromER[0]);
reulermultinom(3, X_IR, &rates[66], dt, &fromIR[0]);
reulermultinom(3, X_TR, &rates[69], dt, &fromTR[0]);
reulermultinom(3, X_RR, &rates[72], dt, &fromRR[0]);

// Drawing number of births from a poisson distribution 
double births;
births = rpois(mu*N*dt);
//Rprintf("births=%.1f, mu=%.1f, dt=%.1f\n", births, mu, dt);

// balance equations

// row 1 of schematic
X_SS += births - fromSS[0] - fromSS[1] + fromRS[0] + fromSR[1] - fromSS[2];
X_ES += fromSS[0] - fromES[0] - fromES[1] + fromER[1] - fromES[2];
X_IS += fromES[0] - fromIS[0] - fromIS[1] + fromIR[1] - fromIS[2];
X_TS += fromIS[0] - fromTS[0] - fromTS[1] + fromTR[1] - fromTS[2];
X_RS += fromTS[0] - fromRS[0] - fromRS[1] + fromRR[1] - fromRS[2];

// row 2
X_SE += fromSS[1] - fromSE[0] - fromSE[1] + fromRE[0] - fromSE[2];
X_EE += fromES[1] + fromSE[0] - fromEE[0] - fromEE[1] - fromEE[2];
X_IE += fromIS[1] + fromEE[0] - fromIE[0] - fromIE[1] - fromIE[2];
X_TE += fromTS[1] + fromIE[0] - fromTE[0] - fromTE[1] - fromTE[2];
X_RE += fromRS[1] + fromTE[0] - fromRE[0] - fromRE[1] - fromRE[2];

// row 3
X_SI += fromSE[1] - fromSI[0] - fromSI[1] + fromRI[0] - fromSI[2];
X_EI += fromEE[1] + fromSI[0] - fromEI[0] - fromEI[1] - fromEI[2];
X_II += fromIE[1] + fromEI[0] - fromII[0] - fromII[1] - fromII[2];
X_TI += fromTE[1] + fromII[0] - fromTI[0] - fromTI[1] - fromTI[2];
X_RI += fromRE[1] + fromTI[0] - fromRI[0] - fromRI[1] - fromRI[2];

// row 4
X_ST += fromSI[1] - fromST[0] - fromST[1] + fromRT[0] - fromST[2];
X_ET += fromEI[1] + fromST[0] - fromET[0] - fromET[1] - fromET[2];
X_IT += fromII[1] + fromET[0] - fromIT[0] - fromIT[1] - fromIT[2];
X_TT += fromTI[1] + fromIT[0] - fromTT[0] - fromTT[1] - fromTT[2];
X_RT += fromRI[1] + fromTT[0] - fromRT[0] - fromRT[1] - fromRT[2];

// row 5
X_SR += fromST[1] - fromSR[0] - fromSR[1] + fromRR[0] - fromSR[2];
X_ER += fromET[1] + fromSR[0] - fromER[0] - fromER[1] - fromER[2];
X_IR += fromIT[1] + fromER[0] - fromIR[0] - fromIR[1] - fromIR[2];
X_TR += fromTT[1] + fromIR[0] - fromTR[0] - fromTR[1] - fromTR[2];
X_RR += fromRT[1] + fromTR[0] - fromRR[0] - fromRR[1] - fromRR[2];

// Total incidence of each virus in the population
v1_T += fromIS[0] + fromIE[0] + fromII[0] + fromIT[0] + fromIR[0];
v2_T += fromSI[1] + fromEI[1] + fromII[1] + fromTI[1] + fromRI[1];

//Rprintf("w1_s=%.4f, t=%.1f, t_si_1=%.1f\n", w1_s, t, t_si_1);

// outputting additional variables to monitor
// w = w1_s;
// delta = delta_i_1;
// lambda_1 = lambda1;
// lambda_2 = lambda2;
// gamma_1 = gamma1;
// beta_1 = beta1;
// s_1 = s1;

//end_rsim
