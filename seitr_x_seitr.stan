// STAN program for an SEITR x SEITR interaction model 
// 
// Useful tutorial: https://mc-stan.org/users/documentation/case-studies/boarding_school_case_study.html
//
// created by: Sarah Pirikahu
// creation date: 4 Sept 2023

  
  // ODEs for the deterministic skeleton 
functions {
  real[] sir(real t, real[] y, real[] theta, 
             real[] x_r, int[] x_i) {
    
    // initialising states y
    // row 1
    real X_SS = y[1];
    real X_ES = y[2];
    real X_IS = y[3];
    real X_TS = y[4];
    real X_RS = y[5];
    
    // row 2
    real X_SE = y[6]; 
    real X_EE = y[7];
    real X_IE = y[8];
    real X_TE = y[9];
    real X_RE = y[10];
    
    // row 3 
    real X_SI = y[11]; 
    real X_EI = y[12];
    real X_II = y[13];
    real X_TI = y[14];
    real X_RI = y[15];
    
    // row 4
    real X_ST = y[16]; 
    real X_ET = y[17];
    real X_IT = y[18];
    real X_TT = y[19];
    real X_RT = y[20];
    
    // row 5
    real X_SR = y[21];
    real X_ER = y[22];
    real X_IR = y[23];
    real X_TR = y[24];
    real X_RR = y[25];
    
    // calculating the incidence of each infection 
    real p1 = X_IS + X_IE + X_II + X_IT + X_IR;
    real p2 = X_SI + X_EI + X_II + X_TI + X_RI;
    
    // specifying the model parameters that will be estimated, i.e. \theta
    real delta1 = theta[1];
    real delta2 = theta[2];
    real theta_lambda1 = theta[3];
    real theta_lambda2 = theta[4];
    real w1 = theta[5];
    real w2 = theta[6];
    
    real gamma1 = theta[7]; // will fix these parameters later
    real gamma2 = theta[8];
    real sigma1 = theta[9];
    real sigma2 = theta[10];
    
    real Ri1 = theta[11];
    real Ri2 = theta[12];
    real R01 = theta[13];
    real R02 = theta[14];
    real R12 = theta[15];
    real A1 = theta[16];
    real A2 = theta[17];
    real phi1 = theta[18];
    real phi2 = theta[19];
    
    // specifying model parameters that are fixed 
    real mu = x_r[1]; // real values used to evaluate f; weekly birth rate
    real nu = x_r[2]; // weekly death rate
    real N = x_i[1]; // integer values used to evaluate f; population size 
    
    // transmission rates
    real beta1 = Ri1 / (1.0 - (R01 + R12)) * gamma1; // virus 1 
    real beta2 = Ri2 / (1.0 - (R02 + R12)) * gamma2; // virus 2
    
    // incorportion of seasonality 
    real omega = (2 * pi())/52;
    real s1 = 1 + A1 * cos(omega * (t - phi1));
    real s2 = 1 + A2 * cos(omega * (t - phi2));
    
    // force of infection rates
    real lambda1 = beta1 * (p1/N) * s1; // virus 1 
    real lambda2 = beta2 * (p2/N) * s2; // virus 2 
    
    // ODEs
    // row 1 of schematic
    real dX_SS = -(lambda1 + lambda2) * X_SS + w2 * X_SR + w1 * X_RS + mu * N - nu * X_SS;
    real dX_ES = lambda1 * X_SS - (sigma1 + lambda2) * X_ES + w2 * X_ER - nu * X_ES;
    real dX_IS = sigma1 * X_ES - (gamma1 + lambda2 * theta_lambda1) * X_IS + w2 * X_IR - nu * X_IS;
    real dX_TS = gamma1 * X_IS - (delta1 + lambda2 * theta_lambda1) * X_TS + w2 * X_TR - nu * X_TS;
    real dX_RS = delta1 * X_TS - lambda2 * X_RS + w2 * X_RR - w1 * X_RS - nu * X_RS;
    
    // row 2 of schematic 
    real dX_SE = lambda2 * X_SS - (lambda1 + sigma2) * X_SE + w1 * X_RE - nu * X_SE;
    real dX_EE = lambda2 * X_ES + lambda1 * X_SE - (sigma1 + sigma2) * X_EE - nu * X_EE;
    real dX_IE = lambda2 * theta_lambda1 * X_IS + sigma1 * X_EE - (gamma1 + sigma2) * X_IE - nu * X_IE;
    real dX_TE = lambda2 * theta_lambda1 * X_TS + gamma1 * X_IE - (delta1 + sigma2) * X_TE - nu * X_TE;
    real dX_RE = lambda2 * X_RS + delta1 * X_TE - sigma2 * X_RE - w1 * X_RE - nu * X_RE;
    
    // row 3 of schematic 
    real dX_SI = sigma2 * X_SE - (lambda1 * theta_lambda2 + gamma2) * X_SI + w1 * X_RI - nu * X_SI;
    real dX_EI = lambda1 * theta_lambda2 * X_SI + sigma2 * X_EE - (sigma1 + gamma2) * X_EI - nu * X_EI;
    real dX_II = sigma1 * X_EI + sigma2 * X_IE - (gamma1 + gamma2) * X_II - nu * X_II;
    real dX_TI = sigma2 * X_TE + gamma1 * X_II - (delta1 + gamma2) * X_TI - nu * X_TI;
    real dX_RI = sigma2 * X_RE + delta1 * X_TI - gamma2 * X_RI - w1 * X_RI - nu * X_RI;
    
    // row 4 of schematic 
    real dX_ST = gamma2 * X_SI - (lambda1 * theta_lambda2 + delta2) * X_ST + w1 * X_RT - nu * X_ST;
    real dX_ET = gamma2 * X_EI + lambda1 * theta_lambda2 * X_ST - (sigma1 + delta2) * X_ET - nu * X_ET;
    real dX_IT = sigma1 * X_ET + gamma2 * X_II - (gamma1 + delta2) * X_IT - nu * X_IT;
    real dX_TT = gamma1 * X_IT + gamma2 * X_TI - (delta1 + delta2)* X_TT - nu * X_TT;
    real dX_RT = gamma2 * X_RI + delta1 * X_TT - delta2* X_RT - w1 * X_RT - nu * X_RT;
    
    // row 5 of schematic  
    real dX_SR = delta2 * X_ST - lambda1 * X_SR + w1 * X_RR - w2 * X_SR - nu * X_SR;
    real dX_ER = lambda1 * X_SR + delta2 * X_ET - sigma1 * X_ER - w2 * X_ER - nu * X_ER;
    real dX_IR = delta2 * X_IT + sigma1 * X_ER - gamma1 * X_IR -w2 * X_IR - nu * X_IR;
    real dX_TR = gamma1 * X_IR + delta2 * X_TT - delta1 * X_TR - w2 * X_TR - nu * X_TR;
    real dX_RR = delta2 * X_RT + delta1  * X_TR - w1 * X_RR - w2 * X_RR - nu * X_RR;
    
    return {dX_SS, dX_ES, dX_IS, dX_TS, dX_RS,
      dX_SE, dX_EE, dX_IE, dX_TE, dX_RE,
      dX_SI, dX_EI, dX_II, dX_TI, dX_RI,
      dX_ST, dX_ET, dX_IT, dX_TT, dX_RT,
      dX_SR, dX_ER, dX_IR, dX_TR, dX_RR};
  }
}

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=1> n_weeks;
  real y0[25]; // inital condition 
  real t0; // the time of the inital condition 
  real ts[n_weeks]; // the time sereis at which we require the solution to be evaluated
  real mu;
  real nu;
  int N;
  int v1_obs[n_weeks];
  int v2_obs[n_weeks];
}

// specifying x_r and x_i 
transformed data {
  real x_r[2] = {mu, nu};
  int  x_i[1] = {N};
}

// model parameters
// this call also defines the space of which to run the Marov chain
parameters {
  real<lower=0, upper=5> delta1;
  real<lower=0, upper=5> delta2;
  real<lower=0, upper=6> theta_lambda1;
  real<lower=0, upper=6> theta_lambda2;
  real<lower=0> w1;
  real<lower=0> w2;
  real<lower=0> gamma1;
  real<lower=0> gamma2;
  real<lower=0> sigma1;
  real<lower=0> sigma2;
  real<lower=0> Ri1;
  real<lower=0> Ri2;
  real<lower=0, upper=1> R01;
  real<lower=0, upper=1> R02;
  real<lower=0, upper=0.3> R12;
  real<lower=0, upper=1> A1;
  real<lower=0, upper=1> A2;
  real<lower=20, upper=40> phi1;
  real<lower=20, upper=40> phi2;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
  real<lower=0, upper=1> rho1;
  real<lower=0, upper=1> rho2;
}

// specifying the parameters 
transformed parameters{
  real y[n_weeks, 25];
  real<lower=0> incidence_v1[n_weeks - 1];
  real<lower=0> incidence_v2[n_weeks - 1];
  {
    // specifying theta so that we can solve the SIR function 
    real theta[19];
    theta[1] = delta1;
    theta[2] = delta2;
    theta[3] = theta_lambda1;
    theta[4] = theta_lambda2;
    theta[5] = w1;
    theta[6] = w2;
    theta[7] = gamma1;
    theta[8] = gamma2;
    theta[9] = sigma1;
    theta[10] = sigma2;
    theta[11] = Ri1;
    theta[12] = Ri2;
    theta[13] = R01;
    theta[14] = R02;
    theta[15] = R12;
    theta[16] = A1;
    theta[17] = A2;
    theta[18] = phi1;
    theta[19] = phi2;

    
    // numerical optimisation to solve for y  
    y = integrate_ode_rk45(sir, y0, t0, ts, theta, x_r, x_i);
    
     // calculating weekly incidence
    for (i in 1:n_weeks-1) {
       // calculating the incidence of each infection 
        real p1 = y[i,3] + y[i,8] + y[i,13] + y[i,18] + y[i,23];
        real p2 = y[i,11] + y[i,12] + y[i,13] + y[i,14] + y[i,15];
      incidence_v1[i] = p1*rho1;
      incidence_v2[i] = p2*rho2;
    }
    
  }
}

// The model to be estimated.
model {
  //priors
  delta1 ~ uniform(0,5); // making the interaction parameters flat distibutions for now
  delta2 ~ uniform(0,5);
  theta_lambda1 ~ uniform(0,6);
  theta_lambda2 ~ uniform(0,6);
  w1 ~ normal(52,8);
  w2 ~ normal(27,8);
  gamma1 ~ gamma(3,7); //puts more weight around 3-8 days  (0.4-1.1 weeks)
  gamma2 ~ gamma(6,7); // more eight around 4-8 days (0.6-1.1 weeks)
  sigma1 ~ gamma(3,14); // puts more weight around 1-3 days (0.1-0.4 weeks)
  sigma2 ~ gamma(3,7); // puts more weight around 3-8 days (0.4-1.1 weeks)
  
  Ri1 ~ gamma(4,3); // puts more weight on Ri1 between 1-2 
  Ri2 ~ gamma(4,1.5); // puts more weight on Ri2 betewen 1-3
  R01 ~ uniform(0,1); // flat prior anywhere between 0-60% can be initially immune to v1
  R02 ~ uniform(0,1); // flat prior anywhere between 0-60% can be initially immune to v2
  R12 ~ uniform(0,0.3); // flat prior anywhere between 0-30% of the population can be initall immune to both  
  A1 ~ uniform(0,1); // amplitude to describe peak of outbreak for v1. Using flat prior
  A2 ~ uniform(0,1); 
  phi1 ~ uniform(20,40); // time till v1 outbreak in weeks from start of season. Outbreak will occur over winter between Nov - Apr this is about a 5 month period 
  phi2 ~ uniform(20,40); 
  
  rho1 ~ uniform(0,1);
  rho2 ~ uniform(0,1);

  //sampling distribution - i.e. likelihood function 
  v1_obs[1:(n_weeks-1)] ~ poisson(incidence_v1);
  v2_obs[1:(n_weeks-1)] ~ poisson(incidence_v2);
}
