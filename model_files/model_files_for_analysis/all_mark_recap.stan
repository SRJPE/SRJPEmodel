data {
  int Nmr;               // Number of MR experiments
  int Ntribs;            // Number of tributaries
  int ind_trib[Nmr]; // Indices for tributaries in MR experiments
  int Recaptures[Nmr];   // Recapture counts in MR experiments
  int Releases[Nmr];     // Release counts in MR experiments
  int Nstrata;           // Number of strata
  int Nstrata_wc;        // Number of strata with unmarked catch observations
  int u[Nstrata_wc];     // Unmarked catch observations
  int K;                 // Number of columns in the bspline basis matrix
  matrix[Nstrata, K] ZP;          // Design matrix for spline-based estimation
  int ind_pCap[Nmr]; // Strata indices for pCap
  vector[Nmr] mr_flow;            // Flow values for MR experiments
  int Uwc_ind[Nstrata_wc]; // Indices for unmarked catch strata
  vector[Nstrata] lgN_max;        // Upper bounds for lg_N
}

parameters {
  real trib_mu_P;
  real<lower=0> trib_tau_P;
  real flow_mu_P;
  real<lower=0> flow_tau_P;
  real<lower=0> pro_tau_P;
  real<lower=0> tau_N;
  real<lower=0> tau_Ne;
  vector[Ntribs] b0_pCap;
  vector[Ntribs] b_flow;
}

transformed parameters {
  // abundance
  vector[Nstrata] Usp;
  vector[Nmr] logit_pCap;
  vector[Nstrata] N;
  vector[Nstrata] pCap_U;
  vector[Nstrata] lg_N;
  // pCap
  vector<lower=0, upper=1>[Nmr] pCap;
  // priors
  real<lower=0> trib_sd_P = 1 / sqrt(trib_tau_P);
  real<lower=0> flow_sd_P = 1 / sqrt(flow_tau_P);
  vector[Nmr] pro_dev_P;
  real<lower=0> pro_sd_P = 1 / sqrt(pro_tau_P);
  real<lower=0> sd_N = 1 / sqrt(tau_N);
  real<lower=0> sd_Ne = 1 / sqrt(tau_Ne);
  vector[K] b_sp;

  // Trib and week effects only
  // Transform logit pCap (parameter) to 0:1 bounded
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap[ind_trib[i]] + b_flow[ind_trib[i]] * mr_flow[i] + pro_dev_P[i];
    pCap[i] = inv_logit(logit_pCap[i]);
  }
  // Spline-based estimate of log U
  Usp = ZP * b_sp; // real Usp = dot_product(ZP[i], b_sp);

  // unmarked abundance
  for (i in 1:Nstrata) {
    lg_N[i] += normal_lcdf(lgN_max[i] | Usp[i], sd_Ne);
    N[i] = round(exp(lg_N[i]) * 1000);
    pCap_U[i] = pCap[ind_pCap[i]];
  }
}

model {
  // Priors for hyperparameters
  trib_mu_P ~ normal(0, 1.0E-06);
  trib_tau_P ~ gamma(0.001, 0.001);
  flow_mu_P ~ normal(0, 1.0E-06);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);
  tau_N ~ gamma(1, 0.05);
  tau_Ne ~ gamma(1, 0.05);

  // Priors for spline coefficients
  b_sp[1] ~ normal(0, 1.0);
  b_sp[2] ~ normal(0, 10); // flat prior can be approximated with a wide normal

  // Likelihood for spline
  for (i in 3:K) {
    b_sp[i] ~ normal(2 * b_sp[i - 1] - b_sp[i - 2], tau_N);
  }

  // Likelihood for flow covariates
  for (j in 1:Ntribs) {
    b0_pCap[j] ~ normal(trib_mu_P, trib_sd_P);
    b_flow[j] ~ normal(flow_mu_P, flow_sd_P);
  }

  // Likelihood for Recaptures
  for (i in 1:Nmr) {
    pro_dev_P[i] ~ normal(0, pro_sd_P);
    Recaptures[i] ~ binomial(Releases[i], pCap[i]);
  }

  // Likelihood for unmarked catch observations
  //for (i in 1:Nstrata_wc) {
   // u[Uwc_ind[i]] ~ binomial(N[Uwc_ind[i]], pCap_U[Uwc_ind[i]]);
  //}
}

generated quantities {
  real Ntot = sum(N);
}
