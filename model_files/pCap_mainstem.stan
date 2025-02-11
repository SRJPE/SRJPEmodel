// model for mainstem - no hierarchy
data {
  int Nmr; // number of MR experiments
  array[Nmr] int Releases; // number of releases in MR experiments
  array[Nmr] int Recaptures; // number of recaptures in MR experiments
  array[Nmr] real mr_flow; // flow data for MR experiments
}

parameters {
  real flow_mu_P; //
  real <lower = 0.01> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real <lower = 0.01> pro_tau_P; // Process error precision
  real b0_pCap; // pCap intercept
  real b_flow; // flow effects
  array[Nmr] real pro_dev_P; // process deviations for MR experiments
}

transformed parameters {
  real flow_sd_P;
  real pro_sd_P;
  array[Nmr] real logit_pCap;// Logit of pCap for MR data

  // Compute derived quantities
  flow_sd_P = 1/sqrt(flow_tau_P);
  pro_sd_P = 1 / sqrt(pro_tau_P);

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap + b_flow * mr_flow[i] + pro_dev_P[i];
  }
}

model {
  // Priors
  flow_mu_P ~ normal(0, 1000);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);

  // Likelihood for MR data
  for (i in 1:Nmr) {
    pro_dev_P[i] ~ normal(0, pro_sd_P);
    Recaptures[i] ~ binomial_logit(Releases[i], logit_pCap[i]);
  }
}

generated quantities {
  array[Nmr] real log_lik; //for loo

  // for loo analysis
  for (i in 1:Nmr) {
    log_lik[i] = binomial_logit_lpmf(Recaptures[i] | Releases[i], logit_pCap[i]);
  }
}
