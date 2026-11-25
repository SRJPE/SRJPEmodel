// model for mainstem - no hierarchy
data {
  int Nmr; // number of MR experiments
  array[Nmr] int Releases; // number of releases in MR experiments
  array[Nmr] int Recaptures; // number of recaptures in MR experiments
  array[Nmr] real mr_flow; // flow data for MR experiments
  array[Nmr] int ind_yr;
  int Nyr_re;
}

parameters {
  real <lower = 0.01> pro_tau_P; // Process error precision
  real b0_pCap; // pCap intercept
  real b_flow; // flow effects
  array[Nmr] real pro_dev_P; // process deviations for MR experiments
  real <lower=0.25> yr_tau_P;//max sd of 2.0 for sites with limited years of MR
  array[Nyr_re] real<lower=-3,upper=3> yr_re;
}

transformed parameters {
  real pro_sd_P;
  real yr_sd_P;
  array[Nmr] real logit_pCap;// Logit of pCap for MR data

  // Compute derived quantities
  pro_sd_P = 1 / sqrt(pro_tau_P);
  yr_sd_P = 1 / sqrt(yr_tau_P);

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap + b_flow * mr_flow[i] +  +yr_re[ind_yr[i]]  + pro_dev_P[i];
  }
}

model {
  // Priors
  //b0_pCap ~ normal(0, 1.0E-06);
  //b_flow ~ normal(0, 1.0E-06);

  pro_tau_P ~ gamma(0.001, 0.001);
  yr_tau_P ~ gamma(0.001, 0.001);

  yr_re~normal(0,yr_sd_P);

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
