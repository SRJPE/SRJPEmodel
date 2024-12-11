data {
  int Ntribs;                 // number of tributaries
  int Nmr;                    // number of MR experiments
  //int Nwmr;   // Number of unmarked catch strata with MR data
  //int Nstrata; // Total number of strata
  //int Nstrata_wc;             // strata with unmarked catch observations
  array[Nmr] int Releases; // number of releases in MR experiments
  array[Nmr] int Recaptures;// number of recaptures in MR experiments
  array[Nmr] real mr_flow;// flow data for MR experiments
  array[Nmr] int ind_trib;// tributary index for each MR experiment
  //array[Nwmr] int ind_pCap;
  //array[Nstrata_wc] int Uwc_ind; //
}

parameters {
  real trib_mu_P;                      // mean pCap across tributaries
  real<lower=0.01> trib_tau_P;            // precision of tributary effects
  real flow_mu_P; //
  real<lower=0.01> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real<lower=0.01> pro_tau_P; // Process error precision
  array[Ntribs] real b0_pCap;// tributary-specific pCap intercepts
  array[Ntribs] real b_flow;// tributary-specific flow effects
  array[Nmr] real pro_dev_P;               // process deviations for MR experiments
}

transformed parameters {
  real trib_sd_P;
  real flow_sd_P;
  real pro_sd_P;
  array[Nmr] real logit_pCap;// Logit of pCap for MR data

  // Compute derived quantities
  trib_sd_P = 1/sqrt(trib_tau_P);
  flow_sd_P = 1/sqrt(flow_tau_P);
  pro_sd_P = 1 / sqrt(pro_tau_P);

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap[ind_trib[i]] + b_flow[ind_trib[i]] * mr_flow[i] + pro_dev_P[i];
  }
}

model {
  // Priors
  trib_mu_P ~ normal(0, 1000);
  trib_tau_P ~ gamma(0.001, 0.001);
  flow_mu_P ~ normal(0, 1000);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);

  // Tributary-specific parameters drawn from across-trib hyper-parameters normal distribution
  for(i in 1:Ntribs){
    b0_pCap[i] ~ normal(trib_mu_P, trib_sd_P);
    b_flow[i] ~ normal(flow_mu_P, flow_sd_P);
  }

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
