data {
  int Ntribs; // Number of tributaries
  int Nmr;    // Number of mark-recapture experiments
  int Nwomr;  // Number of unmarked catch strata without MR data
  int Nstrata; // Total number of strata
  int Nstrata_wc; // Number of strata with unmarked catch observations
  array[Nmr] int ind_trib; // Index for tributaries in MR data
  array[Nmr] int Releases;// Number of releases in MR data
  array[Nmr] int Recaptures; // Number of recaptures in MR data
  array[Nmr] real mr_flow;// Flow values for MR data
  array[Nstrata_wc] real catch_flow;// Flow values for unmarked catch strata without MR data
  array[Nwomr] int Uind_woMR;// Indices for strata without MR data
}

parameters {
  real trib_mu_P; // Hyper-parameter for tributary-specific mean pCap
  real<lower=0.01> trib_tau_P; // Hyper-parameter precision for tributary-specific mean pCap
  real flow_mu_P; //
  real<lower=0.01> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real<lower=0.01> pro_tau_P; // Process error precision
  array[Ntribs] real b_flow;// Flow effect on pCap
  array[Ntribs] real b0_pCap;// Tributary-specific mean pCap
  array[Nmr] real pro_dev_P; // Process deviations for MR data
  array[Nwomr] real pro_dev;// Process deviations for unmarked catch strata without MR data
}

transformed parameters {
  real pro_sd_P; // Standard deviation for process error
  real trib_sd_P; // Standard deviation for trib-specific mean pCap
  real flow_sd_P; // Standard deviation for flow effect on pCap
  array[Nmr] real logit_pCap;// Logit of pCap for MR data
  real logit_b0_pCap; // Logit of b0pCap for simulated data
  real logit_b0_flow;

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
  // Priors on variance terms (in units of precision)
  trib_mu_P ~ normal(0, 1000);
  trib_tau_P ~ gamma(0.001, 0.001);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);

  logit_b0_pCap ~ normal(trib_mu_P, trib_sd_P);
  logit_b0_flow ~ normal(flow_mu_P, flow_sd_P);

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
  array[Nstrata] real lt_pCap_U; // estimate for all weeks (Nstrata), abundance model will use Nstrata_wc
  array[Nwomr] real sim_pro_dev;
  array[Nmr] real log_lik; //for loo

  for (i in 1:Nwomr) {
    //for weeks without efficiency trials
    sim_pro_dev[i] = normal_rng(0, pro_sd_P);
    lt_pCap_U[Uind_woMR[i]] = logit_b0_pCap + logit_b0_flow * catch_flow[Uind_woMR[i]] + sim_pro_dev[i];
  }

  // for loo analysis
  for (i in 1:Nmr) {
    log_lik[i] = binomial_logit_lpmf(Recaptures[i] | Releases[i], logit_pCap[i]);
  }
}