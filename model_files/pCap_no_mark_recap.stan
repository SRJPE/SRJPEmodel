data {
  int Ntribs;          // number of tributaries
  int Nmr;             // number of MR experiments
  int Nwomr;           // number of strata without MR data
  int Nstrata;         // total number of strata
  int Nstrata_wc;      // strata with unmarked catch observations
  int Releases[Nmr];  // number of releases in MR experiments
  int Recaptures[Nmr]; // number of recaptures in MR experiments
  vector[Nmr] mr_flow;          // flow data for MR experiments
  vector[Nwomr] catch_flow;     // flow data for strata without MR data
  array[Nmr] int ind_trib; // tributary index for each MR experiment
  int use_trib; // tributary being modeled
  array[Nwomr] int Uind_woMR; // Indices for strata without MR data
}

parameters {
  real trib_mu_P;               // mean pCap across tributaries
  real<lower=0.01> trib_tau_P;     // precision of tributary effects
  real flow_mu_P;
  real <lower=0.01> flow_tau_P; // This could be modified as needed
  real<lower=0.01> pro_tau_P;      // precision of process error
  array[Ntribs] real b_flow; // tributary-specific flow effects
  array[Ntribs] real b0_pCap; // tributary-specific pCap intercepts
  vector[Nmr] pro_dev_P;        // process deviations for MR experiments
  vector[Nwomr] pro_dev;        // process deviations for unmarked data
}

transformed parameters {
  real pro_sd_P;
  real trib_sd_P;
  real flow_sd_P;
  array[Nmr] real logit_pCap; // Logit of pCap for MR data

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
  trib_mu_P ~ normal(0, 1e-6);
  trib_tau_P ~ gamma(0.001, 0.001);
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

generated quantities{
  array[Nstrata] real lt_pCap_U; // estimate for all weeks (Nstrata), abundance model will use Nstrata_wc
  array[Nwomr] real sim_pro_dev;

  for (i in 1:Nwomr) {
    //for weeks without efficiency trials
    sim_pro_dev[i] = normal_rng(0, pro_sd_P);
    lt_pCap_U[Uind_woMR[i]] = b0_pCap[use_trib] + b_flow[use_trib] * catch_flow[Uind_woMR[i]] + sim_pro_dev[i];
  }
}
