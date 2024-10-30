data {
  int Ntribs; // Number of tributaries
  int Nmr;    // Number of mark-recapture experiments
  int Nwmr;   // Number of unmarked catch strata with MR data
  int Nwomr;  // Number of unmarked catch strata without MR data
  int Nstrata; // Total number of strata
  int Nstrata_wc; // Number of strata with unmarked catch observations
  int use_trib; // index of site being estimated
  int ind_trib[Nmr]; // Index for tributaries in MR data
  int ind_pCap[Nwmr];
  real mr_flow[Nmr]; // Flow values for MR data
  int Releases[Nmr]; // Number of releases in MR data
  int Recaptures[Nmr]; // Number of recaptures in MR data
  real catch_flow[Nstrata_wc]; // Flow values for unmarked catch strata without MR data
  int u[Nstrata_wc]; // Unmarked catch observations
  int Uwc_ind[Nstrata_wc]; // Indices for unmarked catch strata
  int Uind_wMR[Nwmr]; // Indices for strata with MR data
  int Uind_woMR[Nwomr]; // Indices for strata without MR data
  int K;                 // Number of columns in the bspline basis matrix
  matrix[Nstrata, K] ZP; // Design matrix for spline estimation
  real lgN_max[Nstrata]; // Maximum values for log abundances
}

parameters {
  real trib_mu_P; // Hyper-parameter for tributary-specific mean pCap
  real<lower=0.01> trib_tau_P; // Hyper-parameter precision for tributary-specific mean pCap
  real flow_mu_P; //
  real<lower=0.01> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real<lower=0.01> pro_tau_P; // Process error precision
  real b_flow[Ntribs]; // Flow effect on pCap
  real b0_pCap[Ntribs]; // Tributary-specific mean pCap
  real pro_dev_P[Nmr]; // Process deviations for MR data
}

transformed parameters {
  real pro_sd_P; // Standard deviation for process error
  real trib_sd_P; // Standard deviation for trib-specific mean pCap
  real flow_sd_P; // Standard deviation for flow effect on pCap
  real logit_pCap[Nmr]; // Logit of pCap for MR data

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
    Recaptures[i] ~ binomial(Releases[i], inv_logit(logit_pCap[i]));
  }

}
