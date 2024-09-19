data {
  int Ntribs; // Number of tributaries
  int Nmr;    // Number of mark-recapture experiments
  int Nwomr;  // Number of unmarked catch strata without MR data
  int Nstrata; // Total number of strata
  int Nstrata_wc; // Number of strata with unmarked catch observations\
  int ind_trib[Nmr]; // Index for tributaries in MR data
  real mr_flow[Nmr]; // Flow values for MR data
  int Releases[Nmr]; // Number of releases in MR data
  int Recaptures[Nmr]; // Number of recaptures in MR data
  real catch_flow[Nstrata_wc]; // Flow values for unmarked catch strata without MR data
  int u[Nstrata_wc]; // Unmarked catch observations
  int Uwc_ind[Nstrata_wc]; // Indices for unmarked catch strata
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
  real<lower=0.01> tau_N; // Precision for spline coefficients
  real<lower=0.01> tau_Ne; // Extra-spline variation precision
  real b_flow[Ntribs]; // Flow effect on pCap
  real b0_pCap[Ntribs]; // Tributary-specific mean pCap
  vector[K] b_sp; // Spline coefficients
  real pro_dev_P[Nmr]; // Process deviations for MR data
  real pro_dev[Nwomr]; // Process deviations for unmarked catch strata without MR data
  real lg_N[Nstrata]; // Log abundance estimates
}

transformed parameters {
  real pro_sd_P; // Standard deviation for process error
  real trib_sd_P; // Standard deviation for trib-specific mean pCap
  real flow_sd_P; // Standard deviation for flow effect on pCap
  real sd_N; // Standard deviation for spline coefficients
  real sd_Ne; // Standard deviation for extra-spline variation
  real logit_pCap[Nmr]; // Logit of pCap for MR data
  real logit_pCap_Sim[Nwomr]; // Logit of pCap for simulated data
  real logit_b0_pCap; // Logit of b0pCap for simulated data
  real logit_b0_flow;
  real Usp[Nstrata]; // Spline-based estimate of log U
  real pCap_U[Nstrata]; // Estimated pCaps for all strata
  vector[Nstrata] N; // Abundance estimates

  // Calculate logit of pCap for MR data - can't do this

  // Calculate logit of pCap for simulated unmarked catch strata without MR data
  for (i in 1:Nwomr) {
    logit_pCap_Sim[i] = logit_b0_pCap + logit_b0_flow * catch_flow[Uind_woMR[i]] + pro_dev[i];
    pCap_U[Uind_woMR[i]] = inv_logit(logit_pCap_Sim[i]); // assigns simulated pCap for weeks with no MRdata
  }

  // estimate weekly pCap for weeks in the site and year - can't do this

  // Spline-based estimate of log U
  for (i in 1:Nstrata) {
    Usp[i] = dot_product(ZP[i,], b_sp);
    N[i] = exp(lg_N[i]) * 1000.0;
  }
}

model {
  // Priors on variance terms (in units of precision)
  trib_mu_P ~ normal(0, 1000);
  trib_tau_P ~ gamma(0.001, 0.001);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);
  tau_N ~ gamma(1, 0.05);
  tau_Ne ~ gamma(1, 0.05);

  logit_b0_pCap ~ normal(trib_mu_P, trib_sd_P);
  logit_b0_flow ~ normal(flow_mu_P, flow_sd_P);

  // Spline coefficient priors
  for (i in 3:K) {
    real xi = 2 * b_sp[i-1] - b_sp[i-2];
    b_sp[i] ~ normal(xi, sd_N);//replace tau_N with sd_N
  }

  // simulated process error deviations for weeks without MR observations
  for(i in 1:Nwomr){
    pro_dev[i] ~ normal(0, pro_sd_P);
  }

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

  // Simulating log abundance estimates for each week (with bounds), mean is defined by spline, sd_Ne is extra-spline variation
  for(i in 1:Nstrata){
    lg_N[i] ~ normal(Usp[i], sd_Ne) T[,lgN_max[i]]; // temporary
  }

  // Likelihood for unmarked catch observations
  for (i in 1:Nstrata_wc) {
      real bcl=lgamma(N[Uwc_ind[i]]+1)-lgamma(u[Uwc_ind[i]]+1)-lgamma(N[Uwc_ind[i]]-u[Uwc_ind[i]]);//log of binomial coefficent
      real kern=u[Uwc_ind[i]]*log(pCap_U[Uwc_ind[i]]) + (N[i]-u[i])*log(1-pCap_U[Uwc_ind[i]]); //log of binomial kernal
      target += bcl + kern;//log of binomial probability
  }
}

generated quantities {
  real Ntot = sum(N);
}
