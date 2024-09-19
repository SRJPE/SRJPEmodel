data {
  int Ntribs;                 // number of tributaries
  int Nmr;                    // number of MR experiments
  int Nstrata;                // total number of strata
  int Nstrata_wc;             // strata with unmarked catch observations
  int ind_pCap[Nstrata];      // indices for pCap corresponding to strata
  int Releases[Nmr];          // number of releases in MR experiments
  int Recaptures[Nmr];        // number of recaptures in MR experiments
  int u[Nstrata_wc];
  real mr_flow [Nmr];                 // flow data for MR experiments
  int K;                 // Number of columns in the bspline basis matrix
  matrix[Nstrata, K] ZP;               // design matrix for splines
  real lgN_max[Nstrata];               // upper bound for logN in observations
  int ind_trib[Nmr]; // tributary index for each MR experiment
  int Uwc_ind[Nstrata_wc]; // indices for unmarked catch observations
}

parameters {
  real trib_mu_P;                      // mean pCap across tributaries
  real<lower=0.01> trib_tau_P;            // precision of tributary effects
  real flow_mu_P; //
  real<lower=0.01> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real<lower=0.01> pro_tau_P; // Process error precision
  real<lower=0.01> tau_N; // Precision for spline coefficients
  real<lower=0.01> tau_Ne; // Extra-spline variation precision
  real b0_pCap[Ntribs];              // tributary-specific pCap intercepts
  real b_flow[Ntribs];               // tributary-specific flow effects
  vector[K] b_sp;                      // spline coefficients
  vector[Nmr] pro_dev_P;               // process deviations for MR experiments
  vector[Nstrata] lg_N;                // log abundance estimates
}

transformed parameters {
  real trib_sd_P;
  real flow_sd_P;
  real pro_sd_P;
  real sd_N; // Standard deviation for spline coefficients
  real sd_Ne; // Standard deviation for extra-spline variation
  real logit_pCap[Nmr]; // Logit of pCap for MR data
  real Usp[Nstrata]; // Spline-based estimate of log U
  real pCap_U[Nstrata]; // Estimated pCaps for all strata
  vector[Nstrata] N; // Abundance estimates

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap[ind_trib[i]] + b_flow[ind_trib[i]] * mr_flow[i] + pro_dev_P[i];
  }

  // estimate weekly pCap for weeks in the site and year
  for(i in 1:Nstrata){
    pCap_U[i] = inv_logit(logit_pCap[ind_pCap[i]]); // Assign estimated pCaps to U estimation strata (weeks with catch)
  }

  // Don't need to calculate logit of pCap for simulated unmarked catch strata without MR data

  // Spline-based estimate of log U
  for (i in 1:Nstrata) {
    Usp[i] = dot_product(ZP[i,], b_sp);
    N[i] = exp(lg_N[i]) * 1000.0;
  }
}

model {
  // Priors
  trib_mu_P ~ normal(0, 1000);
  trib_tau_P ~ gamma(0.001, 0.001);
  flow_tau_P ~ gamma(0.001, 0.001);
  pro_tau_P ~ gamma(0.001, 0.001);
  tau_N ~ gamma(1, 0.05);
  tau_Ne ~ gamma(1, 0.05);

  // Spline coefficient priors
  for (i in 3:K) {
    real xi = 2 * b_sp[i-1] - b_sp[i-2];
    b_sp[i] ~ normal(xi, sd_N);//replace tau_N with sd_N
  }

  // Don't need to simulate process error deviations for weeks without MR observations

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
      real bcl = lgamma(N[Uwc_ind[i]]+1) - lgamma(u[Uwc_ind[i]]+1) - lgamma(N[Uwc_ind[i]] - u[Uwc_ind[i]]);//log of binomial coefficent
      real kern = u[Uwc_ind[i]] * log(pCap_U[Uwc_ind[i]]) + (N[i]-u[i]) * log(1 - pCap_U[Uwc_ind[i]]); //log of binomial kernal
      target += bcl + kern;//log of binomial probability
  }
}

generated quantities {
  real Ntot = sum(N);
}
