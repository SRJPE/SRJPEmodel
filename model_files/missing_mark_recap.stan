data {
  int Ntribs; // Number of tributaries
  int Nmr;    // Number of mark-recapture experiments
  int Nwmr;   // Number of unmarked catch strata with MR data
  int Nwomr;  // Number of unmarked catch strata without MR data
  int Nstrata; // Total number of strata
  int Nstrata_wc; // Number of strata with unmarked catch observations
  int use_trib; // index of site being estimated
  array[Nmr] int ind_trib; // Index for tributaries in MR data
  array[Nwmr] int ind_pCap;
  array[Nmr] real mr_flow;// Flow values for MR data
  array[Nmr] int Releases; // Number of releases in MR data
  array[Nmr] int Recaptures; // Number of recaptures in MR data
  array[Nstrata_wc] real catch_flow; // Flow values for unmarked catch strata without MR data
  array[Nstrata_wc] int u; // Unmarked catch observations
  array[Nstrata_wc] int Uwc_ind; // Indices for unmarked catch strata
  array[Nwmr] int Uind_wMR;// Indices for strata with MR data
  array[Nwomr] int Uind_woMR; // Indices for strata without MR data
  int K;                 // Number of columns in the bspline basis matrix
  array[Nstrata, K] real ZP; // design matrix for splines
  array[Nstrata] real lgN_max;// Maximum values for log abundances
}

parameters {
  real trib_mu_P; // Hyper-parameter for tributary-specific mean pCap
  real<lower=0.001> trib_tau_P; // Hyper-parameter precision for tributary-specific mean pCap
  real flow_mu_P; //
  real<lower=0.001> flow_tau_P; // Hyper-parameter precision for flow effect on pCap
  real<lower=0.001> pro_tau_P; // Process error precision
  real<lower=0.0> tau_N; // Precision for spline coefficients
  real<lower=0.0> tau_Ne; // Extra-spline variation precision
  array[Ntribs] real b_flow; // Flow effect on pCap
  array[Ntribs] real b0_pCap; // Tributary-specific mean pCap
  array[K] real b_sp; // Spline coefficients
  array[Nmr] real pro_dev_P; // Process deviations for MR data
  array[Nwomr] real pro_dev; // Process deviations for unmarked catch strata without MR data
  array[Nstrata] real lg_N; // Log abundance estimates
}

transformed parameters {
  real pro_sd_P; // Standard deviation for process error
  real trib_sd_P; // Standard deviation for trib-specific mean pCap
  real flow_sd_P; // Standard deviation for flow effect on pCap
  real sd_N; // Standard deviation for spline coefficients
  real sd_Ne; // Standard deviation for extra-spline variation
  array[Nmr] real logit_pCap; // Logit of pCap for MR data
  array[Nwomr] real logit_pCap_Sim; // Logit of pCap for simulated data
  array[Nstrata] real Usp; // Spline-based estimate of log U
  array[Nstrata] real pCap_U; // Estimated pCaps for all strata
  vector[Nstrata] N; // Abundance estimates

  // Compute derived quantities
  trib_sd_P = 1/sqrt(trib_tau_P);
  flow_sd_P = 1/sqrt(flow_tau_P);
  pro_sd_P = 1 / sqrt(pro_tau_P);
  sd_N = 1 / sqrt(tau_N);
  sd_Ne = 1 / sqrt(tau_Ne);

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap[ind_trib[i]] + b_flow[ind_trib[i]] * mr_flow[i] + pro_dev_P[i];
  }
  // estimate weekly pCap for weeks in the site and year
  for(i in 1:Nwmr){
    pCap_U[Uind_wMR[i]] = inv_logit(logit_pCap[ind_pCap[i]]); // Assign estimated pCaps to U estimation strata (weeks with catch)
  }

  // Calculate logit of pCap for simulated unmarked catch strata without MR data
  for (i in 1:Nwomr) {
    logit_pCap_Sim[i] = b0_pCap[use_trib] + b_flow[use_trib] * catch_flow[Uind_woMR[i]] + pro_dev[i];
    pCap_U[Uind_woMR[i]] = inv_logit(logit_pCap_Sim[i]); // assigns simulated pCap for weeks with no MRdata
  }

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

  // Spline coefficient priors
  real xi;
  for (i in 3:K) {
    xi = 2 * b_sp[i-1] - b_sp[i-2];
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
  real bcl;
  real kern;
  for (i in 1:Nstrata_wc) {
      bcl=lgamma(N[Uwc_ind[i]]+1)-lgamma(u[Uwc_ind[i]]+1)-lgamma(N[Uwc_ind[i]]-u[Uwc_ind[i]]);//log of binomial coefficent
      kern=u[Uwc_ind[i]]*log(pCap_U[Uwc_ind[i]]) + (N[i]-u[i])*log(1-pCap_U[Uwc_ind[i]]); //log of binomial kernal
      target += bcl + kern;//log of binomial probability
  }
}

generated quantities{
  real Ntot = sum(N[]);
  real lg_Ntot = sum(lg_N[]);
}
