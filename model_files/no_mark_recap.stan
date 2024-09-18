data {
  int<lower=0> Ntribs;          // number of tributaries
  int<lower=0> Nmr;             // number of MR experiments
  int<lower=0> Nwomr;           // number of strata without MR data
  int<lower=0> Nstrata;         // total number of strata
  int<lower=0> Nstrata_wc;      // strata with unmarked catch observations
  int<lower=0> Releases[Nmr];  // number of releases in MR experiments
  int<lower=0> Recaptures[Nmr]; // number of recaptures in MR experiments
  int u[Nstrata_wc]; // Unmarked catch observations
  int Uind_woMR[Nwomr]; // Indices for strata without MR data
  vector[Nmr] mr_flow;          // flow data for MR experiments
  vector[Nwomr] catch_flow;     // flow data for strata without MR data
  int<lower=1, upper=Ntribs> ind_trib[Nmr]; // tributary index for each MR experiment
  int<lower=1, upper=Ntribs> use_trib; // tributary being modeled
  int K;                 // Number of columns in the bspline basis matrix
  matrix[Nstrata, K] ZP;        // design matrix for splines
  real lgN_max[Nstrata_wc];     // upper bound for logN in observations
  int<lower=0> Uwc_ind[Nstrata_wc]; // indices for unmarked catch observations
}

parameters {
  real trib_mu_P;               // mean pCap across tributaries
  real<lower=0.01> trib_tau_P;     // precision of tributary effects
  real flow_mu_P;
  real <lower=0.01> flow_tau_P; // This could be modified as needed
  real<lower=0.01> pro_tau_P;      // precision of process error
  real<lower=0.01> tau_N;          // precision for spline parameters
  real<lower=0.01> tau_Ne;         // precision for extra-spline variation
  real b_flow[Ntribs];        // tributary-specific flow effects
  real b0_pCap[Ntribs];       // tributary-specific pCap intercepts
  vector[K] b_sp;               // spline coefficients
  vector[Nmr] pro_dev_P;        // process deviations for MR experiments
  vector[Nwomr] pro_dev;        // process deviations for unmarked data
  vector[Nstrata] lg_N;         // log abundance estimates
}

transformed parameters {
  real pro_sd_P;
  real trib_sd_P;
  real flow_sd_P;
  real sd_N; // Standard deviation for spline coefficients
  real sd_Ne; // Standard deviation for extra-spline variation
  vector[Nstrata] N;            // abundance estimates
  real pCap_U[Nstrata];       // pCap estimates for strata without MR data
  real Usp[Nstrata]; // Spline-based estimate of log U
  real logit_pCap[Nmr]; // Logit of pCap for MR data
  real logit_pCap_Sim[Nwomr]; // Logit of pCap for simulated data

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

  // no weeks with mark recap, so move on to simulation
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
  // Priors
  trib_mu_P ~ normal(0, 1e-6);
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

generated quantities{
  real Ntot=sum(N[]);
}
