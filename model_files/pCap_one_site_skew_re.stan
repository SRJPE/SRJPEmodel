// model for mainstem - no hierarchy
data {
  int Nmr; // number of MR experiments
  array[Nmr] int Releases; // number of releases in MR experiments
  array[Nmr] int Recaptures; // number of recaptures in MR experiments
  array[Nmr] real mr_flow; // flow data for MR experiments
  array[Nmr] int ind_yr;
  array[Nmr] real effort;
  int Nyr_re;
}

parameters {
  real b0_pCap; // pCap intercept
  real b_flow; // flow effects
  real b_eff; // effort effects
  real<lower=0> pro_sd_P;
  real<lower=0> yr_sd_P;

  array[Nyr_re] real<lower=-3,upper=3> yr_re;
  array[Nmr] real pro_dev_P; // process deviations for MR experiments

  real <lower=0>alpha;     // Shape/skewness parameter for process error random effects
  //real alpha; didn't coverge
}

transformed parameters {
  array[Nmr] real logit_pCap;// Logit of pCap for MR data

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap + b_flow * mr_flow[i] + yr_re[ind_yr[i]]  + pro_dev_P[i] + b_eff * effort[i];
  }
}

model {
  alpha ~ normal(0, 10);
  pro_sd_P ~ student_t(3, 0, 2.5);
  yr_sd_P ~ student_t(3, 0, 2.5);

  yr_re~normal(0,yr_sd_P);//year effect deviates
  pro_dev_P[] ~ skew_normal(0, pro_sd_P, alpha);//process error deviation

  // Likelihood for MR data
  for (i in 1:Nmr) {
    Recaptures[i] ~ binomial_logit(Releases[i], logit_pCap[i]);
  }
}

generated quantities {
  array[Nmr] real log_lik; //for loo

  for (i in 1:Nmr) {// for loo analysis
    log_lik[i] = binomial_logit_lpmf(Recaptures[i] | Releases[i], logit_pCap[i]);
  }
}
