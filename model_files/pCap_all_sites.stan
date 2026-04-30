data {
  int Ntribs;                 // number of tributaries
  int Nmr;                    // number of MR experiments
  array[Ntribs] int use_in_hyper; // allows user to remove trib from hyper for b0_pCap and b_flow
  array[Nmr] int Releases; // number of releases in MR experiments
  array[Nmr] int Recaptures;// number of recaptures in MR experiments
  array[Nmr] real effort;
  array[Nmr] real mr_flow;// flow data for MR experiments
  array[Nmr] int ind_trib;// tributary index for each MR experiment

  array[Nmr] int ind_yr;
  int Nyr_re;
  array[Nyr_re] int sd_yr_ind;
}

parameters {
  real trib_mu_P;                      // mean pCap across tributaries
  real<lower=0> trib_sd_P;            // precision of tributary effects
  real flow_mu_P; //
  real<lower=0> flow_sd_P; // Hyper-parameter precision for flow effect on pCap
  array[Ntribs] real <lower=0> pro_sd_P;
  array[Ntribs] real<lower = -10, upper = 10> b0_pCap;// tributary-specific pCap intercepts
  array[Ntribs] real b_flow;// tributary-specific flow effects
  array[Ntribs] real b_eff; // tributary-specific effort effects
  array[Nmr] real pro_dev_P;               // process deviations for MR experiments

  array[Ntribs] real <lower=0, upper=2> yr_sd_P;//max sd of 2.0 for sites with limited years of MR
  array[Nyr_re] real<lower=-3, upper=3> yr_re;
}

transformed parameters {
   array[Nmr] real logit_pCap;// Logit of pCap for MR data

  // Calculate logit of pCap for MR data
  for (i in 1:Nmr) {
    logit_pCap[i] = b0_pCap[ind_trib[i]] + b_flow[ind_trib[i]] * mr_flow[i] + yr_re[ind_yr[i]] + pro_dev_P[i] + b_eff[ind_trib[i]] * effort[i];
  }
}

model {
  // Priors
  trib_mu_P ~ normal(0, 1000);
  flow_mu_P ~ normal(0, 1000);

  // Tributary-specific parameters drawn from across-trib hyper-parameters normal distribution
  for(i in 1:Ntribs){
    // remove selected tributaries from estimation of b0_pCap and b_flow
    if(use_in_hyper[i] == 1) {
            b0_pCap[i] ~ normal(trib_mu_P, trib_sd_P); // should not include okie dam and steep riffle
            b_flow[i] ~ normal(flow_mu_P, flow_sd_P);
      }
       pro_sd_P[i] ~ student_t(3, 0, 2.5);
       yr_sd_P[i] ~ student_t(3, 0, 2.5);
    }

  // Likelihood for MR data
  for (i in 1:Nmr) {
    pro_dev_P[i] ~ normal(0, pro_sd_P[ind_trib[i]]);
    Recaptures[i] ~ binomial_logit(Releases[i], logit_pCap[i]);
  }

  for (i in 1:Nyr_re) {
     yr_re[i] ~ normal(0, yr_sd_P[sd_yr_ind[i]]);
  }
}

generated quantities {
  array[Nmr] real log_lik; //for loo

  // for loo analysis
  for (i in 1:Nmr) {
    log_lik[i] = binomial_logit_lpmf(Recaptures[i] | Releases[i], logit_pCap[i]);
  }
}
