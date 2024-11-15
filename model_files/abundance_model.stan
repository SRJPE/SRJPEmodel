data {
  int Nstrata; // Total number of strata
  int Nstrata_wc; // Number of strata with unmarked catch observations
  int use_trib; // index of site being estimated
  real lt_pCap_mu;
  real lt_pCap_sd;
  array[Nstrata_wc] int u; // Unmarked catch observations
  array[Nstrata_wc] int Uwc_ind; // Indices for unmarked catch strata
  int K;                 // Number of columns in the bspline basis matrix
  array[Nstrata, K] real ZP; // design matrix for splines
  array[Nstrata] real lgN_max;// Maximum values for log abundances
}

parameters {
  real<lower=0.0> tau_N; // Precision for spline coefficients
  real<lower=0.0> tau_Ne; // Extra-spline variation precision
  array[K] real b_sp; // Spline coefficients
  array[Nstrata] real lg_N; // Log abundance estimates
}

transformed parameters {
  real sd_N; // Standard deviation for spline coefficients
  real sd_Ne; // Standard deviation for extra-spline variation
  array[Nstrata] real Usp; // Spline-based estimate of log U
  vector[Nstrata] N; // Abundance estimates
  array[Nstrata_wc] real lt_pCap_U;
  array[Nstrata_wc] real pCap_U;

  // Compute derived quantities
  sd_N = 1 / sqrt(tau_N);
  sd_Ne = 1 / sqrt(tau_Ne);

  for(i in 1:Nstrata_wc) {
      pCap_U[i] = inv_logit(lt_pCap_U[i]);
  }

  // Spline-based estimate of log U
  for (i in 1:Nstrata) {
    Usp[i] = dot_product(ZP[i,], b_sp);
    N[i] = exp(lg_N[i]) * 1000.0;
  }
}

model {
  // Priors on variance terms (in units of precision)
  tau_N ~ gamma(1, 0.05);
  tau_Ne ~ gamma(1, 0.05);

  // Spline coefficient priors
  // we do not need to specify priors for b_sp[1] and b_sp[2] as STAN by default assumes an improper uniform
  real xi;
  for (i in 3:K) {
    xi = 2 * b_sp[i-1] - b_sp[i-2];
    b_sp[i] ~ normal(xi, sd_N);//replace tau_N with sd_N
  }

  // Simulating log abundance estimates for each week (with bounds), mean is defined by spline, sd_Ne is extra-spline variation
  for(i in 1:Nstrata){
    lg_N[i] ~ normal(Usp[i], sd_Ne) T[,lgN_max[i]]; // temporary
  }

  // Likelihood for unmarked catch observations
  real bcl;
  real kern;
  for (i in 1:Nstrata_wc) {
    // simulate pCap_U from mean and sd (simulate in logit space and transform using inv_logit)
    lt_pCap_U[i] ~ normal(lt_pCap_mu, lt_pCap_sd);

    bcl=lgamma(N[Uwc_ind[i]]+1)-lgamma(u[Uwc_ind[i]]+1)-lgamma(N[Uwc_ind[i]]-u[Uwc_ind[i]]);//log of binomial coefficent
    kern=u[Uwc_ind[i]]*log(pCap_U[Uwc_ind[i]]) + (N[i]-u[i])*log(1-pCap_U[Uwc_ind[i]]); //log of binomial kernal
    target += bcl + kern;//log of binomial probability
  }
}

generated quantities{
  real Ntot = sum(N[]);
  real lg_Ntot = sum(lg_N[]);
}
