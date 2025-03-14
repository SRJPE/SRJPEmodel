data {
  int N;
  int input_years[N];
  int observed_passage[N];
  int observed_spawners[N]; // this is redd count but also holding count
  real percent_female;
  real environmental_covar[N];
  real ss_total;
  real average_upstream_passage;
}

parameters {
  real log_mean_redds_per_spawner; //no lower constraint of 0 since in log space
  real <lower = 0> sigma_redds_per_spawner;
  real b1_survival;
  real log_redds_per_spawner[N];
}

transformed parameters{

  real conversion_rate[N];
  real redds_per_spawner[N];
  real predicted_spawners[N];
  real mean_redds_per_spawner = exp(log_mean_redds_per_spawner) * percent_female;

  for(i in 1:N) {
    conversion_rate[i] = exp(log_redds_per_spawner[i] + b1_survival * environmental_covar[i]);

    // predicted redds is product of observed passage * conversion rate
    predicted_spawners[i] = observed_passage[i] * conversion_rate[i];
    redds_per_spawner[i] = exp(log_redds_per_spawner[i]); //for ouptut only
  }
}

model {
  log_redds_per_spawner ~ normal(log_mean_redds_per_spawner, sigma_redds_per_spawner);
  observed_spawners ~ poisson(predicted_spawners);
}

generated quantities {

  // report years for ease later
  int years[N];

  // calculate R2 fixed effects
  real y_pred[N];
  real RE_draws[N]; // vector for storing random draws for RE
  real var_fit = 0;
  real var_res = 0.0;
  real ss_res = 0.0;

  for(i in 1:N) {
    years[i] = input_years[i];
    RE_draws[i] = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);
    y_pred[i] = observed_passage[i] * exp(RE_draws[i]);
  }

  real y_pred_mu = mean(y_pred[]);
  real R2_data;
  real R2_fixed;
  array[N] real log_lik; //for look

  // calculate R2_data and R2_fixed
  for(i in 1:N) {
    var_fit = var_fit + (y_pred[i] - y_pred_mu)^2;
    // TODO var_res is actually based off y_pred with full env effects?
    // TODO and then subtract y_pred from that?
    var_res = var_res + (observed_spawners[i] - y_pred[i])^2;
    ss_res = ss_res + (observed_spawners[i] - predicted_spawners[i])^2;

    // for loo analysis
    log_lik[i] = poisson_lpmf(observed_spawners[i] | predicted_spawners[i]);
  }

  R2_data = 1.0 - (ss_res / ss_total);
  R2_fixed = 1.0 - (var_fit / (var_fit + var_res));

  // forecast error in spawner abundance at average upstream passage abundance
  real spawner_abundance_forecast[2];
  real rps_forecast = normal_rng(log_mean_redds_per_spawner, sigma_redds_per_spawner);// random draw from hyper

  spawner_abundance_forecast[1] = average_upstream_passage * exp(rps_forecast); //for dummy variable = 0 condition, or at average covariate conditions
  spawner_abundance_forecast[2] = average_upstream_passage * exp(rps_forecast + b1_survival); //for dummy variable=1 condition, or at 1 SD from mean covariate condition
}
