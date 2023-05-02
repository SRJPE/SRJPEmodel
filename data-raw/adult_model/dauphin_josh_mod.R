battle_data_new <- survival_model_data_raw |>
  filter(stream == "battle creek") |>
  select(year, upstream_count, redd_count, gdd_total) |>
  rename(environmental_covar = gdd_total,
         spawner_count = redd_count) |>
  drop_na(environmental_covar) |>
  glimpse()

battle_data_list <- list("N" = length(unique(battle_data_new$year)),
                         "upstream_count" = battle_data_new$upstream_count,
                         "spawner_count" = battle_data_new$spawner_count,
                         "log_redds_per_spawner" = log(battle_data_new$spawner_count / battle_data_new$upstream_count),
                         "environmental_covar" = battle_data_new$environmental_covar,
                         "tau_obsv" = 5)

single_stream_dauphin_josh <- "
  data {
    int N;
    int upstream_count[N];
    int spawner_count[N]; // this is redd count but also holding count
    real log_redds_per_spawner[N]; // calculated, not sure how to generate
    real environmental_covar[N];
    real tau_obsv;
  }
  parameters {
    real <lower = 0> mean_redds_per_spawner;
    real <lower = 0> tau_redds_per_spawner;
    real b0_survival;
    real b1_survival;

  } model {
    // priors
    mean_redds_per_spawner ~ normal(-0.69, 0.1);
    tau_redds_per_spawner ~ gamma(0.01, 0.01);
    b0_survival ~ normal(0, 0.001);
    b1_survival ~ normal(0, 0.001);

    // transform parameter
    real sigma_redds_per_spawner = pow(tau_redds_per_spawner, -0.5);

    // calculated values
    vector[N] redds_per_spawner;
    vector[N] survival_rate;
    vector[N] predicted_redds;

    // calibration between upstream adults and redd count
    for(i in 1:N) {

      // redds per spawner for current year
      log_redds_per_spawner[i] ~ normal(mean_redds_per_spawner, tau_redds_per_spawner);
      redds_per_spawner[i] = exp(log_redds_per_spawner[i]);

      // predicted logit survival rate
      survival_rate[i] = inv_logit(-b0_survival + b1_survival * environmental_covar[i]);

      // predicted redds
      predicted_redds[i] = upstream_count[i] * survival_rate[i] * redds_per_spawner[i];

      // likelihood
      spawner_count[i] ~ normal(predicted_redds[i], tau_obsv);

    }

  }"


battle_dauphin_josh <- stan(model_code = single_stream_dauphin_josh,
                            data = battle_data_list,
                            chains = 4, iter = 5000*2, seed = 84735)


mcmc_trace(battle_dauphin_josh, pars = c("mean_redds_per_spawner", "tau_redds_per_spawner", "b0_survival", "b1_survival"))
mcmc_areas(battle_dauphin_josh, pars = c("mean_redds_per_spawner", "b0_survival", "b1_survival"))
mcmc_dens_overlay(battle_dauphin_josh, pars = c("mean_redds_per_spawner", "tau_redds_per_spawner", "b0_survival", "b1_survival")) # should be indistinguishable
neff_ratio(battle_dauphin_josh, pars = c("mean_redds_per_spawner", "tau_redds_per_spawner", "b0_survival", "b1_survival")) # should be >0.1
mcmc_acf(battle_dauphin_josh, pars = c("mean_redds_per_spawner", "tau_redds_per_spawner", "b0_survival", "b1_survival")) # should drop to be low
rhat(battle_dauphin_josh, c("mean_redds_per_spawner", "tau_redds_per_spawner", "b0_survival", "b1_survival")) # should be close to 1

model{

  #Data are: Nyrs, us_passage[iyr] (upstream passage estimate),
  # X[iyr] (covariate), tauObs (observation precision in redd counts)

  #Prior for mean redds per spawner. Initialize at log(0.5) (1 redd/2 spawners)
  muRsP~dnorm(-0.69,1.0E-03)	#

  #Prior for precision in redds/spawner (variation across years)
  tauRsP~dgamma(0.01,0.01)	#precision of RsP[iyr] (initialize at tauRsP=1)
  sdRsP<-power(tauRsp,-0.5)	#convert precition to standard deviation for output

  #priors for paramters predicting passage-redd survival rate in logit space
  b0_surv~dnorm(0,1.0E-03);b1_surv~dnorm(0,1.0E-03)

  for(iyr in 1:Nyrs){

    #Redds per Spawner for current year (note can be > 1)
    lgRsP[iyr]~dnorm(muRsp,tauRsP) #in log space
    RsP[iyr]<-exp(lgRsP[iyr])	#transformed so can't be negative

    #Predicted logit survival rate, convereted to linear space
    logit(surv[iyr])>-b0_surv + b1_surv*X[iyr]

    #predicted ress is product of # of upstream passage fish * survival rate * ReddsPerSpawner
    pred_redd[iyr]<-us_passage[iyr]*surv[iyr]*RsP[iyr]

    #likelihood of # of observed redds. tauObs is hardwired observation error in redd count (set to high #
    obs_redd[iyr]~dnorm(pred_redd[iyr],tauObs)

  }

}



# first attempt -----------------------------------------------------------

single_stream_dauphin_josh <- "
  data {
    int N;
    real b0; // output of environmental lm
    real b1; // output of environmental lm
    real b0_rps; // model this intercept? lm of existing data?
    int spawner_count[N];
    int upstream_count[N];
    real ratio_k[N];
    real environmental_covar[N];
  }
  parameters {
    real <lower = 0> survival_rate[N];
    real <lower = 0> redds_per_spawner_devs[N];
    real <lower = 0> sigma_rps;
  }
  model {
    // priors

    // calculated values
    vector[N] lambda;
    vector[N] redds_per_spawner;

    // calibration between upstream adults and redd count
    for(i in 1:N) {
      survival_rate[i] = b0 + b1 * environmental_covar[i];
      survival_rate[i] = logit(survival_rate[i]);

      redds_per_spawner_devs[i] ~ normal(0, sigma_rps);
      redds_per_spawner[i] = b0_rps + redds_per_spawner_devs[i];

      ratio_k[i] ~ redds_per_spawner[i] * survival_rate[i];

      lambda[i] = upstream_count[i] * ratio_k[i];
      spawner_count[i] ~ poisson(lambda[i]);
    }
  }"
