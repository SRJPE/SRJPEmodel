# get_survival_forecast <- function() {
#
# }

read_plad_csv <- function(site) {
  PLAD_results <- suppressWarnings(read_csv(paste0(here::here("data-raw", "PLAD", "PLAD_results"),
                                                   "/SRforklength_", site, ".csv"))) |>
    mutate(site = site) |>
    tidyr::separate_wider_delim("0.975", names = c("0.975", "extra"), delim = ",") |>
    mutate(across(`0.5`:extra, as.numeric))

  colnames(PLAD_results) <- c("week", "plad_0.5", "plad_0.025", "plad_0.25",
                              "plad_0.75", "plad_0.975", "site")

  return(PLAD_results)
}

# args
set.seed(SRJPEmodel::forecast_seed)
water_year_forecast <- 2 # TODO confirm this with Josh, will this change? this represents the element # to use for forecast year (e.g., iwy=1 for critical, iwy=2 for D/BN/AN, iwy=3 for W)
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)

#1) Get the median size of spring run outmigrants by model week
#2) Get the multi-year proportion of spring run outmigrants leaving by model week from inseason model
#3) Get survival from RST to Delta for each model week given fish size for that week (1)
#4) Calculate a weighted average survival across run, with 2) providing the weights.

# load CJS survival model fit object
# TODO if survival covariate is ever null, we will have to use an if/else statement (see DS_Surv.R)
survival_fit <- readRDS(here::here("data-raw", "survival_model", "survival_model_fit_wy3_fl.rds")) # TODO this will be updated to get_most_recent_model_object(trib, model type)
survival_posterior <- as.data.frame(survival_fit,
                                       pars = c("SurvForecastSz","TribSurvForecastSz")) |>
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") |>
  arrange(parameter) |>
  # extract indices from STAN par estimates so we can filter. These will be
  # size class and water year
  separate_wider_delim(parameter, delim = "[", names = c("parameter", "indices"),
                       too_few = "align_start") |>
  mutate(indices = str_remove_all(indices, "\\]")) |>
  separate_wider_delim(indices, ",", names = c("size_class", "wy_class", "reach"),
                       too_few = "align_start") |>
  mutate(across(size_class:reach, as.numeric))


n_sims <- dim(as.data.frame(survival_fit,
                            pars = c("SurvForecastSz","TribSurvForecastSz")))[1]

n_survival_locations <- 3 # Types of survival predictions from CJS model (mainstem, butte, feather)
trib_type_lookup <- c(1, 1, 1, 1, 2, 3) # ubc, ucc, mill, deer = 1, butte = 2, and eventually feather = 3

# Read in CJS forecasted survival posteriors
n_size_classes <- 25 # TODO confirm if this will change
size_covariate  <- tibble("size_class_bin" = seq(from = 10, to = 130, length.out = n_size_classes),
                          "size_class_index" = 1:n_size_classes)

survival_statistics <- survival_posterior |>
  mutate(location = ifelse(parameter == "SurvForecastSz", "mainstem", "tributary")) |>
  group_by(location, wy_class, size_class, reach) |>
  summarise(mean = mean(estimate),
            `10` = quantile(estimate, 0.1),
            `50` = quantile(estimate, 0.5),
            `90` = quantile(estimate, 0.9)) |>
  ungroup() |>
  left_join(size_covariate, by = c("size_class" = "size_class_index")) |>
  relocate(size_class_bin, .after = size_class)

# set up weights
Fl_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year
Surv_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year from weighted posterior sample
Surv_mu_check=vector(length=Ntribs)#pRun_weighted average survival calculated intuitive way for checking

#mean and sd of logit tranformed posterior samples of survival and the pRun-weighted mean forklength for JPE.stan model
DS_surv_mu=vector(length=Ntribs);
DS_surv_sd=DS_surv_mu

# extract PLAD results
plad_fl_results <- purrr::pmap(list(site = SRJPEmodel::forecast_sites),
                            read_plad_csv) |>
  reduce(bind_rows)

# TODO update this function call so you can filter by site, model name, etc. in call
cumulative_proportions <- SRJPEmodel::get_most_recent_model_results(con) |>
  filter(model_name == "inseason",
         parameter == "For_cp",
         statistic %in% c("mean")) |>
  select(site, week_fit, value, statistic) |>
  mutate(statistic = paste0("cp_", statistic)) |>
  pivot_wider(names_from = "statistic", values_from = "value") |>
  mutate(sort = ifelse(week_fit >= 36, 1, 2)) |>
  arrange(site, sort, week_fit) |>
  select(-sort)

# link forecasted cumulative proportion by week to PLAD results by week
fl_cp_results <- cumulative_proportions |>
  left_join(plad_fl_results, by = c("week_fit" = "week", "site"))

  # 1) Get across-year mean size of spring run by week from PLAD  #########

  # 2) Get proportion of run passing trap for every week in year (Sep-Aug) ########

  cp_proportions <- cumulative_proportions |>
    right_join(plad_fl_results, by = c("week_fit" = "week", "site")) |>
    select(-c(plad_0.5:plad_0.975)) |>
    group_by(site) |>
    mutate(weekly_added = cp_mean - lag(cp_mean),
           weekly_added = ifelse(is.na(weekly_added), cp_mean, weekly_added),
           total = sum(weekly_added),
           weekly_added_std = weekly_added / total)

  #3) Get pRun-weighted size and survival rate for outmigrant run ################

  #Set the jtrib index needed for SurvForecastSz. TribSurvForecastSz, and survival_statistics
  if(DoTrib[itrib]=="ubc"| DoTrib[itrib]=="ucc"| DoTrib[itrib]=="mill creek"| DoTrib[itrib]=="deer creek"){
    jtrib=1
  } else if (DoTrib[itrib]=="okie dam"){
    jtrib=2
  } else {#Feather eventually
    jtrib=3
  }

  Surv_mu_check[itrib]=0;
  Fl_mu[itrib=0]
  #Set of samples to grab from posterior sample of CJS survival for size class associated with following weeks
  if(itrib==1) isamps=sample(1:n_sims,size=500,replace=F) #to keep comp_post to reasonable size, sample 100 posterior survival values from each size class (wk)
  for(iwk in 1:Nwks){

    #find the closest size in size_covariate given median size in this week from plad (Sz) to get the size_class index for size_covariate used in CJS model
    irecs=which(size_covariate<=Sz[itrib,iwk])
    if(length(irecs)==0){
      size_class=1 #Sz is smaller than lowest value of size_covariate so set to lowest index
    } else{
      size_class=max(irecs) #For Sz>max(size_covariate) this will set index to last element for size_covariate
    }

    #Create a composite posterior of survival rates which accumulates posterior samples across all weeks/size classes.
    #Posterior sample each wk/size class is weighted based on pRun relative to lowest pRun aross weeks (where reps=1)
    preps=round(pRun[itrib,iwk]/min(pRun[itrib,])) # # of replicates for this iwk
    if(iwk==1){
      comp_post=rep(survival_posterior[isamps,size_class,jtrib],times=preps) #repeat the random sample of posterior for this wk preps times
    } else {
      comp_post=c(comp_post,rep(survival_posterior[isamps,size_class,jtrib],times=preps))#each set of repeated samples is added to the vectory
    }

    Fl_mu[itrib]=Fl_mu[itrib]+pRun[itrib,iwk]*Sz[itrib,iwk]#pRun-weighted mean size

    #Most obvious way to compute a pRun-weighted average survival.
    Surv_mu_check[itrib]=Surv_mu_check[itrib]+pRun[itrib,iwk]*survival_statistics[jtrib,size_class,2]
  }
  Surv_mu[itrib]=mean(comp_post)#To compare against Surv_mu_check

  #mean and sd of survival rate in logit space for outmigrant run for JPE.stan model
  DS_surv_mu[itrib]=mean(logit(comp_post))
  DS_surv_sd[itrib]=sd(logit(comp_post))

  print(c(DoTrib[itrib],round(Surv_mu_check[itrib],digits=4),round(Surv_mu[itrib],digits=4)))

  #Simulated survival rates as done in JPE.stan (but with lower constraint of -6.5 ~ 0.15% survival)
  #sim_surv used to see if composite posterior (comp_post) is accurately modelled by DS_surv_mu and DS_surv_sd and JPE.stan approap
  #For plotting (if PlotType==2) or text output on different plot type (if PlotType==1)
  sim_surv=inv_logit(rnorm(n=5000,mean=DS_surv_mu[itrib],sd=DS_surv_sd[itrib]))
