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

extract_survival_posteriors <- function(con) {
  model_fits <- get_most_recent_model_objects(con, model_component = "model_fit",
                                              model_name = "survival")
  name_lookup <- lapply(names(model_fits), function(x) {
      str_split(x, "_")[[1]][4]
    }) |>
    unlist()
  names(model_fits) <- name_lookup

  posteriors <- lapply(model_fits, function(x) {
    df <- as.data.frame(x,
                        pars = c("SurvForecastSz","TribSurvForecastSz")) |>
      mutate(stream = stream_name)
    return(df)
  })

  # TODO bind rows and mutate the site name into a column
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

# load fit objects
# TODO if survival covariate is ever null, we will have to use an if/else statement (see DS_Surv.R)

# survival
survival_fit <- readRDS(here::here("data-raw", "survival_model", "survival_model_fit_wy3_fl.rds")) # TODO this will be updated to get_most_recent_model_object(trib, model type)
# survival_fit <- get_most_recent_model_objects(con, model_component = "model_fit",
#                                               model_name = "survival")
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

# indexing
n_sims <- dim(as.data.frame(survival_fit,
                            pars = c("SurvForecastSz","TribSurvForecastSz")))[1]

n_survival_locations <- 3 # Types of survival predictions from CJS model (mainstem, butte, feather)
trib_type_lookup <- c(1, 1, 1, 1, 2, 3) # ubc, ucc, mill, deer = 1, butte = 2, and eventually feather = 3

# Read in CJS forecasted survival posteriors
n_size_classes <- 25 # TODO confirm if this will change
size_covariate  <- tibble("size_class_bin" = seq(from = 10, to = 130, length.out = n_size_classes),
                          "size_class_index" = 1:n_size_classes)

# inseason model
cp_fits <- get_most_recent_model_objects(con, model_component = "model_fit",
                                         model_name = "inseason")
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

# 1) Get across-year mean size of spring run by week from PLAD  #########
# PLAD fork length results
plad_fl_results <- purrr::pmap(list(site = SRJPEmodel::forecast_sites),
                               read_plad_csv) |>
  reduce(bind_rows) |>
  select(site, week_fit = week, plad_0.5) |>
  mutate(nearest_fl_bin = round(plad_0.5 / 5) * 5) |>
  # link fork lengths in PLAD to nearest size class bin (which is currently at level of 5)
  left_join(size_covariate, by = c("nearest_fl_bin" = "size_class_bin")) |>
  # catch sizes that are below or above the size class bins in the covariate
  mutate(size_class_index = case_when(is.na(size_class_index) & plad_0.5 < min(size_covariate$size_class_bin) ~ min(size_covariate$size_class_index),
                                      is.na(size_class_index) & plad_0.5 > max(size_covariate$size_class_bin) ~ max(size_covariate$size_class_index),
                                      TRUE ~ size_class_index))

# 2) Get proportion of run passing trap for every week in year (Sep-Aug) ########
cp_proportions <- cumulative_proportions |>
  right_join(plad_fl_results, by = c("week_fit", "site")) |>
  select(-plad_0.5) |>
  group_by(site) |>
  mutate(weekly_added = cp_mean - lag(cp_mean),
         weekly_added = ifelse(is.na(weekly_added), cp_mean, weekly_added),
         total = sum(weekly_added),
         weekly_added_std = weekly_added / total) |>
  ungroup() |>
  select(site, week_fit, cp_mean, weekly_added_std)

#3) Get pRun-weighted size and survival rate for outmigrant run ################

# Create a composite posterior of survival rates which accumulates
# posterior samples across all weeks/size classes.
# Posterior sample each wk/size class is weighted based on pRun
# relative to lowest pRun aross weeks (where reps=1)

# Set of samples to grab from posterior sample of CJS survival for size class associated with following weeks
# this is only for ubc, ucc, mill creek, and deer creek
trial_samples <- sample(1:n_sims, size = 500, replace = F) # to keep comp_post to reasonable size, sample 100 posterior survival values from each size class (wk)

weighted_fl <- cp_proportions |>
  group_by(site) |>
  mutate(weekly_replicates = round(weekly_added_std / min(weekly_added_std))) |>
  ungroup()

# by week

test <- survival_posterior |>
  filter(size_class == 1,
         is.na(reach))

for(i in 1:31) {
  if(i == 1) {
    replicates <- rep(test$estimate[trial_samples],
                      times = weighted_fl$weekly_replicates[1])
  } else {
    replicates <- c(replicates,
                    rep(test$estimate[trial_samples],
                        times = weighted_fl$weekly_replicates[1]))
  }
}

weighted_fl_results <- cp_proportions |>
  # link mean plad by site and week
  left_join(plad_fl_results, by = c("site", "week_fit")) |>
  mutate(weekly_fl_mean_weighted = weekly_added_std * plad_0.5) |> # pRun weighted mean size
  group_by(site) |>
  mutate(weekly_fl_mean_weighted_cumulative = cumsum(weekly_fl_mean_weighted)) |>
  ungroup() |> glimpse()

weighted_mean_survival <- cp_proportions |>
  mutate(reach_lookup = case_when(site == "okie dam" ~ 1,
                                  site == "herringer riffle" ~ 2,
                                  TRUE ~ NA)) |>
  left_join(survival_statistics |>
              select(size_class, reach, mean_survival = mean),
            by = c("reach" = "reach_lookup",))
  # now link size classes by reach from survival model
  left_join(survival_statistics, by = c)

# TODO take mx of weekly_fl_men_weighted_cumulative by site for FL_mu

    # left_join(survival_statistics, by = c("reach_lookup" = "reach")) |> glimpse()
    # mutate(reach_lookup = case_when(site == "okie dam" ~ 1,
    #                                 site == "herringer riffle" ~ 2, # TODO update for feather site
    #                                 TRUE ~ NA)) |>
         weekly_mean_survival_weighted_check = weekly_added_std * )


    #Most obvious way to compute a pRun-weighted average survival.
    Surv_mu_check[itrib]=Surv_mu_check[itrib]+pRun[itrib,iwk]*survival_statistics[jtrib,size_class,2]
  }
  Surv_mu[itrib]=mean(comp_post)#To compare against Surv_mu_check

  #mean and sd of survival rate in logit space for outmigrant run for JPE.stan model
  DS_surv_mu[itrib]=mean(logit(comp_post))
  DS_surv_sd[itrib]=sd(logit(comp_post))

  print(c(DoTrib[itrib],round(Surv_mu_check[itrib],digits=4),round(Surv_mu[itrib],digits=4)))


  # set up weights
  Fl_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year
  Surv_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year from weighted posterior sample
  Surv_mu_check=vector(length=Ntribs)#pRun_weighted average survival calculated intuitive way for checking

  #mean and sd of logit tranformed posterior samples of survival and the pRun-weighted mean forklength for JPE.stan model
  DS_surv_mu=vector(length=Ntribs);
  DS_surv_sd=DS_surv_mu

  # link forecasted cumulative proportion by week to PLAD results by week
  fl_cp_results <- cumulative_proportions |>
    left_join(plad_fl_results, by = c("week_fit", "site"))
