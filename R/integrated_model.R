#' Read in PLAD csv results
#' @description helper function to read in PLAD results as currently provided in a .csv format.
#' @keywords internal
#' @export
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

#' Sample inseason model posterior
#' @description Uses weights to set replicates and then samples based on those
#' replicates for use in calculating / comparing weighted survival
#' @keywords internal
#' @export
sample_composite_posterior <- function(size_class_arg, reach_arg,
                                       week_fit_arg, weekly_replicates_arg) {
  trial_samples <- sample(1:n_sims, size = 500, replace = F) # to keep comp_post to reasonable size, sample 100 posterior survival values from each size class (wk)
  survival_posterior_samples <- survival_posterior |>
    filter(size_class == !!size_class_arg,
           reach == !!reach_arg) |>
    mutate(row_number = row_number()) |>
    filter(row_number %in% trial_samples) |>
    pull(estimate)

  if(week_fit_arg == 45) {
    sample <- rep(survival_posterior_samples, times = weekly_replicates_arg)
  } else {
    sample <- rep(survival_posterior_samples, times = weekly_replicates_arg)
  }
  return(sample)
}


#' Get survival forecasts weighted by FL and run size
#'
#' @description
#' Pulls together estimates from the inseason, PLAD, and survival models to generate a
#' forecasted survival rate from RST to Delta that is weighted by fork length and run size.
#'
#' @details
#' TODO
#'
#' @param con A valid DBI database connection to the SRJPE model run database.

#' @return A tibble containing:
#' \describe{
#'   \item{`site`}{The site for which the survival rate is produced.}
#'   \item{`mean_weighted_survival`}{Mean of forecasted survival rate in logit space for outmigrant run for JPE.stan model.}
#'   \item{`sd_weighted_survival`}{Standard deviation of forecasted survival rate in logit space for outmigrant run for JPE.stan model.}
#' }
#'
#' @family Forecast
#' @export
#' @md
get_survival_forecast <- function(con) {

  # set up
  set.seed(SRJPEmodel::forecast_seed)
  # water_year_forecast <- 2 # TODO confirm this with Josh, will this change? this represents the element # to use for forecast year (e.g., iwy=1 for critical, iwy=2 for D/BN/AN, iwy=3 for W)

  # covariates
  n_size_classes <- 25 # TODO confirm if this will change
  size_covariate  <- tibble("size_class_bin" = seq(from = 10, to = 130, length.out = n_size_classes),
                            "size_class_index" = 1:n_size_classes)

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
    mutate(across(size_class:reach, as.numeric),
           reach = ifelse(is.na(reach), 0, reach))

  survival_statistics <- survival_posterior |>
    mutate(location = ifelse(parameter == "SurvForecastSz", "mainstem", "tributary")) |>
    group_by(location, size_class, reach) |> # TODO grouping by water year types here - confirm with josh what we do if we have 3 WYs ?
    summarise(mean = mean(estimate),
              `10` = quantile(estimate, 0.1),
              `50` = quantile(estimate, 0.5),
              `90` = quantile(estimate, 0.9)) |>
    ungroup() |>
    left_join(size_covariate, by = c("size_class" = "size_class_index")) |>
    relocate(size_class_bin, .after = size_class) |>
    # link mainstem to reach index = 0 for clarity
    mutate(reach = ifelse(is.na(reach), 0, reach))

  # indexing
  n_sims <- dim(as.data.frame(survival_fit,
                              pars = c("SurvForecastSz","TribSurvForecastSz")))[1]

  n_survival_locations <- 3 # Types of survival predictions from CJS model (mainstem, butte, feather)

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
                                 SRJPEmodel::read_plad_csv) |>
    reduce(bind_rows) |>
    select(site, week_fit = week, plad_0.5) |>
    mutate(nearest_fl_bin = round(plad_0.5 / 5) * 5) |>
    # link fork lengths in PLAD to nearest size class bin (which is currently at level of 5)
    left_join(size_covariate, by = c("nearest_fl_bin" = "size_class_bin")) |>
    # catch sizes that are below or above the size class bins in the covariate
    mutate(size_class_index = case_when(is.na(size_class_index) & plad_0.5 < min(size_covariate$size_class_bin) ~ min(size_covariate$size_class_index),
                                        is.na(size_class_index) & plad_0.5 > max(size_covariate$size_class_bin) ~ max(size_covariate$size_class_index),
                                        TRUE ~ size_class_index))

  #2) Get the multi-year proportion of spring run outmigrants leaving by model week from inseason model
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

  #3) Get survival from RST to Delta for each model week given fish size for that week (pRun weighted)

  # Create a composite posterior of survival rates which accumulates
  # posterior samples across all weeks/size classes.
  # Posterior sample each wk/size class is weighted based on pRun
  # relative to lowest pRun aross weeks (where reps=1)

  weighted_fl_results <- cp_proportions |>
    # link mean plad by site and week
    left_join(plad_fl_results, by = c("site", "week_fit")) |>
    mutate(weekly_fl_mean_weighted = weekly_added_std * plad_0.5) |> # pRun weighted mean size
    group_by(site) |>
    mutate(weekly_fl_mean_weighted_cumulative = cumsum(weekly_fl_mean_weighted)) |>
    ungroup()

  weighted_mean_survival_results <- cp_proportions |>
    mutate(reach_lookup = case_when(site == "okie dam" ~ 1,
                                    site == "herringer riffle" ~ 2,
                                    TRUE ~ 0)) |>
    # link mean plad by site and week
    left_join(plad_fl_results |>
                select(site, week_fit, size_class_index),
              by = c("site", "week_fit")) |>
    # now link size classes by reach from survival model
    left_join(survival_statistics |>
                select(size_class, reach, mean_survival = mean),
              by = c("reach_lookup" = "reach", "size_class_index" = "size_class")) |>
    mutate(weighted_survival_weekly = weekly_added_std * mean_survival) |>
    group_by(site) |>
    mutate(weighted_survival_by_site = cumsum(weighted_survival_weekly)) |>
    ungroup()

  # Set of samples to grab from posterior sample of CJS survival for size class associated with following weeks
  # this is only for ubc, ucc, mill creek, and deer creek
  replicates <- cp_proportions |>
    group_by(site) |>
    mutate(weekly_replicates = round(weekly_added_std / min(weekly_added_std))) |>
    ungroup() |>
    mutate(reach = case_when(site == "herringer riffle" ~ 2,
                             site == "okie dam" ~ 1,
                             TRUE ~ 0)) |>
    left_join(plad_fl_results |>
                select(site, week_fit, size_class_index),
              by = c("site", "week_fit")) |>
    select(site, reach, week_fit, size_class_index, weekly_replicates)

  #4) Calculate a weighted average survival across run, with 2) providing the weights.
  mean_weighted_survival <- list()
  mean_weighted_survival_logit <- list()
  sd_weighted_survival_logit <- list()
  for(site_iter in c("lcc", "ubc")) { # unique(replicates$site)) {
    new_replicates <- replicates |>
      filter(site == site_iter)
    # produces samples by week
    weighted_survival_weekly <- purrr::pmap(list(size_class_arg = new_replicates$size_class_index,
                                                 reach_arg = new_replicates$reach,
                                                 week_fit_arg = new_replicates$week_fit,
                                                 weekly_replicates_arg = new_replicates$weekly_replicates),
                                                 SRJPEmodel::sample_composite_posterior) |>
      # concatenate all into one vector
      unlist()
    index <- which(unique(replicates$site) == site_iter)
    mean_weighted_survival[[index]] <- mean(weighted_survival_weekly)
    mean_weighted_survival_logit[[index]] <- mean(qlogis(weighted_survival_weekly))
    sd_weighted_survival_logit[[index]] <- sd(qlogis(weighted_survival_weekly))

    names(mean_weighted_survival)[index] <- site_iter
    names(mean_weighted_survival_logit)[index] <- site_iter
    names(sd_weighted_survival_logit)[index] <- site_iter
  }

  # then compare mean_weighted_survival against weighted_mean_survival_results$weighted_survival_by_site

  # mean and sd of survival rate in logit space for outmigrant run for JPE.stan model
  weighted_survival_results <- tibble("site" = c("lcc", "ubc"), # unique(replicates$site)
                                      "mean_weighted_survival" = mean_weighted_survival_logit,
                                      "sd_weighted_survival" = sd_weighted_survival_logit)


  return(weighted_survival_results)
}



# generate_inseason_forecast <- function(con) {
#
# }

# predIn_stats=array(dim=c(Ntribs,Nfor,3))
# cp_mu=matrix(nrow=Ntribs,ncol=Nfor)
# cp_sd=cp_mu
# obs_Nx_mu=cp_mu;
# obs_Nx_cv=cp_mu;
# obs_Nx_sd=cp_mu
#
# for(itrib in 1:Ntribs){
#   for(ifor in 1:Nfor){
#
#     For_ewk=which(calwk==For_CalWk[ifor])
#
#     #Currently setup to use null model based on proportion of run forecasts caclculated in stan model.
#     #Will eventually replace this with code to use covariate model specified in For_Input.csv
#     #This will require reading in parameters and calculating For_cp within this code
#     fnm=paste0(InseasDir,DoTrib[itrib],"_null.Rdata")
#     load(file=fnm)
#
#     dis=as.data.frame(fit, pars = c("For_cp"))
#     Nsims=dim(dis)[1]
#
#     cp_all=dis[,For_ewk];grecs=which(is.finite(cp_all)==T & cp_all>0 & cp_all<1)
#     Ngsims=length(grecs)
#     cp=cp_all[grecs]
#
#     cp_mu[itrib,ifor]=mean(logit(cp));cp_sd[itrib,ifor]=sd(logit(cp))
#
#     obs_Nx_mu[itrib,ifor]=log(For_Nx_mu[itrib,ifor])#total abundance through this week in log space
#     obs_Nx_cv[itrib,ifor]=For_Nx_sd[itrib,ifor]/For_Nx_mu[itrib,ifor]#the cv in untranformed space
#     obs_Nx_sd[itrib,ifor]=sqrt(log(obs_Nx_cv[itrib,ifor]^2+1))#convert from cv in untransformed space to sd in log space
#
#     #Uncertainty in cummulative abundance on forecast week in log space
#     pNx=rnorm(n=Ngsims,mean=obs_Nx_mu[itrib,ifor],sd=obs_Nx_sd[itrib,ifor])
#
#     #estimate of annual abundance (cum abundance on forecast wk expaned by proportion of run by that wk)
#     predIn=exp(pNx)/cp
#     #predIn=mean(exp(pNx))/cp #no uncertainty
#
#     predIn_stats[itrib,ifor,]=as.double(quantile(predIn,probs=CI))#;predIn_stats[itrib,ifor,2]=mean(predIn)
#     print(c(DoTrib[itrib],For_ewk,round(predIn_stats[itrib,ifor,],digits=0)))
#   }#next ifor
# }
