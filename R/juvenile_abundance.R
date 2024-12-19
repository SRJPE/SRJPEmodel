#' Fit pCap model
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param input A list containing the inputs for the pCap model.
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
fit_pCap_model <- function(input) {

  stan_model <- eval(parse(text = "SRJPEmodel::bt_spas_x_model_code$pCap_all"))

  options(mc.cores=parallel::detectCores())
  cli::cli_alert("running pCap model")
  pcap <- rstan::stan(model_code = stan_model,
                      data = input$data,
                      init = input$inits,
                      # do not save logit_pCap or pro_dev_P (way too big)
                      pars = c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "trib_mu_P", "trib_sd_P",
                               "flow_mu_P", "flow_sd_P"), # for new efficient model
                      chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                      iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                      seed = 84735)

  return(pcap)
}

#' Generate lt_pCap_U values from the pCap model object
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param pCap_inputs A list containing the inputs for the pCap model including `model_name`, `Nstrata`,
#' `catch_flow`, `use_trib`, `Ntribs`, `Nmr`, `Nwmr`, `Nwomr`, `Uind_wMR`, `Uind_woMR`,
#' `ind_pCap`.
#' @param pCap_model_object A STANfit object with the output of the pCap hierarchical model.
#' @returns a named list with the required elements for that model run.
#' @export
#' @md
generate_lt_pCap_Us <- function(abundance_inputs, pCap_model_object){

  # set up objects
  ModelName <- abundance_inputs$model_name
  Nstrata <- abundance_inputs$inputs$data$Nstrata
  catch_flow <- abundance_inputs$lt_pCap_U_data$catch_flow
  use_trib <- abundance_inputs$lt_pCap_U_data$use_trib
  Nwmr <- abundance_inputs$lt_pCap_U_data$Nwmr
  Nwomr <- abundance_inputs$lt_pCap_U_data$Nwomr
  Uind_wMR <- abundance_inputs$lt_pCap_U_data$Uind_wMR
  Uind_woMR <- abundance_inputs$lt_pCap_U_data$Uind_woMR
  ind_pCap <- abundance_inputs$lt_pCap_U_data$ind_pCap

  # extract from pCap model fit object

  samples <- rstan::extract(pCap_model_object, pars = c("logit_pCap", "b0_pCap", "b_flow", "pro_sd_P", "trib_mu_P", "trib_sd_P",
                                           "flow_mu_P", "flow_sd_P"),
                            permuted = TRUE)
  Ntrials <- dim(samples$logit_pCap)[1] # of saved posterior samples from pCap model in stan
  logit_pCap <- samples$logit_pCap # logit_pCap[1:Ntrials,1:Nmr] # The estimated logit pCap posterior for each efficiency trial
  b0_pCap <- samples$b0_pCap # b0_pCap[1:Ntrials,1:Ntribs] #mean logit pCap for each site (at mean discharge)
  b_flow <- samples$b_flow # b_flow[1:Ntrials,1:Ntribs] #flow effect for each site
  pro_sd_P <- samples$pro_sd_P # pro_sd_P[1:Ntrials]        #process error (sd)
  trib_mu_P <- samples$trib_mu_P # trib_mu_P[1:Ntrials] #hyper mean for b0_pCap
  trib_sd_P <- samples$trib_sd_P # trib_sd_P[1:Ntrials] #hyper sd for b0_pCap
  flow_mu_P <- samples$flow_mu_P # flow_mu_P[1:Ntrials] #hyper mean for b_flow
  flow_sd_P <- samples$flow_sd_P # flow_sd_P[1:Ntrials] #hyper sd for b_flow

  # calculations
  lt_pCap_U=matrix(nrow=Ntrials,ncol=Nstrata)
  sim_pro_dev=vector(length=Ntrials)
  lt_pCap_mu=matrix(nrow=Nstrata,ncol=Ntrials) #function needs to return this
  lt_pCap_sd=lt_pCap_mu                        #function needs to return this

  if(ModelName=="all_mark_recap" ){

    for(i in 1:Nstrata){
      lt_pCap_U[,i] = logit_pCap[,ind_pCap[i]];
    }

  } else if (ModelName=="missing_mark_recap"){

    for(i in 1:Nwmr){
      #Assign estimated pCaps for strata with efficiency data
      lt_pCap_U[,Uind_wMR[i]] = logit_pCap[,ind_pCap[i]];
    }
    for (i in 1:Nwomr) {
      #for weeks without efficiency trials
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0,sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials]
    }

  } else if (ModelName=="no_mark_recap"){

    for (i in 1:Nwomr) {
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0, sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = b0_pCap[,use_trib] + b_flow[,use_trib] * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials]
    }

  } else if (ModelName=="no_mark_recap_no_trib"){

    logit_b0=vector(length=Ntrials); logit_bflow=logit_b0
    for(itrial in 1:Ntrials){
      logit_b0[itrial] = rnorm(n=1, mean=trib_mu_P[itrial], sd=trib_sd_P[itrial])
      logit_bflow[itrial] = rnorm(n=1, mean=flow_mu_P[itrial], sd=flow_sd_P[itrial])
    }

    for (i in 1:Nwomr) {
      for(itrial in 1:Ntrials) sim_pro_dev[itrial] = rnorm(n=1, mean=0, sd=pro_sd_P[itrial]);
      lt_pCap_U[,Uind_woMR[i]] = logit_b0 + logit_bflow * catch_flow[Uind_woMR[i]] + sim_pro_dev[1:Ntrials];
    }

  }#end if on ModelName


  #Calculate mean and sd for each lt_pCap_U and return these from function
  for(i in 1:Nstrata){
    lt_pCap_mu[i,]=mean(lt_pCap_U[,i])
    lt_pCap_sd[i,]=sd(lt_pCap_U[,i])
    # lt_pCap_mu[,i]=mean(lt_pCap_U[i])
    # lt_pCap_sd[,i]=sd(lt_pCap_U[i])
  }

  return(list("lt_pCap_mu" = lt_pCap_mu |> rowMeans(),
              "lt_pCap_sd" = lt_pCap_sd |> rowMeans()))
}


#' Run BUGS abundance model.
#' @details This function runs the BUGS abundance model.
#' @param abundance_inputs the object produced by `prepare_abundance_inputs()`
#' @param lt_pCap_Us the object produced by `generate_lt_pCap_Us()`
#' @param bugs_model_file the filepath pointing to where your BUGS abundance model file is
#' @param bugs_directory the filepath pointing to where your WinBUGS14/ directory is
#' @returns a BUGS object.
#' @export
#' @md
fit_abundance_model_BUGS <- function(abundance_inputs, lt_pCap_Us,
                                     bugs_model_file,
                                     bugs_directory) {

  parameters <- c("lt_pCap_U", "pCap_U", "N", "Ntot", "sd.N", "sd.Ne", "lg_CumN")
  Nmcmc = 2000
  Nburnin = 500
  Nthin = 2
  Nchains = 3

  # clean up later
  data <- abundance_inputs$inputs$data
  data$lt_pCap_mu <- lt_pCap_Us$lt_pCap_mu
  data$lt_pCap_tau <- 1/lt_pCap_Us$lt_pCap_sd^2

  inits_with_lt_pCap_U <- abundance_inputs$inputs$inits[[1]]

  inits_with_lt_pCap_U$lt_pCap_U <- data$lt_pCap_mu
  inits_with_lt_pCap_U$tau.N <- 1
  inits_with_lt_pCap_U$tau.Ne <- 1

  new_inits <- list(inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U,
                    inits_with_lt_pCap_U)

  abundance <- bugs(data,
                    new_inits,
                    parameters,
                    bugs_model_file, debug = F,
                    n.chains = Nchains, n.burnin = Nburnin, n.thin = Nthin, n.iter = Nmcmc,
                    codaPkg = F, DIC = T, clearWD = T,
                    bugs.directory = bugs_directory)
                    #bugs.directory = here::here("data-raw", "WinBUGS14"))

  return(abundance)
}



#' Fit abundance model in STAN
#' @details This function prepares the data list for the abundance STAN model based on what the model name is.
#' @param input A list containing the inputs for the abundance model.
#' @param pCap_fit A STANfit object resulting from running `fit_pCap_model()`
#' @returns a named list with the required elements for that model run.
#' @keywords internal
#' @export
#' @md
fit_abundance_model_STAN <- function(input, pCap_fit, lt_pCap_Us) {

  # TODO this is in progress

  input$data$lt_pCap_mu <- lt_pCap_Us$lt_pCap_mu
  input$data$lt_pCap_sd <- lt_pCap_Us$lt_pCap_sd

  # generate inits for lt_pCap_U
  inits_with_lt_pCap_U <- input$inits[[1]]
  inits_with_lt_pCap_U$lt_pCap_U <- input$data$lt_pCap_mu
  new_inits <- list(inits = inits_with_lt_pCap_U,
                    inits = inits_with_lt_pCap_U,
                    inits = inits_with_lt_pCap_U)

  stan_model <- eval(parse(text = "SRJPEmodel::bt_spas_x_model_code$abundance"))

  options(mc.cores=parallel::detectCores())

  cli::cli_alert("running abundance model")
  abundance <- rstan::stan(model_code = stan_model,
                           data = input$data,
                           # algorithm = "HMC",
                           init = new_inits,
                           chains = SRJPEmodel::bt_spas_x_bayes_params$number_chains,
                           iter = SRJPEmodel::bt_spas_x_bayes_params$number_mcmc,
                           seed = 84735)

  return(abundance)
}

#' BT SPAS X diagnostic plots
#' @details This function produces a plot with data and results of fitting the pCap and abundance models.
#' @param site_arg The site being fit
#' @param run_year_arg The run year being fit
#' @param abundance_model The STANfit model object produced by running `fit_abundance_model()`
#' @returns A plot
#' @export
#' @md
generate_diagnostic_plot_juv <- function(site_arg, run_year_arg,
                                         abundance_model) {


  julian_week_to_date_lookup <- SRJPEmodel::julian_week_to_date_lookup

  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(life_stage != "yearling") |>
    filter(run_year == run_year_arg,
           site == site_arg,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = if(all(is.na(catch_standardized_by_hours_fished))) NA_real_ else sum(catch_standardized_by_hours_fished, na.rm = TRUE),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # change all NaNs to NAs
           across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"),
           week_index = row_number())

  pCap_estimates <- rstan::summary(abundance_model,pars=c("lt_pCap_U"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "lt_pCap_U")

  N_estimates <- rstan::summary(abundance_model,pars=c("N"))$summary |>
    data.frame() |>
    tibble::rownames_to_column("parameter") |>
    mutate(week_index = readr::parse_number(parameter)) |>
    select(week_index, mean, sd, `50` = X50., `2.5` = X2.5., `97.5` = X97.5.) |>
    mutate(parameter = "N")


  data |>
    ggplot(aes(x = final_date, y = catch_standardized_by_hours_fished)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    #scale_x_date(date_breaks = "1 week", date_labels = "%V") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    labs(x = "Date",
         y = "Raw catch",
         title = paste0("Catch at ", site_arg, " for run year ",
                        run_year_arg))

  # abundance plot
  abundance_plot <- data |>
    left_join(N_estimates) |>
    #mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_abundance),
               shape = 1, color = "blue") +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    scale_color_manual(values = c("TRUE" = "white", "FALSE" = "red")) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(site_arg, run_year_arg)) +
    #theme(axis.text.x=element_blank()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # efficiency
  efficiency_plot <- data |>
    left_join(pCap_estimates) |>
    mutate(across(c(mean, `50`, `2.5`, `97.5`), plogis),
           number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = final_date, y = `50`)) +
    geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_errorbar(aes(x = final_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
    geom_point(aes(x = final_date, y = lincoln_peterson_efficiency),
               shape = 1, color = "blue") +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    # geom_point(aes(x = fake_date, y = Inf, color = sampled),
    #            size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)
}


#' BT SPAS X raw data plots
#' @details This function produces a plot with data used to fit BT SPAS X.
#' @param site The site being fit
#' @param run_year The run year being fit
#' @returns A plot
#' @export
#' @md
plot_juv_data <- function(site, run_year) {
  data <- SRJPEdata::weekly_juvenile_abundance_catch_data |>
    filter(life_stage != "yearling") |>
    filter(run_year == !!run_year,
           site == !!site,
           week %in% c(seq(45, 53), seq(1, 22))) |>
    group_by(year, week, stream, site, run_year) |>
    # keep NAs in count columns
    summarise(count = if(all(is.na(count))) NA_real_ else sum(count, na.rm = TRUE),
              mean_fork_length = mean(mean_fork_length, na.rm = T),
              hours_fished = mean(hours_fished, na.rm = T),
              catch_flow_cfs = mean(flow_cfs, na.rm = T),
              average_stream_hours_fished = mean(average_stream_hours_fished, na.rm = T),
              standardized_flow = mean(standardized_flow, na.rm = T),
              catch_standardized_by_hours_fished = if(all(is.na(catch_standardized_by_hours_fished))) NA_real_ else sum(catch_standardized_by_hours_fished, na.rm = TRUE),
              lgN_prior = mean(lgN_prior, na.rm = T)) |>
    ungroup() |>
    left_join(SRJPEdata::weekly_juvenile_abundance_efficiency_data,
              by = c("year", "run_year", "week", "stream", "site")) |>
    mutate(count = round(count, 0),
           catch_standardized_by_hours_fished = round(catch_standardized_by_hours_fished, 0),
           # change all NaNs to NAs
           across(mean_fork_length:lgN_prior, ~ifelse(is.nan(.x), NA, .x)),
           # plot things
           lincoln_peterson_abundance = count * (number_released / number_recaptured),
           lincoln_peterson_efficiency = number_recaptured / number_released) |>
    left_join(SRJPEmodel::julian_week_to_date_lookup, by = c("week" = "Jwk")) |>
    mutate(year = ifelse(week >= 43, run_year - 1, run_year),
           fake_date = ymd(paste0(year, "-01-01")),
           final_date = fake_date + weeks(week - 1),
           date = format(final_date, "%b-%d"),
           week_index = row_number())

  abundance_plot <- data |>
    ggplot(aes(x = final_date, y = count)) +
    geom_bar(stat = "identity", fill = "grey", width = 5) +
    # geom_point(aes(x = final_date, y = lincoln_peterson_abundance),
    #            shape = 1, color = "blue") +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(count),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "",
         #x = "Date",
         y = "Abundance",
         title = paste(site, run_year)) +
    #theme(axis.text.x=element_blank()) +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  efficiency_plot <- data |>
    mutate(number_released_label = ifelse(is.na(number_released), "", number_released),
           number_recaptured_label = ifelse(is.na(number_recaptured), "", number_recaptured)) |>
    ggplot(aes(x = final_date, y = lincoln_peterson_efficiency)) +
    geom_point(shape = 1, color = "blue") +
    # geom_bar(stat = "identity", fill = "grey", width = 4) +
    geom_text(aes(x = final_date, y = Inf,
                  label = paste(number_released_label, number_recaptured_label),
                  angle = 90),
              hjust = 1,
              size = 3) +
    theme_minimal() +
    labs(x = "Date", y = "Weekly Efficiency") +
    scale_x_date(date_breaks = "1 week", date_labels = "%b %d") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  # arrange together
  gridExtra::grid.arrange(abundance_plot, efficiency_plot)

}

