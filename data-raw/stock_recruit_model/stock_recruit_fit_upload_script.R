cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)
library(tidyverse)

for(i in c("battle creek", "clear creek", "deer creek",
           "feather river", "mill creek", "yuba river")) {

  print(paste("running for", i))
  data_type = case_when(i %in% c("battle creek", "clear creek", "mill creek") ~ "redd",
                        i == "yuba river" ~ "passage",
                        i == "deer creek" ~ "holding",
                        i %in% c("butte creek", "feather river") ~ "carcass")

  sr_inputs <- prepare_stock_recruit_inputs(con, stream = i, adult_data_type = data_type,
                                            covariate = "spawning_above_13_temp_week")
  # modify inits for beta
  if(i == "feather river") {
    sr_inputs$inits[[1]]$beta <- -0.0001
    sr_inputs$inits[[2]]$beta <- -0.0001
    sr_inputs$inits[[3]]$beta <- -0.0001
  }

  print(paste("fitting for", i))
  fit <- fit_stock_recruit_model(sr_inputs)

  # push to cloud
  print(paste("Uploading to cloud for", i))
  store_model_fit(con,
                  model_fit_object = fit,
                  model_inputs = sr_inputs,
                  results_name = "stock_recruit",
                  description = paste(i, "model fit object from auto-run tests using spawning_above_13_temp_week"))

}

# test <- SRJPEmodel::get_most_recent_model_results(con)


# sql code used in process, not relevant anymore ----------------------------------------------

# to_insert <- tibble("definition" = unique(pars$parameter),
#                     "description" = "stock-recruit parameter") |>
#   filter(!definition %in% c("log_lik", "lp__"))
# unique(pars$parameter)
# query <- glue::glue_sql("INSERT INTO parameter (definition, description)
#                         VALUES (
#                               {to_insert$definition},
#                               {to_insert$description}
#                      );",
#                         .con = con)
# for(i in 1:length(query)) {
#   res <- DBI::dbSendQuery(con, query[i])
#   DBI::dbClearResult(res)
# }
#
#
#
# query <- glue::glue_sql("DELETE FROM model_run WHERE updated_at > '2025-04-12'")
# res <- DBI::dbSendQuery(con, query)
# DBI::dbClearResult(res)
#
# query <- glue::glue_sql("DELETE FROM parameter WHERE id = 109")
# res <- DBI::dbSendQuery(con, query)
# DBI::dbClearResult(res)
#
# # fix NA sites in trap location
# to_insert <- tibble("stream" = c("butte creek", "clear creek", "deer creek",
#                                  "feather river", "mill creek", "yuba river"),
#                     "site" = NA,
#                     "site_group" = c("butte creek", "clear creek", "deer creek",
#                                      "feather river", "mill creek", "yuba river"),
#                     "description" = "site is not recorded",
#                     "id" = 74:79)
# query <- glue::glue_sql("INSERT INTO trap_location (stream, site, site_group, description, id)
#                         OVERRIDING SYSTEM VALUE VALUES (
#                         {to_insert$stream},
#                         {to_insert$site},
#                         {to_insert$site_group},
#                         {to_insert$description},
#                         {to_insert$id}
#                         );",
#                         .con = con)
# for(i in 1:length(query)) {
#   res <- DBI::dbSendQuery(con, query[i])
#   DBI::dbClearResult(res)
# }

