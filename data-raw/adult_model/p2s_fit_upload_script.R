cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)
library(tidyverse)

# what streams can we do for P2S?
for(i in c("battle creek", "clear creek", "deer creek", "mill creek")) {

  print(paste("running for", i))

  P2S_inputs <- prepare_P2S_inputs(i, "wy_type")
  fit <- fit_passage_to_spawner_model(P2S_inputs)

  # push to cloud
  print(paste("Uploading to cloud for", i))
  store_model_fit(con,
                  storage_account = "jpemodelresults",
                  container_name = "model-results",
                  access_key = Sys.getenv("AZ_CONTAINER_ACCESS_KEY"),
                  model_fit_object = fit,
                  results_name = "p2s",
                  stream = i,
                  covariate = "wy_type",
                  description = paste(i, "model fit object from auto-run tests using wy_type"))

}

test <- SRJPEmodel::get_most_recent_model_output(con)

# plot
for(i in c("battle creek", "clear creek")) {
  inputs <- prepare_P2S_inputs(i, "wy_type")
  generate_results_plot_p2s()
}
# sql code used in process, not relevant anymore ----------------------------------------------

query <- glue::glue_sql("DELETE from model_parameters WHERE location_id IS null AND updated_at > '2025-05-12'")
res <- DBI::dbSendQuery(con, query)
DBI::dbClearResult(res)
