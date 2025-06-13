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
                  model_fit_object = fit,
                  model_inputs = P2S_inputs,
                  results_name = "p2s",
                  description = paste(i, "model fit object from auto-run tests using wy_type"))

}


# sql code used in process, not relevant anymore ----------------------------------------------

query <- glue::glue_sql("DELETE from model_parameters WHERE location_id IS null AND updated_at > '2025-05-12'")
query <- glue::glue_sql("DELETE from model_run")
res <- DBI::dbSendQuery(con, query)
DBI::dbClearResult(res)
