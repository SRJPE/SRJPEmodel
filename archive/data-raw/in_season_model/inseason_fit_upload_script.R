cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)
library(tidyverse)

for(i in c("battle creek", "clear creek")) {

  print(paste("running for", i))

  site <- case_when(i == "battle creek" ~ "ubc",
                    i == "clear creek" ~ "lcc")

  inseason_inputs <- prepare_inseason_inputs(con, i, site,
                                             covariate_effect = FALSE,
                                             autocorrelation = FALSE)

  print(paste("fitting for", i))
  inseason_fit <- fit_inseason_model(inseason_inputs)

  # push to cloud
  print(paste("Uploading to cloud for", i))
  store_model_fit(con,
                  model_fit_object = inseason_fit,
                  model_inputs = inseason_inputs,
                  results_name = inseason_inputs$model_name,
                  description = paste(i, "model fit object from auto-run tests no covariate effect or autocorrelation"))

}
