#### Adult modeling

The adult model repo consists of modeling scripts to produce estimates of spawner abundance for seven tributaries and a vignette describing those methods and results. They should be run in the following order: 

* `adult_model_data_prep.R`: pulls in raw adult model data sources; performs regressions to identify environmental covariates; formats adult data objects for input into STAN model. 
* `predicted_spawners_STAN.R`: pulls in prepared adult data objects; runs predicted_spawners STAN model for Battle, Clear, Mill, and Deer Creeks; formats output into a summarized table with parameter estimates of interest. 
* `adult_table_for_spawn_recruit.R`: pulls in prepared adult data objects and output of predicted_spawners STAN model; formats results into a table of `spawner_count` by stream and `data_type` for input into the SR JPE stock-recruit model.

The `predicted_spawners_STAN.R` code takes a while to run and fit the four tributaries.

All results and data files are written to the jpe-dev-bucket and to the `data-raw/adult_model/adult_model_data/` file in the repo.
