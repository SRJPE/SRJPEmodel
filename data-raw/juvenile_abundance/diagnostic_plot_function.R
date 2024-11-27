library(gridExtra)

# plot code for diagnostic plots
options(mc.cores=parallel::detectCores())
ubc_2002 <- run_single_bt_spas_x_stan(SRJPEmodel::bt_spas_x_bayes_params,
                                      SRJPEdata::weekly_juvenile_abundance_catch_data,
                                      SRJPEdata::weekly_juvenile_abundance_efficiency_data,
                                      site = "ubc",
                                      run_year = 2002,
                                      lifestage = NA, # group
                                      effort_adjust = T)

model_fit_summary_object <- ubc_2002$summary_table
site <- "ubc"
run_year <- 2002



