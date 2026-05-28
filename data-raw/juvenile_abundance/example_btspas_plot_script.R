library(tidyverse)

# connect to db
cfg <- config::get()
con <- DBI::dbConnect(
  RPostgres::Postgres(),
  dbname   = cfg$db_name,
  host     = cfg$db_host,
  port     = cfg$db_port,
  user     = cfg$db_user,
  password = cfg$db_password
)
on.exit(DBI::dbDisconnect(con), add = TRUE)

pcap_all_sites <- get_model_fit(results_name = "pcap_all_sites", con = con)

plots <- plot_pCap_all_sites(pcap_all_sites)
plots$obs_vs_pred_plot
plots$site_mean_plot
plots$flow_plot
plots$hist_plot

pcap_kl <- get_model_fit(results_name = "pcap_one_site",
                         con = con,
                         site_selection = "knights landing")

plots_kl <- plot_pCap_main(site_name = "knights landing",
                           pcap = pcap_kl)
plots_kl$hist_plot
plots_kl$flow_plot

pcap_ubc_2010 <- get_model_fit(results_name = "abundance",
                               con = con,
                               site = "ubc",
                               run_year = 2010)
plots_ubc_2010 <- plot_pCap_site_year(site = "ubc",
                                      run_year = 2010,
                                      pcap = pcap_all_sites,
                                      abundance = pcap_ubc_2010)
plots_ubc_2010$pcap_plot
plots_ubc_2010$abundance_plot

plots_re <- plot_year_re_with_effort(pcap = pcap_all_sites,
                                     model_type = "all_sites")
plots_re$yr_re_plot
plots_re$site_year_plot
