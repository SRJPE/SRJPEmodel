# re-fit models and upload

# important - need to be running on a PC
# pCap and BT-SPAS-X
source(here::here("data-raw", "juvenile_abundance", "btspasx_fit_upload_script.R"))

# stock recruit
source(here::here("data-raw", "stock_recruit_model", "stock_recruit_fit_upload_script.R"))

# inseason
source(here::here("data-raw", "in_season_model", "inseason_fit_upload_script.R"))
