# prep data for josh's code
# fix feather sites

library(googleCloudStorageR)
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))
gcs_get_object(object_name = "jpe-model-data/feather_annual_site_selection.csv",
               bucket = gcs_get_global_bucket(),
               saveToDisk = here::here("data-raw", "first_draft", "feather_annual_site_selection.csv"),
               overwrite = TRUE)

feather_site_selection <- read.csv(here::here("data-raw", "first_draft", "feather_annual_site_selection.csv")) |>
  mutate(site_group = case_when(site_group == "feather river hfc" ~ "upper feather hfc",
                                site_group == "feather river lfc" ~ "upper feather lfc",
                                TRUE ~ site_group),
    stream_site = paste0(stream, "_", site_group)) |>
  rename(site_to_use = site)

RST_input <- read.table(here::here("data-raw", "first_draft", "RST_Input.txt"), header = T) |>
  glimpse()

RST_input_new <- RST_input |>
  left_join(feather_site_selection, by = c("CalYr" = "year", "Trib" = "stream_site")) |>
  mutate(Trib = ifelse(Trib %in% c("feather river_upper feather hfc", "feather river_upper feather lfc"),
                       paste0("feather river_", site_to_use),
                       Trib)) |>
  select(-c(stream, site_group, site_to_use)) |>
  glimpse()

write.table(RST_input_new, here::here("data-raw", "first_draft", "RST_input_new.txt"))

# process data from josh's code

read_bugs_output <- function(filepath) {

  total_abundance <- read.table(filepath, header = T) |>
    pull(Ntot)

  mean_total_abundance = mean(total_abundance)
  ucl_total_abundance = quantile(total_abundance, 0.975)
  lcl_total_abundance = quantile(total_abundance, 0.025)

  metadata <- str_remove(filepath, pattern ="/Users/liz/Documents/code/SRJPEmodel/data-raw/first_draft/OutSpecPriors/")
  year <- as.numeric(gsub("\\D", "", metadata))
  stream <- case_when(str_detect(metadata, "battle creek") ~ "battle creek",
                      str_detect(metadata, "clear creek") ~ "clear creek",
                      str_detect(metadata, "feather river") ~ "feather river")
  site <- case_when(str_detect(metadata, "ubc") ~ "ubc",
                    str_detect(metadata, "lcc") ~ "lcc",
                    str_detect(metadata, "ucc") ~ "ucc",
                    str_detect(metadata, "eye riffle") ~ "eye riffle",
                    str_detect(metadata, "gateway riffle") ~ "gateway riffle",
                    str_detect(metadata, "herringer riffle") ~ "herringer riffle",
                    str_detect(metadata, "steep riffle") ~ "steep riffle")

  output <- tibble(filepath = filepath,
                   year = year,
                   stream = stream,
                   site = site,
                   mean_total_abundance = mean_total_abundance,
                   ucl_total_abundance = ucl_total_abundance,
                   lcl_total_abundance = lcl_total_abundance)

  return(output)
}

output_files <- list.files(path = here::here("data-raw", "first_draft", "OutSpecPriors"),
                           full.names = TRUE) |>
  as_tibble() |>
  filter(str_detect(value, "_post.out"))

annual_bt_spas_estimates <- output_files$value |>
  map_dfr(read_bugs_output) |>
  glimpse()
