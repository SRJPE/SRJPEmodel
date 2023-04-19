# create plot for adult data report using butte creek estimates (not raw data)
library(tidyverse)
library(readxl)
library(googleCloudStorageR)

# read in data ------------------------------------------------------------

butte_estimates_raw <- read_xlsx(here::here("data-raw", "qc-markdowns", 
                                            "adult-holding-redd-and-carcass-surveys", "butte-creek",
                                            "butte_creek_historic_estimates.xlsx"),
                                 skip = 1)

butte_estimates <- butte_estimates_raw |> 
  rename(carcass = `Post Spawn Est.`,
         upstream_passage = Vaki) |> 
  select(year = Year, carcass, upstream_passage) |> 
  pivot_longer(cols = c(carcass, upstream_passage),
               names_to = "data_type",
               values_to = "count") |> 
  mutate(stream = "butte creek",
         data_type = if_else(data_type == "upstream_passage", "upstream passage", data_type)) |> 
  glimpse()


# load in holding data ----------------------------------------------------

# upload data
gcs_auth(json_file = Sys.getenv("GCS_AUTH_FILE"))
gcs_global_bucket(bucket = Sys.getenv("GCS_DEFAULT_BUCKET"))

# holding
holding_butte_raw <- read_csv(gcs_get_object(object_name = "standard-format-data/standard_holding.csv",
                                   bucket = gcs_get_global_bucket())) |> 
  filter(stream == "butte creek") |> 
  glimpse()

holding_butte <- holding_butte_raw |> 
  group_by(year) |> 
  summarise(count = sum(count, na.rm = T)) |> 
  mutate(data_type = "holding",
         stream = "butte creek") |> 
  select(year, count, data_type, stream) |> 
  glimpse()


# join estimates with raw holding -----------------------------------------
# add blank rows for holding for later years
holding_blanks <- tibble(year = seq(2018,2022),
                         data_type = "holding",
                         stream = "butte creek",
                         count = 0)
butte <- bind_rows(holding_butte, butte_estimates,
                   holding_blanks) |> 
  glimpse()


# plot --------------------------------------------------------------------

colors <- c("#DD5E66FF","#FEB72DFF", "#5901A5FF", "#F7E225FF") # for data type (upstream, holding, redd annual, carcass)
data_type_palette = c("upstream passage" = colors[1],
                      "holding" = colors[2],
                      "redd" = colors[3],
                      "carcass" = colors[4])

butte |>
  ggplot(aes(x=year, y=count, fill = data_type)) +
  geom_col(position = "dodge", width = 0.75) +
  xlab("Year") + ylab("Count") +
  ggtitle("Butte Creek: Counts by Data Type") +
  scale_fill_manual(values = data_type_palette, drop = T) +
  theme_minimal() +
  labs(fill = "") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position="bottom")
