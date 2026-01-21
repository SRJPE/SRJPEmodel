# Data loading and data frame creation

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management

# Read data -----------------------------------------------------------------------
## Inp data
d_Sac <- read.csv(here('data','SacInp.csv'),  stringsAsFactors = F)
d_FeaBut <-  read.csv(here('data','FeatherButteInp.csv'),  stringsAsFactors = F)

## Cov data for from Flow West 
Envdata_flowwest <- SRJPEdata::forecast_covariates
unique(Envdata_flowwest$name)

maxflow_Sac <- Envdata_flowwest %>% 
              filter(name=="monthly_max_flow" & stream == "sacramento river") %>% 
              mutate(Maxflow = value) %>% 
              select(year,month,Maxflow)

maxflow_Butte <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "butte creek") %>% 
  mutate(Maxflowbut = value) %>% 
  select(year,month,Maxflowbut)

maxflow_Feather <- Envdata_flowwest %>% 
  filter(name=="monthly_max_flow" & stream == "feather river") %>% 
  mutate(Maxflowfea = value) %>% 
  select(year,month,Maxflowfea)

maxflow_ButteFeather <- left_join(maxflow_Butte,maxflow_Feather, by=c("year","month"))

flowexceed_Sac <-  Envdata_flowwest %>% 
                   filter(name== "3_category_flow_exceedance_year_type" & stream == "sacramento river") %>% 
                    mutate(flowexceedtype = text_value,
                           year=water_year) %>% 
                    select(year,flowexceedtype)
  
flowexceed_Butte <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "butte creek") %>% 
  mutate(flowexceedtype_but = text_value,
         year=water_year)%>% 
  select(year,flowexceedtype_but)


flowexceed_Feather <-  Envdata_flowwest %>% 
  filter(name== "3_category_flow_exceedance_year_type" & stream == "feather river") %>% 
  mutate(flowexceedtype_fea = text_value,
         year=water_year)%>% 
  select(year,flowexceedtype_fea) %>% 
  add_row(year=2024,flowexceedtype_fea="Average") # Added this row for now because no value in 2024

flowexceed_ButteFeather <- left_join(flowexceed_Butte,flowexceed_Feather, by="year")

# Clean data -----------------------------------------------------------------------------
#Reorder release groups from order in file so they are ordered by year and release group in upstream-downstream direction
d_Sac_sort <- d_Sac %>% 
  drop_na(fish_length) %>% 
#  filter(StudyID != "DeerCk_Wild_CHK_2020") %>%  # remove this release group because of unusual large sizes compared to other groups
  # Add month number to match with flowwest enviro monthly data
  mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) %>% 
  left_join(maxflow_Sac, by = c('year','month')) %>% 
   arrange(year,rl) %>% 
  left_join(flowexceed_Sac, by = 'year') %>% 
  # Add dummy water year type variable
  mutate(Maxflow.z = scale(Maxflow),
        WY2 = case_when(year %in% c(2013,2015,2016,2018,2020, 2021,2022) ~ 0,
                         year %in% c(2017, 2019, 2023,2024) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
         WY3 = case_when(year %in% c(2015, 2021, 2022) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023,2024) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN, AN-W water year
         Flowexceed_Sac = case_when(flowexceedtype =="Dry" ~ 0,
                                 flowexceedtype == "Average" ~ 1,
                                 flowexceedtype == "Wet" ~ 2),
         firstCap = 1,  # define first capture location, it is always the release location
         length.z = scale(fish_length), #standardized length
         weight.z = scale(fish_weight),#standardized weight
         k.z = scale(fish_k)) %>%  #standardized condition factor
  group_by(FishID) %>% 
  # find last capture location for each fish and each potential capture history ch
  mutate(lastCap = case_when(ch == 10000 ~ 1,
                             ch == 11000 ~ 2,
                             ch %in% c(11100,10100) ~ 3,
                             ch %in% c(11110,10110,10010,11010) ~ 4 ,
                             ch %in% c(11111,10111,11011,11101,11001,10011,10101,10001) ~ 5),
         dist_rlwoodson = case_when(rl == "1S" ~ 91.8, # dist from Battle Creek to Woodson Bridge
                                    rl == "2S" ~ 36.9, # dist from RBDD to Woodson Bridge
                                    rl == "3S" ~ 36.9, # dist from Altube Island to Woodson Bridge, same as RBDD
                                    rl == "4S" ~ 36.9, # dist from RB River Park to Woodson Bridge, same as RBDD
                                    rl == "5S" ~ 36.9, # dist from Irvine Finch to Woodson Bridge, same as RBDD
                                    rl == "6S" ~ 25.7, # dist from Mill Creek to Woodson Bridge
                                    rl == "7S" ~ 16.6), # dist from Deer Creek to Woodson Bridge
         dist_woodsonbutte = 88, # distance in km from Woodson Bridge to Butte Bridge
         dist_buttesac = 170, # distance in km from Butte Bridge to Sac
         dist_sacdelta = 110,
         dist_rlwoodson.z = dist_rlwoodson/100, # standardize distances per 100km
         dist_woodsonbutte.z = dist_woodsonbutte/100,# standardize distances per 100km
         dist_buttesac.z= dist_buttesac/100, # standardize distances per 100km
         dist_sacdelta.z =dist_sacdelta/100) %>% # standardize distances per 100km
  ungroup()

d_FeaBut_sort <- d_FeaBut %>% 
  # Add month number to match with flowwest enviro monthly data
  mutate(month = month(as.Date(fish_release_date,format="%m/%d/%Y"))) %>% 
  left_join(maxflow_ButteFeather, by = c('year','month')) %>% 
  left_join(flowexceed_ButteFeather, by = 'year') %>% 
  arrange(year,rl) %>% 
  # Add dummy water year type variable
  mutate(WY2 = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021) ~ 0,
                         year %in% c(2017, 2019, 2023,2024) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
         WY3 = case_when(year %in% c(2014,2015, 2021) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023,2024) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN, AN-W water year
         firstCap = 1,  # define first capture location, it is always the release location
         length.z = scale(fish_length), #standardized length
         weight.z = scale(fish_weight),#standardized weight
         k.z = scale(fish_k)) %>%  #standardized condition factor)  
  group_by(FishID) %>% 
  # find last capture location for each fish and each potential capture history ch
  mutate(lastCap = case_when(ch == 100 ~ 1,
                             ch == 110 ~ 2,
                             ch == 111 ~ 3,
                             ch == 101 ~ 3),
         trib_ind = case_when(release_location %in% c('FR_Boyds_Rel','FR_Gridley_Rel') ~ 2,
                              TRUE ~ 1),
         dist_rlsac = case_when(rl == "1B" ~ 120, # dist from Butte_Blw_Sanborn to Sac
                                rl == "2B" ~ 103.5, # dist from Laux Road to Sac
                                rl == "3B" ~ 117, # dist from North Weir to Sac
                                rl == "4B" ~ 116.7, # dist from Sanborn Slough to Sac
                                rl == "5B" ~ 78, # dist from Sutter Bypass Weir 2 to Sac
                                rl == "6B" ~ 168.5, # dist from Upper Butte to Sac
                                rl == "1F" ~ 115, # dist from Gridley to Sac
                                rl == "2F" ~ 69), # dist from Boyds to Sac
         dist_sacdelta = 110,
         dist_rlsac.z =dist_rlsac/100, # standardize distances per 100km
         dist_sacdelta.z = dist_sacdelta/100,
         Maxflow = case_when(rl %in% c('1B','2B','3B','4B','5B','6B') ~ Maxflowbut,
                             TRUE ~ Maxflowfea),
         flowexceedtype = case_when(rl %in% c('1B','2B','3B','4B','5B','6B') ~ flowexceedtype_but,
                      TRUE ~ flowexceedtype_fea),
         FlowexceedT = case_when(flowexceedtype =="Dry" ~ 0,
                                flowexceedtype == "Average" ~ 1,
                                flowexceedtype == "Wet" ~ 2)) %>% 
  ungroup() %>% 
  mutate(Maxflow.z = scale(Maxflow))

write.csv(d_Sac_sort,here("outputs", "Sac_data.csv"), row.names=FALSE)
write.csv(d_FeaBut_sort,here("outputs", "FeaBut_data.csv"), row.names=FALSE)
