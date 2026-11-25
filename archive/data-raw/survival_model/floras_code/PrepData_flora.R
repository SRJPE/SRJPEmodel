# Data loading and data frame creation

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management

# Read data -----------------------------------------------------------------------
## Inp data
#d=read.csv(file="SacInp_withcov.csv",header=T)
d_Sac <- read.csv(here('data','SacInp.csv'),  stringsAsFactors = F)
d_FeaBut <-  read.csv(here('data','FeatherButteInp.csv'),  stringsAsFactors = F)

# Clean data -----------------------------------------------------------------------------
#Reorder release groups from order in file so they are ordered by year and release group in upstream-downstream direction
d_Sac_sort <- d_Sac %>% 
  drop_na(fish_length) %>% 
#  filter(StudyID != "DeerCk_Wild_CHK_2020") %>%  # remove this release group because of unusual large sizes compared to other groups
  arrange(year,rl) %>% 
  # Add dummy water year type variable
  mutate(WY2 = case_when(year %in% c(2013,2015,2016,2018,2020, 2021,2022) ~ 0,
                         year %in% c(2017, 2019, 2023) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
         WY3 = case_when(year %in% c(2015, 2021, 2022) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN-AN, W water year
         firstCap = 1, # define first capture location, it is always the release location
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
                                    rl == "3S" ~ 36.9, # same as RBDD
                                    rl == "4S" ~ 25.7, # dist from Mill Creek to Woodson Bridge
                                    rl == "5S" ~ 16.6), # dist from Deer Creek to Woodson Bridge
         dist_woodsonbutte = 88, # distance in km from Woodson Bridge to Butte Bridge
         dist_buttesac = 170, # distance in km from Butte Bridge to Sac
         dist_sacdelta = 110,
         dist_rlwoodson.z = dist_rlwoodson/100, # standardize distances per 100km
         dist_woodsonbutte.z = dist_woodsonbutte/100,# standardize distances per 100km
         dist_buttesac.z= dist_buttesac/100, # standardize distances per 100km
         dist_sacdelta.z =dist_sacdelta/100) %>% # standardize distances per 100km
  ungroup()

d_FeaBut_sort <- d_FeaBut %>% 
  arrange(year,rl) %>% 
  # Add dummy water year type variable
  mutate(WY2 = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021) ~ 0,
                         year %in% c(2017, 2019, 2023) ~ 1), #0 or 1 for 2 water year type categories: dry (C,D,BN) and wet (AN,W) water year 
         WY3 = case_when(year %in% c(2014,2015, 2021) ~ 0,
                         year %in% c(2013, 2016, 2018, 2020) ~ 1,
                         year %in% c(2017, 2019, 2023) ~ 2), #0, 1 or 2 for 3 water year type categories: C, D-BN-AN, W water year
         firstCap = 1, # define first capture location, it is always the release location
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
         dist_sacdelta.z = dist_sacdelta/100) %>% # standardize distances per 100km
  ungroup()

write.csv(d_Sac_sort,here("outputs", "Sac_data.csv"), row.names=FALSE)
write.csv(d_FeaBut_sort,here("outputs", "FeaBut_data.csv"), row.names=FALSE)
