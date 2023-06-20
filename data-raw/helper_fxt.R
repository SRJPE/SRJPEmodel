# Author: Tom Pham modified by Flora Cordoleani
# Description: This script contains all of the datasets and functions required
# to perform the survival analyses (reach per 10km, region raw, region per 10km,
# and cumulative), to plot them as figures, and to format them for use in 
# reports. These act as standardized methods to repeat these analyses for 
# each StudyID in the BOR-EAT 2019-2020 Report. 

library(RMark) # For running program MARK
library(tidyverse) # Data manipulations
library(rerddap) # To retrieve NOAA ERDDAP data
library(lubridate) # Date time manipulations
library(clusterPower) # Need for cumulative survival
library(leaflet) # To visualize receiver locations quickly on a map
library(vroom) # Read CSV's quickly

###### Load TaggedFish and ReceiverDeployments tables through ERDDAP ---------------------------------------------------------------------
# TaggedFish is used to get all FishIDs for a study, and release information
# ReceiverDeployments is used to get Region

my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_taggedfish', url = my_url)
TaggedFish <- tabledap(JSATSinfo, url = my_url)  

JSATSinfo <- info('FED_JSATS_receivers', url = my_url)
ReceiverDeployments <- tabledap(JSATSinfo, url = my_url)

# Establish ERDDAP url and database name
my_url <- "https://oceanview.pfeg.noaa.gov/erddap/"
JSATSinfo <- info('FED_JSATS_detects', url = my_url)

# Retrieve list of all studyIDs on FED_JSATS
studyid_list <- tabledap(JSATSinfo,
                         fields = c('study_id'),
                         url = my_url,
                         distinct = TRUE
) %>% 
  filter(study_id != "2017_BeaconTag") %>% 
  pull(study_id)

# Load functions ----------------------------------------------------------
get_receiver_GEN <- function(all_detections) {
  # Get a list of all receiver sites and metadata for a given detections df
  #
  # Arguments:
  #  all_detections: detections df 
  #
  # Return:
  #  df of receiver sites along with RKM, Lat, Lon, Region
  
  reach.meta <- all_detections %>% 
    bind_rows() %>% 
    distinct(GEN, GenRKM, GenLat, GenLon, Region) %>% 
    # Necessary because detections files shows differing RKM, Lat, Lon for some 
    # GEN sometimes
    group_by(GEN) %>% 
    summarise(
      GenRKM = mean(GenRKM),
      GenLat = mean(GenLat),
      GenLon = mean(GenLon),
      Region = first(Region)
    ) %>% 
    arrange(desc(GenRKM))
}

aggregate_GEN_Sac <- function(detections, 
                          replace_dict = list(replace_with = list(c("Releasepoint"),
                                                                  c("WoodsonBridge"),
                                                                 # c("ButteBridge"),
                                                                  c("Sacramento"),
                                                                  c("Endpoint")),
                                              replace_list = list(c("BattleCk_CNFH_Rel","RBDD_Rel","RBDD_Rel_Rec",
                                                                    "Altube Island","Abv_Altube1",
                                                                    "MillCk_RST_Rel","MillCk2_Rel","DeerCk_RST_Rel"), 
                                                                  c("Mill_Ck_Conf","Abv_WoodsonBr","Blw_Woodson"),
                                                                 # c("ButteBr","BlwButteBr"),
                                                                  c("TowerBridge","I80-50_Br",
                                                                    "ToeDrainBase","Hwy84Ferry"),
                                                                  c("BeniciaE","BeniciaW",
                                                                    "ChippsE","ChippsW"
                                                                        )))) {
  
  # Replace GEN in detections df according to replace_list, basically aggregate
  # sites into one. By default this is done for ChippsE/W, BeniciaE/W, and
  # SacTrawl1/2
  #
  # Arguments:
  #  detections: a detections df
  #  replace_dict: list of receiver locations to aggregate and, the aggregated
  #  name
  #     
  # Return:
  #  a detections df that has replaced list of GEN with the aggregated GEN,
  #  mean RKM, mean Lat, mean Lon. Creates reach.meta.aggregate which is the list
  #  of receiver sites with new aggregated GEN's, along with RKM, Lat, Lon
  
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <<- reach.meta
  
  # Walk through each key/pair value
  for (i in 1:length(replace_dict$replace_with)) {
    # Unlist for easier to use format
    replace_list <- unlist(replace_dict[[2]][i])
    replace_with <- unlist(replace_dict[[1]][i])
  
      replace <- reach.meta %>% 
        select(GEN, GenRKM, GenLat, GenLon, Region) %>% 
        filter(GEN %in% c(replace_list, replace_with)) %>%
        distinct() %>% 
        select(-GEN) %>% 
        group_by(Region) %>% 
        summarise_all(mean)

    # Replace replace_list GENs name with replace_with GEN, and replace all of 
    # their genrkm with the averaged val
    detections <- detections %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)
      )
    
    
    # This new df shows receiver metadata and reflects the aggregation done
    reach.meta.aggregate <<- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)) %>% 
      distinct()
    
    detections <- detections %>%
      filter(GEN %in% reach.meta.aggregate$GEN) %>% 
      group_by(FishID) %>% 
      arrange(FishID,desc(GenRKM)) %>% 
      ungroup()
  }
  detections
}

aggregate_GEN_Feather <- function(detections, 
                          replace_dict = list(replace_with = list(c("Releasepoint"),
                                                                  c("Sacramento"),
                                                                  c("Endpoint")),
                                              replace_list = list(c("FR_Gridley_Rel","FR_Boyds_Rel","FR_Boyds_Rel_Rec"), 
                                                                  c("TowerBridge","I80-50_Br",
                                                                    "ToeDrainBase","Hwy84Ferry"),
                                                                  c("BeniciaE","BeniciaW",
                                                                    "ChippsE","ChippsW"
                                                                  )))) {
  
  # Replace GEN in detections df according to replace_list, basically aggregate
  # sites into one. By default this is done for ChippsE/W, BeniciaE/W, and
  # SacTrawl1/2
  #
  # Arguments:
  #  detections: a detections df
  #  replace_dict: list of receiver locations to aggregate and, the aggregated
  #  name
  #     
  # Return:
  #  a detections df that has replaced list of GEN with the aggregated GEN,
  #  mean RKM, mean Lat, mean Lon. Creates reach.meta.aggregate which is the list
  #  of receiver sites with new aggregated GEN's, along with RKM, Lat, Lon
  
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <<- reach.meta
  
  # Walk through each key/pair value
  for (i in 1:length(replace_dict$replace_with)) {
    # Unlist for easier to use format
    replace_list <- unlist(replace_dict[[2]][i])
    replace_with <- unlist(replace_dict[[1]][i])
    
    replace <- reach.meta %>% 
      select(GEN, GenRKM, GenLat, GenLon, Region) %>% 
      filter(GEN %in% c(replace_list, replace_with)) %>%
      distinct() %>% 
      select(-GEN) %>% 
      group_by(Region) %>% 
      summarise_all(mean)
    
    # Replace replace_list GENs name with replace_with GEN, and replace all of 
    # their genrkm with the averaged val
    detections <- detections %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)
      )
    
    
    # This new df shows receiver metadata and reflects the aggregation done
    reach.meta.aggregate <<- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)) %>% 
      distinct()
    
    detections <- detections %>%
      filter(GEN %in% reach.meta.aggregate$GEN) %>% 
      group_by(FishID) %>% 
      arrange(FishID,desc(GenRKM)) %>% 
      ungroup()
  }
  detections
}

aggregate_GEN_Butte <- function(detections, 
                          replace_dict = list(replace_with = list(c("Releasepoint"),
                                                                  c("Sacramento"),
                                                                  c("Endpoint")),
                                              replace_list = list(c("UpperButte_RST_Rel","UpperButte_RST","UpperButte_SKWY",
                                                                    "SutterBypass_Weir2_RST_Rel","SutterBypass Weir2 RST"), 
                                                                  c("TowerBridge","I80-50_Br",
                                                                    "ToeDrainBase","Hwy84Ferry"),
                                                                  c("BeniciaE","BeniciaW",
                                                                    "ChippsE","ChippsW"
                                                                  )))) {
  
  # Replace GEN in detections df according to replace_list, basically aggregate
  # sites into one. By default this is done for ChippsE/W, BeniciaE/W, and
  # SacTrawl1/2
  #
  # Arguments:
  #  detections: a detections df
  #  replace_dict: list of receiver locations to aggregate and, the aggregated
  #  name
  #     
  # Return:
  #  a detections df that has replaced list of GEN with the aggregated GEN,
  #  mean RKM, mean Lat, mean Lon. Creates reach.meta.aggregate which is the list
  #  of receiver sites with new aggregated GEN's, along with RKM, Lat, Lon
  
  # Make a copy of reach.meta (receiver metadata)
  reach.meta.aggregate <<- reach.meta
  
  # Walk through each key/pair value
  for (i in 1:length(replace_dict$replace_with)) {
    # Unlist for easier to use format
    replace_list <- unlist(replace_dict[[2]][i])
    replace_with <- unlist(replace_dict[[1]][i])
    
    replace <- reach.meta %>% 
      select(GEN, GenRKM, GenLat, GenLon, Region) %>% 
      filter(GEN %in% c(replace_list, replace_with)) %>%
      distinct() %>% 
      select(-GEN) %>% 
      group_by(Region) %>% 
      summarise_all(mean)
    
    # Replace replace_list GENs name with replace_with GEN, and replace all of 
    # their genrkm with the averaged val
    detections <- detections %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)
      )
    
    
    # This new df shows receiver metadata and reflects the aggregation done
    reach.meta.aggregate <<- reach.meta.aggregate %>% 
      mutate(
        GEN = ifelse(GEN %in% replace_list, replace_with, GEN),
        GenRKM = ifelse(GEN %in% c(replace_list, replace_with), replace$GenRKM, GenRKM),
        GenLat = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLat, GenLat),
        GenLon = ifelse(GEN %in% c(replace_list, replace_with), replace$GenLon, GenLon),
        Region = ifelse(GEN %in% c(replace_list, replace_with), replace$Region, Region)) %>% 
      distinct()
    
    detections <- detections %>%
      filter(GEN %in% reach.meta.aggregate$GEN) %>% 
      group_by(FishID) %>% 
      arrange(FishID,desc(GenRKM)) %>% 
      ungroup()
  }
  detections
}


make_DH <- function(detections) {
  # Make a detection history df
  # Get earliest detection for each fish at each GEN
  min_detects <- detections %>% 
    filter(GEN %in% reach.meta.aggregate$GEN) %>% 
    group_by(FishID, GEN, GenRKM,Region) %>% 
    summarise(
      min_time = min(time)
    ) %>% 
    arrange(
      FishID, min_time
    )
  
  min_detects
  
}

make_EH <- function(detections) {
  # Make an encounter history df
  #
  # Arguments:
  #  detections: a detections df
  #     
  # Return:
  #  Encounter history df. A matrix of every fish tagged for a given studyID
  #  at every given receiver site (that is in reach.meta.aggregate) and whether
  #  it was present 1 or absent 0 in the detection df
  
  # Get earliest detection for each fish at each GEN
  min_detects <- detections %>% 
    filter(GEN %in% reach.meta.aggregate$GEN) %>% 
    group_by(FishID, GEN, GenRKM) %>% 
    summarise(
      min_time = min(time)
    ) %>% 
    arrange(
      FishID, min_time
    )
  
  # Get list of all tagged fish for the studyID
  fish <- TaggedFish %>% 
    filter(study_id == detections$StudyID[1]) %>% 
    arrange(fish_id) %>% 
    pull(fish_id)
  
  # Create matrix of all combinations of fish and GEN
  EH <- expand.grid(
    fish,
    reach.meta.aggregate$GEN, stringsAsFactors = FALSE 
  )
  
  names(EH) <- c('FishID', 'GEN')  
  
  # Add col detect to min_detects, these fish get a 1
  min_detects$detect <- 1
  
  # Join in detections to the matrix, fish detected a GEN will be given a 1
  # otherwise it will be given a 0
  EH <- EH %>% 
    left_join(
      min_detects %>% 
        select(
          FishID, GEN, detect
        ), by = c("FishID", "GEN")
    ) %>% 
    # Replace NA with 0 https://stackoverflow.com/questions/28992362/dplyr-join-define-na-values
    mutate_if(
      is.numeric, coalesce, 0
    )
  
  # Reshape the df wide, so that columns are GEN locations, rows are fish, 
  # values are 1 or 0 for presence/absence
  EH <- reshape(EH, idvar = 'FishID', timevar = 'GEN', direction = 'wide')
  colnames(EH) <- gsub('detect.', '', colnames(EH))
  # Manually make the release column a 1 because all fish were released there
  # sometimes detections df does not reflect that accurately
  EH[2] <- 1
  EH
}

create_inp <- function(detections, EH) { 
  # Create an inp df
  #
  # Arguments:
  #  detections: a detections df (make sure to use the aggregated version)
  #  EH: an encounter history df
  #     
  # Return:
  #  inp df i.e. Fish01 | 11101, a record of a fish and it's presence/absence
  #  at each given receiver location. 
  
  EH.inp <- EH %>% 
    # Collapse the encounter columns into a single column of 1's and 0's
    unite("ch", 2:(length(EH)), sep ="") %>% 
    # Use the detections df to get the StudyID assignment
    mutate(StudyID = unique(detections$StudyID))
  EH.inp
}
