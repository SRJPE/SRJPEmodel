# Dataframe loading and setting up model variables

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management

# Load data frame ---------------------------------------------------
d_Sac_sort <- read.csv(here("outputs", "Sac_data.csv"),  stringsAsFactors = F)
d_FeaBut_sort <- read.csv(here("outputs", "FeaBut_data.csv"),  stringsAsFactors = F)

# Set variables for Call_Model function -------------------------------------------------------------------
# For Sac mainstem model
Nind <- dim(d_Sac_sort)[1]
Nreaches <- 4 #1=release-woodson, 2=woodson-butte, 3=butte-sac, 4=sac-delta
Ndetlocs <- 5
Nsz=25 # of size classes to plot size effect over 
rch_covind=c(1,2,3,3)#index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index)
# All sac and tribs years combined
(Year_sac <- sort(unique(d_Sac_sort$year)))
(Year_FeaBut <- sort(unique(d_FeaBut_sort$year)))
(Year_all <- sort(unique(c(Year_sac,Year_FeaBut))))
Nyrs <- length(Year_all)
RelGp <- unique(d_Sac_sort$StudyID)
Nrg <- length(RelGp)
firstCap <- d_Sac_sort$firstCap
lastCap <- d_Sac_sort$lastCap
WY2 <- d_Sac_sort$WY2
WY3 <- d_Sac_sort$WY3
CH <- data.frame(cbind(as.integer(substr(d_Sac_sort$ch, 1, 1)),as.integer(substr(d_Sac_sort$ch,2, 2)),
                       as.integer(substr(d_Sac_sort$ch, 3, 3)),as.integer(substr(d_Sac_sort$ch, 4, 4)),
                       as.integer(substr(d_Sac_sort$ch, 5, 5))))

Rmult <- data.frame(cbind(d_Sac_sort$dist_rlwoodson.z, d_Sac_sort$dist_woodsonbutte.z, # reach distances standardized per 100km
                          d_Sac_sort$dist_buttesac.z,d_Sac_sort$dist_sacdelta.z))
RmultSac <- 0.43 # distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge

yrind <- vector(length=Nind)
rgind <- yrind
for(i in 1:Nind){
  yrind[i] <- which(Year_all == d_Sac_sort$year[i])
  rgind[i] <- which(RelGp == d_Sac_sort$StudyID[i])
}

rgwy3.df <- d_Sac_sort %>% 
            group_by(StudyID) %>% 
            dplyr::summarise(ind = unique(WY3), 
                             year = unique(year),
                             rl= unique(rl)) %>% 
            arrange(year,rl) %>% 
            ungroup()

rgwy3_ind <- rgwy3.df$ind

rgwy2.df <- d_Sac_sort %>% 
  group_by(StudyID) %>% 
  dplyr::summarise(ind = unique(WY2), 
                   year = unique(year),
                   rl= unique(rl)) %>% 
  arrange(year,rl) %>% 
  ungroup()

rgwy2_ind <- rgwy2.df$ind

FL <- d_Sac_sort$fish_length
WGT <- d_Sac_sort$fish_weight
CF <- d_Sac_sort$fish_k

# For Feather/Butte model
Ntribs <- 2 
trib_ind <- d_FeaBut_sort$trib_ind
RelGpT <- unique(d_FeaBut_sort$StudyID)
NrgT <- length(RelGpT)
NindT <- dim(d_FeaBut_sort)[1]
NreachesT <- 2
NdetlocsT <- 3  # of detection locations for tribs = release, sacramento and delta stations
firstCapT <- d_FeaBut_sort$firstCap
lastCapT <-d_FeaBut_sort$lastCap
WY2T <- d_FeaBut_sort$WY2
WY3T <- d_FeaBut_sort$WY3
CHT <- data.frame(cbind(as.integer(substr(d_FeaBut_sort$ch, 1, 1)),as.integer(substr(d_FeaBut_sort$ch,2, 2)),
                        as.integer(substr(d_FeaBut_sort$ch, 3, 3))))

RmultT <- d_FeaBut_sort$dist_rlsac.z # reach distances standardized per 100km
RmultTrib <- c(1.17,0.92) # distances for prediction model correspond to the average distance from all release locations for Butte and Feather respectively

yrindT <- vector(length=NindT)
rgindT <- yrindT
for(i in 1:NindT){
  yrindT[i] <- which(Year_all == d_FeaBut_sort$year[i])
  rgindT[i] <- which(RelGpT == d_FeaBut_sort$StudyID[i])
}

rgwy3_T.df <- d_FeaBut_sort %>% 
             group_by(StudyID) %>% 
              dplyr::summarise(ind = unique(WY3),
                      year = unique(year)) %>% 
              arrange(year) %>% 
              ungroup()

rgwy3_indT <- rgwy3_T.df$ind

rgwy2_T.df <- d_FeaBut_sort %>% 
  group_by(StudyID) %>% 
  dplyr::summarise(ind = unique(WY2),
                   year = unique(year)) %>% 
  arrange(year) %>% 
  ungroup()

rgwy2_indT <- rgwy2_T.df$ind

trib_rg <- vector(length=NrgT) #the tributary index for each release groups
for(irg in 1:NrgT){
  irecs <- which(rgindT == irg) #identify all records with the current release group index irg
  trib_rg[irg] <- d_FeaBut_sort$trib_ind[irecs[1]] #Get the tributary index. Only need first records as all individuals with same irg will be from same trib
}

# Upper Butte 2019 does not have weight information so use average values across all release group instead. 
d_FeaBut_sort <- d_FeaBut_sort %>% 
  mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
         across(fish_k, ~ replace_na(., mean(., na.rm=TRUE))))

FL_T <- d_FeaBut_sort$fish_length
WGT_T <- d_FeaBut_sort$fish_weight 
CF_T <- d_FeaBut_sort$fish_k

# Summary stats -----------------------------------------------------------------------------
d_Sac_summary <- d_Sac_sort %>% 
  group_by(StudyID) %>% 
  rename(studyid = StudyID) %>% 
  dplyr::summarise(Mean_FL = mean(fish_length, na.rm=TRUE),
                   Mean_Wt = mean(fish_weight, na.rm=TRUE),
                   Mean_k = mean(fish_k, na.rm=TRUE),
                   Min_FL = min(fish_length, na.rm=TRUE),
                   Min_Wt = min(fish_weight, na.rm=TRUE),
                   Max_FL = max(fish_length, na.rm=TRUE),
                   Max_Wt = max(fish_weight, na.rm=TRUE),
                   CV_FL = sd(fish_length, na.rm=TRUE) / mean(fish_length, na.rm=TRUE),
                   CV_Wt = sd(fish_weight, na.rm=TRUE) / mean(fish_weight, na.rm=TRUE),
                   CV_k = sd(fish_k, na.rm=TRUE) / mean(fish_k, na.rm=TRUE),
                   ntot = n(),
                   ndetect_Woodson = sum(as.numeric(substr(as.character(ch), 2, 2))),
                   ndetect_Butte = sum(as.numeric(substr(as.character(ch), 3, 3))),
                   ndetect_Sac = sum(as.numeric(substr(as.character(ch), 4, 4))),
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 5, 5))))

d_FeaBut_summary <- d_FeaBut_sort %>% 
  group_by(StudyID) %>% 
  rename(studyid = StudyID) %>% 
  dplyr::summarise(Mean_FL = mean(fish_length, na.rm=TRUE),
                   Mean_Wt = mean(fish_weight, na.rm=TRUE),
                   Mean_k = mean(fish_k, na.rm=TRUE),
                   Min_FL = min(fish_length, na.rm=TRUE),
                   Min_Wt = min(fish_weight, na.rm=TRUE),
                   Max_FL = max(fish_length, na.rm=TRUE),
                   Max_Wt = max(fish_weight, na.rm=TRUE),
                   CV_FL = sd(fish_length, na.rm=TRUE) / mean(fish_length, na.rm=TRUE),
                   CV_Wt = sd(fish_weight, na.rm=TRUE) / mean(fish_weight, na.rm=TRUE),
                   CV_k = sd(fish_k, na.rm=TRUE) / mean(fish_k, na.rm=TRUE),
                   ntot = n(),
                   ndetect_Sac = sum(as.numeric(substr(as.character(ch), 2, 2))),
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 3, 3))))

write.csv(d_Sac_summary,here("outputs", "Sac_summary.csv"), row.names=FALSE)
write.csv(d_FeaBut_summary,here("outputs", "FeaBut_summary.csv"), row.names=FALSE)
