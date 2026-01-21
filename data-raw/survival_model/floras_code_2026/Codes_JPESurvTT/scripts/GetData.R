# Dataframe loading and setting up model variables

#Clear the workspace
rm(list = ls())

# Prepare the environment ------------------------------------------------------------
library(tidyverse)
library(here) #better working directory management
library(ggpubr)
library(lubridate)

# Liz note - edited the paths here to direct to the right subfolders
# TODO confirm that MaxFTPflow.z was renamed to Maxflow.z and MaxFTPflow renamed to Maxflow

# Load data frame ---------------------------------------------------
d_Sac_sort <- read.csv(here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "Sac_data.csv"),  stringsAsFactors = F)
d_FeaBut_sort <- read.csv(here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "FeaBut_data.csv"),  stringsAsFactors = F)

# Set variables for Call_Model function -------------------------------------------------------------------
######## For Sac mainstem model ######################
Nind <- dim(d_Sac_sort)[1]
Nreaches <- 4 #1=release-woodson, 2=woodson-butte, 3=butte-sac, 4=sac-delta
Ndetlocs <- 5
Nsz <- 25 # Number of size classes to plot size effect over
NsX <- 25  # Number of continuous variable values
Nrst <- 4
rch_covind=c(1,2,3,4)
#rch_covind=c(1,2,3,3)#index pointing to covariate effect for each reach (note Butte-Sac and Sac-Delta have same fixed effect index)
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
d_Sac_sort$DateRelease <- as.Date(d_Sac_sort$fish_release_date, format = "%m/%d/%Y %H:%M:%S")
d_Sac_sort$DoR <- as.numeric(strftime(d_Sac_sort$DateRelease, format = "%j"))
DoR <- d_Sac_sort$DoR
MaxflowSac.z <- d_Sac_sort$Maxflow.z
MaxflowSac <- d_Sac_sort$Maxflow
MaxflowDelta.z <- d_Sac_sort$Maxflow.z
MaxflowDelta <- d_Sac_sort$Maxflow
FlowexceedSac <- d_Sac_sort$Flowexceed_Sac
FL <- d_Sac_sort$fish_length
WGT <- d_Sac_sort$fish_weight
CF <- d_Sac_sort$fish_k

CH <- data.frame(cbind(as.integer(substr(d_Sac_sort$ch, 1, 1)),as.integer(substr(d_Sac_sort$ch,2, 2)),
                       as.integer(substr(d_Sac_sort$ch, 3, 3)),as.integer(substr(d_Sac_sort$ch, 4, 4)),
                       as.integer(substr(d_Sac_sort$ch, 5, 5))))

Rmult <- data.frame(cbind(d_Sac_sort$dist_rlwoodson.z, d_Sac_sort$dist_woodsonbutte.z, # reach distances standardized per 100km
                          d_Sac_sort$dist_buttesac.z,d_Sac_sort$dist_sacdelta.z))
RmultSac <- 0.43 # distance for prediction model corresponds to the average distance from all release locations to Woodson Bridge per 100km
RmultTis <- 1.23 # distance for prediction model from Tisdale to Sacramento per 100km
dist_rstwoodson <- c(124.9512976,93.2391355490999,26.1273401272,18.3412666561999) # Distance from Upper Clear C, Battle C, Mill C, Deer C to Woodson Bridge
dist_rstwoodson.z <- dist_rstwoodson/100
Rmultrst <- dist_rstwoodson.z # distance for prediction model corresponds to the distance from RST locations to Woodson Bridge

yrind <- vector(length=Nind)
rgind <- yrind
for(i in 1:Nind){
  yrind[i] <- which(Year_all == d_Sac_sort$year[i])
  rgind[i] <- which(RelGp == d_Sac_sort$StudyID[i])
}

### Add index for covariates
rgwy3.df <- d_Sac_sort %>%
            group_by(StudyID) %>%
            dplyr::reframe(ind = unique(WY3),
                             year = unique(year)) %>% #rl= unique(rl)
            arrange(year) %>%
            ungroup()

rgwy3_ind <- rgwy3.df$ind

rgwy2.df <- d_Sac_sort %>%
  group_by(StudyID) %>%
  dplyr::reframe(ind = unique(WY2),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

rgwy2_ind <- rgwy2.df$ind

flowexceedSac.df <- d_Sac_sort %>%
  group_by(StudyID) %>%
  dplyr::reframe(ind = unique(Flowexceed_Sac),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

flowexceedSac_ind <- flowexceedSac.df$ind

maxflowSac.df <- d_Sac_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(Maxflow.z),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

maxflowSac_ind <- maxflowSac.df$ind

maxflowDelta.df <- d_Sac_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(Maxflow.z),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

maxflowDelta_ind <- maxflowDelta.df$ind

DoR.df <- d_Sac_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(DoR),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

DoR_ind <- DoR.df$ind

### Travel Time component
#relKM <- c(517.344,461.579,461.579,450.703,441.728) # 'BattleCk_CNFH_Rel','RBDD_Rel','Altube Island','Irvine_finch,'DeerCk_RST_Rel','MillCk_RST_Rel'
# Woodson bridge = 423.2790, Butte bridge = 343.19467, Sacramento = 150.0867

ReachKM <- c(45, 88, 170, 110)
ReachKM_ind <- data.frame(cbind(d_Sac_sort$dist_rlwoodson,
                                d_Sac_sort$dist_woodsonbutte,
                                d_Sac_sort$dist_buttesac,
                                d_Sac_sort$dist_sacdelta))
ReachKMrst <- dist_rstwoodson

#identify records with one or more detections after release and get their FishID.
#Fish not seen after release provide no data for travel time
d_DHSac <- read.csv(here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", 'data',"DetectionHistorySac.csv"),stringsAsFactors = F) %>%   #detection history file with times of date/time of detection
          arrange(FishID)

d_Sac_ord <- d_Sac_sort %>%
            arrange(FishID)

missing_ids <- d_DHSac %>%
  anti_join(d_Sac_sort, by = "FishID") %>%
  pull(FishID)

release_times <- d_DHSac %>%
  filter(GEN == "Releasepoint") %>%
  select(FishID, release_time = min_time)

Woodson_times <- d_DHSac %>%
  filter(GEN == "WoodsonBridge") %>%
  select(FishID, woodson_time = min_time)

Butte_times <- d_DHSac %>%
  filter(GEN == "ButteBridge") %>%
  select(FishID, butte_time = min_time)

Sacramento_times <- d_DHSac %>%
  filter(GEN == "Sacramento") %>%
  select(FishID, sacramento_time = min_time)

Endpoint_times <- d_DHSac %>%
  filter(GEN == "Endpoint") %>%
  select(FishID, endpoint_time = min_time)

TTfR <- release_times %>%
  left_join(Butte_times, by = "FishID") %>%
  left_join(Woodson_times, by = "FishID") %>%
  left_join(Sacramento_times, by = "FishID") %>%
  left_join(Endpoint_times, by = "FishID") %>%
  mutate(TTfR1 = as.numeric(difftime(woodson_time, release_time, units = "days")),
         TTfR2 = as.numeric(difftime(butte_time, release_time, units = "days")),
         TTfR3 = as.numeric(difftime(sacramento_time, release_time, units = "days")),
         TTfR4 = as.numeric(difftime(endpoint_time, release_time, units = "days")),
         TTfR4 = if_else(TTfR4 < 0, NA_real_, TTfR4)) %>%# Replace negative TTfR4 with NA, example fish SP2023-810
  arrange(FishID) %>%
  select(TTfR1,TTfR2,TTfR3,TTfR4)

Nobs_df <- TTfR %>%
  mutate(Obs = 4 - rowSums(is.na(across(everything()))))

Nobs <- sum(Nobs_df$Obs)

ObsTT <- TTfR %>%
  rowwise() %>%
  mutate(row_vec = list(c_across(everything()))) %>%
  pull(row_vec) %>%
  unlist() %>%
  na.omit()

vec <- as.vector(t(TTfR))                    # row-wise flattening
counter <- seq_len(sum(!is.na(vec)))       # counter for non-NA values
vec[!is.na(vec)] <- counter                # replace only non-NA

# Rebuild dataframe with same shape and column names
TTind <- matrix(vec, nrow = nrow(TTfR), byrow = TRUE) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# Save travel time info in csv file
TTfinalSac_df <- data.frame(cbind(d_Sac_ord$FishID,d_Sac_ord$year,d_Sac_ord$fish_release_date, d_Sac_ord$StudyID,
                                  rgind,d_Sac_ord$fish_length, d_Sac_ord$fish_weight,d_Sac_ord$WY3,
                                  d_Sac_ord$ch,round(TTfR,digits=2))) %>%
  rename(c(FishID=d_Sac_ord.FishID,Year=d_Sac_ord.year,RelDate=d_Sac_ord.fish_release_date,
           StudyID=d_Sac_ord.StudyID,RelGp=rgind,FL=d_Sac_ord.fish_length,WGT=d_Sac_ord.fish_weight,
           WY = d_Sac_ord.WY3,CH=d_Sac_ord.ch))

write.csv(TTfinalSac_df ,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "TravTime_Sac.csv"), row.names=FALSE)


summary_TT <-TTfinalSac_df %>%
  group_by(Year,StudyID) %>%
  summarise(
    TTfR1 = sum(!is.na(TTfR1)),
    TTfR2 = sum(!is.na(TTfR2)),
    TTfR3 = sum(!is.na(TTfR3)),
    TTfR4 = sum(!is.na(TTfR4))
  ) %>%
  ungroup()

summary_TT_SacDelta <-TTfinalSac_df %>%
                  group_by(FishID,StudyID) %>%
                  summarise(TTfR3_and_TTfR4_ind = sum(!is.na(TTfR3) & !is.na(TTfR4))) %>%
                  group_by(StudyID) %>%
                  summarise(TTfR3_and_TTfR4 = sum(TTfR3_and_TTfR4_ind))

summary_TT_final <- summary_TT %>%
                    left_join(summary_TT_SacDelta,by="StudyID")


write.csv(summary_TT_final ,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "TravTime_Sac_summary.csv"), row.names=FALSE)

# Figure of fish size for winner fish
(Sac_fishsize_fig <- ggplot()+
    geom_density(data= TTfinalSac_df %>% filter(!is.na(TTfR3)), aes(x=FL),alpha=0.3,adjust=1.3,size=1) +
    geom_density(data= TTfinalSac_df %>% filter(!is.na(TTfR4)), aes(x=FL, fill="red",color="red"),
                 alpha=0.3,adjust=1.3, size=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
#    scale_x_continuous(limits=c(70,140), breaks=seq(70,140,20))
)

(Sac_fishspeed_fig <- ggplot()+
    geom_density(data= TTfinalSac_df, aes(x=TTfR2),alpha=0.3,adjust=1.3,size=1) +
    geom_density(data= TTfinalSac_df %>% filter(!is.na(TTfR4)), aes(x=TTfR4, fill="red",color="red"),
                 alpha=0.3,adjust=1.3, size=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  #    scale_x_continuous(limits=c(70,140), breaks=seq(70,140,20))
)


##### For Feather/Butte model ###################
Ntribs <- 2
trib_ind <- d_FeaBut_sort$trib_ind
RelGpT <- unique(d_FeaBut_sort$StudyID)
NrgT <- length(RelGpT)
NindT <- dim(d_FeaBut_sort)[1]
NrstT <- 3
TribForRST <- c(1,2,2)
NreachesT <- 2
NdetlocsT <- 3  # of detection locations for tribs = release, sacramento and delta stations
firstCapT <- d_FeaBut_sort$firstCap
lastCapT <-d_FeaBut_sort$lastCap
WY2T <- d_FeaBut_sort$WY2
WY3T <- d_FeaBut_sort$WY3
MaxflowT.z <- d_FeaBut_sort$Maxflow.z
MaxflowT <- d_FeaBut_sort$Maxflow
MaxflowDeltaT.z <- d_FeaBut_sort$Maxflow.z
MaxflowDeltaT <- d_FeaBut_sort$Maxflow
FlowexceedT <- d_FeaBut_sort$FlowexceedT
d_FeaBut_sort$DateReleaseT <- as.Date(d_FeaBut_sort$fish_release_date, format = "%m/%d/%Y %H:%M:%S")
d_FeaBut_sort$DoRT <- as.numeric(strftime(d_FeaBut_sort$DateRelease, format = "%j"))
DoRT <- d_FeaBut_sort$DoRT

# Upper Butte 2019 does not have weight information so use average values across all release group instead.
d_FeaBut_sort <- d_FeaBut_sort %>%
  mutate(across(fish_weight, ~ replace_na(., mean(., na.rm=TRUE))),
         across(fish_k, ~ replace_na(., mean(., na.rm=TRUE))))

FL_T <- d_FeaBut_sort$fish_length
WGT_T <- d_FeaBut_sort$fish_weight
CF_T <- d_FeaBut_sort$fish_k

CHT <- data.frame(cbind(as.integer(substr(d_FeaBut_sort$ch, 1, 1)),as.integer(substr(d_FeaBut_sort$ch,2, 2)),
                        as.integer(substr(d_FeaBut_sort$ch, 3, 3))))

RmultT <- d_FeaBut_sort$dist_rlsac.z # reach distances standardized per 100km
RmultTrib <- c(1.17,0.92) # distances for prediction model correspond to the average distance from all release locations for Butte and Feather respectively per 100km
dist_rstsac <- c(201.583633469999,89.5817927391999, 131.0731) # Distance from PPDD, Yuba Hallwood, and average distance of eye riffle, steep riffle, and gateway riffle to Delta entry
dist_rstsac.z =dist_rstsac/100 # standardize distances per 100km
RmultTrst <- dist_rstsac.z

yrindT <- vector(length=NindT)
rgindT <- yrindT
for(i in 1:NindT){
  yrindT[i] <- which(Year_all == d_FeaBut_sort$year[i])
  rgindT[i] <- which(RelGpT == d_FeaBut_sort$StudyID[i])
}

trib_rg <- vector(length=NrgT) #the tributary index for each release groups
for(irg in 1:NrgT){
  irecs <- which(rgindT == irg) #identify all records with the current release group index irg
  trib_rg[irg] <- d_FeaBut_sort$trib_ind[irecs[1]] #Get the tributary index. Only need first records as all individuals with same irg will be from same trib
}

#### Add index for covariates
rgwy3_T.df <- d_FeaBut_sort %>%
             group_by(StudyID) %>%
              dplyr::reframe(ind = unique(WY3),
                      year = unique(year)) %>%
              arrange(year) %>%
              ungroup()

rgwy3_indT <- rgwy3_T.df$ind

rgwy2_T.df <- d_FeaBut_sort %>%
  group_by(StudyID) %>%
  dplyr::reframe(ind = unique(WY2),
                   year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

rgwy2_indT <- rgwy2_T.df$ind

flowexceedT.df <- d_FeaBut_sort %>%
  group_by(StudyID) %>%
  dplyr::reframe(ind = unique(FlowexceedT),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

flowexceedT_ind <- flowexceedT.df$ind

maxflowT.df <- d_FeaBut_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(Maxflow.z),
                 year = unique(year),
                 month=unique(month)) %>%
  arrange(year,month) %>%
  ungroup()

maxflowT_ind <- maxflowT.df$ind

maxflowDeltaT.df <- d_FeaBut_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(Maxflow.z),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

maxflowDeltaT_ind <- maxflowDeltaT.df$ind

DoRT.df <- d_FeaBut_sort %>%
  group_by(FishID) %>%
  dplyr::reframe(ind = unique(DoRT),
                 year = unique(year)) %>%
  arrange(year) %>%
  ungroup()

DoRT_ind <- DoRT.df$ind

### Travel Time component
#identify records with one or more detections after release and get their FishID.
#Fish not seen after release provide no data for travel time
d_DHFea <- read.csv(here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", 'data',"DetectionHistoryFea.csv"),stringsAsFactors = F) #detection history file with times of date/time of detection
d_DHBut <- read.csv(here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", 'data',"DetectionHistoryBut.csv"),stringsAsFactors = F) #detection history file with times of date/time of detection
d_DHFeaBut <- data.frame(rbind(d_DHBut,d_DHFea)) %>%
              arrange(FishID)

# relRKM <- 'UpperButte_RST_Rel' =340.854, 'Butte_Blw_Sanborn_Rel' =293.010, 'Laux Rd'=276.000 ,'North_Weir_Rel' =289.490
#'Sanborn_Slough_Rel' =288.650, 'SutterBypass_Weir2_RST_Rel' =249.541
#' 'FR_Gridley_Rel' = 287.387, 'FR_Boyds_Rel' = 240.755
#' # Sacramento = 150.0867

ReachKMT <- data.frame(rbind(cbind(117,110),cbind(92,110)))
ReachKMT_ind <- data.frame(cbind(d_FeaBut_sort$dist_rlsac,
                                 d_FeaBut_sort$dist_sacdelta))
ReachKMTrst <- dist_rstsac

# ReachKMT_rg <- data.frame(cbind(RelGpT,rep(0,NrgT),rep(110,NrgT))) %>%
#   mutate(R1 = case_when(substr(RelGpT, 1, 2) =="FR" ~ 92, # to update with right numbers
#                         TRUE ~ 117),  #to update with right numbers
#          R2 = V3) %>%
#   select(R1,R2)

release_timesT <- d_DHFeaBut %>%
  filter(GEN == "Releasepoint") %>%
  select(FishID, release_time = min_time)

Sacramento_timesT <- d_DHFeaBut %>%
  filter(GEN == "Sacramento") %>%
  select(FishID, sacramento_time = min_time)

Endpoint_timesT <- d_DHFeaBut %>%
  filter(GEN == "Endpoint") %>%
  select(FishID, endpoint_time = min_time)

TTfRT <- release_timesT %>%
  left_join(Sacramento_timesT, by = "FishID") %>%
  left_join(Endpoint_timesT, by = "FishID") %>%
  mutate(TTfR1 = as.numeric(difftime(sacramento_time, release_time, units = "days")),
         TTfR2 = as.numeric(difftime(endpoint_time, release_time, units = "days"))) %>%
  arrange(FishID) %>%
  select(TTfR1,TTfR2)

NobsT_df <- TTfRT %>%
  mutate(Obs = 2 - rowSums(is.na(across(everything()))))
NobsT <- sum(NobsT_df$Obs)

ObsTTT <- TTfRT %>%
  rowwise() %>%
  mutate(row_vec = list(c_across(everything()))) %>%
  pull(row_vec) %>%
  unlist() %>%
  na.omit()

vecT <- as.vector(t(TTfRT))                    # row-wise flattening
counterT <- seq_len(sum(!is.na(vecT)))       # counter for non-NA values
vecT[!is.na(vecT)] <- counterT                # replace only non-NA

# Rebuild dataframe with same shape and column names
TTindT <- matrix(vecT, nrow = nrow(TTfRT), byrow = TRUE) %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ replace_na(., 0)))

# Save travel time info in csv file
d_FeaBut_ord <- d_FeaBut_sort %>%
  arrange(FishID)

TTfinalFeaBut_df <- data.frame(cbind(d_FeaBut_ord$FishID,d_FeaBut_ord$year,d_FeaBut_ord$fish_release_date,
                                     d_FeaBut_ord$StudyID,rgindT,d_FeaBut_ord$fish_length, d_FeaBut_ord$fish_weight,
                                     d_FeaBut_ord$WY3,d_FeaBut_ord$ch,round(TTfRT,digits=2))) %>%
  rename(c(FishID=d_FeaBut_ord.FishID,Year=d_FeaBut_ord.year,RelDate=d_FeaBut_ord.fish_release_date,
           StudyID=d_FeaBut_ord.StudyID,RelGp=rgindT,FL=d_FeaBut_ord.fish_length,WGT=d_FeaBut_ord.fish_weight,
           WY = d_FeaBut_ord.WY3,CH=d_FeaBut_ord.ch))

write.csv(TTfinalFeaBut_df ,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "TravTime_FeaBut.csv"), row.names=FALSE)

summary_TTT <-TTfinalFeaBut_df %>%
  group_by(Year,StudyID) %>%
  summarise(
    TTfR1 = sum(!is.na(TTfR1)),
    TTfR2 = sum(!is.na(TTfR2))
  ) %>%
  ungroup()

summary_TTT_SacDelta <-TTfinalFeaBut_df %>%
  group_by(FishID,StudyID) %>%
  summarise(TTfR1_and_TTfR2_ind = sum(!is.na(TTfR1) & !is.na(TTfR2))) %>%
  group_by(StudyID) %>%
  summarise(TTfR1_and_TTfR2 = sum(TTfR1_and_TTfR2_ind))

summary_TTT_final <- summary_TTT %>%
  left_join(summary_TTT_SacDelta,by="StudyID")


write.csv(summary_TTT_final ,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "TravTime_FeaBut_summary.csv"), row.names=FALSE)

# Summary stats -----------------------------------------------------------------------------
d_Sac_summary <- d_Sac_sort %>%
  group_by(StudyID) %>%
  rename(studyid = StudyID) %>%
  dplyr::reframe(Mean_FL = mean(fish_length, na.rm=TRUE),
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
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 5, 5)))) %>%
  mutate(watershed ="Sac")

d_FeaBut_summary <- d_FeaBut_sort %>%
  group_by(StudyID) %>%
  rename(studyid = StudyID) %>%
  dplyr::reframe(Mean_FL = mean(fish_length, na.rm=TRUE),
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
                   ndetect_Delta = sum(as.numeric(substr(as.character(ch), 3, 3)))) %>%
  mutate(watershed ="FeaBut")

write.csv(d_Sac_summary,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "Sac_summary.csv"), row.names=FALSE)
write.csv(d_FeaBut_summary,here("data-raw", "survival_model", "floras_code_2026", "Codes_JPESurvTT", "outputs", "FeaBut_summary.csv"), row.names=FALSE)

# Summary figures -------------------------------------------------------------
### mainstem vs trib fish size
d_SacFeaBut_summary <- data.frame(rbind(d_Sac_summary[,c(1:12,17)],d_FeaBut_summary[,c(1:12,15)]))

(fishsizedens <- ggplot()+
    geom_density(data=d_SacFeaBut_summary, aes(x=Mean_FL, fill = watershed,color=watershed),
                 alpha=0.3,adjust=1.3, size=1) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    scale_x_continuous(limits=c(70,140), breaks=seq(70,140,20))
)

 ggsave('figures/fishsizdens.jpg',plot=fishsizedens, dpi = 350, height = 4, width = 5)
