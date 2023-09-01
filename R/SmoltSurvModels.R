##%##########################################################################################################%##
#                                                                                                              #
#       Survival Analysis for spring-run Chinook through the Sacramento River, Butte Creek and Feather River   #
#                                                                                                              #
##%##########################################################################################################%##
# Authors: Jeremy Notch and Tom Pham, modified by Flora Cordoleani
# Script and data info: This script performs survival analysis using CJS model
# and outputs detection probabilities using RMark
# Data source: Data extracted from NOAA ERDDAP data server (oceanview)

#source('Script/helper_fxt.R')
source(here::here("data-raw", "helper_fxt.R"))

library(here)
library(reshape2)
library(ggplot2)
library(readxl)
library(locfit)
library(gtools)
library(mapview)
library(ggpubr)


# unzip source data -------------------------------------------------------

unzip(here::here("data-raw", "survival_model_data", "ButteCreekAllYears.zip"),
      exdir = here::here("data-raw", "survival_model_data"))
unzip(here::here("data-raw", "survival_model_data", "FeatherRiverAllYears.zip"),
      exdir = here::here("data-raw", "survival_model_data"))
unzip(here::here("data-raw", "survival_model_data", "SacRiverAllYears.zip"),
      exdir = here::here("data-raw", "survival_model_data"))

#################### Sacramento River model ##############################################

##### Create detection history list for Sac River model  -----------------------------------------------------------------------------------------
name <- "SacRiverAllYears"
#path <- paste0("./Outputs/", name, "/")

path <- paste0(here::here("data-raw", "survival_model_data"), "/", name, "/")
dir.create(path, showWarnings = FALSE)

studyIDs <- c("ColemanFall_2013","ColemanFall_2016","ColemanFall_2017",
              "CNFH_FMR_2019","CNFH_FMR_2020","CNFH_FMR_2021",
              "RBDD_2017","RBDD_2018",
              "DeerCk_Wild_CHK_2018","DeerCk_Wild_CHK_2020",
              "MillCk_Wild_CHK_2013","MillCk_Wild_CHK_2015","MillCk_Wild_CHK_2017")

## Retrieve ERDDAP detection data saved as csv files
#names <- lapply(studyIDs, function(x) paste0("./Outputs/",name, "/", x, ".csv"))
names <- lapply(studyIDs, function(x) paste0(here::here("data-raw", "survival_model_data"),
                                             "/", name, "/", x, ".csv"))
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Manually select receiver locations to use and combine for Sac River study
reach.meta <- reach.meta %>%
  filter(GEN %in% c("BattleCk_CNFH_Rel","RBDD_Rel","RBDD_Rel_Rec","Altube Island","MillCk_RST_Rel",
                    "MillCk2_Rel","DeerCk_RST_Rel","Mill_Ck_Conf",
                    "Abv_WoodsonBr","Blw_Woodson",
                    "I80-50_Br","TowerBridge",
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>%
  mutate(Region = case_when(Region == 'Battle Ck' ~ 'Release',
                            Region == 'DeerCk' ~ 'Release',
                            Region == 'Mill Ck' ~ 'Release',
                            GEN == 'RBDD_Rel'& Region == 'Upper Sac R' ~ 'Release',
                            GEN == 'RBDD_Rel_Rec'& Region == 'Upper Sac R' ~ 'Release',
                            Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))

# Aggregate receiver locations and detections
all_aggregated <- lapply(all_detections, aggregate_GEN_Sac)

# View study reaches in map
reach.map_agg <- reach.meta.aggregate %>%
                filter(GEN != "Releasepoint") %>%
                 rbind(c("Battle Creek",517.344,40.39816,-122.1456,"Release"),
                       c("RBDD",461.579, 40.15444, -122.2025,'Release'),
                       c("Mill Creek",450.703,40.05479, -122.0321,'Release'),
                       c('Deer Creek',441.728,39.99740,-121.9677,'Release')) %>%
                mutate(GenLat = as.numeric(GenLat),
                       GenLon = as.numeric(GenLon))

(map_receiv_agg <- leaflet(data = reach.map_agg) %>% addTiles() %>%
    addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
               label = ~as.character(GEN),
               labelOptions = labelOptions(noHide = T, textOnly = TRUE)) %>%
    addProviderTiles("Esri.WorldTopoMap")
)

mapshot(map_receiv_agg, file = paste0(getwd(), "/Figures/mapreceiv_Sac.png"),
        width = 20, height = 36)

# Create Encounter History list and inp file for Sac River model------------------------------------------------------------------
all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# Add in fish information to inp file
# First add in fish info
all.inp <- all.inp %>%
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date,
                                  release_location),by = c("FishID" = "fish_id"))

# Rename release location
all.inp$rl <- all.inp$release_location

# add year
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp, paste0(here::here("data-raw", "survival_model_data"), "/SacInp.csv"))
#write.csv(all.inp,"./Outputs/SacInp.csv", row.names=FALSE)

# Summarize fish info
fish_summary <- all.inp %>%
                group_by(year) %>%
                summarise(minFL = min(as.numeric(fish_length),na.rm=TRUE),
                          maxFL = max(as.numeric(fish_length),na.rm=TRUE),
                          minweight = min(as.numeric(fish_weight),na.rm=TRUE),
                          maxweight = max(as.numeric(fish_weight),na.rm=TRUE),
                          N=n())

# Run RMark Sac River model ----------------------------------------------------------
Sac.inp <- all.inp

## load the CH file. This is where you specify groups, covariates, etc. time.intervals=reach_length,groups="StudyID")
Sac.process <- process.data(Sac.inp, model="CJS", begin.time=1, groups=c("year"))
Sac.ddl <- make.design.data(Sac.process)

## Find what is the best p and phi model
Phi.dot <- list(formula= ~1)
p.dot <- list(formula= ~1)
Phi.t <- list(formula= ~time)
p.t <- list(formula= ~time)

Phi.t.x.y <- list(formula= ~time* year)
Phi.t.plus.y <- list(formula= ~time+ year)
p.t.x.y <- list(formula= ~time* year)
p.t.plus.y <- list(formula= ~time + year)

cml = create.model.list("CJS")

## Run mark.wrapper for all model structures (the combination of all phi and p models possible).
model.outputs <- mark.wrapper(cml, data=Sac.process, ddl=Sac.ddl) #,adjust=FALSE
model_table <- model.outputs$model.table

write.csv(model_table, paste0(here::here("data-raw", "survival_model_data"), "/SacSurvmodel_AICtable.csv"),
          row.names = FALSE)
#write.csv(model_table,"Outputs/SacSurvmodel_AICtable.csv",row.names = FALSE)

lfc.Phi.t.x.y.p.t.plus.y <- mark(Sac.process, Sac.ddl, model.parameters=list(Phi=Phi.t.x.y, p=p.t.plus.y),
                              realvcv = TRUE)
# NOTE: these will not run unless you have mark v9 installed. current version is v10.
# when installing previous versions of mark, you will have to make sure you have the
# right version of gcc installed as well

# Reach-specific Survival estimates for Sacramento River model ------------------------------------------------------------
Phi.t.x.y.p.t.plus.y.means <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$estimate[1:24],3)
Phi.t.x.y.p.t.plus.y.se <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$se[1:24],3)
Phi.t.x.y.p.t.plus.y.lcl <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$lcl[1:24],3)
Phi.t.x.y.p.t.plus.y.ucl <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$ucl[1:24],3)
Phi.t.x.y.p.t.plus.y.vcv <- lfc.Phi.t.x.y.p.t.plus.y$results$real.vcv[1:24,1:24]

relwood_surv <- lfc.Phi.t.x.y.p.t.plus.y$results$real[c(1,4,7,10,13,16,19,22),] %>%
  mutate(year = as.factor(c(2013,2015,2016,2017,2018,2019,2020,2021)))

woodsac_surv <- lfc.Phi.t.x.y.p.t.plus.y$results$real[c(2,5,8,11,14,17,20,23),] %>%
  mutate(year = as.factor(c(2013,2015,2016,2017,2018,2019,2020,2021)))
write.csv(relwood_surv, paste0(here::here("data-raw", "survival_model_data"), "/relwood_surv.csv"),
          row.names = FALSE)
write.csv(woodsac_surv, paste0(here::here("data-raw", "survival_model_data"), "/woodsac_surv.csv"),
          row.names = FALSE)
# write.csv(relwood_surv,"Outputs/relwood_surv.csv",row.names = FALSE)
# write.csv(woodsac_surv,"Outputs/woodsac_surv.csv",row.names = FALSE)

# Per 10km survival for Sacramento River model -----------------------------------------------------------
# Assign new rounded rkm values that are not average of several receiver locations for reach length calculation
# we used Battle Creek rkm for release point which is the most upstream release location and Tower Bridge for Sacramento location
reach.meta.aggregate <- reach.meta.aggregate %>%
  mutate(GenRKM = case_when(GEN =='Releasepoint' ~ ceiling(517.344),
                            GEN == "WoodsonBridge" ~ ceiling(429.292),
                            GEN == "Sacramento" ~ ceiling(171.374),
                            TRUE ~ ceiling(107.512)))

# Calculate the reach lengths. Here I divided reach lengths by 10 so that my survival estimates later will be survival per 10km
KM <- reach.meta.aggregate[order(reach.meta.aggregate$GenRKM, decreasing= TRUE),2] %>%
  data.frame()
KM <- KM[,1]

reach_length <- abs(diff(KM))/10

# set up the basic model structure. Here we are using the CJS model for live recaptures
# we are setting time intervals to reach lenth, which means our reach survival estimates will be on a per 10km scale
Sac.process10km <- process.data(Sac.inp, model="CJS", begin.time=1,time.intervals=reach_length,groups="year")
Sac.ddl10km <- make.design.data(Sac.process10km)

# pull out the survival results from the best model
lfc.Phi.t.x.y.p.t.plus.y.per.10km <- mark(Sac.process10km, Sac.ddl10km, model.parameters=list(Phi=Phi.t.x.y, p=p.t.plus.y),
                                       realvcv = TRUE)

Phi.t.x.y.p.t.plus.y.per.10km.means <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$estimate[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.se <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$se[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.lcl <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$lcl[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.ucl <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$ucl[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.vcv <- lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real.vcv[1:24,1:24]

# Cumulative survival estimates for Sacramento River model -----------------------------------------------------------------------------
# Reach selection
reach_vec <- c("Release to WoodsonBr","WoodsonBr to Sacramento")
reach_numb <-length(reach_vec)

###### 2013
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2013 <- cumprod(Phi.t.x.y.p.t.plus.y.means[1:2])

# calculate standard errors for the cumulative product.
cum.phi.se_2013 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[1:2],
                                       Phi.t.x.y.p.t.plus.y.vcv[1:2,
                                                             1:2])


# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2013 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2013[,1] <- cum.phi_2013
Cumul_all_2013[,2] <- cum.phi.se_2013

# Also calculate the LCI for UCI for all the reach estimates
LCI_2013 <- expit(logit(cum.phi_2013)-1.96*sqrt(cum.phi.se_2013^2/((exp(logit(cum.phi_2013))/
                                                                      (1+exp(logit(cum.phi_2013)))^2)^2)))
UCI_2013 <- expit(logit(cum.phi_2013)+1.96*sqrt(cum.phi.se_2013^2/((exp(logit(cum.phi_2013))/
                                                                      (1+exp(logit(cum.phi_2013)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2013[,3] <- LCI_2013
Cumul_all_2013[,4] <- UCI_2013

# Now round both matrices to 3 decimal places
Cumul_all_2013 <- round(Cumul_all_2013,4)

###### 2015
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2015 <- cumprod(Phi.t.x.y.p.t.plus.y.means[4:5])

# calculate standard errors for the cumulative product.
cum.phi.se_2015 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[4:5],
                                       Phi.t.x.y.p.t.plus.y.vcv[4:5,
                                                             4:5])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2015 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2015[,1] <- cum.phi_2015
Cumul_all_2015[,2] <- cum.phi.se_2015

# Also calculate the LCI for UCI for all the reach estimates
LCI_2015 <- expit(logit(cum.phi_2015)-1.96*sqrt(cum.phi.se_2015^2/((exp(logit(cum.phi_2015))/
                                                                      (1+exp(logit(cum.phi_2015)))^2)^2)))
UCI_2015 <- expit(logit(cum.phi_2015)+1.96*sqrt(cum.phi.se_2015^2/((exp(logit(cum.phi_2015))/
                                                                      (1+exp(logit(cum.phi_2015)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2015[,3] <- LCI_2015
Cumul_all_2015[,4] <- UCI_2015

# Now round both matrices to 3 decimal places
Cumul_all_2015 <- round(Cumul_all_2015,4)

###### 2016
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2016 <- cumprod(Phi.t.x.y.p.t.plus.y.means[7:8])

# calculate standard errors for the cumulative product.
cum.phi.se_2016 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[7:8],
                                       Phi.t.x.y.p.t.plus.y.vcv[7:8,
                                                             7:8])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2016 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2016[,1] <- cum.phi_2016
Cumul_all_2016[,2] <- cum.phi.se_2016

# Also calculate the LCI for UCI for all the reach estimates
LCI_2016 <- expit(logit(cum.phi_2016)-1.96*sqrt(cum.phi.se_2016^2/((exp(logit(cum.phi_2016))/
                                                                      (1+exp(logit(cum.phi_2016)))^2)^2)))
UCI_2016 <- expit(logit(cum.phi_2016)+1.96*sqrt(cum.phi.se_2016^2/((exp(logit(cum.phi_2016))/
                                                                      (1+exp(logit(cum.phi_2016)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2016[,3] <- LCI_2016
Cumul_all_2016[,4] <- UCI_2016

# Now round both matrices to 3 decimal places
Cumul_all_2016 <- round(Cumul_all_2016,4)

###### 2017
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2017 <- cumprod(Phi.t.x.y.p.t.plus.y.means[10:11])

# calculate standard errors for the cumulative product.
cum.phi.se_2017 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[10:11],
                                       Phi.t.x.y.p.t.plus.y.vcv[10:11,
                                                             10:11])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2017 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2017[,1] <- cum.phi_2017
Cumul_all_2017[,2] <- cum.phi.se_2017

# Also calculate the LCI for UCI for all the reach estimates
LCI_2017 <- expit(logit(cum.phi_2017)-1.96*sqrt(cum.phi.se_2017^2/((exp(logit(cum.phi_2017))/
                                                                      (1+exp(logit(cum.phi_2017)))^2)^2)))
UCI_2017 <- expit(logit(cum.phi_2017)+1.96*sqrt(cum.phi.se_2017^2/((exp(logit(cum.phi_2017))/
                                                                      (1+exp(logit(cum.phi_2017)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2017[,3] <- LCI_2017
Cumul_all_2017[,4] <- UCI_2017

# Now round both matrices to 3 decimal places
Cumul_all_2017 <- round(Cumul_all_2017,4)

###### 2018
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2018 <- cumprod(Phi.t.x.y.p.t.plus.y.means[13:14])

# calculate standard errors for the cumulative product.
cum.phi.se_2018 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[13:14],
                                       Phi.t.x.y.p.t.plus.y.vcv[13:14,
                                                             13:14])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2018 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2018[,1] <- cum.phi_2018
Cumul_all_2018[,2] <- cum.phi.se_2018

# Also calculate the LCI for UCI for all the reach estimates
LCI_2018 <- expit(logit(cum.phi_2018)-1.96*sqrt(cum.phi.se_2018^2/((exp(logit(cum.phi_2018))/
                                                                      (1+exp(logit(cum.phi_2018)))^2)^2)))
UCI_2018 <- expit(logit(cum.phi_2018)+1.96*sqrt(cum.phi.se_2018^2/((exp(logit(cum.phi_2018))/
                                                                      (1+exp(logit(cum.phi_2018)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2018[,3] <- LCI_2018
Cumul_all_2018[,4] <- UCI_2018

# Now round both matrices to 3 decimal places
Cumul_all_2018 <- round(Cumul_all_2018,4)

###### 2019
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2019 <- cumprod(Phi.t.x.y.p.t.plus.y.means[16:17])

# calculate standard errors for the cumulative product.
cum.phi.se_2019 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[16:17],
                                       Phi.t.x.y.p.t.plus.y.vcv[16:17,
                                                             16:17])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2019 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2019[,1] <- cum.phi_2019
Cumul_all_2019[,2] <- cum.phi.se_2019

# Also calculate the LCI for UCI for all the reach estimates
LCI_2019 <- expit(logit(cum.phi_2019)-1.96*sqrt(cum.phi.se_2019^2/((exp(logit(cum.phi_2019))/
                                                                      (1+exp(logit(cum.phi_2019)))^2)^2)))
UCI_2019 <- expit(logit(cum.phi_2019)+1.96*sqrt(cum.phi.se_2019^2/((exp(logit(cum.phi_2019))/
                                                                      (1+exp(logit(cum.phi_2019)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2019[,3] <- LCI_2019
Cumul_all_2019[,4] <- UCI_2019

# Now round both matrices to 3 decimal places
Cumul_all_2019 <- round(Cumul_all_2019,4)

###### 2020
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2020 <- cumprod(Phi.t.x.y.p.t.plus.y.means[19:20])

# calculate standard errors for the cumulative product.
cum.phi.se_2020 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[19:20],
                                       Phi.t.x.y.p.t.plus.y.vcv[19:20,
                                                             19:20])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2020 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2020[,1] <- cum.phi_2020
Cumul_all_2020[,2] <- cum.phi.se_2020

# Also calculate the LCI for UCI for all the reach estimates
LCI_2020<- expit(logit(cum.phi_2020)-1.96*sqrt(cum.phi.se_2020^2/((exp(logit(cum.phi_2020))/
                                                                     (1+exp(logit(cum.phi_2020)))^2)^2)))
UCI_2020 <- expit(logit(cum.phi_2020)+1.96*sqrt(cum.phi.se_2020^2/((exp(logit(cum.phi_2020))/
                                                                      (1+exp(logit(cum.phi_2020)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2020[,3] <- LCI_2020
Cumul_all_2020[,4] <- UCI_2020

# Now round both matrices to 3 decimal places
Cumul_all_2020 <- round(Cumul_all_2020,4)

###### 2021
## Let's run delta method
# Calculate the cumulative survival, do not use last reach as this is just for better detection efficiency estimation upstream.
cum.phi_2021 <- cumprod(Phi.t.x.y.p.t.plus.y.means[22:23])

# calculate standard errors for the cumulative product.
cum.phi.se_2021 <- deltamethod.special("cumprod",Phi.t.x.y.p.t.plus.y.means[22:23],
                                       Phi.t.x.y.p.t.plus.y.vcv[22:23,
                                                             22:23])

# Now make another matrix to fill in with cumprod estimates of survival and propagated SE to all reaches
Cumul_all_2021 <- matrix(0,reach_numb,4, dimnames = list(c(1:2),c("Phi","SE","LCI","UCI")))

# Also put all the phi and SE cumul estimates into the other matrix
Cumul_all_2021[,1] <- cum.phi_2021
Cumul_all_2021[,2] <- cum.phi.se_2021

# Also calculate the LCI for UCI for all the reach estimates
LCI_2021 <- expit(logit(cum.phi_2021)-1.96*sqrt(cum.phi.se_2021^2/((exp(logit(cum.phi_2021))/
                                                                      (1+exp(logit(cum.phi_2021)))^2)^2)))
UCI_2021 <- expit(logit(cum.phi_2021)+1.96*sqrt(cum.phi.se_2021^2/((exp(logit(cum.phi_2021))/
                                                                      (1+exp(logit(cum.phi_2021)))^2)^2)))

# Also plug in the LCI and UCI into the all reach matrix
Cumul_all_2021[,3] <- LCI_2021
Cumul_all_2021[,4] <- UCI_2021

# Now round both matrices to 3 decimal places
Cumul_all_2021 <- round(Cumul_all_2021,4)

########## All years combined
Cumul_all_years <- rbind(Cumul_all_2013,Cumul_all_2015,Cumul_all_2016,Cumul_all_2017,Cumul_all_2018,
                         Cumul_all_2019,Cumul_all_2020,Cumul_all_2021)

Cumul_Sac <- Cumul_all_years[c(2,4,6,8,10,12,14,16),]

#################### Feather River model ##############################################

##### Create detection history list for Feather River model -----------------------------------------------------------------------------------------
name <- "FeatherRiverAllYears"
path <- paste0(here::here("data-raw", "survival_model_data"), "/", name, "/")
dir.create(path, showWarnings = FALSE)

studyIDs <- c("FR_Spring_2013","FR_Spring_2015","FR_Spring_2019","FR_Spring_2020","FR_Spring_2021")

## Retrieve ERDDAP detection data saved as csv files
names <- lapply(studyIDs, function(x) paste0(here::here("data-raw", "survival_model_data"),                                              "/", name, "/", x, ".csv"))
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Manually select receiver locations to use and combine for Sac River study
reach.meta <- reach.meta %>%
  filter(GEN %in% c("FR_Gridley_Rel","FR_Boyds_Rel","FR_Boyds_Rel_Rec",
                    "I80-50_Br","TowerBridge",
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>%
  mutate(Region = case_when(Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))


# Aggregate receiver locations and detections
all_aggregated <- lapply(all_detections, aggregate_GEN_Feather)

# View study reaches in map
reach.map_agg <- reach.meta.aggregate %>%
  filter(GEN != "Releasepoint") %>%
  rbind(c("Feather_Gridley",287.387,39.35788,-121.6360,"Feather_R"),
        c("Feather_Boyds",240.755,39.05734, -121.6107,"Feather_R")) %>%
  mutate(GenLat = as.numeric(GenLat),
         GenLon = as.numeric(GenLon))

(map_receiv_agg <- leaflet(data = reach.map_agg) %>% addTiles() %>%
    addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
               label = ~as.character(GEN),
               labelOptions = labelOptions(noHide = T, textOnly = TRUE)) %>%
    addProviderTiles("Esri.WorldTopoMap")
)

mapshot(map_receiv_agg, file = paste0(getwd(), "/Figures/mapreceiv_Feather.png"))

# Create Encounter History list and inp file for Feather River model------------------------------------------------------------------
all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# Add in fish information to inp file
# First add in fish info
all.inp <- all.inp %>%
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date,
                                  release_location),by = c("FishID" = "fish_id"))

# Rename release location
all.inp$rl <- all.inp$release_location

# add year
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))
write.csv(all.inp, paste0(here::here("data-raw", "survival_model_data"), "/FeatherInp.csv"),
          row.names = FALSE)
#write.csv(all.inp,"./Outputs/FeatherInp.csv", row.names=FALSE)

# Summarize fish info
fish_summary <- all.inp %>%
  group_by(year) %>%
  summarise(minFL = min(as.numeric(fish_length),na.rm=TRUE),
            maxFL = max(as.numeric(fish_length),na.rm=TRUE),
            minweight = min(as.numeric(fish_weight),na.rm=TRUE),
            maxweight = max(as.numeric(fish_weight),na.rm=TRUE),
            N=n())

# Run RMark Feather River model ----------------------------------------------------------
Feather.inp <- all.inp

#remove previously saved lists
rm(list = ls()[grep('Phi.', ls())])
rm(list = ls()[grep('p.', ls())])

## load the CH file. This is where you specify groups, covariates, etc. time.intervals=reach_length,groups="StudyID")
Feather.process <- process.data(Feather.inp, model="CJS", begin.time=1, groups=c("year"))
Feather.ddl <- make.design.data(Feather.process)

## Find what is the best p and phi model
Phi.dot <- list(formula= ~1)
p.dot <- list(formula= ~1)
Phi.t <- list(formula= ~time)
p.t <- list(formula= ~time)

Phi.t.x.y <- list(formula= ~time* year)
Phi.t.plus.y <- list(formula= ~time+ year)
p.t.x.y <- list(formula= ~time* year)
p.t.plus.y <- list(formula= ~time + year)

cml = create.model.list("CJS")

## Run mark.wrapper for all model structures (the combination of all phi and p models possible).
model.outputs <- mark.wrapper(cml, data=Feather.process, ddl=Feather.ddl) #,adjust=FALSE
model_table <- model.outputs$model.table

write.csv(model_table, paste0(here::here("data-raw", "survival_model_data"),
                              "/FeatherSurvmodel_AICtable.csv"),
          row.names = FALSE)
#write.csv(model_table,"Outputs/FeatherSurvmodel_AICtable.csv",row.names = FALSE)

lfc.Phi.t.x.y.p.dot <- mark(Feather.process,Feather.ddl, model.parameters=list(Phi=Phi.t.x.y, p=p.dot),
                                 realvcv = TRUE)

# Reach-specific Survival estimates for Feather River model ------------------------------------------------------------
lfc.Phi.t.x.y.p.dot.means <- round(lfc.Phi.t.x.y.p.dot$results$real$estimate[1:10],3)
lfc.Phi.t.x.y.p.dot.se <- round(lfc.Phi.t.x.y.p.dot$results$real$se[1:10],3)
lfc.Phi.t.x.y.p.dot.lcl <- round(lfc.Phi.t.x.y.p.dot$results$real$lcl[1:10],3)
lfc.Phi.t.x.y.p.dot.ucl <- round(lfc.Phi.t.x.y.p.dot$results$real$ucl[1:10],3)

feather_surv <- lfc.Phi.t.x.y.p.dot$results$real[c(1,3,5,7,9),]
feather_surv <- feather_surv %>%
  mutate(year = as.factor(c(2013,2015,2019,2020,2021)))

write.csv(feather_surv, paste0(here::here("data-raw", "survival_model_data"),
                               "/feather_surv.csv"),
          row.names = FALSE)
#write.csv(feather_surv,"Outputs/feather_surv.csv",row.names = FALSE)

# Per 10km survival for Sacramento River model -----------------------------------------------------------
# Assign new rounded rkm values that are not average of several receiver locations for reach length calculation
# we used Battle Creek rkm for release point which is the most upstream release location and Tower Bridge for Sacramento location
reach.meta.aggregate <- reach.meta.aggregate %>%
  mutate(GenRKM = case_when(GEN =='Releasepoint' ~ ceiling(517.344),
                            GEN == "WoodsonBridge" ~ ceiling(429.292),
                            GEN == "Sacramento" ~ ceiling(171.374),
                            TRUE ~ ceiling(107.512)))

# Calculate the reach lengths. Here I divided reach lengths by 10 so that my survival estimates later will be survival per 10km
KM <- reach.meta.aggregate[order(reach.meta.aggregate$GenRKM, decreasing= TRUE),2] %>%
  data.frame()
KM <- KM[,1]

reach_length <- abs(diff(KM))/10

# TODO problem with time.intervals = reach_length
# set up the basic model structure. Here we are using the CJS model for live recaptures
# we are setting time intervals to reach lenth, which means our reach survival estimates will be on a per 10km scale
Sac.process10km <- process.data(Sac.inp, model="CJS", begin.time=1,time.intervals=reach_length,groups="year")
Sac.ddl10km <- make.design.data(Sac.process10km)

# pull out the survival results from the best model
lfc.Phi.t.x.y.p.t.plus.y.per.10km <- mark(Sac.process10km, Sac.ddl10km, model.parameters=list(Phi=Phi.t.x.y, p=p.t.plus.y),
                                          realvcv = TRUE)

Phi.t.x.y.p.t.plus.y.per.10km.means <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$estimate[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.se <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$se[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.lcl <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$lcl[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.ucl <- round(lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real$ucl[1:24],3)
Phi.t.x.y.p.t.plus.y.per.10km.vcv <- lfc.Phi.t.x.y.p.t.plus.y.per.10km$results$real.vcv[1:24,1:24]


#################### Butte Creek model ##############################################

##### Create detection history list for Butte Creek model  -----------------------------------------------------------------------------------------
name <- "ButteCreekAllYears"
path <- paste0(here::here("data-raw", "survival_model_data"), "/", name, "/")
dir.create(path, showWarnings = FALSE)

studyIDs <- c("SB_Spring_2015","SB_Spring_2016","SB_Spring_2017","SB_Spring_2018",
              "SB_Spring_2019","Upper_Butte_2019")


## Retrieve ERDDAP detection data saved as csv files
names <- lapply(studyIDs, function(x) paste0(here::here("data-raw", "survival_model_data"),                                              "/", name, "/", x, ".csv"))
all_detections <- lapply(names, vroom)

# Get list of all receiver GEN
reach.meta <- get_receiver_GEN(all_detections)

# Manually select receiver locations to use and combine for Sac River study
reach.meta <- reach.meta %>%
  filter(GEN %in% c("UpperButte_RST_Rel","UpperButte_RST","UpperButte_SKWY",
                    "SutterBypass_Weir2_RST_Rel","SutterBypass Weir2 RST",
                    "I80-50_Br","TowerBridge",
                    "ToeDrainBase","Hwy84Ferry",
                    "BeniciaE","BeniciaW","ChippsE","ChippsW"))%>%
  mutate(Region = case_when(Region == 'Sutter Bypass' ~ 'Butte Creek',
                            Region == 'Yolo Bypass' ~ 'Lower Sac R',
                            Region == 'North Delta' ~ 'Lower Sac R',
                            Region == 'West Delta' ~ 'End',
                            Region == 'Carquinez Strait' ~ 'End',
                            TRUE ~ Region))


# Aggregate receiver locations and detections
all_aggregated <- lapply(all_detections, aggregate_GEN_Butte)

# View study reaches in map
reach.map_agg <- reach.meta.aggregate %>%
  filter(GEN != "Releasepoint") %>%
  rbind(c("UpperButte", 340.854,39.70936, -121.7506,'Butte Creek'),
        c('SutterBypass',249.540,39.10239, -121.7588,'Butte Creek')) %>%
  mutate(GenLat = as.numeric(GenLat),
         GenLon = as.numeric(GenLon))

(map_receiv_agg <- leaflet(data = reach.map_agg) %>% addTiles() %>%
    addMarkers(~GenLon, ~GenLat, popup = ~as.character(GEN),
               label = ~as.character(GEN),
               labelOptions = labelOptions(noHide = T, textOnly = TRUE)) %>%
    addProviderTiles("Esri.WorldTopoMap")
)

mapshot(map_receiv_agg, file = paste0(getwd(), "/Figures/mapreceiv_Butte.png"))

# Create Encounter History list and inp file for Butte Creek model------------------------------------------------------------------
all_EH <- lapply(all_aggregated, make_EH)

all.inp <- pmap(list(all_detections,all_EH), create_inp) %>%
  bind_rows()

# Add in fish information to inp file
# First add in fish info
all.inp <- all.inp %>%
  left_join(TaggedFish %>% select(fish_id, fish_length, fish_weight, fish_type,fish_release_date,
                                  release_location),by = c("FishID" = "fish_id"))

# add year
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))
all.inp$rl <- all.inp$release_location

write.csv(all.inp, paste0(here::here("data-raw", "survival_model_data"),
                               "/ButteInp.csv"),
          row.names = FALSE)
#write.csv(all.inp,"./Outputs/ButteInp.csv", row.names=FALSE)

# Summarize fish info
fish_summary <- all.inp %>%
  group_by(year) %>%
  summarise(minFL = min(as.numeric(fish_length),na.rm=TRUE),
            maxFL = max(as.numeric(fish_length),na.rm=TRUE),
            minweight = min(as.numeric(fish_weight),na.rm=TRUE),
            maxweight = max(as.numeric(fish_weight),na.rm=TRUE),
            N=n())

# Run RMark Butte Creek model ----------------------------------------------------------
Butte.inp <- all.inp

## load the CH file. This is where you specify groups, covariates, etc. time.intervals=reach_length,groups="StudyID")
Butte.process <- process.data(Butte.inp, model="CJS", begin.time=1, groups=c("year"))
Butte.ddl <- make.design.data(Butte.process)

## Find what is the best p and phi model
Phi.dot <- list(formula= ~1)
p.dot <- list(formula= ~1)
Phi.t <- list(formula= ~time)
p.t <- list(formula= ~time)

Phi.t.x.y <- list(formula= ~time* year)
Phi.t.plus.y <- list(formula= ~time+ year)
p.t.x.y <- list(formula= ~time* year)
p.t.plus.y <- list(formula= ~time + year)

cml = create.model.list("CJS")

## Run mark.wrapper for all model structures (the combination of all phi and p models possible).
model.outputs <- mark.wrapper(cml, data=Butte.process, ddl=Butte.ddl) #,adjust=FALSE
model_table <- model.outputs$model.table

write.csv(model_table, paste0(here::here("data-raw", "survival_model_data"),
                               "/ButteSurvmodel_AICtable.csv"),
          row.names = FALSE)
#write.csv(model_table,"Outputs/ButteSurvmodel_AICtable.csv",row.names = FALSE)

lfc.Phi.t.x.y.p.dot <- mark(Butte.process,Butte.ddl, model.parameters=list(Phi=Phi.t.x.y, p=p.dot),
                                 realvcv = TRUE)

# Reach-specific Survival estimates for Butte Creek model ------------------------------------------------------------
Phi.t.x.y.p.dot.means <- round(lfc.Phi.t.x.y.p.dot$results$real$estimate[1:10],3)
Phi.t.x.y.p.dot.se <- round(lfc.Phi.t.x.y.p.dot$results$real$se[1:10],3)
Phi.t.x.y.p.dot.lcl <- round(lfc.Phi.t.x.y.p.dot$results$real$lcl[1:10],3)
Phi.t.x.y.p.dot.ucl <- round(lfc.Phi.t.x.y.p.dot$results$real$ucl[1:10],3)

butte_surv <- lfc.Phi.t.x.y.p.dot$results$real[c(1,3,5,7,9),] %>%
  mutate(year = as.factor(c(2015,2016,2017,2018,2019)))

write.csv(butte_surv, paste0(here::here("data-raw", "survival_model_data"),
                               "/butte_surv.csv"),
          row.names = FALSE)
#write.csv(butte_surv,"Outputs/butte_surv.csv",row.names = FALSE)

# Figures -----------------------------------------------------------------------
### Plot survival from release to Woodson Bridge and Woodson to Sacramento
(relwoodFig <- ggplot(relwood_surv)+
    geom_errorbar(aes(x=year, ymin=ucl, ymax=lcl),
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=year,y=estimate),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    #  scale_y_continuous(limits=c(0.05,0.35), breaks=seq(0.05,0.35,0.05))+
    xlab("")+ylab("Survival rate")
)

(woodsacFig <- ggplot(woodsac_surv)+
    geom_errorbar(aes(x=year, ymin=ucl, ymax=lcl),
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=year,y=estimate),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+ylab("Survival rate")
)

(SacSurvFig <- ggarrange(relwoodFig,woodsacFig,
                        ncol=2,
                        labels =c('A','B'))
)

ggsave(plot=SacSurvFig, "Figures/ReachSacSurvFig.png", width =10, height = 5)


## Plot cumulative survival from release to Sacramento
# TODO these weren't working
CumulSacdf <- data.frame(Cumul_Sac) %>%
  mutate(year = as.factor(c(2013,2015,2016,2017,2018,2019,2020,2021)))

(CumSacSurv <- ggplot(CumulSacdf)+
    geom_errorbar(aes(x=year, ymin=LCI, ymax=UCI),
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=year,y=Phi),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    #  scale_y_continuous(limits=c(0.05,0.35), breaks=seq(0.05,0.35,0.05))+
    xlab("")+ylab("Survival rate")
)

ggsave(plot=CumSacSurv, "Figures/CumSacSurvFig.png", width = 8, height = 6)


### Plot survival from Feather River release to Sacramento
(feather_surv <- ggplot(feather_surv)+
    geom_errorbar(aes(x=year, ymin=ucl, ymax=lcl),
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=year,y=estimate),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+ylab("Survival rate")
)

ggsave(plot=feather_surv, "Figures/FeatherSurvFig.png", width =5, height = 4)

### Plot survival from Butte Creek release to Sacramento
(butte_surv <- ggplot(butte_surv)+
    geom_errorbar(aes(x=year, ymin=ucl, ymax=lcl),
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=year,y=estimate),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    xlab("")+ylab("Survival rate")
)

ggsave(plot=butte_surv , "Figures/ButteSurvFig.png", width =5, height = 4)
