##%##########################################################################################################%##
#                                                                                                              #
#       Survival Analysis for spring-run Chinook through the Sacramento River, Butte Creek and Feather River   #
#                                                                                                              #  
##%##########################################################################################################%##
# Authors: Jeremy Notch and Tom Pham, modified by Flora Cordoleani
# Script and data info: This script performs survival analysis using CJS model
# and outputs detection probabilities using RMark
# Data source: Data extracted from NOAA ERDDAP data server (oceanview)

source('Script/helper_fxt.R')

library(here)
library(reshape2)
library(ggplot2)
library(readxl)
library(locfit)
library(gtools)
library(mapview)
library(ggpubr)

#################### Sacramento River model ##############################################

##### Create detection history list for Sac River model  -----------------------------------------------------------------------------------------
name <- "SacRiverAllYears"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("ColemanFall_2013","ColemanFall_2016","ColemanFall_2017",
              "CNFH_FMR_2019","CNFH_FMR_2020","CNFH_FMR_2021",
              "RBDD_2017","RBDD_2018",
              "DeerCk_Wild_CHK_2018","DeerCk_Wild_CHK_2020",
              "MillCk_Wild_CHK_2013","MillCk_Wild_CHK_2015","MillCk_Wild_CHK_2017")

## Retrieve ERDDAP detection data saved as csv files
names <- lapply(studyIDs, function(x) paste0("./Outputs/",name, "/", x, ".csv"))
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

# add year 
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp,"./Outputs/SacInp.csv", row.names=FALSE)

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
Sac.process <- process.data(Sac.inp, model="CJS", begin.time=1, groups="year")
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

write.csv(model_table,"Outputs/SacSurvmodel_AICtable.csv",row.names = FALSE)

lfc.Phi.t.x.y.p.t.plus.y <- mark(Sac.process, Sac.ddl, model.parameters=list(Phi=Phi.t.x.y, p=p.t.plus.y),
                              realvcv = TRUE) 

# Reach-specific Survival estimates for Sac River model ------------------------------------------------------------
Phi.t.x.y.p.t.x.y.means <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$estimate[1:24],3)
Phi.t.x.y.p.t.x.y.se <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$se[1:24],3)
Phi.t.x.y.p.t.x.y.lcl <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$lcl[1:24],3)
Phi.t.x.y.p.t.x.y.ucl <- round(lfc.Phi.t.x.y.p.t.plus.y$results$real$ucl[1:24],3)

relwood_surv <- lfc.Phi.t.x.y.p.t.plus.y$results$real[c(1,4,7,10,13,16,19,22),] %>% 
  mutate(year = as.factor(c(2013,2015,2016,2017,2018,2019,2020,2021)))

woodsac_surv <- lfc.Phi.t.x.y.p.t.plus.y$results$real[c(2,5,8,11,14,17,20,23),] %>% 
  mutate(year = as.factor(c(2013,2015,2016,2017,2018,2019,2020,2021)))

write.csv(relwood_surv,"Outputs/relwood_surv.csv",row.names = FALSE)
write.csv(woodsac_surv,"Outputs/woodsac_surv.csv",row.names = FALSE)


#################### Feather River model ##############################################

##### Create detection history list  for Feather River model -----------------------------------------------------------------------------------------
name <- "FeatherRiverAllYears"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("FR_Spring_2013","FR_Spring_2015","FR_Spring_2019","FR_Spring_2020","FR_Spring_2021") 

## Retrieve ERDDAP detection data saved as csv files
names <- lapply(studyIDs, function(x) paste0("./Outputs/",name, "/", x, ".csv"))
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

# add year 
all.inp$year <- as.factor(year(as.Date(all.inp$fish_release_date,format="%m/%d/%Y")))

write.csv(all.inp,"./Outputs/FeatherInp.csv", row.names=FALSE)

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

## load the CH file. This is where you specify groups, covariates, etc. time.intervals=reach_length,groups="StudyID")
Feather.process <- process.data(Feather.inp, model="CJS", begin.time=1, groups="year")
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

write.csv(model_table,"Outputs/FeatherSurvmodel_AICtable.csv",row.names = FALSE)

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

write.csv(feather_surv,"Outputs/feather_surv.csv",row.names = FALSE)


#################### Butte Creek model ##############################################

##### Create detection history list for Butte Creek model  -----------------------------------------------------------------------------------------
name <- "ButteCreekAllYears"
path <- paste0("./Outputs/", name, "/")
dir.create(path, showWarnings = FALSE) 

studyIDs <- c("SB_Spring_2015","SB_Spring_2016","SB_Spring_2017","SB_Spring_2018",
              "SB_Spring_2019","Upper_Butte_2019")


## Retrieve ERDDAP detection data saved as csv files
names <- lapply(studyIDs, function(x) paste0("./Outputs/",name, "/", x, ".csv"))
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

write.csv(all.inp,"./Outputs/ButteInp.csv", row.names=FALSE)

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
Butte.process <- process.data(Butte.inp, model="CJS", begin.time=1, groups="year")
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

write.csv(model_table,"Outputs/ButteSurvmodel_AICtable.csv",row.names = FALSE)

lfc.Phi.t.x.y.p.dot <- mark(Butte.process,Butte.ddl, model.parameters=list(Phi=Phi.t.x.y, p=p.dot),
                                 realvcv = TRUE) 

# Reach-specific Survival estimates for Butte Creek model ------------------------------------------------------------
Phi.t.x.y.p.dot.means <- round(lfc.Phi.t.x.y.p.dot$results$real$estimate[1:10],3)
Phi.t.x.y.p.dot.se <- round(lfc.Phi.t.x.y.p.dot$results$real$se[1:10],3)
Phi.t.x.y.p.dot.lcl <- round(lfc.Phi.t.x.y.p.dot$results$real$lcl[1:10],3)
Phi.t.x.y.p.dot.ucl <- round(lfc.Phi.t.x.y.p.dot$results$real$ucl[1:10],3)

butte_surv <- lfc.Phi.t.x.y.p.dot$results$real[c(1,3,5,7,9),] %>% 
  mutate(year = as.factor(c(2015,2016,2017,2018,2019)))

write.csv(butte_surv,"Outputs/butte_surv.csv",row.names = FALSE)

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

ggsave(plot=SacSurvFig, "Figures/SacSurvFig.png", width =10, height = 5)

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
