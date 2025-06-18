library(tidyverse)
library(DBI)
library(SRJPEdata)

# args
stream <- "battle creek"
# get connection
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)


# TODO our dates are slightly off (3 days?) from josh's weeks
# calwk <- SRJPEmodel::julian_week_to_date_lookup
#ddate=read.table(file="../../RST/Data/Jwk_Dates.txt",header=T) #if you want to see calendar date and model index for each week
# calwk=c('Sep-07','Sep-14','Sep-21','Sep-28','Oct-05','Oct-12','Oct-19','Oct-26','Nov-02','Nov-09',
#         'Nov-16','Nov-23','Nov-30','Dec-07','Dec-14','Dec-21','Dec-28','Jan-04','Jan-11', 'Jan-18',
#         'Jan-25','Feb-01','Feb-08','Feb-15','Feb-22','Mar-01','Mar-08','Mar-15','Mar-22','Mar-29',
#         'Apr-05','Apr-12','Apr-19','Apr-26','May-03','May-10','May-17','May-24','May-31','Jun-07',
#         'Jun-14','Jun-21','Jun-28','Jul-05','Jul-12','Jul-19','Jul-26','Aug-02','Aug-09','Aug-16',
#         'Aug-23','Aug-30','Aug-31')

#Define period that will be modelled
# fwk=36;lwk=35 #All year Sep-01 - Aug-31
# mwk=c(fwk:53,1:lwk)#the calendar week associated with each modelled wk
# Nwks=length(mwk) #the number of modelled weeks

# set up weeks processing
number_weeks <- nrow(SRJPEmodel::julian_week_to_date_lookup)
weeks_ordered <- c(36:53, 1:35)
theta <- (1:number_weeks) / number_weeks # set proportional week (proportion of total year through week iwk)
theta[number_weeks] <- 0.999 # set last value to 0.999 instead of 1

#theta=(1:Nwks)/Nwks;theta[Nwks]=0.999#the proportional week (proportion of total year through week iwk)

# if(RunMainstem==F){
#   RSTpath="../../RST/RunSize/StanVersion/Output_Trib_SR"
#   RSTpathA="../../RST/RunSize/StanVersion/Output_Trib"
# } else {
#   RSTpath="../../RST/RunSize/StanVersion/Output_Main_SR"
#   RSTpathA="../../RST/RunSize/StanVersion/Output_Main"
# }

# set up data processing of RSTs
# TODO this right now is pulling lbc AND ubc
juv_results_from_db <- tbl(con, "model_parameters") |>
  left_join(tbl(con, "trap_location") |>
              select(location_id = id, site, stream),
            by = "location_id") |>
  left_join(tbl(con, "parameter") |>
              select(parameter_id = id, parameter = definition),
            by = "parameter_id") |>
  left_join(tbl(con, "statistic") |>
              select(statistic_id = id, statistic = definition),
            by = "statistic_id") |>
  filter(stream == !!stream,
         parameter %in% c("Ntot", "N"),
         statistic %in% c("mean", "sd")) |>
  collect() |>
  group_by(updated_at, year) |>
  mutate(recent = cur_group_id()) |>
  ungroup() |>
  group_by(year) |>
  mutate(most_recent_by_year = max(recent)) |>
  ungroup() |>
  # filter to only keep the most recent run for each year
  filter(recent == most_recent_by_year) |>
  select(-c(id, location_id, parameter_id, statistic_id,
            updated_at, model_run_id, location_fit_id,
            recent, most_recent_by_year)) |>
  distinct_all() # TODO this is an error in the model parameter upload. we are joining on location_id but there are multiple matches for knights landing, so we are getting duplicates of all the parameter esitimates.

number_years <- length(unique(juv_results_from_db$year)) # number of years fit for that stream
estimated_weeks_matrix <- matrix(nrow = number_weeks, ncol = number_years, 0) # number of weeks fit within a year
proportional_weeks_matrix <- matrix(nrow = number_weeks, ncol = number_years, 0) # proportional week for each week-year combination. Varies by year based on when they sampled at RST

annual_abundance <- juv_results_from_db |>
  filter(parameter == "Ntot") |>
  pivot_wider(names_from = "statistic",
              values_from = "value") |>
  mutate(cv = sd / mean,
         log_N = log(mean))

weekly_abundance <- juv_results_from_db |>
  filter(parameter == "N") |>
  pivot_wider(names_from = "statistic",
              values_from = "value") |>
  mutate(cv = sd / mean) |>
  arrange(year, week_fit)


get_input_variables <- function(site, run_year) {

  inputs <- prepare_abundance_inputs(site, run_year, effort_adjust = T)

  # TODO include "trim run" functionality here, which strips out some of initial and last weeks with 0 catch
  # if(TrimRun==T){
  #   wkbuff=3
  #   fuwk=min(which(u>0))-wkbuff;  if(fuwk<=0) fuwk=1
  #   luwk=max(which(u>0))+wkbuff;if(luwk>Nstrata) luwk=Nstrata
  #   Nestwks[iyr]=luwk-fuwk+1
  # } else { #Use all weeks where abundance was estimated
  #   fuwk=1
  #   Nestwks[iyr]=Nstrata
  # }

  return(list(max_flow = max(inputs$catch_flow_raw, na.rm = T),
              catch = inputs$inputs$data$u,
              nstrata = inputs$inputs$data$Nstrata,
              weeks_fit = inputs$weeks_fit))
}

inputs <- purrr::pmap(list(site = "ubc",
                           run_year = annual_abundance$year),
                      get_input_variables)


for(iyr in 1:Nyrs){

  Jwk=abundance_inputs$weeks_fit
  Flow=abundance_inputs$catch_flow_raw
  Nstrata=abundance_inputs$inputs$data$Nstrata #u is in here (now called dd$u), as is Nstrata
  u=abundance_inputs$inputs$data$u

  TrimRun=F

  #strip out some of initial and last weeks with 0 catch


  for(iwk in 1:Nestwks[iyr]){
    ewk[iwk,iyr]=which(mwk==Jwk[iwk+fuwk-1])
    pwk[iwk,iyr]=theta[ewk[iwk,iyr]]

    icols=fuwk:(iwk+fuwk-1)#the

    #N=abundance$sims.list$N[,icols]
    N=abundance[,icols]#Assumes weekly estimates come first followed by Ntot

    if(iwk==1){
      cumN=N
    } else{
      cumN=rowSums(N) #abundance from start of sampling through week iwk
    }
    mu=mean(cumN)#the mean across posterior samples
    if(mu==0)  mu=.01 #can be 0 in cases with multiple weeks of zero catch at start of sampling
    Nx_mu[iwk,iyr]=log(mu)#total abundance through this week in log space
    sdx=sd(cumN)#sd across posterior samples of cummulative abundance
    if(sdx==0) sdx=0.01 #can also happen if multiple weeks of 0 catch at start of sampling
    Nx_cv[iwk,iyr]=sdx/mu#the cv in untranformed space
    Nx_sd[iwk,iyr]=sqrt(log(Nx_cv[iwk,iyr]^2+1))#convert from cv in untransformed space to sd in log space
    #Nx_sd[iwk,iyr]=sqrt(log(0.05^2+1)) #check to see if high cv is what is messing up fit
  }#iwk

  #Get annual flows statistic

  irecs=which(Jwk<=53 |(Jwk>=1 & Jwk<=5)) #Flows prior to Feb in yr 't'
  CovX0[iyr,1]=max(Flow[irecs])
  CovX0[iyr,2]=CovX0[iyr,1] #use same covariate for prediction of phi and lambda but this structure allows one to use different covariates

}#iyr

#Standardize annual covariate values
CovLabel=rep("Peak flow prior to February (cfs)",2)
for(j in 1:2){
  muX=mean(CovX0[,j]);sdX=sd(CovX0[,j])
  CovX[,j]=(CovX0[,j]-muX)/sdX
}

dat=as.data.frame(cbind(year,Nestwks,Ntot_cv,t(Nx_cv)))
names(dat)=c("Year","NestWks","lgNtot_cv",paste0("lgNx_cv",1:Nwks))
fnout=paste0(OutDir,doTrib,".out")
write.table(file=fnout,x=dat,row.names=F,col.names=T)


data <- list(year = annual_abundance$year,
             number_estimated_weeks,
             annual_abudnance_cv,
             covariate_matrix)

