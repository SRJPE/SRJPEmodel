library("R2WinBUGS")
library(splines2)
library(tidyverse)

EffortAdjust = T
MultiRun_Mode = T

doTrib <- "battle creek"
doYr <- 2009

# TODO initiate multirun logic

if (MultiRun_Mode == F) {
  #doTrib="deer creek_deer creek";doYr=2010
  doTrib = "battle creek_ubc"
  doYr = 2009
  #doTrib="clear creek_lcc";doYr=2005
  #doTrib="feather river_eye riffle";doYr=2020

  doWks = c(seq(45, 53), seq(1, 22)) #can't have a wider range than range used to build RST_input.txt (as specified in BuildData.R)
  #doWks=seq(1,19)

  fnprefix = paste0(doTrib, "_", doYr)
  #fnprefix=paste0("output/",doTrib,"_",doYr)
} else {
  fnprefix = paste0("data-raw/first_draft/OutSpecPriors/", doTrib, "_", doYr)
}

# TODO add to cache_params.R or define as model parameter at model run
number_mcmc <- 10000
number_burnin <- 5000
number_thin <- 5
number_chains <- 3


### Setup Data and initial values for pCap component of model ########
# TODO migrate to data repo
weekly_tributary_RST = as_tibble(read.table(file = "data-raw/first_draft/RST_Input.txt", header = T)) |>
  rename(stream = Trib, calendar_year = CalYr, run_year = RunYr, week = Week,
         catch_all = ua, catch_spring = us, fork_length = Spr_Sz,
         lifestage = Spr_Stage, release = Rel1, recapture = Recap1,
         release_size = Rel_Sz, effort = Effort, flow_at_release = mr_flow,
         flow_at_catch = catch_flow, origin_released = rel_origin) |>
  separate(stream, c("stream", "site"), "_", remove = T) |>
  mutate(catch_spring = ifelse(stream %in% c("mill creek", "deer creek"), catch_all, catch_spring)) |>
  select(-c(catch_all, rel_daynight)) |>
  filter(stream != "sacramento river") |> # exclude sac from trib model
  glimpse()


selected_stream <- "battle creek" # added to test model
all_streams = unique(weekly_tributary_RST$stream)
number_streams = length(all_streams)

#Trib index to use to predict efficiency for missing strata.
#Note length=0 if doTrib is not part of trib set that has MR data. In this case model that samples from trib hyper will be called
use_stream = which(all_streams == selected_stream)
selected_run_years <- 2009

select_stream_data <- function(all_data, selected_stream, selected_run_years) {
  all_data |>
  filter(stream == selected_stream,
         run_year %in% selected_run_years)
}

selected_stream_data <- select_stream_data(weekly_tributary_RST, selected_stream, selected_run_years)

select_efficiency_data <- function(all_data) {
  all_data |>
    filter(!is.na(release), !is.na(recapture), !is.na(flow_at_release))
}

efficiency_data <- select_efficiency_data(weekly_tributary_RST)


# Unique MR experiments across tribs, years and weeks
efficiency_trials = unique(efficiency_data[c("stream", "run_year", "week")])
number_efficiency_trials = nrow(efficiency_trials)
efficiency_trials_years = unique(efficiency_trials$run_year)

selected_weeks <- 1:53
number_weeks <- length(selected_weeks)

indexes_stream <- tibble(stream = all_streams,
                         stream_idx = row_number(all_streams))
indexes_week <- tibble(week = selected_weeks,
                       week_idx = row_number(selected_weeks))
indexes_year <- tibble(run_year = selected_run_years,
                       year_idx = row_number(selected_run_years))

stream_flows <- efficiency_data |>
  group_by(stream, site) |>
  summarise(stream_flow_mean = mean(flow_at_catch, na.rm = T),
            stream_flow_sd = sd(flow_at_catch, na.rm = T)) |>
  glimpse()

efficiency_data_with_indices <- efficiency_trials |>
  left_join(efficiency_data |> select(stream, run_year, week, release, recapture, flow_at_release),
            by = c("stream", "run_year", "week")) |>
  left_join(indexes_stream, by = "stream") |>
  left_join(indexes_week, by = "week") |>
  left_join(indexes_year, by = "run_year") |>
  left_join(stream_flows, by = "stream") |>
  filter(run_year %in% selected_run_years,
         week %in% selected_weeks,
         stream %in% selected_stream) |>
  mutate(standardized_flow = (flow_at_release - stream_flow_mean) / stream_flow_sd) |>
  select(-c(stream_flow_mean, stream_flow_sd, flow_at_release)) |>
  glimpse()

# TODO check this logic
if(EffortAdjust == FALSE) {
  abundance_data <- efficiency_data |>
    filter(stream %in% selected_stream,
           week %in% selected_weeks,
           !is.na(week),
           run_year %in% selected_run_years,
           !is.na(effort)) |>
    mutate(abundance = catch_spring)
} else {
  abundance_data <- efficiency_data |>
    filter(stream %in% selected_stream,
           week %in% selected_weeks,
           !is.na(week),
           run_year %in% selected_run_years,
           !is.na(effort)) |>
    mutate(effort = ifelse(stream %in% c("deer creek", "mill creek"), 1, effort),
           mean_effort = ifelse(stream %in% c("deer creek", "mill creek"), 1, mean(effort, na.rm = T)),
           abundance = round((catch_spring * mean_effort)/effort, 0))
}


#Prior for upper limit on log abundance for any week. Note lgN estimated in units of thousands to avoid values in millions which occurs becasue of huge weekly catches in Butte
#Likely results in reduced precision problems and better mising for other systems with big numbers as well

#BT-SPAS constraint is 20 in log space but units are in 1's, not '000s.
#So for BT-SPAS-X convert 20, divide by 1000, then put back into log units.
#lgN_max=rep(log(exp(20)/1000),times=Nstrata)
u <- abundance_data$abundance

#Default prior fornow
lgN_max = log((u / 1000 + 1) / 0.025)		#maxmim possible value for log N across strata
imiss = which(is.na(u) == T)
lgN_max[imiss] = mean(lgN_max, na.rm = T) #for missing strata set to average max across strata

#lgN_max for special cases
dsp = read.csv(file = here::here("data-raw", "first_draft", "Special_Priors.csv"), header = T)
dsp1 = subset(dsp, Stream_Site == doTrib & RunYr == doYr)
istrata = which(is.na(match(doWks, dsp1$Jweek)) == F)
lgN_max[istrata] = dsp1$lgN_max


#Calculate standardized flow for each weekly strata (can be different than mr_flow which only averages over days of recovery
#Note the mean and sd used for standardization has to be the same used for pCap model (above) for this trip
dX = subset(abundance_data, stream == doTrib & is.na(match(week, doWks)) == F)
cf = abundance_data$flow_at_catch
irecs = which(is.na(abundance_data$flow_at_catch) == T)
if (length(irecs) > 0)
  cf[irecs] = mean(cf, na.rm = T) #if missing values for one or more strata set to mean (but no missing values as of Aug 02,2022.
catch_flow = (cf - mean(dX$flow_at_catch, na.rm = T)) / sd(dX$flow_at_catch, na.rm =
                                                          T)

dC <- abundance_data
#Identify the elements in 1:Nstrata (unmarked catch set) without and with pCap and corresponding flow data
#Alternate is to identify records without mr_flow but with catch_flow, and then sub catch_flow into mr_flow for these cases. Complicated and a bit manky.
Uind_woMR = which(is.na(dC$release) == T |
                    is.na(dC$flow_at_release) == T)
Nwomr = length(Uind_woMR)
Jwks_womr = dC$week[Uind_woMR]
Uind_wMR = which(is.na(dC$release) == F &
                   is.na(dC$flow_at_release) == F)
Nwmr = length(Uind_wMR)
Jwks_wmr = dC$week[Uind_wMR]
if (Nwomr == 1)
  Uind_woMR = c(Uind_woMR, -99)
if (Nwmr == 1)
  Uind_wMR = c(Uind_wMR, -99)#so bugs doesn't bomb if only one strata with missing MR

MR <- read.table(here::here("data-raw", "first_draft", "MRdata.txt"), header = T)

#identify the elements in pCap from full MR dataset for U Strata being estimated
ind_pCap = which(MR$Trib == "battle creek_ubc" &
                   MR$RunYr == doYr & is.na(match(MR$Wk, Jwks_wmr)) == F)
Uwc_ind = which(is.na(u) == F)
Nstrata = length(Uwc_ind) #The elements of 1:Nstrata that have catch data (RST fished)

#Setup a data file for current run which lines up with output. Used for Plotmodel
trel = rep(NA, Nstrata)
trel[Uind_wMR] = MR$Rel1[ind_pCap]
trec = rep(NA, Nstrata)
trec[Uind_wMR] = MR$Recap1[ind_pCap]
write.table(
  file = here::here(paste0(fnprefix, "_data.out")),
  cbind(
    dC$week,
    u,
    trel,
    trec,
    round(abundance_data$effort, digits = 1),
    round(catch_flow, digits = 1),
    round(dC$flow_at_catch, digits = 1),
    round(exp(lgN_max), digits = 1)
  ),
  quote = F,
  row.names = F,
  col.names = c(
    "Jwk",
    "u",
    "Releases",
    "Recaptures",
    "Effort",
    "Std_Flow",
    "Flow",
    "N_max_000s"
  )
)

#### Setup B-spline basis matrix
k_int = 4	#Rule of thumb is 1 knot for every 4 data points for a cubic spline (which has 4 parameters)
Nknots = round(Nstrata / k_int, digits = 0)
fkp = 2
lkp = Nstrata - 1 #Keep first and/or last knot positions away from tails if there are intervals with no sampling on the tails
kknots = seq(fkp, lkp, length.out = Nknots) #Define position of b-spline knots using even interval if no missing data
ZP <-
  bSpline(
    x = 1:Nstrata,
    knots = kknots,
    deg = 3,
    intercept = T
  )	 #bspline basis matrix. One row for each data point (1:Nstrata), and one column for each term in the cubic polynomial function (4) + number of knots
K = dim(ZP)[2]
#######################################################################################################

### Pass data and initial values to lists, define parameters to save, and run model ##################

#Select correct model given data for this run. Note slightly different data inputs across models
if (length(doTrib) == 1) {
  #Some MR was done in this tributary

  if (Nwomr == 0) {
    #all strata have corresponding MR data
    ModName = "estN_allMR"
    data <-
      list(
        "Nmr",
        "Ntribs",
        "ind_trib",
        "Releases",
        "Recaptures",
        "Nstrata",
        "u",
        "K",
        "ZP",
        "ind_pCap",
        "Nstrata_wc",
        "Uwc_ind",
        "mr_flow",
        "lgN_max"
      )

    #ModName="gamma_estN_allMR"
    #data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nstrata_wc","Uwc_ind","mr_flow","lgN_max")

  } else {
    #some or all strata don't have MR data

    if (Nwmr > 0) {
      #some strata have MR data
      ModName = "estN_missMR"		#MR data for some strata and tributary was included in the MR analysis

      data <-
        list(
          "Nmr",
          "Ntribs",
          "ind_trib",
          "Releases",
          "Recaptures",
          "Nstrata",
          "u",
          "K",
          "ZP",
          "ind_pCap",
          "Nwmr",
          "Nwomr",
          "Uind_wMR",
          "Uind_woMR",
          "use_trib",
          "Nstrata_wc",
          "Uwc_ind",
          "mr_flow",
          "catch_flow",
          "lgN_max"
        )

      #ModName="gamma_estN_missMR"		#MR data for some strata and tributary was included in the MR analysis
      #data<-list("Nmr","Ntribs","ind_trib","Releases","Recaptures","Nstrata","u","K","ZP","ind_pCap","Nwmr","Nwomr","Uind_wMR","Uind_woMR","use_trib","Nstrata_wc","Uwc_ind","mr_flow","catch_flow","lgN_max")

    } else if (Nwmr == 0) {
      #No strata have MR data
      ModName = "estN_noMR"		#No MR data for any strata and trib was included in MR analysis
      data <-
        list(
          "Nmr",
          "Ntribs",
          "ind_trib",
          "Releases",
          "Recaptures",
          "Nstrata",
          "u",
          "K",
          "ZP",
          "Nwomr",
          "Uind_woMR",
          "use_trib",
          "Nstrata_wc",
          "Uwc_ind",
          "mr_flow",
          "catch_flow",
          "lgN_max"
        )
    }
  }

} else if (length(doTrib) == 0) {
  #No MR was done in this trib
  ModName = "estN_noMR_notrib"	#No MR data for any strata and trib was not included in MR analysis
  data <-
    list(
      "Nmr",
      "Ntribs",
      "ind_trib",
      "Releases",
      "Recaptures",
      "Nstrata",
      "u",
      "K",
      "ZP",
      "Nwomr",
      "Uind_woMR",
      "Nstrata_wc",
      "Uwc_ind",
      "mr_flow",
      "catch_flow",
      "lgN_max"
    )#
}


#"beta_dev","tau.eta",
#parameters<-c("trib_mu.P","trib_sd.P","flow_mu.P","flow_sd.P","pro_sd.P","b0_pCap","b_flow","pCap_U","N","Ntot","sd.N","sd.Ne","sd.Ntot","pSplineVar")
parameters <-
  c(
    "trib_mu.P",
    "trib_sd.P",
    "flow_mu.P",
    "flow_sd.P",
    "pro_sd.P",
    "b0_pCap",
    "b_flow",
    "pCap_U",
    "N",
    "Ntot",
    "sd.N",
    "sd.Ne"
  )

Ntribs <- number_streams
ini_b0_pCap = vector(length = Ntribs)
for (itrib in 1:Ntribs) {
  irows = which(indexes_stream == itrib)
  ini_b0_pCap[itrib] = logit(sum(MR$Recap1[irows]) / sum(MR$Rel1[irows]))
}
#N is estimated in units of thousands in log space
ini_lgN = log(u / 1000 + 2)
irecs = which(is.na(ini_lgN) == T)
ini_lgN[irecs] = log(2 / 1000)

inits1 <-
  list(
    trib_mu.P = logit(sum(MR$Recap1) / sum(MR$Rel1)),
    b0_pCap = ini_b0_pCap,
    flow_mu.P = 0,
    b_flow = rep(0, Ntribs),
    trib_tau.P = 1,
    flow_tau.P = 1,
    pro_tau.P = 1,
    b_sp = rep(1, K),
    lg_N = ini_lgN
  )
inits2 <- inits1
inits3 <- inits1
inits <- list(inits1, inits2, inits3)


ModName2 = paste0(ModName, ".bug")
post <-
  bugs(
    data,
    inits,
    parameters,
    ModName2,
    n.chains = number_chains,
    n.burnin = number_burnin,
    n.thin = number_thin,
    n.iter = number_mcmc,
    debug = F,
    codaPkg = F,
    DIC = T,
    clearWD = T,
    bugs.directory = "c:/WinBUGS14/"
  )	#run bugs model

fnstats = paste0(fnprefix, "_post.out")
write.table(
  post$sims.list,
  file = fnstats,
  col.names = T,
  row.names = T
)
fnstats = paste0(fnprefix, "_sum.out")
write.table(
  round(post$summary, digits = 3),
  file = fnstats,
  col.names = T,
  row.names = T
)
fn_dic = paste0(fnprefix, "_dic.out")
write(file = fn_dic, c(post$pD, post$DIC), ncolumns = 2)
fn_knots = paste0(fnprefix, "_knots.out")
write(file = fn_knots, kknots, ncolumns = Nknots)

if (MultiRun_Mode == T) {
  print(c(doTrib, doYr, ModName))
  Sys.sleep(0.01)
  flush.console()
}
