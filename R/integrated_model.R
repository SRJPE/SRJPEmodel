# get_survival_forecast <- function() {
#
# }

# args
set.seed(SRJPEmodel::forecast_seed)
water_year_forecast <- 2 # TODO confirm this with Josh, will this change? this represents the element # to use for forecast year (e.g., iwy=1 for critical, iwy=2 for D/BN/AN, iwy=3 for W)
cfg <- config::get()

con <- DBI::dbConnect(RPostgres::Postgres(),
                      dbname = cfg$db_name,
                      host = cfg$db_host,
                      port = cfg$db_port,
                      user = cfg$db_user,
                      password = cfg$db_password)

#1) Get the median size of spring run outmigrants by model week
#2) Get the multi-year proportion of spring run outmigrants leaving by model week from inseason model
#3) Get survival from RST to Delta for each model week given fish size for that week (1)
#4) Calculate a weighted average survival across run, with 2) providing the weights.

# load CJS survival model fit object
# TODO if survival covariate is ever null, we will have to use an if/else statement (see DS_Surv.R)
survival_fit <- readRDS(here::here("data-raw", "survival_model", "survival_model_fit_wy3_fl.rds")) # TODO this will be updated to get_most_recent_model_object(trib, model type)
survival_posterior <- as.data.frame(survival_fit,
                                       pars = c("SurvForecastSz","TribSurvForecastSz")) |>
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "estimate") |>
  arrange(parameter) |>
  # extract indices from STAN par estimates so we can filter. These will be
  # size class and water year
  separate_wider_delim(parameter, delim = "[", names = c("parameter", "indices"),
                       too_few = "align_start") |>
  mutate(indices = str_remove_all(indices, "\\]")) |>
  separate_wider_delim(indices, ",", names = c("size_class", "wy_class", "reach"),
                       too_few = "align_start") |>
  mutate(across(size_class:reach, as.numeric))


n_sims <- dim(survival_param_matrix)[1]

n_survival_locations <- 3 # Types of survival predictions from CJS model (mainstem, butte, feather)
trib_type_lookup <- c(1, 1, 1, 1, 2, 3) # ubc, ucc, mill, deer = 1, butte = 2, and eventually feather = 3

# Read in CJS forecasted survival posteriors
n_size_classes <- 25 # TODO confirm if this will change
size_covariate  <- seq(from = 10, to = 130, length.out = n_size_classes)
survival_posterior <- array(dim = c(n_sims, n_size_classes, n_survival_locations))

survival_statistics <- survival_posterior |>
  mutate(location = ifelse(parameter == "SurvForecastSz", "mainstem", "tributary")) |>
  group_by(location, wy_class, size_class, reach) |>
  summarise(mean = mean(estimate),
            `10` = quantile(estimate, 0.1),
            `50` = quantile(estimate, 0.5),
            `90` = quantile(estimate, 0.9)) |>
  ungroup()

# set up weights
Fl_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year
Surv_mu=vector(length=Ntribs)#pRun-weighted average forklength in forecast year from weighted posterior sample
Surv_mu_check=vector(length=Ntribs)#pRun_weighted average survival calculated intuitive way for checking

#mean and sd of logit tranformed posterior samples of survival and the pRun-weighted mean forklength for JPE.stan model
DS_surv_mu=vector(length=Ntribs);
DS_surv_sd=DS_surv_mu

# extract PLAD results
plad_results <- purrr::pmap(list(site = SRJPEmodel::forecast_sites),
                            read_plad_csv) |>
  reduce(bind_rows)
# TODO we need julian_week,




for(itrib in 1:Ntribs){

  # 1) Get across-year mean size of spring run by week from PLAD  #########
  #First step is to load the PLAD sizes for each trib and week
  fn_sz=paste0(PLADpath,"/SRforklength_",DoTrib[itrib],".csv")
  dsz=read.csv(file=fn_sz,header=T)
  Jwk=as.integer(rownames(dsz))
  Nwks=length(Jwk)# # of wks where sizes are available from PLAT

  # 2) Get proportion of run passing trap for every week in year (Sep-Aug) ########
  if(itrib==1){
    pRun=matrix(nrow=Ntribs,ncol=Nwks)
    Sz=pRun;  Surv=pRun

    labwk=character(length=Nwks)#for plotting
    for(iwk in 1:Nwks)labwk[iwk]=calwk[which(mwk==Jwk[iwk])]
  }
  Sz[itrib,1:Nwks]=dsz$X0.5# put size in a matrix for later plotting

  #Load inseason outmigrant timing model to calculate proportion of run for each week
  fnm=paste0(InseasDir,DoTrib[itrib],"_null.Rdata")
  load(file=fnm)
  dcp=as.data.frame(fit, pars = c("For_cp"))

  iwk=1;imwk=which(mwk==Jwk[iwk])#identify the inseaon model week for the Julian week designation in PLAD
  pRun[itrib,iwk]=mean(dcp[,imwk],na.rm=T)
  for(iwk in 2:Nwks){
    imwk=which(mwk==Jwk[iwk])
    pRun[itrib,iwk]=mean(dcp[,imwk],na.rm=T)-sum(pRun[itrib,1:(iwk-1)])
  }
  pRun[itrib,]=pRun[itrib,]/sum(pRun[itrib,])#Make sure it sums to one to avoid any discretization error


  #3) Get pRun-weighted size and survival rate for outmigrant run ################

  #Set the jtrib index needed for SurvForecastSz. TribSurvForecastSz, and survival_statistics
  if(DoTrib[itrib]=="ubc"| DoTrib[itrib]=="ucc"| DoTrib[itrib]=="mill creek"| DoTrib[itrib]=="deer creek"){
    jtrib=1
  } else if (DoTrib[itrib]=="okie dam"){
    jtrib=2
  } else {#Feather eventually
    jtrib=3
  }

  Surv_mu_check[itrib]=0;Fl_mu[itrib=0]
  #Set of samples to grab from posterior sample of CJS survival for size class associated with following weeks
  if(itrib==1) isamps=sample(1:n_sims,size=500,replace=F) #to keep comp_post to reasonable size, sample 100 posterior survival values from each size class (wk)
  for(iwk in 1:Nwks){

    #find the closest size in size_covariate given median size in this week from plad (Sz) to get the size_class index for size_covariate used in CJS model
    irecs=which(size_covariate<=Sz[itrib,iwk])
    if(length(irecs)==0){
      size_class=1 #Sz is smaller than lowest value of size_covariate so set to lowest index
    } else{
      size_class=max(irecs) #For Sz>max(size_covariate) this will set index to last element for size_covariate
    }

    #Create a composite posterior of survival rates which accumulates posterior samples across all weeks/size classes.
    #Posterior sample each wk/size class is weighted based on pRun relative to lowest pRun aross weeks (where reps=1)
    preps=round(pRun[itrib,iwk]/min(pRun[itrib,])) # # of replicates for this iwk
    if(iwk==1){
      comp_post=rep(survival_posterior[isamps,size_class,jtrib],times=preps) #repeat the random sample of posterior for this wk preps times
    } else {
      comp_post=c(comp_post,rep(survival_posterior[isamps,size_class,jtrib],times=preps))#each set of repeated samples is added to the vectory
    }

    Fl_mu[itrib]=Fl_mu[itrib]+pRun[itrib,iwk]*Sz[itrib,iwk]#pRun-weighted mean size

    #Most obvious way to compute a pRun-weighted average survival.
    Surv_mu_check[itrib]=Surv_mu_check[itrib]+pRun[itrib,iwk]*survival_statistics[jtrib,size_class,2]
  }
  Surv_mu[itrib]=mean(comp_post)#To compare against Surv_mu_check

  #mean and sd of survival rate in logit space for outmigrant run for JPE.stan model
  DS_surv_mu[itrib]=mean(logit(comp_post))
  DS_surv_sd[itrib]=sd(logit(comp_post))

  print(c(DoTrib[itrib],round(Surv_mu_check[itrib],digits=4),round(Surv_mu[itrib],digits=4)))

  #Simulated survival rates as done in JPE.stan (but with lower constraint of -6.5 ~ 0.15% survival)
  #sim_surv used to see if composite posterior (comp_post) is accurately modelled by DS_surv_mu and DS_surv_sd and JPE.stan approap
  #For plotting (if PlotType==2) or text output on different plot type (if PlotType==1)
  sim_surv=inv_logit(rnorm(n=5000,mean=DS_surv_mu[itrib],sd=DS_surv_sd[itrib]))

  if(PlotType==1){
    #Plot proportion of run and size by week with weighted survival mean and other stats printed on each panel
    plot(1:Nwks,pRun[itrib,],type='l',bty='n',main=DoTrib[itrib],xlab="",ylab="",axes=F)
    axis(2);axis(1,at=1:Nwks,labels=labwk,las=3)
    mu1=mean(sim_surv); sd1=sd(sim_surv)
    #text(x=1,y=max(pRun[itrib,])*0.9, labels=paste0("Wgt'd mean FL ",round(Fl_mu[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.85, labels=paste0("Wgt'd mean Surv (Surv_mu_check) ",round(Surv_mu_check[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.8, labels=paste0("Mean Surv from comp_post (Surv_mu) ",round(Surv_mu[itrib],digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.75, labels=paste0("mean of tranformed DSsurv ",round(mu1,digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.7, labels=paste0("sd Surv from comp_post ",round(sd(comp_post),digits=2)),pos=4)
    #text(x=1,y=max(pRun[itrib,])*0.65, labels=paste0("sd of sim surv's from DSsurv_mu and _sd ",round(sd1,digits=2)),pos=4)


    par(new=T)
    plot(1:Nwks,Sz[itrib,],type='l',bty='n',col="red",xlab="",ylab="",axes=F)
    axis(4)
    if(itrib==Ntribs){
      legend("bottom",legend=c("proportion","forklength"),lty=c(1,1),col=c("black","red"),bty='n')
      mtext("Date",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
      mtext("Proportion of Outmigrant Run Passing Trap",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
      mtext("Median Forklength of Spring Run Outmigrants (mm)",side=4, las=3,line=-1,outer=T,cex=1.3,font=2,col="red")
    }

  } else if (PlotType==2){
    Pr=density(sim_surv)
    Po=hist(comp_post,plot=F)
    ymax=max(c(max(Po$density),max(Pr$y)))

    hist(comp_post,ylim=c(0,ymax),xlab="",ylab="",main=DoTrib[itrib],freq=F,axes=F)
    axis(1);lines(density(sim_surv),lty=2)

    if(itrib==Ntribs){
      legend("topright",legend=c("pRun-weighted composite posterior","JPE.stan simulated posterior"),lty=c(NA,2),pch=c(22,NA),pt.bg=c("gray","NA"),bty='n')
      mtext("Downstream Survival Rate",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
      mtext("Frequency",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)
    }

  }
}#next itrib


#Plot survival stats from CJS model and overlay stats that will be used in JPE.stan based on pRun-weighted mean size. Should be close
par(mfrow=c(2,2),xaxs="i",mai=c(.8,.7,.5,.5), omi=c(0.25,0.25,0.5,0.25),cex.axis=0.9)
xcol=c("black","blue","red","green")#multiple tribs use the same mainstem survival relationship
for(jtrib in 1:3){
  xmain=switch(jtrib,"ubc/ucc/mill/deer","Butte","Feather")
  plot(size_covariate,survival_statistics[jtrib,,2],type='l',bty='n',ylim=c(0,max(survival_statistics[jtrib,,3])),xlab="",ylab="",main=xmain)
  polygon(x = c(rev(size_covariate), size_covariate), y = c(rev(survival_statistics[jtrib,,1]), survival_statistics[jtrib,,3]), col = alpha("grey", 0.25), border = NA)

  #check to see if simulated survival going into stan model (based on DS_surv_mu and DS_surv_sd
  #is reproducing what is coming out of CJS model
  if(jtrib<=2){ #exclude feather for now
    trib_type_lookup=c(1,1,1,1,2,3)
    itribs=which(trib_type_lookup==jtrib)

    for(j in 1:length(itribs)){

      y1=inv_logit(rnorm(n=5000,mean=DS_surv_mu[itribs[j]],sd=DS_surv_sd[itribs[j]]))
      dstats=as.double(quantile(y1,probs=CI));dstats[2]=mean(y1)

      #For x-axis position, find the size_covariate index closes to FL_mu, the pRun-weighted mean size
      irecs=which(size_covariate<=Fl_mu[itribs[j]])
      if(length(irecs)==0){size_class_mu=1} else{size_class_mu=max(irecs)}
      #points(size_covariate[size_class_mu],dstats[2],pch=19,cex=1.2,col=xcol[j])
      #lines(x=c(size_covariate[size_class_mu],size_covariate[size_class_mu]),y=c(dstats[1],dstats[3]),col=xcol[j])
    }
  }
  #if(jtrib==1) legend("topright",legend=c(DoTrib[1:4]),col=xcol,pch=rep(19,4),bty='n')
}
mtext("Forklength (mm)",side=1, las=1,line=-1,outer=T,cex=1.3,font=2)
mtext("Survival Rate to Delta",side=2, las=3,line=-1,outer=T,cex=1.3,font=2)



read_plad_csv <- function(site) {
  PLAD_results <- read_csv(paste0(here::here("data-raw", "PLAD", "PLAD_results"),
                                  "SR_forklength_", site, ".csv")) |>
    mutate(site = site)

  return(PLAD_results)
}
