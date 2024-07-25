library(scales)
library(tidyverse)

results <- readRDS("data-raw/juvenile_abundance/ubc_2004_2024-07-23.rds")
# deer_YOY <- readRDS("data-raw/juvenile_abundance/deer_2023_YOY_results.rds")
site <- "ubc"
run_year <- 2004
lcl <- 0.025
ucl <- 0.975

# TODO scale results up to 1000s of fish

# TODO use lubridate to create this
julian_week_to_date_lookup <- read.table(file = "data-raw/juvenile_abundance/btspas_model_code/Jwk_Dates.txt", header = F) |>
  tibble() |>
  filter(V1 != "Jwk") |>
  mutate(V1 = as.numeric(V1)) |>
  select(Jwk = V1, date = V2)
# Jdt=dx$Date

# TODO if reading in _sum.out, use $summary in model object

# read in input data
input_data <- SRJPEdata::weekly_juvenile_abundance_model_data |>
  filter(site == "ubc",
         run_year == 2004,
         week %in% c(seq(45, 53), seq(1, 22))) |> # TODO update with arguments
  mutate(lincoln_peterson_abundance = count * (number_released / number_recaptured)) # TODO sd calculation is what
         # lincoln_peterson_abundance_sd = (number_released * number_recaptured * (number_released - number_recaptured) * count) /
         #   (number_recaptured^2 * number_recaptured))

summary_output <- results$final_results |>
  select(-c(model_name, srjpedata_version)) |>
  pivot_wider(id_cols = site:parameter,
              values_from = value,
              names_from = statistic) |>
  mutate(cv = round(100 * (sd / mean), digits = 0)) # TODO confirm we use cv for more than just Ntot plot
# TODO what do we need from "input data" besides Nstrata (which we can get from result object?)
# TODO output from larger model results object
sims_list_output <- results$full_object$sims.list # TODO what do we use this for? in separate object
n_sims <- results$full_object$n.sims

total_abundance <- summary_output |>
  filter(parameter == "Ntot")

# plot abundance only

# plot all model output
# don't need to use the "sims list" because summary table in BUGS calculates quantiles

# set up cv, lcl, and ucl and make sure to scale to 1000s

plot_data <- summary_output |>
  filter(str_detect(parameter, "N\\[")) |>
  left_join(julian_week_to_date_lookup, by = c("week_fit" = "Jwk")) |>
  mutate(fake_date = ifelse(week_fit > 35, paste0(run_year - 1, "-", date),
                            paste0(run_year, "-", date)),
         fake_date = as.Date(fake_date, format = "%Y-%b-%d")) |>
  left_join(input_data |>
              select(week, site, run_year, lincoln_peterson_abundance),
            by = c("site", "run_year", "week_fit" = "week"))

# abundance bar plot
# TODO scale?
plot_title <- paste0(str_to_title(site), " ", run_year, " predicted annual abundance = ",
                     prettyNum(round(total_abundance$`50`, 0), big.mark = ","), " (",
                     prettyNum(round(total_abundance$`2.5`, 0), big.mark = ","), "-",
                     prettyNum(round(total_abundance$`97.5`, 0), big.mark = ","), ")")
plot_data |>
  ggplot(aes(x = fake_date, y = `50`)) +
  geom_bar(stat = "identity", fill = "grey") +
  geom_errorbar(aes(x = fake_date, ymin = `2.5`, ymax = `97.5`), width = 0.2) +
  geom_point(aes(x = fake_date, y = lincoln_peterson_abundance),
             shape = 1, color = "blue") +
  theme_minimal() +
  labs(x = "Date", y = "Abundance",
       title = plot_title)



xmain=paste0(RunNm," Ntot=",NtotLab)
ylims=c(-0.1,max(c(N,pN))*1.1)
irecs=which(is.na(d$u)==T);#points(irecs,rep(ylims[1],length(irecs)),pch=19,cex=1.1,col='red') #identify strata that were not sampled
bcols=rep("grey",Nstrata);bcols[irecs]="red"
x=barplot(height=N[,2],col=bcols,names.arg=Jdt[d$Jwk],cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l',main=xmain, xlab="",ylab="Abundance ('000s)")

  for(i in 1:Nstrata){#error bars and catch
    arrows(x0=x[i],x1=x[i],y0=N[i,1],y1=N[i,3],angle=90,code=3,length=0.025,col="black")
    par(srt=90);text(x=x[i],y=ylims[2]*0.9,labels=d$u[i],cex=0.7);par(srt=0)
    if(is.na(match(i,irecs))==F) points(x[i],ylims[2]*0.9,pch=19,col="red",cex=0.8)#identify strata with no sampling (red bars may not show up if estimated abundance is very low)
  }

  #Plot Peterson estimates and CL's
  coff=0.1
  nn=length(irecs_N)
  if(nn>0){
    for(i in 1:nn){
      points(x[irecs_N[i]]-coff,pN[i],pch=21,cex=1.5,col='blue')	#data-driven predictions of N
      if(sdN[i]>0){
        lc=pN[i]-1.96*sdN[i];uc=pN[i]+1.96*sdN[i]
        arrows(x0=x[irecs_N[i]]-coff,x1=x[irecs_N[i]]-coff,y0=lc,y1=uc,angle=90,code=3,length=0.025,col="blue")
      }
    }
  }
  #legend("left",legend=c("Bayesian","Peterson","No sampling","Spline knot"),pch=c(22,21,19,25),pt.cex=rep(1.2,4),pt.bg=c("grey",NA,NA,"green"),col=c("grey","blue","red","green"),cex=0.8,bty='n')
  legend("right",legend=c("Bayesian","Peterson","No sampling or no releases","Discharge"),pch=c(22,21,19,NA),pt.cex=c(rep(1.2,3),NA),pt.bg=c("grey",NA,NA,NA),col=c("grey","blue","red","darkgreen"),lwd=c(NA,NA,NA,1.2),cex=0.8,bty='n')

  if(PlotAbOnly==T){
    par(new=T)
    plot(x-coff,d$Flow/1000,ylim=c(min(d$Flow/1000),max(d$Flow/1000)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)
    mtext("Discharge ('000s cfs)",side=4, las=3,adj=0.5,outer=T,cex=1.0,font=2,line=-1,col="darkgreen")
  }

  if(PlotAbOnly==F){
    ###Plot weekly capture probability #########################
    pCap=matrix(nrow=Nstrata,ncol=3)
    for(i in 1:Nstrata){
      vn=paste0("pCap_U.",i);icol=which(names(dp)==vn)
      pCap[i,]=as.double(quantile(dp[,icol],prob=c(lb,0.5,ub)))
    }
    ylims=c(0,max(pCap)*1.1)
    irecs=which(is.na(d$Releases)==T)#;points(irecs,pCap[irecs,2],pch=19,cex=1.1,col='red')
    bcols=rep("grey",Nstrata);bcols[irecs]="red"
    x=barplot(height=pCap[,2],col=bcols,names.arg=Jdt[d$Jwk],cex.names=0.7, las=3,axis.lty=1,ylim=ylims,  bty='l', xlab="First Date of Week",ylab="Capture Probability")

    irecs=which(is.na(d$Releases)==F)
    points(x[irecs]-coff,d$Recaptures[irecs]/d$Releases[irecs],pch=21,cex=1.5,col='blue')

    for(i in 1:Nstrata) {
      arrows(x0=x[i],x1=x[i],y0=pCap[i,1],y1=pCap[i,3],angle=90,code=3,length=0.025,col="black")
      par(srt=90)
      text(x=x[i],y=ylims[2]*0.95,labels=d$Recaptures[i],cex=0.7)
      text(x=x[i],y=ylims[2]*0.82,labels=d$Releases[i],cex=0.7)
      par(srt=0)
    }

    par(new=T)
    plot(x-coff,d$Flow/1000,ylim=c(min(d$Flow/1000),max(d$Flow/1000)*1.2), xlim=c(0,trunc(max(x))+1),type='o',pch=19,cex=0.7,lwd=2,col="darkgreen",axes=F,xlab="",ylab="");axis(4,cex.axis=0.9)
    mtext("                           Discharge ('000s cfs)",side=4, las=3,adj=0,outer=T,,cex=1.0,font=2,line=-1,col="darkgreen")
  }
  print(c(doTrib,doYr,Ntot_Rhat))

  ### Plot pCap distributions ####
  if(MultiRun_Mode==F){

    #pCap distributions for each tributary in one panel and overlay hyper distribution
    par(mfcol=c(1,1),mai=c(0.85,2.25,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
    xmax=0.3
    pstat=matrix(nrow=Ntribs,ncol=3)
    for(i in 1:Ntribs){
      vn=paste0("b0_pCap.",i);icol=which(names(dp)==vn)
      pCap=exp(dp[,icol])/(1+exp(dp[,icol]))
      pstat[i,]=as.double(quantile(pCap,prob=c(0.025,0.5,0.975)))
    }
    plot(pstat[,2],1:Ntribs,pch=19,bty='l',xlim=c(0,xmax),axes=F,ylab="",xlab="Mean Trap Efficiency")

    tnm=character(length=Ntribs);mup=vector(length=Ntribs)
    xTrib=Trib#make a copy since subset gets confused when dataframe variable name is the same as the one it is checking against
    mup=vector(length=Ntribs)
    for(i in 1:Ntribs){
      dmr1=subset(dmr,Trib==xTrib[i] & is.na(Rel1)==F)
      mup[i]=sum(dmr1$Recap1)/sum(dmr1$Rel1)	#mup[i]=mean(dmr1$Recap1/dmr1$Rel1)
      nsamps=dim(dmr1)[1]
      tnm[i]=paste0(gsub(pattern="riffle",replacement="",x=Trib[i])," (",nsamps,")")
    }
    axis(1);axis(2,at=1:Ntribs,labels=tnm,cex.axis=0.75)

    for(i in 1:Ntribs){#plot CI's and point estimate for each trib
      arrows(x0=pstat[i,1],x1=pstat[i,3],y0=i,y1=i,angle=90,code=3,length=0.025,col="black")
      points(x=mup[i],y=i,pch=21,cex=1.5,col='blue')
    }
    hmu=mean(dp$trib_mu.P);	hmu=exp(hmu)/(1+exp(hmu)); abline(v=hmu,lty=2)

    #hyper distribution
    pvec=seq(0,xmax,length.out=50);logit_pvec=log(pvec/(1-pvec))
    prob=dnorm(x=logit_pvec,mean=mean(dp$trib_mu.P),sd=mean(dp$trib_sd.P))
    par(new=T); plot(pvec,prob,lty=1,type='l',xlim=c(0,xmax),axes=F,xlab="",ylab="")
    legend("topright",pch=c(19,21,NA),lty=c(NA,NA,1),col=c("black","blue","black"),legend=c("Model Estimate","Point Estimate","Hyper-distribution"),bty='n')


    par(mfcol=c(1,2),mai=c(1,0.75,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
    #hyper-distributions from which pCaps will be drawn for tribs with no efficiency estimates
    b0=vector(length=Nsims);b1=b0
    for(i in 1:Nsims){
      b0[i]=rnorm(n=1,dp$trib_mu.P[i],dp$trib_sd.P[i])
      b1[i]=rnorm(n=1,dp$flow_mu.P[i],dp$flow_sd.P[i])
    }
    pCap3=exp(b0)/(1+exp(b0))
    hist(pCap3,col="grey",xlab="Trap Efficiency for Tribs Without MR at Average Flow",breaks=seq(0,1,by=0.005),xlim=c(0,xmax),main="",yaxt='n',ylab="")


    #Plot discharge-trap efficiency relationship for a trib without MR data
    Q=seq(-2,2,by=0.01)
    pstat=matrix(nrow=length(Q),ncol=3)
    for(i in 1:length(Q)){
      pCap4=exp(b0+b1*Q[i])/(1+exp(b0+b1*Q[i]))
      pstat[i,]=as.double(quantile(pCap4,probs=c(lb,0.5,ub)))
      pstat[i,2]=mean(pCap4)
    }

    plot(Q,pstat[,2],type='l',lty=1,bty='l',ylim=c(0,max(pstat)),ylab="Trap Efficiency",xlab="Standardized Discharge")
    polygon(x = c(rev(Q), Q), y = c(rev(pstat[,1]), pstat[,3]), col = alpha("grey", 0.3), border = NA)

  }

}
