#Create data file for RunSize and pCap models. 
#Assume that if trap was fished and no chinook caught, a value of u=0 is in the input file
#and that if trap was not fished there is no record in file

rm(list=ls(all=TRUE)) 

fnout="RST_Input.txt"
#fnout="RST_Input_June.txt"

Fyr=2003;Lyr=2020	#RunYr = 45-53 from CalYr-1 and 1:22 from CalYr. So if last CalYr in unmarked catch is 2021, last year to read is 2020
Year=seq(Fyr,Lyr);Nyrs=length(Year)

Fwk=45; Lwk=22	#Nov 01 - end of May
#Fwk=40;Lwk=26	#Sep 30 - end of June

Week=c(seq(1,Lwk),seq(Fwk,53));Nwks=length(Week)

ddate=read.table(file="../Data/Jwk_Dates.txt",header=T)
Jmon=ddate$Month;Jday=ddate$Day

dC0=read.csv(file="../data/jpe_weekly_unmarked_catch.csv",header=T)
dMR0=read.csv(file="../data/Weekly_Release_Recaptures.csv",header=T)
dQ0=read.csv(file="../data/Daily/standard_flow.csv",header=T)

Stream_Site_Q=dQ0$stream_site
Yr_Q=as.numeric(strftime(dQ0$date,format="%Y"))
Jwk_Q=as.numeric(strftime(dQ0$date,format="%V"))
dQ=cbind(dQ0,Stream_Site_Q,Yr_Q,Jwk_Q)

dMR=subset(dMR0,is.na(site)==F & is.na(yr)==F & Rel1>0)# & is.na(recaptured_fl)==F

Tribs=sort(unique(dC0$stream_site));Ntribs=length(Tribs)

#data frames for chinook of any kind and for spring run only (exclude hatchery in both cases)
dCa=subset(dC0,  origin!="hatchery" & is.na(match(stream_site,Tribs))==F & (Yr>=Fyr & Yr<=Lyr) & ((Jwk>=1 & Jwk<=Lwk) | (Jwk>=Fwk & Jwk<=53)) )
dCs=subset(dCa,run=="spring")

nrows=0#determine how many trib-year combinations with at least one weekly catch observation and only save these to RST_Input.txt
Process_TribYr=matrix(data=F,nrow=Ntribs,ncol=Nyrs)

for(itrib in 1:Ntribs){
	for(iyr in 1:Nyrs){
		dCax=subset(dCa, is.na(match(stream_site,Tribs[itrib]))==F & Yr==Year[iyr])
		if(dim(dCax)[1]>0){
			nrows=nrows+Nwks
			Process_TribYr[itrib,iyr]=T
		}
	}
}


Trib=character(length=nrows);Spr_Stage=Trib;rel_origin=Trib;rel_daynight=Trib
RunYr=vector(length=nrows);CalYr=RunYr;Weeks=RunYr
ua=RunYr;us=RunYr; Rel1=ua;Recap1=ua;Spr_Size=ua;Effort=ua;mr_flow=ua;catch_flow=ua

print("Funk_Record-Record_Num-Trib-CalYr-RunYr-Jwk where this is no catch of unmarked chinook but this is an efficiency trial (uh oh)")
jj=0
k=0
for(itrib in 1:Ntribs){
	dQ1=subset(dQ,Stream_Site_Q==Tribs[itrib])
	
	for(iyr in 1:Nyrs){
		
		if(Process_TribYr[itrib,iyr]==T){
						
			dCa_yr=subset(dCa, is.na(match(stream_site,Tribs[itrib]))==F & Yr==Year[iyr])
			dCs_yr=subset(dCs, is.na(match(stream_site,Tribs[itrib]))==F & Yr==Year[iyr])

			for(iwk in 1:Nwks){
				k=k+1
				
				dCa1=subset(dCa_yr, Jwk==Week[iwk])
				dCs1=subset(dCs_yr, Jwk==Week[iwk])

				#As only a subset of Tribs is in MR file (not all have MR data), this will grab up any tribs in common
				#But note if you have excluded some stream_site from MR (e.g., sacramento) they won't be grabbed up
				dMR1=subset(dMR,is.na(match(stream_site,Tribs[itrib]))==F & yr==Year[iyr] & Jwk==Week[iwk])
				
				Trib[k]=Tribs[itrib]

				#RunYr is the calendar year for Jwks 1:Lwk (e.g., 1:25). Note catch from 
				#Jwks Fwk:53 (e.g., 45:53) in yr t-1 are assigned to RunYr=t 
				#e.g. catch from Jwk 45:53 from calendar yr 2003 are assigned to RunYr 2004 along with Jwk 1:25 from calendar year 2004
				CalYr[k]=Year[iyr]
				RunYr[k]=Year[iyr]
				if(Week[iwk]>=Fwk & Week[iwk]<=53) RunYr[k]=Year[iyr]+1

				Weeks[k]=Week[iwk]
				ua[k]=NA; us[k]=NA; Spr_Stage[k]=NA; Spr_Size[k]=NA; Effort[k]=NA; catch_flow[k]=NA
				
				#Set average flow for week
				dQ2=subset(dQ1,Yr_Q==Year[iyr] & Jwk_Q==iwk)
				if(dim(dQ2)[1]>0) catch_flow[k]=mean(dQ2$flow_cfs,na.rm=T)

				NotSampled=F
				if(dim(dCa1)[1]==0) NotSampled=T #no record for this trib-site-yr-wk. If they didn't sample there is no record in the catch table. If they sampled but caught no fish u=0

				if(NotSampled==F){
					ua[k]=sum(dCa1$u)	#sum across potentially multiple run types or originsaces
					Effort[k]=dCa1$Eff[1]

					if(dim(dCs1)[1]==0){ #it was sampled but no spring run records, so 0 spring run catch
						us[k]=0
					} else {
						us[k]=sum(dCs1$u) #if more than one record usually origin=natural and origin=not recorded

						#Stage classificaiton based on forklength. Function provided by Flora C on May 23, 2022 
						sz=mean(dCs1$u_sz);Spr_Size[k]=sz
						if(is.na(sz)==F){
							if((Jmon[Weeks[iwk]]>9 & sz>50) | 
								(Jmon[Weeks[iwk]]<=2 & sz>60)|
								(Jmon[Weeks[iwk]]==3 & sz>76)|
								(Jmon[Weeks[iwk]]==4 & Jday[Weeks[iwk]]>=1 & Jday[Weeks[iwk]]<15 & sz>85)|
								(Jmon[Weeks[iwk]]==4 & Jday[Weeks[iwk]]>=15 & Jday[Weeks[iwk]]<30 & sz>95)|
								(Jmon[Weeks[iwk]]>=4 & Jmon[Weeks[iwk]]<7 & sz>100)){
								Spr_Stage[k]="Yearling"
							} else if ((Jmon[Weeks[iwk]]>=1 & Jmon[Weeks[iwk]]<3 & sz>45 & sz<=60) |
								 (Jmon[Weeks[iwk]]>=3 & Jmon[Weeks[iwk]]<7 & sz>45 & sz<=100)){
								Spr_Stage[k]="Smolt"
							} else if ((Jmon[Weeks[iwk]]>10 & sz<=45) | (Jmon[Weeks[iwk]]<=6 & sz<=45)){
								Spr_Stage[k]="Fry"
							}
						}
					}
				}

				Rel1[k]=NA;Recap1[k]=NA;mr_flow[k]=NA;rel_origin[k]=NA;rel_daynight[k]=NA
				if(dim(dMR1)[1]>0){
					Rel1[k]=dMR1$Rel1
					Recap1[k]=dMR1$Recap1
					if(is.nan(dMR1$flow_cfs)==F) mr_flow[k]=dMR1$flow_cfs
					rel_origin[k]=dMR1$rel_origin
					rel_daynight[k]=dMR1$rel_daynight
				}

				if(is.na(ua[k])==T & is.na(Rel1[k])==F){
					jj=jj+1
					print(paste(jj,k,Trib[k],CalYr[k],RunYr[k],Weeks[k],sep="-"))
				}

			}#wk
			
		}#if to see if any sampling for this year
	}#yr
}#trib
out=data.frame(cbind(Trib,CalYr,RunYr,Weeks,ua,us,Spr_Stage,Spr_Size,Rel1,Recap1,Effort,mr_flow,catch_flow,rel_origin,rel_daynight));names(out)=c("Trib","CalYr","RunYr","Week","ua","us","Spr_Stage","Spr_Sz","Rel1","Recap1","Effort","mr_flow","catch_flow","rel_origin","rel_daynight")
write.table(file=fnout,out,quote=1,row.names=F)
