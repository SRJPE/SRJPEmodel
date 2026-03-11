#Temporary routine to calculate mean and sd of cummulative abundance through forecast week For_ewk across years
#This is used to populate Nx_mu and Nx_sd in For_Inputs.csv temporarily

RSTpath="c:/projects/BayDelta/SAC_JPE/RST/RunSize/StanVersion/Output_Trib_SR"
RSTpathA="c:/projects/BayDelta/SAC_JPE/RST/RunSize/StanVersion/Output_Trib"

#wks for inseason model (entire year)
calwk=c('Sep-07','Sep-14','Sep-21','Sep-28','Oct-05','Oct-12','Oct-19','Oct-26','Nov-02','Nov-09',
        'Nov-16','Nov-23','Nov-30','Dec-07','Dec-14','Dec-21','Dec-28','Jan-04','Jan-11', 'Jan-18',
        'Jan-25','Feb-01','Feb-08','Feb-15','Feb-22','Mar-01','Mar-08','Mar-15','Mar-22','Mar-29',
        'Apr-05','Apr-12','Apr-19','Apr-26','May-03','May-10','May-17','May-24','May-31','Jun-07',
        'Jun-14','Jun-21','Jun-28','Jul-05','Jul-12','Jul-19','Jul-26','Aug-02','Aug-09','Aug-16',
        'Aug-23','Aug-30','Aug-31')


mwk=c(36:53,1:35)#the calendar week associated with each inseason wk (Sep-Aug)
#bt_wk=c(45:53,1:22) #the weeks btspas abundance estimates were computed for (Nov-May)

dfor=read.csv(file="For_Input.csv",header=T)
DoTrib=dfor$Trib
For_CalWk='Jan-04'
#For_CalWk='Feb-01'
#For_CalWk='Mar-01'
#For_CalWk='Apr-05'


Ntribs=dim(dfor)[1]
Nx_mu=vector(length=Ntribs);Nx_sd=Nx_mu;Nx_cv=Nx_mu

for(itrib in 1:Ntribs){
  For_ewk=which(calwk==For_CalWk)#the index for forecast week for inseason model period
  
  fn1=list.files(path=RSTpath,paste0("N_",DoTrib[itrib]))
  Nyrs=length(fn1)
  year=character(length=Nyrs)
  
  print(DoTrib[itrib])
  sum_mu=0;sum_cv=0
  for(iyr in 1:Nyrs){
    year[iyr]=substr(fn1[iyr],start=nchar(fn1[iyr])-9,stop=nchar(fn1[iyr])-6)
    fn_fit=paste0(RSTpath,"/N_",DoTrib[itrib],"_",year[iyr],".rdata")  
    load(fn_fit)
    
    fn_dat=paste0(RSTpathA,"/Inputs_",DoTrib[itrib],"_",year[iyr],".rdata")  
    load(fn_dat)
    Jwk=abundance_inputs$weeks_fit 
    
    
    #Cummulative abundance through forecast week
    btwk=which(Jwk==mwk[For_ewk])#the index for week in btspas timestep for the forecast week
    N=rowSums(abundance[,1:btwk])#for each posterior sample (row) calculate the sum of abundance from wk 1 through forecast wk
    sum_mu=sum_mu+mean(N)#the mean of cum abundandance across posterior samples
    sum_cv=sum_cv+sd(N)/mean(N)
  }
  
  Nx_mu[itrib]=sum_mu/Nyrs
  Nx_cv[itrib]=sum_cv/Nyrs
  Nx_sd[itrib]=Nx_mu[itrib]*Nx_cv[itrib]
  print(c(Nx_mu[itrib],Nx_cv[itrib],Nx_sd[itrib])) #manually entere these in For_Input.csv
}
