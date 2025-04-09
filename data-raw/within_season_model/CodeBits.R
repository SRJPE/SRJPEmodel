for(erun in 2:2){#do eruns without and with error in Nx
  NxError=switch(erun,0,1)
  niter=switch(erun,1000,5000)
  
  for(irun in 1:2){ #without and with covariate effect on phi
    
    if(irun==1){
      if(erun==1){
        OutDir="Output/NoNxError/NoCovEffect/"
      } else if (erun==2){
        OutDir="Output/NxError/NoCovEffect/"
      }
      UseCovX=c(0,0)
      
    } else if (irun==2){
      if(erun==1){
        OutDir="Output/NoNxError/WithCovEffect/"
      } else if (erun==2){
        OutDir="Output/NxError/WithCovEffect/"
      }
      UseCovX=c(1,0)
    }