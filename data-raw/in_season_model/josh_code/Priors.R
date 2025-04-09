#From R, demonstrates it produces thesame distribution as stan simulation

graphics.off()
if(exists(".SavedPlots",where=1)==T){rm(.SavedPlots,pos=1)}
windows(record=T,width=12,height=12);par(las=1)

par(mfrow=c(2,2),mai=c(0.5,0.2,0.5,0.2), omi=c(0.5,0.5,0.5,0.5),cex.main=1.1,cex.lab=1.2,font.lab=2)
#par(mfrow=c(1,1))
for(i in 1:4){
  xshape=switch(i,5,10,20,10)
  xrate=switch(i,2.5,5,10,20)
  
  
  y=rgamma(n=1000,shape=xshape,rate=xrate)
  hist(y,main=c(xshape,xrate))
  
  #https://en.wikipedia.org/wiki/Gamma_distribution
  #mean=shape/rate, variance=shape/rate^2
  print(c(xshape,xrate))
  print("mean, sd, cv");print(c(xshape/xrate,sqrt(xshape/xrate^2),1/sqrt(xshape)))
  print(quantile(y,probs=c(0.025,0.5,0.975)))
  print("")
}