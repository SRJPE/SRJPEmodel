
fit <- battle_sr
params <- battle_sr_params
inputs <- sr_inputs
CoVarNm <- "spawning_above_13_temp_week"

# calculate covariate effects
# TODO can we add this code into the mutate?
pred_covar <- c()
pred_covar_w_gamma <- c()
for(i in 1:sr_inputs$data$Nyrs) {
  pred_covar[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i])) * 0.001
  pred_covar_w_gamma[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i] +
                                                                   posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
}

obsv_data <- obsv_data |>
  mutate(pred_covar = pred_covar,
         pred_covar_w_gamma = pred_covar_w_gamma,
         covar_effect_color = ifelse(pred_covar_w_gamma > pred_covar & (obsv_recruits > pred_covar) | pred_covar_w_gamma < pred_covar & (obsv_recruits < pred_covar),
                                     "covariate effect (consistent)",
                                     "covariate effect (inconsistent)"))


# clean code --------------------------------------------------------------

# pull raw covariate
# TODO confirm these are being updated, can we have them added to df?
raw_covar <- SRJPEdata::stock_recruit_covariates |>
  ungroup() |>
  mutate(covar_nm = ifelse(lifestage == "spawning and incubation",
                           paste0("spawning_", covariate_structure),
                           paste0(lifestage, "_", covariate_structure))) |>
  filter(covar_nm == CoVarNm,
         stream == sr_inputs$stream,
         year %in% sr_inputs$year_lookup$brood_year) |>
  pull(value)

# create covariate vector across range of values
covar_vector <- seq(min(raw_covar), max(raw_covar), length.out = 50)
std_covar_vector <- (covar_vector - mean(covar_vector)) / sd(covar_vector)
posteriors <- as.data.frame(fit, pars = c("alpha", "beta", "gamma", "sd_pro"))
pred_recruitment_matrix <- matrix(nrow = 50, ncol = 3)
mean_spawners <- mean(sr_inputs$data$SP)

for(i in 1:50) {
  pred_recruitment <- mean_spawners * exp(posteriors$alpha + posteriors$beta * mean_spawners +
                                            posteriors$gamma * std_covar_vector[i]) * 0.001 # scale
  pred_recruitment_matrix[i, ] <- quantile(pred_recruitment, probs = c(0.025, 0.5, 0.975))
  pred_recruitment_matrix[i, 2] <- mean(pred_recruitment) # replace median with mean
}


# create observed data frame that also includes effect of covariate size
# estimated using the posteriors
obsv_data <- tibble("obsv_spawners" = sr_inputs$data$SP,
                    "obsv_recruits" = sr_inputs$data$R * 0.001,
                    "obsv_raw_covar" = raw_covar,
                    "obsv_covar" = sr_inputs$data$X,
                    "brood_year" = paste0("\'", substr(sr_inputs$year_lookup$brood_year, 3, 4)))

# calculate covariate effects
pred_covar_1 <- c()
pred_covar_2 <- c()
for(i in 1:sr_inputs$data$Nyrs) {
  pred_covar_1[i] <- mean(mean_spawners * exp(posteriors$alpha + posteriors$beta *
                                                mean_spawners + posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
  pred_covar_2[i] <- mean(obsv_data$obsv_spawners[i] * exp(posteriors$alpha + posteriors$beta * obsv_data$obsv_spawners[i] +
                                                                   posteriors$gamma * obsv_data$obsv_covar[i])) * 0.001
}

obsv_data <- obsv_data |>
  mutate(pred_covar_1 = pred_covar_1,
         pred_covar_2 = pred_covar_2)



# create data frame with a range of spawning stock sizes and
# associated predicted recruitment across that range
pred_data <- tibble("covar_vector" = covar_vector,
                    "pred_recruit_025" = pred_recruitment_matrix[, 1],
                    "pred_recruit_mean" = pred_recruitment_matrix[, 2],
                    "pred_recruit_975" = pred_recruitment_matrix[, 3]) |>
  mutate(color = "@ average cov value")

ggplot() +
  # range of covar values
  geom_ribbon(data = pred_data,
              aes(x = covar_vector, ymin = pred_recruit_025, ymax = pred_recruit_975),
              alpha = 0.6, fill = "grey") +
  geom_line(data = pred_data,
            aes(x = covar_vector, y = pred_recruit_mean,
                color = "@ average cov value")) +
  # observed data
  geom_errorbar(data = obsv_data,
                aes(x = obsv_raw_covar,
                    ymin = pred_covar_1,
                    ymax = pred_covar_2,
                    color = "with spawner effect"),
                width = 0) +
  geom_point(data = obsv_data,
             aes(x = obsv_raw_covar, y = obsv_recruits))  +
  geom_text(data = obsv_data,
            aes(x = obsv_raw_covar + 0.5, y = obsv_recruits + 0.5, # TODO this "jigger" will be different for different covariates
                label = brood_year),
            size = 3) +
  labs(x = "", # TODO what is this?
       y = "Outmigrant abundance ('000s)") +
  theme_minimal() +
  scale_color_manual(name = "",
                     breaks = c("@ average cov value", "with spawner effect"),
                     values = c("@ average cov value" = "black",
                                "with spawner effect" = "red")) +
  theme(legend.position = "bottom")




# to refactor ---------------------------------------------------------------

one_plot <- FALSE
# TODO can we get the raw covariate, confirm this table is being updated


# if(!one_plot) {
#   # plot covariate relationship
#   if(!CoVarNM %in% c("mu", "null")) {
#
#
#     plot <- ggplot()
#
#   }
# }

  if(OnePlot==F){
    #Plot covariate relationshp
    if(CoVarNm!="mu" & CoVarNm!="null"){

      #ymax=max(c(stat[,2],R*Nscale))*1.25
      #plot(X0,R*Nscale,pch=19,cex=0.75,xlim=c(min(x),max(x)),ylim=c(0,ymax),bty='l',xlab=CoVarNm,ylab="")
      #polygon(x = c(rev(x), x), y = c(rev(stat[,1]), stat[,3]), col = alpha("grey", 0.25), border = NA)
      #lines(x,stat[,2],lty=1,lwd=2)
      #xoff=(max(X0)-min(X0))*0.025
      for(iyr in 1:Nyrs){
        #text(x=X0[iyr]+xoff,y=R[iyr]*Nscale,labels=yr[iyr],cex=0.75)
        p1=mean(muS*exp(alpha+beta*muS + gamma*X[iyr]))*Nscale
        p2=mean(SP[iyr]*exp(alpha+beta*SP[iyr] + gamma*X[iyr]))*Nscale
        lines(x=c(X0[iyr],X0[iyr]),y=c(p1,p2),col="red")
      }
      legend("topright",legend=c("@ avg spawners", "with spawner effect"),lty=c(1,1),lwd=c(2,2),col=c("black","red"),bty='n')
      graph_let(2)
      mtext("Outmigrant abundance ('000s)",side=2, line=-1,las=3,outer=T,cex=1.3,font=2)
      mtext(paste0(DoSite,"_",Am),side=3, line=-1,outer=T,cex=1.3,font=2)

      #Stock-recruit at lowest and highest covariate values
      stat=array(dim=c(2,NxS,3))
      for(i in 1:NxS){
        pred=xS[i]*exp(alpha + beta*xS[i] + gamma*min(X))*Nscale
        stat[1,i,]=quantile(pred,probs=c(0.025,0.5,0.975));  stat[1,i,2]=mean(pred)
        pred=xS[i]*exp(alpha + beta*xS[i] + gamma*max(X))*Nscale
        stat[2,i,]=quantile(pred,probs=c(0.025,0.5,0.975));stat[2,i,2]=mean(pred)
      }
      plot(xS,stat[1,,2],type='l',lwd=2,col="blue",bty='n',ylim=c(0,max(stat[,,3])),xlab="Spawner Abundance",ylab="")
      polygon(x = c(rev(xS), xS), y = c(rev(stat[1,,1]), stat[1,,3]), col = alpha("blue", 0.25), border = NA)
      lines(xS,stat[2,,2],col="red")
      polygon(x = c(rev(xS), xS), y = c(rev(stat[2,,1]), stat[2,,3]), col = alpha("red", 0.25), border = NA)
      legend("topleft",legend=c("min cov value", "max cov value"),lty=c(1,1),lwd=c(2,2),col=c("blue","red"),bty='n')
      graph_let(3)
    }

    #par(mfcol=c(2,1),mai=c(0.85,1,0.25,1),omi=c(0.1,0.1,0.1,0.1))

    #Juveniles with error time series and covariates on second axis
    statR=matrix(nrow=Nyrs,ncol=2)
    for(iyr in 1:Nyrs){
      simR=exp(rnorm(n=1000,mean=mu_obslgR[iyr],sd=sd_obslgR[iyr]))*Nscale
      statR[iyr,]=quantile(simR,probs=c(0.025,0.975))
    }

    x=barplot(R*Nscale,space=0.1,names.arg="",col="grey",bty='l',xlab="Brood Year",ylab="",axes=F,ylim=c(0,max(statR)),xlim=c(0,Nyrs+1))
    axis(2);axis(1,at=x,labels=yr,las=3,cex.axis=0.9)
    abline(h=mean(R)*0.001,lty=2)
    for(iyr in 1:Nyrs){
      arrows(x0=x[iyr],x1=x[iyr],y0=statR[iyr,1],y1=statR[iyr,2],angle=90,code=3,length=0.025)
    }

    if(CoVarNm!="mu" & CoVarNm!="null"){
      par(new=T)
      plot(x,X0,pch=19,cex=1.3,axes=F,ylab="",xlab="",ylim=c(min(X0)*0.95,max(X0)*1.15),xlim=c(0,Nyrs+1)); axis(4)
      mtext(CoVarNm, side=4, line=3,las=3)
      legend("top",legend=c("outmigrants","covariate"),pch=c(15,19),col=c("grey","black"),bty='y',ncol=2)
    }
    graph_let(4)

    mtext("Outmigrant abundance ('000s)",side=2, line=-1,las=3,outer=F,cex=1.3,font=2)
    mtext(paste0(DoSite,"_",Am),side=3, line=-1,outer=F,cex=1.3,font=2)
  }
  #Do a forecast assuming spawning stock is hitorical average spawner abundance, and covariate is at average historical conditon
  #Differentiate between parameter uncertainty and process error in plot
  par(mfcol=c(1,1),mai=c(1,1,0.25,0.1),omi=c(0.1,0.1,0.1,0.1))
  muSP=mean(SP)
  pred1=muSP*exp(alpha + beta*muSP)*Nscale
  nsims=length(pred1)
  dev=rnorm(n=nsims,mean=0,sd=mean(sd_pro))
  pred2=pred1*exp(dev)
  print("")
  print(c(DoSite," ",Am))
  print(quantile(pred1,prob=c(0.025,0.5,0.975)));print(sd(pred1)/mean(pred1))
  print(quantile(pred2,prob=c(0.025,0.5,0.975)));print(sd(pred2)/mean(pred2))
}
