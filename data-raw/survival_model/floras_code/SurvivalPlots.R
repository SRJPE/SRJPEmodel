# Figures for Bayesian models
library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(dplyr)
library(ggpubr)
library(viridis)#Color scale

source("scripts/GetData_flora.R")

# Figures for NoCov model --------------------------------------------
load(file="outputs/fit_NoCov.Rdata")
posterior_draws_NoCov <- rstan::extract(fit)
posterior.df_NoCov <-as.data.frame(posterior_draws_NoCov)
names(posterior.df_NoCov)

fit_summary_NoCov <- summary(fit,probs=c(0.025,0.975))# shows only 10% and 90% quantiles
data_summary_NoCov <- as.data.frame(fit_summary_NoCov$summary)
print(fit_summary_NoCov$summary)

### Detection probability figure
data_pcap_NoCov <- data.frame(matrix(NA, nrow = 4*length(Year_all) , ncol = 9))
pars <- c(paste0("pred_pcap[",1,",1]"),paste0("pred_pcap[",1,",2]"),
          paste0("pred_pcap[",1,",3]"),paste0("pred_pcap[",1,",4]"))
loc <- c('Woodson','Butte','Sacramento','Delta')
loc <- factor(loc, levels = unique(loc))

data_pcap_NoCov <- data.frame(cbind(loc,Year_all[1],data_summary_NoCov[pars,]))
colnames(data_pcap_NoCov) <- c('loc','year','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

for(i in 2:length(Year_all)){
pars <- c(paste0("pred_pcap[",i,",1]"),paste0("pred_pcap[",i,",2]"),
        paste0("pred_pcap[",i,",3]"),paste0("pred_pcap[",i,",4]"))

data_temp <- data.frame(cbind(loc,Year_all[i],data_summary_NoCov[pars,]))
names(data_temp) <- names(data_pcap_NoCov) # to be able to bind data_temp and data make sure the column names at the same

data_pcap_NoCov <- data.frame(rbind(data_pcap_NoCov,data_temp))

}


(pcap_nocov_fig <- ggplot(data_pcap_NoCov)+
    geom_errorbar(aes(x=loc, ymin=q2.5, ymax=q97.5, color=loc), 
                  width=0.2, size=1)+
    geom_point(aes(x=loc,y=mean, color=loc),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1))+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    xlab("")+ylab("Detection probability")+
    facet_grid(~year)
)

ggsave(plot=pcap_nocov_fig , "figures/pcap.png", width =12, height = 6)

### Survival rates figures
# Release-Sacramento Survival
data_surv_relsac_novcov <- data.frame(matrix(NA, nrow = dim(rgwy2.df)[1] , ncol = 8))

for(i in 1: dim(rgwy2.df)[1]){
  data_surv_relsac_novcov [i,] <- data.frame(cbind(rgwy2.df$StudyID[i],
                                                    data_summary_NoCov[paste0('SurvRelSac[',i,"]"),]))
  
}

colnames(data_surv_relsac_novcov ) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_relsac_novcov$studyid <- factor(data_surv_relsac_novcov$studyid, 
                                           levels = unique(data_surv_relsac_novcov$studyid))

(surv_relsac_nocov_Fig <- ggplot(data_surv_relsac_novcov)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5), 
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=studyid,y=mean),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.6), breaks=seq(0,0.6,0.1))+
    xlab("")+ylab("Survival rate")
)

ggsave(plot=surv_relsac_nocov_Fig, "figures/SurvRelSac_Nocov.png", width =15, height = 8)

# Butte and Feather-Sacramento Survival
data_surv_tribsac_novcov <- data.frame(matrix(NA, nrow = dim(rgwy2_T.df)[1] , ncol = 8))

for(i in 1: dim(rgwy2_T.df)[1]){
  data_surv_tribsac_novcov[i,] <- data.frame(cbind(rgwy2_T.df$StudyID[i],
                                                   data_summary[paste0('pred_survT[',i,",1]"),]))
  
}

colnames(data_surv_tribsac_novcov) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_tribsac_novcov$studyid <- factor(data_surv_tribsac_novcov$studyid, 
                                           levels = unique(data_surv_tribsac_novcov$studyid))

(surv_tribsac_nocov_Fig <- ggplot(data_surv_tribsac_novcov)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5), 
                  width=0.2, size=1, color="blue")+
    geom_point(aes(x=studyid,y=mean),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.7,0.1))+
    xlab("")+ylab("Survival rate")
)

ggsave(plot=surv_tribsac_nocov_Fig, "figures/SurvTribSac_Nocov.png", width =15, height = 8)


# Figures for WY2 model --------------------------------------------
load(file="outputs/fit_CovWY2.Rdata")
posterior_draws_CovWY2 <- rstan::extract(fit)
posterior.df_CovWY2 <-as.data.frame(posterior_draws_CovWY2)
names(posterior.df_CovWY2)

fit_summary_CovWY2 <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY2 <- as.data.frame(fit_summary_CovWY2$summary)
print(fit_summary_CovWY2$summary)

### Survival rates figures
# Release-Sacramento Survival
data_surv_relsac_wy2 <- data.frame(matrix(NA, nrow = dim(rgwy2.df)[1] , ncol = 8))

for(i in 1: dim(rgwy2.df)[1]){
  data_surv_relsac_wy2[i,] <- data.frame(cbind(rgwy2.df$StudyID[i],data_summary_CovWY2[paste0('SurvRelSac[',i,"]"),]))
  
}

colnames(data_surv_relsac_wy2) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_relsac_wy2$studyid <- factor(data_surv_relsac_wy2$studyid, 
                                       levels = unique(data_surv_relsac_wy2$studyid))

data_surv_relsac_wy2$wytype <- c('C/D/BN','C/D/BN','C/D/BN','C/D/BN','AN/W','AN/W','AN/W','C/D/BN','C/D/BN',
                                 'AN/W','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN',
                                 'AN/W','AN/W')
data_surv_relsac_wy2$wytype <- factor(data_surv_relsac_wy2$wytype, 
                                      levels = unique(data_surv_relsac_wy2$wytype))

(surv_relsac_wy2_Fig <- ggplot(data= data_surv_relsac_wy2)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_relsac_wy2_Fig, "figures/SurvRelSac_wy2.png", width =8, height = 5)

# Woodson-Sacramento Survival
data_surv_woodsac_wy2 <- data.frame(matrix(NA, nrow = dim(rgwy2.df)[1] , ncol = 8))

for(i in 1: dim(rgwy2.df)[1]){
  data_surv_woodsac_wy2[i,] <- data.frame(cbind(rgwy2.df$StudyID[i],data_summary_CovWY2[paste0('SurvWoodSac[',i,"]"),]))
  
}

colnames(data_surv_woodsac_wy2) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_woodsac_wy2$studyid <- factor(data_surv_woodsac_wy2$studyid, 
                                        levels = unique(data_surv_woodsac_wy2$studyid))

data_surv_woodsac_wy2$wytype <- c('C/D/BN','C/D/BN','C/D/BN','C/D/BN','AN/W','AN/W','AN/W','C/D/BN','C/D/BN',
                                 'AN/W','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN',
                                 'AN/W','AN/W')
data_surv_woodsac_wy2$wytype <- factor(data_surv_woodsac_wy2$wytype, 
                                       levels = unique(data_surv_woodsac_wy2$wytype))

(surv_woodsac_wy2_Fig <- ggplot(data= data_surv_woodsac_wy2)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.6), breaks=seq(0,0.6,0.1))+
    xlab("")+ylab("Woodson - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_woodsac_wy2_Fig, "figures/SurvWoodSac_wy2.png", width =8, height = 5)

# Butte and Feather-Sacramento Survival
data_surv_tribsac_wy2 <- data.frame(matrix(NA, nrow = dim(rgwy2_T.df)[1] , ncol = 8))

for(i in 1: dim(rgwy2_T.df)[1]){
  data_surv_tribsac_wy2[i,] <- data.frame(cbind(rgwy2_T.df$StudyID[i],
                                                data_summary[paste0('pred_survT[',i,",1]"),]))
  
}

colnames(data_surv_tribsac_wy2) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_tribsac_wy2$studyid <- factor(data_surv_tribsac_wy2$studyid, 
                                        levels = unique(data_surv_tribsac_wy2$studyid))

data_surv_tribsac_wy2$wytype <- c('C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','AN/W','C/D/BN','AN/W',
                                  'AN/W','AN/W','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','AN/W','AN/W',
                                  'AN/W')
data_surv_tribsac_wy2$wytype <- factor(data_surv_tribsac_wy2$wytype, 
                                       levels = unique(data_surv_tribsac_wy2$wytype))

(surv_tribsac_wy2_Fig <- ggplot(data_surv_tribsac_wy2)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.7,0.1))+
    xlab("")+ylab("Release -Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_tribsac_wy2_Fig, "figures/SurvTribSac_wy2.png", width =8, height = 5)

# forecast figures
data_survforecast_wy2 <- data.frame(rbind(data_summary_CovWY2['SurvForecast[1]',],
                                          data_summary_CovWY2['SurvForecast[2]',],
                                          data_summary_CovWY2['SurvForecast[3]',],
                                          data_summary_CovWY2['TribSurvForecast[1,1]',],
                                          data_summary_CovWY2['TribSurvForecast[1,2]',],
                                          data_summary_CovWY2['TribSurvForecast[2,1]',],
                                          data_summary_CovWY2['TribSurvForecast[2,2]',],
                                          data_summary_CovWY2['TribSurvForecast[3,1]',],
                                          data_summary_CovWY2['TribSurvForecast[3,2]',]))
colnames(data_survforecast_wy2) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_survforecast_wy2$loc <- c('Sacramento','Sacramento','Sacramento','Butte','Feather','Butte','Feather',
                               'Butte','Feather')
data_survforecast_wy2$loc <- factor(data_survforecast_wy2$loc,
                                    levels = unique(data_survforecast_wy3$loc))


data_survforecast_wy2$wytype <- c('Critical','D/BN/AN','Wet','Critical','Critical','D/BN/AN','D/BN/AN',
                                  'Wet','Wet')
data_survforecast_wy2$wytype <- factor(data_survforecast_wy2$wytype,
                                       levels = unique(data_survforecast_wy2$wytype))

(survforecast_wy2_Fig <- ggplot(data_survforecast_wy2)+
    geom_errorbar(aes(x=wytype, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=wytype,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.9), breaks=seq(0,0.9,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")+
    facet_wrap(~loc)
)

ggsave(plot=survforecast_wy2_Fig, "figures/SurvForecast_wy2.png", width =10, height = 6)


# Figures for WY3 model --------------------------------------------
load(file="outputs/fit_CovWY3.Rdata")
posterior_draws_CovWY3 <- rstan::extract(fit)
posterior.df_CovWY3 <-as.data.frame(posterior_draws_CovWY3)
names(posterior.df_CovWY3)

fit_summary_CovWY3 <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY3 <- as.data.frame(fit_summary_CovWY3$summary)
print(fit_summary_CovWY3$summary)

### Release-Sacramento Survival figure
# Collect survival estimates and stats from model posterior
data_surv_relsac_wy3 <- data.frame(matrix(NA, nrow = dim(rgwy3.df)[1] , ncol = 8))
for(i in 1: dim(rgwy3.df)[1]){
  data_surv_relsac_wy3[i,] <- data.frame(cbind(rgwy3.df$StudyID[i],data_summary_CovWY3[paste0('SurvRelSac[',i,"]"),]))
  
}

colnames(data_surv_relsac_wy3) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_relsac_wy3$studyid <- factor(data_surv_relsac_wy3$studyid, 
                                        levels = unique(data_surv_relsac_wy3$studyid))
# Add water year type info
data_surv_relsac_wy3$wytype <- c('D/BN/AN','D/BN/AN','Critical','D/BN/AN','Wet','Wet','Wet','D/BN/AN','D/BN/AN',
                                  'Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Critical',
                                  'Wet','Wet')
data_surv_relsac_wy3$wytype <- factor(data_surv_relsac_wy3$wytype, 
                                       levels = unique(data_surv_relsac_wy3$wytype))

# Add fish size and detection info
data_surv_relsac_wy3 <- data_surv_relsac_wy3 %>% 
                        left_join(d_Sac_summary, by='studyid')

(surv_relsac_wy3_Fig <- ggplot(data= data_surv_relsac_wy3)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_relsac_wy3_Fig, "figures/SurvRelSac_wy3.png", width =9, height =7)

### Butte and Feather release to Sacramento Survival figure
data_surv_tribsac_wy3 <- data.frame(matrix(NA, nrow = dim(rgwy3_T.df)[1] , ncol = 8))

for(i in 1: dim(rgwy3_T.df)[1]){
  data_surv_tribsac_wy3[i,] <- data.frame(cbind(rgwy3_T.df$StudyID[i],
                                                   data_summary_CovWY3[paste0('pred_survT[',i,",1]"),]))
  
}

colnames(data_surv_tribsac_wy3) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_tribsac_wy3$studyid <- factor(data_surv_tribsac_wy3$studyid, 
                                           levels = unique(data_surv_tribsac_wy3$studyid))

data_surv_tribsac_wy3$wytype <- c('D/BN/AN','Critical','Critical','Critical','D/BN/AN','Wet','D/BN/AN','Wet',
                                  'Wet','Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Wet','Wet',
                                  'Wet')
data_surv_tribsac_wy3$wytype <- factor(data_surv_tribsac_wy3$wytype, 
                                       levels = unique(data_surv_tribsac_wy3$wytype))

data_surv_tribsac_wy3 <- data_surv_tribsac_wy3 %>% 
                        left_join(d_FeaBut_summary, by='studyid')

(surv_tribsac_wy3_Fig <- ggplot(data_surv_tribsac_wy3)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.7,0.1))+
    xlab("")+ylab("Release -Sac survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_tribsac_wy3_Fig, "figures/SurvTribSac_wy3.png", width =8, height = 6)

### Survival vs fish size figures
(sacsurvvssize_wy3_Fig <- ggplot(data= data_surv_relsac_wy3)+
    geom_point(aes(x=Mean_Wt,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
    scale_y_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))+
    xlab("Mean FL")+ylab("Release - Sac survival rate")+
    ggtitle("upper Sac")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

(tribsurvvssize_wy3_Fig <- ggplot(data= data_surv_tribsac_wy3)+
    geom_point(aes(x=Mean_Wt,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+
  # scale_y_continuous(limits=c(0,0.4), breaks=seq(0,0.4,0.1))+
    xlab("Mean FL")+ylab("")+ ggtitle("Butte-Feather")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)


(SurvSize_CovWY3_fig <- ggarrange(sacsurvvssize_wy3_Fig ,tribsurvvssize_wy3_Fig,
                                 ncol=2, common.legend=TRUE,legend="right",
                                # labels=c('A','B'),
                                 align="h")
)

ggsave(plot=SurvSize_CovWY3_fig, "figures/SurvSize_CovWY3.png", width =11, height =4)


### forecast survival figure
data_survforecast_wy3 <- data.frame(rbind(data_summary_CovWY3['SurvForecast[1]',],
                                          data_summary_CovWY3['SurvForecast[2]',],
                                          data_summary_CovWY3['SurvForecast[3]',],
                                          data_summary_CovWY3['TribSurvForecast[1,1]',],
                                          data_summary_CovWY3['TribSurvForecast[1,2]',],
                                          data_summary_CovWY3['TribSurvForecast[2,1]',],
                                          data_summary_CovWY3['TribSurvForecast[2,2]',],
                                          data_summary_CovWY3['TribSurvForecast[3,1]',],
                                          data_summary_CovWY3['TribSurvForecast[3,2]',]))
colnames(data_survforecast_wy3) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_survforecast_wy3$loc <- c('Sacramento','Sacramento','Sacramento','Butte','Feather','Butte','Feather',
                               'Butte','Feather')
data_survforecast_wy3$loc <- factor(data_survforecast_wy3$loc,
                                    levels = unique(data_survforecast_wy3$loc))


data_survforecast_wy3$wytype <- c('Critical','D/BN/AN','Wet','Critical','Critical','D/BN/AN','D/BN/AN',
                                  'Wet','Wet')
data_survforecast_wy3$wytype <- factor(data_survforecast_wy3$wytype,
                                       levels = unique(data_survforecast_wy3$wytype))

(survforecast_wy3_Fig <- ggplot(data_survforecast_wy3)+
    geom_errorbar(aes(x=wytype, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=wytype,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.9), breaks=seq(0,0.9,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")+
    facet_wrap(~loc)
)

ggsave(plot=survforecast_wy3_Fig, "figures/SurvForecast_wy3.png", width =10, height = 6)


# Figures for WY2 + FL model --------------------------------------------
load(file="outputs/fit_CovWY2_FL.Rdata")
posterior_draws_CovWY2_FL <- rstan::extract(fit)
posterior.df_CovWY2_FL <-as.data.frame(posterior_draws_CovWY2_FL)
names(posterior.df_CovWY2_FL)

fit_summary_CovWY2_FL <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY2_FL <- as.data.frame(fit_summary_CovWY2_FL$summary)
print(fit_summary_CovWY2_FL$summary)

# Figures for WY2 + Wgt model --------------------------------------------
load(file="outputs/fit_CovWY2_Wgt.Rdata")
posterior_draws_CovWY2_Wgt <- rstan::extract(fit)
posterior.df_CovWY2_Wgt <-as.data.frame(posterior_draws_CovWY2_Wgt)
names(posterior.df_CovWY2_Wgt)

fit_summary_CovWY2_Wgt <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY2_Wgt <- as.data.frame(fit_summary_CovWY2_Wgt$summary)
print(fit_summary_CovWY2_Wgt$summary)

# Figures for WY2 + CF model --------------------------------------------
load(file="outputs/fit_CovWY2_CF.Rdata")
posterior_draws_CovWY2_CF <- rstan::extract(fit)
posterior.df_CovWY2_CF <-as.data.frame(posterior_draws_CovWY2_CF)
names(posterior.df_CovWY2_CF)

fit_summary_CovWY2_CF <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY2_CF <- as.data.frame(fit_summary_CovWY2_CF$summary)
print(fit_summary_CovWY2_CF$summary)

# Figures for WY3 + FL model --------------------------------------------
load(file="outputs/fit_CovWY3_FL.Rdata")
posterior_draws_CovWY3_FL <- rstan::extract(fit)
posterior.df_CovWY3_FL <-as.data.frame(posterior_draws_CovWY3_FL)
names(posterior.df_CovWY3_FL)

fit_summary_CovWY3_FL <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY3_FL <- as.data.frame(fit_summary_CovWY3_FL$summary)
print(fit_summary_CovWY3_FL$summary)

### Detection probability figures
data_pcap_CovWY3_FL <- data.frame(matrix(NA, nrow = 4*length(Year_all) , ncol = 9))
pars <- c(paste0("pred_pcap[",1,",1]"),paste0("pred_pcap[",1,",2]"),
          paste0("pred_pcap[",1,",3]"),paste0("pred_pcap[",1,",4]"))
loc <- c('Woodson','Butte','Sacramento','Delta')
loc <- factor(loc, levels = unique(loc))

data_pcap_CovWY3_FL <- data.frame(cbind(loc,Year_all[1],data_summary_CovWY3_FL[pars,]))
colnames(data_pcap_CovWY3_FL) <- c('loc','year','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

for(i in 2:length(Year_all)){
  pars <- c(paste0("pred_pcap[",i,",1]"),paste0("pred_pcap[",i,",2]"),
            paste0("pred_pcap[",i,",3]"),paste0("pred_pcap[",i,",4]"))
  
  data_temp <- data.frame(cbind(loc,Year_all[i],data_summary_CovWY3_FL[pars,]))
  names(data_temp) <- names(data_pcap_CovWY3_FL) # to be able to bind data_temp and data make sure the column names at the same
  
  data_pcap_CovWY3_FL <- data.frame(rbind(data_pcap_CovWY3_FL,data_temp))
  
}


(pcap_CovWY3_FL_fig <- ggplot(data_pcap_CovWY3_FL)+
    geom_errorbar(aes(x=loc, ymin=q2.5, ymax=q97.5, color=loc), 
                  width=0.2, size=1)+
    geom_point(aes(x=loc,y=mean, color=loc),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1))+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    xlab("")+ylab("Detection probability")+
    facet_grid(~year)
)

ggsave(plot=pcap_CovWY3_FL_fig , "figures/pcap_CovWY3_FL.png", width =10, height = 4)

### Release-Sacramento Survival figure
data_surv_relsac_CovWY3_FL <- data.frame(matrix(NA, nrow = dim(rgwy3.df)[1] , ncol = 8))

for(i in 1: dim(rgwy3.df)[1]){
  data_surv_relsac_CovWY3_FL[i,] <- data.frame(cbind(rgwy3.df$StudyID[i],
                                                      data_summary_CovWY3_FL[paste0('SurvRelSac[',i,"]"),]))
  
}

colnames(data_surv_relsac_CovWY3_FL) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_relsac_CovWY3_FL$studyid <- factor(data_surv_relsac_CovWY3_FL$studyid, 
                                              levels = unique(data_surv_relsac_CovWY3_FL$studyid))

data_surv_relsac_CovWY3_FL$wytype <- c('D/BN/AN','D/BN/AN','Critical','D/BN/AN','Wet','Wet','Wet','D/BN/AN','D/BN/AN',
                                        'Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Critical',
                                        'Wet','Wet')
data_surv_relsac_CovWY3_FL$wytype <- factor(data_surv_relsac_CovWY3_FL$wytype, 
                                             levels = unique(data_surv_relsac_CovWY3_FL$wytype))

(surv_relsac_CovWY3_FL_fig <- ggplot(data= data_surv_relsac_CovWY3_FL)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_relsac_CovWY3_FL_fig, "figures/SurvRelSac_CovWY3_FL.png", width =7, height = 6)

### Butte and Feather-Sacramento Survival figure
data_surv_tribsac_CovWY3_FL <- data.frame(matrix(NA, nrow = dim(rgwy3_T.df)[1] , ncol = 8))

for(i in 1: dim(rgwy3_T.df)[1]){
  data_surv_tribsac_CovWY3_FL[i,] <- data.frame(cbind(rgwy3_T.df$StudyID[i],
                                                       data_summary_CovWY3_FL[paste0('pred_survT[',i,",1]"),]))
  
}

colnames(data_surv_tribsac_CovWY3_FL) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_tribsac_CovWY3_FL$studyid <- factor(data_surv_tribsac_CovWY3_FL$studyid, 
                                               levels = unique(data_surv_tribsac_CovWY3_FL$studyid))

data_surv_tribsac_CovWY3_FL$wytype <- c('D/BN/AN','Critical','Critical','Critical','D/BN/AN','Wet','D/BN/AN','Wet',
                                         'Wet','Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Wet','Wet',
                                         'Wet')
data_surv_tribsac_CovWY3_FL$wytype <- factor(data_surv_tribsac_CovWY3_FL$wytype, 
                                              levels = unique(data_surv_tribsac_CovWY3_FL$wytype))

(surv_tribsac_CovWY3_FL_fig <- ggplot(data_surv_tribsac_CovWY3_FL)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.6), breaks=seq(0,0.6,0.1))+
    xlab("")+ylab("")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_tribsac_CovWY3_FL_fig, "figures/SurvTribSac_CovWY3_FL.png", width =7, height = 5)


(Surv_CovWY3_FL_fig <- ggarrange(surv_relsac_CovWY3_FL_fig ,surv_tribsac_CovWY3_FL_fig,
                                  ncol=2, common.legend=TRUE,
                                  labels=c('A','B'),
                                  align="h")
)

ggsave(plot=Surv_CovWY3_FL_fig, "figures/Surv_CovWY3_FL.png", width =9, height =6)


### WY type forecast figures
data_survforecast_wy3_FL <- data.frame(rbind(data_summary_CovWY3_FL['SurvForecast[1]',],
                                          data_summary_CovWY3_FL['SurvForecast[2]',],
                                          data_summary_CovWY3_FL['SurvForecast[3]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[1,1]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[1,2]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[2,1]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[2,2]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[3,1]',],
                                          data_summary_CovWY3_FL['TribSurvForecast[3,2]',]))
colnames(data_survforecast_wy3_FL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_survforecast_wy3_FL$loc <- c('Sacramento','Sacramento','Sacramento','Butte','Feather','Butte','Feather',
                               'Butte','Feather')
data_survforecast_wy3_FL$loc <- factor(data_survforecast_wy3_FL$loc,
                                    levels = unique(data_survforecast_wy3_FL$loc))

data_survforecast_wy3_FL$wytype <- c('Critical','D/BN/AN','Wet','Critical','Critical','D/BN/AN','D/BN/AN',
                                  'Wet','Wet')
data_survforecast_wy3_FL$wytype <- factor(data_survforecast_wy3_FL$wytype,
                                       levels = unique(data_survforecast_wy3_FL$wytype))

(survforecast_wy3_FL_Fig <- ggplot(data_survforecast_wy3_FL)+
    geom_errorbar(aes(x=wytype, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=wytype,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.9), breaks=seq(0,0.9,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")+
    facet_wrap(~loc)
)

ggsave(plot=survforecast_wy3_FL_Fig, "figures/SurvForecast_wy3_FL.png", width =10, height = 6)

### FL forecast figures
FL_vec  <- seq(from=round(min(c(FL,FL_T),na.rm=TRUE)),to=round(max(c(FL,FL_T),na.rm=TRUE)),length.out=Nsz)

data_survforecast_CovWY3_FL_Sac <- array(data = NA, dim = c(Nsz,3,7),
                                          dimnames =list(Size = FL_vec, WY = c('C','D/BN/AN','W'),var= c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')))

data_survforecast_CovWY3_FL_Trib <- array(data = NA, dim = c(2,Nsz,3,7),
                                           dimnames =list(Trib = c('Butte','Feather'),Size = FL_vec, WY = c('C','D/BN/AN','W'),var= c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')))

for(i in 1:Nsz){ # number of class size
  for (j in 1:3){ # Number of water year type 
    for (k in 1:7){ # number of variables saved in data frame
      data_survforecast_CovWY3_FL_Sac[i,j,k] <- data_summary_CovWY3_FL[paste0('SurvForecastSz[',i,",",j,"]"),k]
      for (l in 1:2){ # tributary index: 1 = Buttte, 2 = Feather
        data_survforecast_CovWY3_FL_Trib[l,i,j,k] <- data_summary_CovWY3_FL[paste0('TribSurvForecastSz[',i,",",j,",",l,"]"),k]
      }
    }
  }
}

# Sacramento fish
(survforecast_CovWY3_C_FL_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_FL_Sac[,1,],aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_FL_Sac[,1,],
                aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
     scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('upper Sac - C Year')
)

(survforecast_CovWY3_D_FL_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_FL_Sac[,2,],aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_FL_Sac[,2,],
                aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("")+
    ggtitle('upper Sac - D/BN/AN Year')
)

(survforecast_CovWY3_W_FL_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_FL_Sac[,3,],aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_FL_Sac[,3,],
                aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
     scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("")+
    ggtitle('upper Sac - W Year')
)

(survforecast_CovWY3_FL_Sac_fig <- ggarrange(survforecast_CovWY3_C_FL_Sac_fig,survforecast_CovWY3_D_FL_Sac_fig,
                                             survforecast_CovWY3_W_FL_Sac_fig,
                                  ncol=3, common.legend=TRUE,
                                # labels=c('A','B','C'),
                                 align="h")
)

ggsave(plot=survforecast_CovWY3_FL_Sac_fig, "figures/SurvForecast_CovWY3_FL_Sac.png", width =6, height = 4)

# Butte Creek fish
(survforecast_CovWY3_C_FL_But_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[1,,1,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - C Year')
)

(survforecast_CovWY3_D_FL_But_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[1,,2,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - D/BN/AN Year')
)

(survforecast_CovWY3_W_FL_But_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[1,,3,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - W Year')
)

(survforecast_CovWY3_FL_But_fig <- ggarrange(survforecast_CovWY3_C_FL_But_fig,survforecast_CovWY3_D_FL_But_fig,
                                             survforecast_CovWY3_W_FL_But_fig,
                                             ncol=3, common.legend=TRUE,
                                             # labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_FL_But_fig, "figures/SurvForecast_CovWY3_FL_But.png", width =6, height = 4)

# Feather River fish
(survforecast_CovWY3_C_FL_Fea_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[2,,1,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - C Year')
)

(survforecast_CovWY3_D_FL_Fea_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[2,,2,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - D/BN/AN Year')
)

(survforecast_CovWY3_W_FL_Fea_fig <- ggplot(data=data_survforecast_CovWY3_FL_Trib[2,,3,])+
    geom_line(aes(x=FL_vec,y=mean),size=1)+
    geom_ribbon(aes(x=FL_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("FL (mm)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - W Year')
)

(survforecast_CovWY3_FL_Fea_fig <- ggarrange(survforecast_CovWY3_C_FL_Fea_fig,survforecast_CovWY3_D_FL_Fea_fig,
                                             survforecast_CovWY3_W_FL_Fea_fig,
                                             ncol=3, common.legend=TRUE,
                                             # labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_FL_Fea_fig, "figures/SurvForecast_CovWY3_FL_Fea.png", width =6, height = 4)

# All fish combined
(survforecast_CovWY3_FL_SacButFea_fig <- ggarrange(survforecast_CovWY3_FL_Sac_fig,
                                                   survforecast_CovWY3_FL_But_fig,
                                                    survforecast_CovWY3_FL_Fea_fig,
                                             nrow=3, 
                                             labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_FL_SacButFea_fig, "figures/SurvForecast_CovWY3_FL_SacButFea.png", 
       width =14, height =11)

# Figures for WY3 + Wgt model --------------------------------------------
load(file="outputs/fit_CovWY3_Wgt.Rdata")
posterior_draws_CovWY3_Wgt <- rstan::extract(fit)
posterior.df_CovWY3_Wgt <-as.data.frame(posterior_draws_CovWY3_Wgt)
names(posterior.df_CovWY3_Wgt)

parlist <- c("S_bTrib","S_bCovT","RE_sdT","S_bCov","S_bSz","S_bReach","RE_sd")
print(summary(fit,pars=parlist)$summary) 

fit_summary_CovWY3_Wgt <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY3_Wgt <- as.data.frame(fit_summary_CovWY3_Wgt$summary)
print(fit_summary_CovWY3_Wgt$summary)

### Detection probability figures
data_pcap_CovWY3_Wgt<- data.frame(matrix(NA, nrow = 4*length(Year_all) , ncol = 9))
pars <- c(paste0("pred_pcap[",1,",1]"),paste0("pred_pcap[",1,",2]"),
          paste0("pred_pcap[",1,",3]"),paste0("pred_pcap[",1,",4]"))
loc <- c('Woodson','Butte','Sacramento','Delta')
loc <- factor(loc, levels = unique(loc))

data_pcap_CovWY3_Wgt <- data.frame(cbind(loc,Year_all[1],data_summary_CovWY3_Wgt[pars,]))
colnames(data_pcap_CovWY3_Wgt) <- c('loc','year','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

for(i in 2:length(Year_all)){
  pars <- c(paste0("pred_pcap[",i,",1]"),paste0("pred_pcap[",i,",2]"),
            paste0("pred_pcap[",i,",3]"),paste0("pred_pcap[",i,",4]"))
  
  data_temp <- data.frame(cbind(loc,Year_all[i],data_summary_CovWY3_Wgt[pars,]))
  names(data_temp) <- names(data_pcap_CovWY3_Wgt) # to be able to bind data_temp and data make sure the column names at the same
  
  data_pcap_CovWY3_Wgt <- data.frame(rbind(data_pcap_CovWY3_Wgt,data_temp))
  
}


(pcap_CovWY3_Wgt_fig <- ggplot(data_pcap_CovWY3_Wgt)+
    geom_errorbar(aes(x=loc, ymin=q2.5, ymax=q97.5, color=loc), 
                  width=0.2, size=1)+
    geom_point(aes(x=loc,y=mean, color=loc),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.1))+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    xlab("")+ylab("Detection probability")+
    facet_grid(~year)
)

ggsave(plot=pcap_CovWY3_Wgt_fig , "figures/pcap_CovWY3_Wgt.png", width =10, height = 4)

### Release-Sacramento Survival figure
data_surv_relsac_CovWY3_Wgt <- data.frame(matrix(NA, nrow = dim(rgwy3.df)[1] , ncol = 8))

for(i in 1: dim(rgwy3.df)[1]){
  data_surv_relsac_CovWY3_Wgt[i,] <- data.frame(cbind(rgwy3.df$StudyID[i],
                                                      data_summary_CovWY3_Wgt[paste0('SurvRelSac[',i,"]"),]))
  
}

colnames(data_surv_relsac_CovWY3_Wgt) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_relsac_CovWY3_Wgt$studyid <- factor(data_surv_relsac_CovWY3_Wgt$studyid, 
                                       levels = unique(data_surv_relsac_CovWY3_Wgt$studyid))

data_surv_relsac_CovWY3_Wgt$wytype <- c('D/BN/AN','D/BN/AN','Critical','D/BN/AN','Wet','Wet','Wet','D/BN/AN','D/BN/AN',
                                 'Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Critical',
                                 'Wet','Wet')
data_surv_relsac_CovWY3_Wgt$wytype <- factor(data_surv_relsac_CovWY3_Wgt$wytype, 
                                      levels = unique(data_surv_relsac_CovWY3_Wgt$wytype))

(surv_relsac_CovWY3_Wgt_fig <- ggplot(data= data_surv_relsac_CovWY3_Wgt)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean,color=wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_relsac_CovWY3_Wgt_fig, "figures/SurvRelSac_CovWY3_Wgt.png", width =7, height = 6)

### Butte and Feather-Sacramento Survival figure
data_surv_tribsac_CovWY3_Wgt <- data.frame(matrix(NA, nrow = dim(rgwy3_T.df)[1] , ncol = 8))

for(i in 1: dim(rgwy3_T.df)[1]){
  data_surv_tribsac_CovWY3_Wgt[i,] <- data.frame(cbind(rgwy3_T.df$StudyID[i],
                                                data_summary_CovWY3_Wgt[paste0('pred_survT[',i,",1]"),]))
  
}

colnames(data_surv_tribsac_CovWY3_Wgt) <- c('studyid','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_surv_tribsac_CovWY3_Wgt$studyid <- factor(data_surv_tribsac_CovWY3_Wgt$studyid, 
                                        levels = unique(data_surv_tribsac_CovWY3_Wgt$studyid))

data_surv_tribsac_CovWY3_Wgt$wytype <- c('D/BN/AN','Critical','Critical','Critical','D/BN/AN','Wet','D/BN/AN','Wet',
                                  'Wet','Wet','D/BN/AN','D/BN/AN','Critical','Critical','Critical','Wet','Wet',
                                  'Wet')
data_surv_tribsac_CovWY3_Wgt$wytype <- factor(data_surv_tribsac_CovWY3_Wgt$wytype, 
                                       levels = unique(data_surv_tribsac_CovWY3_Wgt$wytype))

(surv_tribsac_CovWY3_Wgt_fig <- ggplot(data_surv_tribsac_CovWY3_Wgt)+
    geom_errorbar(aes(x=studyid, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=studyid,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    scale_y_continuous(limits=c(0,0.6), breaks=seq(0,0.6,0.1))+
    xlab("")+ylab("")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")
)

ggsave(plot=surv_tribsac_CovWY3_Wgt_fig, "figures/SurvTribSac_CovWY3_Wgt.png", width =7, height = 5)


(Surv_CovWY3_Wgt_fig <- ggarrange(surv_relsac_CovWY3_Wgt_fig ,surv_tribsac_CovWY3_Wgt_fig,
          ncol=2, common.legend=TRUE,
          labels=c('A','B'),
          align="h")
)

ggsave(plot=Surv_CovWY3_Wgt_fig, "figures/Surv_CovWY3_Wgt.png", width =9, height =6)


### WY type forecast figure
data_survforecast_wy3 <- data.frame(rbind(data_summary_CovWY3_Wgt['SurvForecast[1]',],
                                          data_summary_CovWY3_Wgt['SurvForecast[2]',],
                                          data_summary_CovWY3_Wgt['SurvForecast[3]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[1,1]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[1,2]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[2,1]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[2,2]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[3,1]',],
                                          data_summary_CovWY3_Wgt['TribSurvForecast[3,2]',]))
colnames(data_survforecast_wy3) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_survforecast_wy3$loc <- c('Sacramento','Sacramento','Sacramento','Butte','Feather','Butte','Feather',
                               'Butte','Feather')
data_survforecast_wy3$loc <- factor(data_survforecast_wy3$loc,
                                    levels = unique(data_survforecast_wy3$loc))

data_survforecast_wy3$wytype <- c('Critical','D/BN/AN','Wet','Critical','Critical','D/BN/AN','D/BN/AN',
                                  'Wet','Wet')
data_survforecast_wy3$wytype <- factor(data_survforecast_wy3$wytype,
                                       levels = unique(data_survforecast_wy3$wytype))

(survforecast_wy3_Fig <- ggplot(data_survforecast_wy3)+
    geom_errorbar(aes(x=wytype, ymin=q2.5, ymax=q97.5, color= wytype), 
                  width=0.2, size=1)+
    geom_point(aes(x=wytype,y=mean, color= wytype),size=3)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(angle = 45, hjust = 1))+
    scale_y_continuous(limits=c(0,0.9), breaks=seq(0,0.9,0.1))+
    xlab("")+ylab("Release - Sacramento survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "WY")+
    facet_wrap(~loc)
)

ggsave(plot=survforecast_wy3_Fig, "figures/SurvForecast_wy3.png", width =10, height = 6)

### Weight forecast figures
Wgt_vec  <- seq(from=round(min(c(WGT,WGT_T),na.rm=TRUE)),to=round(max(c(WGT,WGT_T),na.rm=TRUE)),length.out=Nsz)

data_survforecast_CovWY3_Wgt_Sac <- array(data = NA, dim = c(Nsz,3,7),
                                      dimnames =list(Size = Wgt_vec, WY = c('C','D/BN/AN','W'),var= c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')))

data_survforecast_CovWY3_Wgt_Trib <- array(data = NA, dim = c(2,Nsz,3,7),
                                      dimnames =list(Trib = c('Butte','Feather'),Size = Wgt_vec, WY = c('C','D/BN/AN','W'),var= c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')))

for(i in 1:Nsz){ # number of class size
  for (j in 1:3){ # Number of water year type 
    for (k in 1:7){ # number of variables saved in data frame
    data_survforecast_CovWY3_Wgt_Sac[i,j,k] <- data_summary_CovWY3_Wgt[paste0('SurvForecastSz[',i,",",j,"]"),k]
      for (l in 1:2){ # tributary index: 1 = Buttte, 2 = Feather
      data_survforecast_CovWY3_Wgt_Trib[l,i,j,k] <- data_summary_CovWY3_Wgt[paste0('TribSurvForecastSz[',i,",",j,",",l,"]"),k]
      }
    }
  }
}

# Sacramento fish
(survforecast_CovWY3_C_Wgt_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_Wgt_Sac[,1,],aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_Wgt_Sac[,1,],
                aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('upper Sac - C Year')
)

(survforecast_CovWY3_D_Wgt_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_Wgt_Sac[,2,],aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_Wgt_Sac[,2,],
                aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("")+
    ggtitle('upper Sac - D/BN/AN Year')
)

(survforecast_CovWY3_W_Wgt_Sac_fig <- ggplot()+
    geom_line(data=data_survforecast_CovWY3_Wgt_Sac[,3,],aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(data=data_survforecast_CovWY3_Wgt_Sac[,3,],
                aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
   scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("")+
    ggtitle('upper Sac - W Year')
)

(survforecast_CovWY3_Wgt_Sac_fig <- ggarrange(survforecast_CovWY3_C_Wgt_Sac_fig,survforecast_CovWY3_D_Wgt_Sac_fig,
                                             survforecast_CovWY3_W_Wgt_Sac_fig,
                                             ncol=3, common.legend=TRUE,
                                             # labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_Wgt_Sac_fig, "figures/SurvForecast_CovWY3_Wgt_Sac.png", width =6, height = 4)

# Butte Creek fish
(survforecast_CovWY3_C_Wgt_But_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[1,,1,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
   scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - C Year')
)

(survforecast_CovWY3_D_Wgt_But_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[1,,2,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
   scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - D/BN/AN Year')
)

(survforecast_CovWY3_W_Wgt_But_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[1,,3,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Butte Creek - W Year')
)

(survforecast_CovWY3_Wgt_But_fig <- ggarrange(survforecast_CovWY3_C_Wgt_But_fig,survforecast_CovWY3_D_Wgt_But_fig,
                                             survforecast_CovWY3_W_Wgt_But_fig,
                                             ncol=3, common.legend=TRUE,
                                             # labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_Wgt_But_fig, "figures/SurvForecast_CovWY3_Wgt_But.png", width =6, height = 4)

# Feather River fish
(survforecast_CovWY3_C_Wgt_Fea_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[2,,1,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - C Year')
)

(survforecast_CovWY3_D_Wgt_Fea_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[2,,2,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - D/BN/AN Year')
)

(survforecast_CovWY3_W_Wgt_Fea_fig <- ggplot(data=data_survforecast_CovWY3_Wgt_Trib[2,,3,])+
    geom_line(aes(x=Wgt_vec,y=mean),size=1)+
    geom_ribbon(aes(x=Wgt_vec,y=mean,ymin = q2.5, ymax = q97.5), alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+        
    scale_y_continuous(limits=c(0,0.95), breaks=seq(0,0.95,0.1))+
    xlab("Weight (g)")+ylab("Rel - Sac survival rate")+
    ggtitle('Feather River - W Year')
)

(survforecast_CovWY3_Wgt_Fea_fig <- ggarrange(survforecast_CovWY3_C_Wgt_Fea_fig,survforecast_CovWY3_D_Wgt_Fea_fig,
                                             survforecast_CovWY3_W_Wgt_Fea_fig,
                                             ncol=3, common.legend=TRUE,
                                             # labels=c('A','B','C'),
                                             align="h")
)

ggsave(plot=survforecast_CovWY3_Wgt_Fea_fig, "figures/SurvForecast_CovWY3_Wgt_Fea.png", width =6, height = 4)

# All fish combined
(survforecast_CovWY3_Wgt_SacButFea_fig <- ggarrange(survforecast_CovWY3_Wgt_Sac_fig,
                                                   survforecast_CovWY3_Wgt_But_fig,
                                                   survforecast_CovWY3_Wgt_Fea_fig,
                                                   nrow=3, 
                                                   labels=c('A','B','C'),
                                                   align="h")
)

ggsave(plot=survforecast_CovWY3_Wgt_SacButFea_fig, "figures/SurvForecast_CovWY3_Wgt_SacButFea.png", 
       width =14, height =11)

# Figures for WY3 + CF model --------------------------------------------
load(file="outputs/fit_CovWY3_CF.Rdata")
posterior_draws_CovWY3_CF <- rstan::extract(fit)
posterior.df_CovWY3_CF <-as.data.frame(posterior_draws_CovWY3_CF)
names(posterior.df_CovWY3_CF)

fit_summary_CovWY3_CF <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovWY3_CF <- as.data.frame(fit_summary_CovWY3_CF$summary)
print(fit_summary_CovWY3_CF$summary)

