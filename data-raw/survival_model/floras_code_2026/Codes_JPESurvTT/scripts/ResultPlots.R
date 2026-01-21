# Figures for Bayesian models
library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(dplyr)
library(ggpubr)
library(viridis)#Color scale
library(loo)
library(reshape2)
library(purrr)
library(glue)
library(here)
library(boot)

source("scripts/GetData.R")

# Extract LooIC and RE_SD values for model comparison -------------------------------------------
analyze_model <- function(file_path, parlist = c("RE_sd[1]","RE_sd[2]","RE_sd[3]","RE_sd[4]","RE_sdT[1]","RE_sdT[2]"), 
                          log_lik_name = "log_lik") {
  # Load the RData file (assumes object named `fit` is in the file)
  load(file_path)  # This will create `fit` in your environment
  
  if (!exists("fit")) stop("Model object 'fit' not found in the Rdata file.")
  
  # Extract posterior draws (if needed later)
  posterior_draws <- rstan::extract(fit)
  
  # Summarize parameters
  summary_df <- data.frame(summary(fit, pars = parlist)$summary)
  
  # Extract log_lik and compute LOOIC
  log_lik <- extract_log_lik(fit, parameter_name = log_lik_name)
  loo_result <- loo(log_lik)
  
  # Build performance metrics: transposed means + LOOIC estimate
  perf_df <- data.frame(cbind(t(summary_df$mean), loo_result$estimates["looic", "Estimate"]))
  
  # Optionally assign column names for clarity
  colnames(perf_df) <- c(parlist, "looic")
  
  return(perf_df)
}

# For multiple files
file_list <- c("outputs/fit_NoCov.Rdata",
               "outputs/fit_CovIndWY2.Rdata",
               "outputs/fit_CovIndWY3.Rdata",
               "outputs/fit_CovIndWY_FL.Rdata",
               "outputs/fit_CovIndWY_Wgt.Rdata",
               "outputs/fit_CovIndWY_CF.Rdata",
               "outputs/fit_CovIndWY_Fexceed.Rdata",
               "outputs/fit_CovIndContTT_MaxFlow.Rdata",
               "outputs/fit_CovIndContTT_MaxFlow_FL.Rdata",
)

all_results <- lapply(file_list, analyze_model)

# Combine into one dataframe (with file labels)
names(all_results) <- basename(file_list)
combined_results <- do.call(rbind, all_results)
combined_results$model <- names(all_results)

combined_results_ord <- combined_results %>% 
                        arrange(looic)

write.csv(combined_results_ord, "outputs/model_performance_summary.csv", row.names = FALSE)

# Figures for each model --------------------------------------------------
process_fit_model <- function(fit_file, output_prefix, rg_data, rg_t_data) {
  load(fit_file)  # loads 'fit'
  
  # Extract posterior draws
  posterior_draws <- rstan::extract(fit)
  posterior.df <- as.data.frame(posterior_draws)
  
  # Fit summary
  fit_summary <- summary(fit, probs=c(0.025, 0.975))
  data_summary <- as.data.frame(fit_summary$summary)
  
  # Save full summary
  write.csv(data_summary, paste0("outputs/data_summary_", output_prefix, ".csv"), row.names = TRUE)
  
  #### Random effects figures
  pars <- c(paste0("S_RE[", 1:length(unique(d_FeaBut_sort$StudyID)), "]"), paste0("S_REt[", 
                                                                       1:length(unique(d_FeaBut_sort$StudyID)), "]"))
  trib_vec <- c(rep("Sac",length(unique(d_Sac_sort$StudyID))), rep("ButFea", length(unique(d_FeaBut_sort$StudyID))))
  data_RE <- data_summary[pars,] %>%
    tibble::rownames_to_column("Cov") %>%
    mutate(trib = trib_vec)
  
  sac_plot <- ggplot(data = data_RE %>% filter(trib == "Sac"), aes(x = Cov, y = sd)) +
    geom_point() + theme_bw()
  
  butfea_plot <- ggplot(data = data_RE %>% filter(trib == "ButFea"), aes(x = Cov, y = sd)) +
    geom_point() + theme_bw()
  
  ggsave(paste0("figures/RE_", output_prefix, "_Sac.jpg"), plot = sac_plot, dpi = 350, height = 4, width = 12)
  ggsave(paste0("figures/RE_", output_prefix, "_ButFea.jpg"), plot = butfea_plot, dpi = 350, height = 4, width = 12)

  
#### Release-Sacramento survival figure
  n_studies <- nrow(rg_data)
  surv_relsac_rows <- paste0("SurvRelSac[", 1:n_studies, "]")
  
  data_surv_relsac <- data.frame(matrix(NA, nrow = n_studies, ncol = 8))
  for (i in 1:n_studies) {
    data_surv_relsac[i, ] <- data_summary[surv_relsac_rows[i], ]
  }
  
  colnames(data_surv_relsac) <- c("mean", "se_mean", "sd", "q2.5", "q97.5", "n_eff", "rhat")
  data_surv_relsac$studyid <- factor(rg_data$StudyID, levels = unique(rg_data$StudyID))
  
  # Assign WY types (should be updated as needed)
  wytypes <- c('C/D/BN','C/D/BN','C/D/BN','C/D/BN','AN/W','AN/W','AN/W','C/D/BN','C/D/BN',
               'AN/W','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN','C/D/BN',
               'AN/W','AN/W','AN/W','AN/W','AN/W','AN/W')
  
  data_surv_relsac$wytype <- factor(wytypes[1:n_studies], levels = unique(wytypes))
  
  surv_relsac_plot <- ggplot(data_surv_relsac) +
    geom_errorbar(aes(x = studyid, ymin = q2.5, ymax = q97.5, color = wytype),
                  width = 0.2, size = 1) +
    geom_point(aes(x = studyid, y = mean, color = wytype), size = 3) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size = 18),
          axis.text.x = element_text(size = 10, angle = 90, hjust = 1)) +
    scale_y_continuous(limits = c(0, 0.5), breaks = seq(0, 0.5, 0.1)) +
    xlab("") + ylab("Release - Sacramento survival rate") +
    scale_colour_viridis(discrete = TRUE, direction = 1, name = "WY")
  
  ggsave(plot = surv_relsac_plot, filename = paste0("figures/SurvRelSac_", output_prefix, ".png"),
         width = 8, height = 5)

  
#### Forecast survival figure
surv_forecast_pars <- c("SurvForecast[1]", "SurvForecast[2]", "SurvForecast[3]",
                        "TribSurvForecast[1,1]", "TribSurvForecast[1,2]",
                        "TribSurvForecast[2,1]", "TribSurvForecast[2,2]",
                        "TribSurvForecast[3,1]", "TribSurvForecast[3,2]")

data_survforecast <- data.frame(data_summary[surv_forecast_pars, ])
data_survforecast$loc <- factor(c("Sacramento", "Sacramento", "Sacramento",
                                  "Butte", "Feather", "Butte", "Feather", "Butte", "Feather"))
data_survforecast$wytype <- factor(c("Critical", "D/BN/AN", "Wet",
                                     "Critical", "Critical", "D/BN/AN", "D/BN/AN", "Wet", "Wet"))

survforecast_fig <- ggplot(data_survforecast) +
  geom_errorbar(aes(x = wytype, ymin = q2.5, ymax = q97.5, color = wytype), width = 0.2, size = 1) +
  geom_point(aes(x = wytype, y = mean, color = wytype), size = 3) +
  theme_bw() + facet_wrap(~loc) +
  scale_colour_viridis(discrete = TRUE, direction = 1, name = "WY") +
  labs(y = "Release - Sacramento survival rate", x = "")

ggsave(plot = survforecast_fig, filename = paste0("figures/SurvForecast_", output_prefix, ".png"), width = 10, height = 6)

}

process_fit_model("outputs/fit_CovIndWY2.Rdata", output_prefix = "CovIndWY2", rg_data = rgwy2.df, rg_t_data = rgwy2_T.df)

# Figures for peak monthly flow + FL model --------------------------------------------
load(file="outputs/fit_CovIndContTT_MaxFlow_FL.Rdata")
posterior_draws_CovInd_MaxFlowFL <- rstan::extract(fit)
posterior.df_CovInd_MaxFlowFL <-as.data.frame(posterior_draws_CovInd_MaxFlowFL)
names(posterior.df_CovInd_MaxFlowFL)

fit_summary_CovInd_MaxFlowFL <- summary(fit,probs=c(0.025,0.975))# shows only 2.5% and 97.5% quantiles
data_summary_CovInd_MaxFlowFL <- as.data.frame(fit_summary_CovInd_MaxFlowFL$summary)
print(fit_summary_CovInd_MaxFlowFL$summary)

write.csv(data_summary_CovInd_MaxFlowFL ,"outputs/data_summary_MaxFlowFL.csv", row.names = TRUE)

#Do LOO calculations
log_lik <- extract_log_lik(fit, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik))
(loo <- loo(log_lik, r_eff = r_eff))
(pointwise_elpd <- loo$pointwise[ , "elpd_loo"])
loopointwise_stats <- loo$pointwise 
write.csv(loopointwise_stats ,"outputs/loopointwise_stats_MaxFlowFL.csv", row.names = FALSE)

### Detection probability figures
loc <- factor(c('Woodson', 'Butte', 'Sacramento', 'Delta'), levels = c('Woodson', 'Butte', 'Sacramento', 'Delta'))

data_pcap_CovInd_MaxFlowFL <- map_dfr(seq_along(Year_all), function(i) {
  pars <- paste0("pred_pcap[", i, ",", 1:4, "]")
  data.frame(
    loc = loc,
    year = Year_all[i],
    data_summary_CovInd_MaxFlowFL[pars, ],
    row.names = NULL
  )
}, .id = NULL)

colnames(data_pcap_CovInd_MaxFlowFL) <- c('loc','year','mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

write.csv(data_pcap_CovInd_MaxFlowFL,"./outputs/pcap_MaxFlowFL.csv", row.names=FALSE)

(pcap_CovInd_MaxFlowFL_fig <- ggplot(data_pcap_CovInd_MaxFlowFL)+
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
    scale_colour_viridis(discrete=T, direction = 1,name = "Location")+
    xlab("")+ylab("Detection probability")+
    facet_grid(~year)
)

ggsave(plot=pcap_CovInd_MaxFlowFL_fig ,filename = "figures/pcap_CovInd_MaxFlowFL.png", width =10, height = 4)

### Sacramento River fish survival figures
# Identify study years
study_years_Sac <- d_Sac_sort %>%
  select(StudyID, year) %>%
  distinct()

## plot reach specific survival for each release group
# Extract reach specific survival data and attach StudyID
Reach_vec <-c('Rel-Wood','Wood-But','But-Sac','Sac-Delta')
Reach_mat <- matrix(rep(Reach_vec, times = nrow(d_Sac_sort)), nrow = 4, ncol = nrow(d_Sac_sort))

data_surv_perreach_CovInd_MaxFlowFL <- map_dfr(rep(1:nrow(d_Sac_sort),4), function(i) {
                                      data_summary_CovInd_MaxFlowFL[paste0("pred_surv[", i, ",", 1:4, "]"),] %>%
                                      as.data.frame() %>%
                                      mutate(fishid = d_Sac_sort$FishID[i],
                                             StudyID = d_Sac_sort$StudyID[i],
                                             Reach= Reach_mat[1:4,i])
                                      })

colnames(data_surv_perreach_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

data_surv_perreach_CovInd_MaxFlowFL <- data_surv_perreach_CovInd_MaxFlowFL %>% 
                                        mutate(Reach = factor(Reach, levels = unique(Reach)))

surv_summary_perreach <- data_surv_perreach_CovInd_MaxFlowFL %>%
  left_join(study_years_Sac, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
 arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)))

write.csv(surv_summary_perreach,"./outputs/surv_perreach_summary_MaxFlowFL.csv", row.names=FALSE)

(surv_perreach_CovInd_MaxFlowFL_fig <- ggplot(surv_summary_perreach)+
    geom_point(aes(x = Reach, y = mean_surv,color=Reach), size = 3) +
    geom_errorbar(aes(x = Reach, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    #   scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
   facet_wrap(year~StudyID,nrow=6)
)

ggsave(plot=surv_perreach_CovInd_MaxFlowFL_fig, filename = "figures/surv_perreach_CovInd_MaxFlowFL.png", width =13, height = 10)

## Plot reach-specific survival for each reach standardized per 100km
data_surv_perreach100_CovInd_MaxFlowFL <- map_dfr(rep(1:nrow(d_Sac_sort),4), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("pred_surv_per100[", i, ",", 1:4, "]"),] %>%
    as.data.frame() %>%
    mutate(fishid = d_Sac_sort$FishID[i],
           StudyID = d_Sac_sort$StudyID[i],
           Reach= Reach_mat[1:4,i])
})

colnames(data_surv_perreach100_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

data_surv_perreach100_CovInd_MaxFlowFL <- data_surv_perreach100_CovInd_MaxFlowFL %>% 
  mutate(Reach = factor(Reach, levels = unique(Reach)))

surv_summary_perreach100 <- data_surv_perreach100_CovInd_MaxFlowFL %>%
  left_join(study_years_Sac, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)))

write.csv(surv_summary_perreach100,"./outputs/surv_perreach100_summary_MaxFlowFL.csv", row.names=FALSE)

(surv_perreach100_CovInd_MaxFlowFL_fig <- ggplot(surv_summary_perreach100)+
    geom_point(aes(x = Reach, y = mean_surv,color=Reach), size = 3) +
    geom_errorbar(aes(x = Reach, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    facet_wrap(year~StudyID,nrow=6)
)

ggsave(plot=surv_perreach100_CovInd_MaxFlowFL_fig, filename = "figures/surv_perreach100_CovInd_MaxFlowFL.png", width =15, height = 12)


## plot survival from release to Sac for each release group
# Extract release to Sac survival data and attach StudyID
data_surv_relsac_CovInd_MaxFlowFL <- map_dfr(1:nrow(d_Sac_sort), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("SurvRelSac[", i, "]"), ] %>%
    as.data.frame() %>%
    mutate(
      fishid = d_Sac_sort$FishID[i],
      StudyID = d_Sac_sort$StudyID[i]
    )
})
colnames(data_surv_relsac_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID')

surv_relsac_summary_study <- data_surv_relsac_CovInd_MaxFlowFL %>%
  left_join(study_years_Sac, by = "StudyID") %>%
  group_by(StudyID, year) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)),
         WY = case_when(year %in% c(2013,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W')) %>% 
  left_join(d_Sac_sort %>% select(StudyID,Maxflow,fish_length),by="StudyID")

write.csv(surv_relsac_summary_study,"./outputs/surv_relsac_summary_MaxFlowFL.csv", row.names=FALSE)

(surv_relsac_CovInd_MaxFlowFL_fig <- ggplot(surv_relsac_summary_study, aes(x = StudyID, y = mean_surv,color=WY))+
    geom_point( size = 3) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 45, hjust = 1))+
    #   scale_y_continuous(limits=c(0,0.5), breaks=seq(0,0.5,0.1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=surv_relsac_CovInd_MaxFlowFL_fig, "figures/surv_relsac_CovInd_MaxFlowFL.png", width =7, height = 6)

#### Butte and Feather fish Survival figures
## Reach specific survival figure for each release group

# Include year in the summary
study_yearsT <- d_FeaBut_sort %>%
  select(StudyID, year) %>%
  distinct()

ReachT_vec <-c('Rel-Sac','Sac-Delta')
ReachT_mat <- matrix(rep(ReachT_vec, times = nrow(d_FeaBut_sort)), nrow = 2, ncol = nrow(d_FeaBut_sort))

data_surv_tribperreach_CovInd_MaxFlowFL <-  map_dfr(rep(1:nrow(d_FeaBut_sort),2), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("pred_survT[", i, ",", 1:2, "]"), ]  %>%
    as.data.frame() %>%
    mutate(
      fishid = d_FeaBut_sort$FishID[i],
      StudyID = d_FeaBut_sort$StudyID[i],
      Reach= ReachT_mat[1:2,i])
})
colnames(data_surv_tribperreach_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

surv_summary_tribperreach <- data_surv_tribperreach_CovInd_MaxFlowFL %>%
  left_join(study_yearsT, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)))

(surv_tribperreach_CovInd_MaxFlowFL_fig <- ggplot(surv_summary_tribperreach)+
    geom_point(aes(x = Reach, y = mean_surv,color=Reach), size = 3) +
    geom_errorbar(aes(x = Reach, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    facet_wrap(year~StudyID,nrow=6)
)

ggsave(plot=surv_tribperreach_CovInd_MaxFlowFL_fig, "figures/surv_tribperreach_CovInd_MaxFlowFL.png", width =13, height = 10)

## Reach specific survival per 100km figure for each release group
data_surv_tribperreach100_CovInd_MaxFlowFL <-  map_dfr(rep(1:nrow(d_FeaBut_sort),2), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("pred_survT_per100[", i, ",", 1:2, "]"), ]  %>%
    as.data.frame() %>%
    mutate(
      fishid = d_FeaBut_sort$FishID[i],
      StudyID = d_FeaBut_sort$StudyID[i],
      Reach= ReachT_mat[1:2,i])
})
colnames(data_surv_tribperreach100_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

surv_summary_tribperreach100 <- data_surv_tribperreach100_CovInd_MaxFlowFL %>%
  left_join(study_yearsT, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)))

(surv_tribperreach100_CovInd_MaxFlowFL_fig <- ggplot(surv_summary_tribperreach100)+
    geom_point(aes(x = Reach, y = mean_surv,color=Reach), size = 3) +
    geom_errorbar(aes(x = Reach, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Survival rate")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    facet_wrap(year~StudyID,nrow=6)
)

ggsave(plot=surv_tribperreach100_CovInd_MaxFlowFL_fig, "figures/surv_tribperreach100_CovInd_MaxFlowFL.png", width =15, height = 12)

## Release to Sacramento survival figure
data_surv_tribsac_CovInd_MaxFlowFL <- map_dfr(1:nrow(d_FeaBut_sort), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0('pred_survT[',i,",1]"), ] %>%
    as.data.frame() %>%
    mutate(
      fishid = d_FeaBut_sort$FishID[i],
      StudyID = d_FeaBut_sort$StudyID[i]
    )
})
colnames(data_surv_tribsac_CovInd_MaxFlowFL) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID')

surv_summary_studyT <- data_surv_tribsac_CovInd_MaxFlowFL %>%
  left_join(study_yearsT, by = "StudyID") %>%
  group_by(StudyID, year) %>%
  summarise(
    mean_surv = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)),
         WY = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W')) %>% 
  left_join(d_FeaBut_sort %>% select(StudyID,Maxflow),by="StudyID") %>% 
  mutate(Location = case_when(str_starts(StudyID, "FR") ~ 'Feather',
                              TRUE ~ 'Butte'))

(surv_tribsac_CovInd_MaxFlowFL_fig <- ggplot(surv_summary_studyT, aes(x = StudyID, y = mean_surv,color=WY))+
    geom_point( size = 3) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 45, hjust = 1))+
    xlab("")+ylab("Rel - Sac survival rate")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=surv_tribsac_CovInd_MaxFlowFL_fig, "figures/surv_tribsac_CovInd_MaxFlowFL.png", width =7, height = 6)

(Surv_allrelsac_CovIndMaxFlowFL_fig <- ggarrange(surv_relsac_CovInd_MaxFlowFL_fig ,surv_tribsac_CovInd_MaxFlowFL_fig,
                                       ncol=2, common.legend=TRUE,
                                       labels=c('A','B'),
                                       align="h")
)

ggsave(plot=Surv_allrelsac_CovIndMaxFlowFL_fig, "figures/Surv_allrelsac_CovIndMaxFlowFL.png", width =11, height =6)

######## forecast figures
#### Flow based figures
# Define index ranges
i_vals <- 1:25   #flow range
j_vals <- 1:2   # trib locations
n_cov <- 25

## Sacramento River fish figures
# flow values
mux <- mean(MaxflowSac,na.rm=TRUE)
sdx <- sd(MaxflowSac,na.rm=TRUE)
Xvec <- seq(from =4000,to=40000,length.out=NsX)
flow_vec <- Xvec

locations <- rep("Sacramento", n_cov)

# Build the dataframes
data_survforecast_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("SurvForecast[", i, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_MaxFlowFL <- data_survforecast_MaxFlowFL %>% 
  mutate(flow_vec = flow_vec,
         Location =locations,
         Scenario= "With RE")

data_survforecast_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("SurvForecast_nore[", i, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_nore_MaxFlowFL <- data_survforecast_nore_MaxFlowFL %>% 
  mutate(flow_vec = flow_vec,
         Location =locations,
         Scenario= "No RE")

(survforecast_MaxFlowFL_Sac_Fig <- ggplot() +
    geom_line(data=data_survforecast_MaxFlowFL, aes(x = flow_vec, y = mean),size = 1, color = "black") +
    geom_ribbon(data=data_survforecast_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5), fill="black",alpha = .25)+
    geom_rug(data=d_Sac_sort,aes(x = Maxflow), sides = "b", alpha = 0.5) + 
    geom_line(data=data_survforecast_nore_MaxFlowFL,aes(x = flow_vec, y = mean),size = 1, color = "darkblue") +
    geom_ribbon(data=data_survforecast_nore_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5), 
                fill="darkblue",alpha = .25)+
    geom_point(data=surv_summary_study, aes(x = Maxflow, y = mean_surv), size = 2) +
    geom_errorbar(data=surv_summary_study,aes(x = Maxflow, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Sacramento')+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(5000, 40000, by = 10000),limits = c(4000, 40000))
)

ggsave(plot=survforecast_MaxFlowFL_Sac_Fig, "figures/survforecast_MaxFlowFL_Sac.png", width =10, height = 6)

## Butte and Feather fish figures
locationBut <- rep("Butte", n_cov)
locationFea <- rep("Feather", n_cov)
muxT <- mean(MaxflowT,na.rm=TRUE)
sdxT <- sd(MaxflowT,na.rm=TRUE)
XvecT <- seq(from =150,to=30500,length.out=NsX)
flow_vecT <- XvecT

data_survforecast_but_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecast[", i, ",", 1, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_but_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_but_MaxFlowFL <- data_survforecast_but_MaxFlowFL %>% 
  mutate(flow_vec = flow_vecT,
         Location =locationBut)

data_survforecast_but_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecast_nore[", i, ",", 1, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_but_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_but_nore_MaxFlowFL <- data_survforecast_but_nore_MaxFlowFL %>% 
  mutate(flow_vec = flow_vecT,
         Location =locationBut)

data_survforecast_fea_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecast[", i, ",", 2, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_fea_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_fea_MaxFlowFL <- data_survforecast_fea_MaxFlowFL %>% 
  mutate(flow_vec = flow_vecT,
         Location =locationFea)

data_survforecast_fea_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecast_nore[", i, ",", 2, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecast_fea_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecast_fea_nore_MaxFlowFL <- data_survforecast_fea_nore_MaxFlowFL %>% 
  mutate(flow_vec = flow_vecT,
         Location =locationFea)

(survforecast_MaxFlowFL_Fea_Fig <- ggplot() +
    geom_line(data=data_survforecast_fea_MaxFlowFL, aes(x = flow_vec, y = mean),size = 1, color = "black") +
    geom_ribbon(data=data_survforecast_fea_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5),
               fill="black", alpha = .25)+
    geom_rug(data=d_FeaBut_sort %>% filter(rl %in% c("1F","2F")),aes(x = Maxflow), sides = "b", alpha = 0.5) + 
    geom_line(data=data_survforecast_fea_nore_MaxFlowFL, aes(x = flow_vec, y = mean),size = 1, 
              color = "darkblue") +
    geom_ribbon(data=data_survforecast_fea_nore_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="darkblue", alpha = .25)+
    geom_point(data=surv_summary_studyT %>% filter(Location =="Feather"), aes(x = Maxflow, y = mean_surv), size = 2) +
    geom_errorbar(data=surv_summary_studyT %>% filter(Location =="Feather"),aes(x = Maxflow, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Feather')+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(500,30500, by = 5000),limits = c(150,30500))
)

ggsave(plot=survforecast_MaxFlowFL_Fea_Fig, "figures/survforecast_MaxFlowFL_Fea.png", width =10, height = 6)

(survforecast_MaxFlowFL_But_Fig <- ggplot() +
    geom_line(data=data_survforecast_but_MaxFlowFL,aes(x = flow_vec, y = mean),size = 1, color = "black") +
    geom_ribbon(data=data_survforecast_but_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5), 
                fill="black",alpha = .25)+
    geom_rug(data=d_FeaBut_sort %>% filter(rl != c("1F","2F")),aes(x = Maxflow), sides = "b", alpha = 0.5) + 
    theme_bw()+
    geom_line(data=data_survforecast_but_nore_MaxFlowFL,aes(x = flow_vec, y = mean),size = 1, color = "darkblue") +
    geom_ribbon(data=data_survforecast_but_nore_MaxFlowFL,aes(x = flow_vec, y = mean,ymin = q2.5, ymax = q97.5), 
                fill="darkblue",alpha = .25)+
    geom_point(data=surv_summary_studyT %>% filter(Location =="Butte"), aes(x = Maxflow, y = mean_surv), size = 2) +
    geom_errorbar(data=surv_summary_studyT %>% filter(Location =="Butte"),aes(x = Maxflow, y = mean_surv,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Flow",
          y = "Release - Sacramento survival rate")+
    ggtitle('Butte')+
    scale_y_continuous(limits = c(0, 1))+    
    scale_x_continuous(breaks = seq(500,30500, by = 5000),limits = c(150,30500))
)

ggsave(plot=survforecast_MaxFlowFL_But_Fig, "figures/survforecast_MaxFlowFL_But.png", width =10, height = 6)

(survforecast_alltrib_MaxFlow_Fig <- ggarrange(survforecast_MaxFlowFL_Sac_Fig,
                                               survforecast_MaxFlowFL_But_Fig,
                                               survforecast_MaxFlowFL_Fea_Fig,
                                               nrow=3, 
                                               labels=c('A','B','C'),
                                               align="h")
)

ggsave(plot=survforecast_alltrib_MaxFlow_Fig , "figures/survforecast_alltrib_MaxFlow.png", width =8, height =14)

#### Size based figures
k_vals <- 1:25  # sizes range
size_vec <- seq(from=10,to=130,length.out=Nsz)

##  Sacramento River fish figures
locations_Sac <- rep("Sacramento", n_cov)

# Build the dataframes
data_survforecastSz_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("SurvForecastSz[", i, ",", 12, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_MaxFlowFL <- data_survforecastSz_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locations)

data_survforecastSz_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("SurvForecastSz_nore[", i, ",", 12, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_nore_MaxFlowFL <- data_survforecastSz_nore_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locations)

(survforecastSz_MaxFlowFL_Sac_Fig <- ggplot() +
    geom_line(data=data_survforecastSz_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "black") +
    geom_ribbon(data=data_survforecastSz_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="black",alpha = .25)+
    geom_rug(data=d_Sac_sort,aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=data_survforecastSz_nore_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=data_survforecastSz_nore_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Fish FL (mm)",
          y = "Rel - Sac survival rate")+
    ggtitle("Sacramento")+
   scale_y_continuous(limits = c(0, 1))+
  scale_x_continuous(breaks = seq(10,130, by = 40),limits = c(10,130))
)

ggsave(plot=survforecastSz_MaxFlowFL_Sac_Fig, "figures/survforecastSz_MaxFlowFL_Sac.png", width =10, height = 6)

## Butte Creek fish figure
data_survforecastSz_but_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecastSz[", i, ",", 12, ",", 1, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_but_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_but_MaxFlowFL <- data_survforecastSz_but_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locationBut)

data_survforecastSz_but_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecastSz_nore[", i, ",", 12, ",", 1, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_but_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_but_nore_MaxFlowFL <- data_survforecastSz_but_nore_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locationBut)

(survforecastSz_MaxFlowFL_But_Fig <- ggplot() +
    geom_line(data=data_survforecastSz_but_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "black") +
    geom_ribbon(data=data_survforecastSz_but_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="black",alpha = .25)+
    geom_rug(data=d_Sac_sort,aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=data_survforecastSz_but_nore_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=data_survforecastSz_but_nore_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Fish FL (mm)",
          y = "Rel - Sac survival rate")+
    ggtitle("Butte")+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(10,130, by = 40),limits = c(10,130))
)

ggsave(plot=survforecastSz_MaxFlowFL_But_Fig, "figures/survforecastSz_MaxFlowFL_But.png", width =10, height = 6)

## Feather River fish figure
data_survforecastSz_fea_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecastSz[", i, ",", 12, ",", 2, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_fea_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_fea_MaxFlowFL <- data_survforecastSz_fea_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locationFea)

data_survforecastSz_fea_nore_MaxFlowFL <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TribSurvForecastSz_nore[", i, ",", 12, ",", 2, "]"), ] %>%
    as.data.frame()
})
colnames(data_survforecastSz_fea_nore_MaxFlowFL)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')
data_survforecastSz_fea_nore_MaxFlowFL <- data_survforecastSz_fea_nore_MaxFlowFL %>% 
  mutate(size_vec = size_vec,
         Location =locationFea)

(survforecastSz_MaxFlowFL_Fea_Fig <- ggplot() +
    geom_line(data=data_survforecastSz_fea_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "black") +
    geom_ribbon(data=data_survforecastSz_fea_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="black",alpha = .25)+
    geom_rug(data=d_Sac_sort,aes(x = fish_length), sides = "b", alpha = 0.5) + 
    geom_line(data=data_survforecastSz_fea_nore_MaxFlowFL,aes(x = size_vec, y = mean), size = 1, color = "darkblue") +
    geom_ribbon(data=data_survforecastSz_fea_nore_MaxFlowFL,aes(x = size_vec, y = mean,ymin = q2.5, ymax = q97.5),
                fill="darkblue",alpha = .25)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18))+   
    labs( x = "Fish FL (mm)",
          y = "Rel - Sac survival rate")+
    ggtitle("Feather")+
    scale_y_continuous(limits = c(0, 1))+
    scale_x_continuous(breaks = seq(10,130, by = 40),limits = c(10,130))
)

ggsave(plot=survforecastSz_MaxFlowFL_Fea_Fig, "figures/survforecastSz_MaxFlowFL_Fea.png", width =10, height = 6)

## All fish combined
(survforecastSz_MaxFlowFL_alltrib_Fig <- ggarrange(survforecastSz_MaxFlowFL_Sac_Fig,
                                           survforecastSz_MaxFlowFL_But_Fig,
                                           survforecastSz_MaxFlowFL_Fea_Fig,
                                                      nrow=3, 
                                                      labels=c('A','B','C'),
                                                      align="h")
)

ggsave(plot=survforecastSz_MaxFlowFL_alltrib_Fig , "figures/survforecastSz_MaxFlow_alltrib.png", width =8, height =14)

#### Random effects figures
study_ids_sac <- rep(unique(d_Sac_sort$StudyID),each=4)
Location <- rep('Sacramento', length(unique(d_Sac_sort$StudyID))*4)
reach_sac <- factor(c('Woodson', 'Butte', 'Sacramento', 'Delta'), levels = c('Woodson', 'Butte', 'Sacramento', 'Delta'))
Reach <- rep(reach_sac,length(unique(d_Sac_sort$StudyID)))

data_RE <- map_dfr(1:length(unique(d_Sac_sort$StudyID)), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("S_RE[", i, ",", 1:4, "]"), ] %>%
    as.data.frame() 
})
colnames(data_RE)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_RE <- data_RE %>% 
  mutate(Location=Location,
         Reach=Reach,
         StudyID = factor(study_ids_sac,levels=unique(study_ids_sac)))

study_ids_butfea <- rep(unique(d_FeaBut_sort$StudyID),each=2)
location_butfea <- rep("Tributaries",length(unique(d_FeaBut_sort$StudyID))*2)
reach_butfea <- factor(c('Sacramento', 'Delta'), levels = c('Sacramento', 'Delta'))
Reacht <- rep(reach_butfea,length(unique(d_FeaBut_sort$StudyID)))

data_REt <- map_dfr(1:length(unique(d_FeaBut_sort$StudyID)), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("S_REt[", i, ",", 1:2, "]"), ] %>%
    as.data.frame() 
})
colnames(data_REt)[1:7] <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

data_REt <- data_REt %>% 
  mutate(Location=location_butfea,
         Reach=Reacht,
         StudyID = factor(study_ids_butfea,levels=unique(study_ids_butfea)))

(Re_Fig <- ggplot(data = data_RE ,aes(x = Reach, y = mean,color=Reach)) +
    geom_point(aes(color=Reach), size = 2) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2,size=1) +
    xlab('')+ ylab("S_RE")+
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        legend.title = element_blank(),text = element_text(size=16))+
    facet_wrap(Location~StudyID)
)

ggsave(plot=Re_Fig , "figures/Re.png", width =16, height =12)

(Ret_Fig <- ggplot(data = data_REt ,aes(x = Reach, y = mean,color=Reach)) +
    geom_point(aes(color=Reach), size = 2) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2,size=1) +
    xlab('')+ ylab("S_RE")+
    theme_bw() +
    theme(axis.text.x = element_blank(), 
          legend.title = element_blank(),text = element_text(size=16))+
    facet_wrap(Location~StudyID)
)

ggsave(plot=Ret_Fig , "figures/Ret.png", width =16, height =12)

data_re_MaxFlowFL <- data.frame(rbind(data_RE,data_REt))
write.csv(data_re_MaxFlowFL ,here("outputs", "data_re_MaxFlowFL.csv"), row.names=TRUE)

##### Travel Time figures 
## Sacramento River fish figures
data_tt_perreach <- map_dfr(1:nrow(d_Sac_sort), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TT_reach[", i, ",", 1:4, "]"),] %>%
    as.data.frame() %>% 
    mutate(fishid = d_Sac_sort$FishID[i],
           StudyID = d_Sac_sort$StudyID[i],
           Reach= Reach_mat[1:4,i])
})

colnames(data_tt_perreach) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

data_tt_perreach <- data_tt_perreach %>% 
  mutate(Reach = factor(Reach, levels = unique(Reach)))

tt_summary_perreach <- data_tt_perreach %>%
  left_join(study_years_Sac, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_tt = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)))


data_tt_relsac <- map_dfr(1:nrow(d_Sac_sort), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TT_RelSac[", i, "]"),] %>%
    as.data.frame() %>% 
    mutate(fishid = d_Sac_sort$FishID[i],
           StudyID = d_Sac_sort$StudyID[i])
})

colnames(data_tt_relsac) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID')

tt_summary_relsac <- data_tt_relsac %>%
  left_join(study_years_Sac, by = "StudyID") %>%
  group_by(StudyID, year) %>%
  summarise(
    mean_tt = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)),
         WY = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W'))

data_tt_forecast <- map_dfr(1:n_cov, function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TTForecast[", i, "]"),] %>%
    as.data.frame()
})

colnames(data_tt_forecast) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat')

tt_forecast_summary <- data_tt_forecast %>%
  summarise(mean_tt = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) 

(tt_perreach_fig <- ggplot(tt_summary_perreach)+
    geom_point(aes(x = Reach, y = mean_tt,color=Reach), size = 2) +
    geom_errorbar(aes(x = Reach, y = mean_tt,ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Travel time")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    facet_wrap(year~StudyID)
  
)

ggsave(plot=tt_perreach_fig, "figures/tt_perreach.png", width =14, height = 9)

(tt_relsac_fig <- ggplot(tt_summary_relsac, aes(x = StudyID, y = mean_tt,color=WY))+
    geom_point( size = 3) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    xlab("")+ylab("Travel time")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=tt_relsac_fig, "figures/tt_relsac.png", width =7, height = 6)

## Butte and Feather fish figures
data_tt_perreachT <- map_dfr(1:nrow(d_FeaBut_sort), function(i) {
  data_summary_CovInd_MaxFlowFL[paste0("TT_reachT[", i, ",", 1:2, "]"),] %>%
    as.data.frame() %>% 
    mutate(fishid = d_FeaBut_sort$FishID[i],
           StudyID = d_FeaBut_sort$StudyID[i],
           Reach= ReachT_mat[1:2,i])
})

colnames(data_tt_perreachT) <- c('mean','se_mean','sd','q2.5','q97.5','n_eff','rhat','FishID','StudyID','Reach')

data_tt_perreachT <- data_tt_perreachT %>% 
  mutate(Reach = factor(Reach, levels = unique(Reach)))

tt_summary_perreachT <- data_tt_perreachT %>%
  left_join(study_yearsT, by = "StudyID") %>%
  group_by(StudyID, year,Reach) %>%
  summarise(
    mean_tt = mean(mean, na.rm = TRUE),
    q2.5 = mean(q2.5, na.rm = TRUE),
    q97.5 = mean(q97.5, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(year) %>%
  mutate(StudyID = factor(StudyID, levels = unique(StudyID)),
         WY = case_when(year %in% c(2013,2014,2015,2016,2018,2020, 2021, 2022) ~ 'D',
                        TRUE ~ 'W'))

(tt_perreachT_fig <- ggplot(tt_summary_perreachT, aes(x = Reach, y = mean_tt, color=Reach))+
    geom_point( size = 2) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_blank())+
    xlab("")+ylab("Travel time")+
    scale_colour_viridis(discrete=T, direction = 1,name = "Reach")+
    facet_wrap(year~StudyID)
)

ggsave(plot=tt_perreachT_fig, "figures/tt_perreachT.png", width =12, height =7)

(tt_relsacT_fig <- ggplot(tt_summary_perreachT %>% filter(Reach =="Rel-Sac") %>% select(-Reach),
                          aes(x = StudyID, y = mean_tt,color=WY))+
    geom_point( size = 3) +
    geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          text = element_text(size=18),axis.text.x = element_text(size=10,angle = 90, hjust = 1))+
    xlab("")+ylab("Travel time")+
    scale_color_manual(values = c("D"="red","W"="blue"),name = NULL)
)

ggsave(plot=tt_relsacT_fig, "figures/tt_tt_relsacT.png", width =12, height =7)

## All tribs combined
tt_summary_relsac_allloc <- data.frame(rbind(tt_summary_relsac,
                                             tt_summary_perreachT %>% filter(Reach =="Rel-Sac") %>% select(-Reach)))
tt_summary_relsac_allloc$Location <- c(rep('Sacramento',23),rep('ButteFeather',20))


(tt_relsac_allloc_fig <- ggarrange(tt_relsac_fig ,tt_relsacT_fig,
                                       ncol=2, common.legend=TRUE,
                                       labels=c('A','B'),
                                       align="h")
)

ggsave(plot=tt_relsac_allloc_fig, "figures/tt_RelSac_allloc.png", width =9, height = 5)


