# LOO analysis for P2S
library(loo)
library(SRJPEdata)
library(SRJPEmodel)
library(tidyverse)

# Battle
battle <- run_passage_to_spawner_model("battle creek", "wy_type", FALSE)

log_lik_1 <- extract_log_lik(battle$full_object, merge_chains = FALSE)
r_eff <- relative_eff(exp(log_lik_1), cores = 2)

loo_1 <- loo(log_lik_1, r_eff = r_eff, cores = 2)
print(loo_1) # Pareto k diagnostics are bad, but not very bad

# Clear
clear <- run_passage_to_spawner_model("clear creek", "wy_type", FALSE)

log_lik_clear <- extract_log_lik(clear$full_object, merge_chains = FALSE)
r_eff_clear <- relative_eff(exp(log_lik_clear), cores = 2)

loo_clear <- loo(log_lik_clear, r_eff = r_eff_clear, cores = 2)
print(loo_clear) # Pareto k diagnostics are bad, but not very bad

