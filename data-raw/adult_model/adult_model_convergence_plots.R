

# https://mc-stan.org/bayesplot/articles/visual-mcmc-diagnostics.html
# https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
# https://cran.r-project.org/web/packages/bayesplot/vignettes/plotting-mcmc-draws.html

# if you return the stan fit object
test <- extract(battle_results$full_object) # see all chains

mcmc_pairs(battle_results$full_object, pars = c("b1_survival", "mean_redds_per_spawner"))
mcmc_trace(battle_results, pars = c("b1_survival", "mean_redds_per_spawner"))
mcmc_areas(battle_results, pars = c("b1_survival", "mean_redds_per_spawner"))

mcmc_rhat(rhat(mill_results$full_object))
