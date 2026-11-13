# Archive and migration plan

Draft, 2026-11-12. Updated with decisions from review.

## Goal

Get SRJPEmodel down to the code that is live today (BT-SPAS-X, model storage, canonical runs) so the stock-recruit, survival/travel-time, in-season and PLAD functions can be added back as they reach a natural break (likely Jan-Apr 2027). Archive everything else without losing it.

## Decisions so far

- P2S is retired from the JPE. Its code is archived and won't be migrated.
- WinBUGS stays in both this repo and runSRJPE.
- Each migrated submodel gets one small regression test.
- Vignettes stay where they are, marked superseded until they're rewritten.
- Plot functions are named `plot_*`. Use `if (`, not `if(`.
- `storage.R` and `canonical.R` header comments have been shortened (done, comments only).

## What stays and what gets archived

| Keep | Archive |
|---|---|
| `R/juvenile_abundance.R` | `R/stock_recruit.R` |
| `R/plot.R` (BT-SPAS-X plots only) | `R/survival_model.R` |
| `R/storage.R`, `R/canonical.R`, `db/` | `R/in_season_model.R` |
| `R/data.R` (BT-SPAS-X and shared entries only) | `R/passage_to_spawner_model.R` |
| `model_files/pCap_all_sites.stan`, `pCap_one_site.stan`, `pCap_one_site_skew_re.stan`, `abundance_model.bug` | `R/integrated_model.R` (runSRJPE owns this now) |
| `data/bt_spas_x_bayes_params.rda`, `bt_spas_x_model_code.rda`, `site_order_north_south.rda`, `julian_week_to_date_lookup.rda` | `model_files/` for SR, survival, in-season, P2S, plus `pCap_trib_skew_re.stan` (nothing references it) |
| `data-raw/cache-*.R` (trimmed), `data-raw/seed_data.R`, `data-raw/upload_models/` | `data/` objects: `p2s_model_code`, `P2S_model_fits`, `P2S_comparison_results`, `comparison_results`, `stock_recruit_model_code`, `survival_model_code`, `in_season_model_code`, `forecast_seed` |
| `data-raw/juvenile_abundance/` (current scripts) | `data-raw/adult_model/`, `in_season_model/`, `integrated_model/`, `stock_recruit_model/`, `survival_model/` |
| `data-raw/WinBUGS14/` | `data-raw/juvenile_abundance/archive/` (move it to the top-level archive so there's one place to look) |
| All six vignettes and `vignettes/images/`; the five non-BT-SPAS-X ones marked superseded | `MRdata.txt` (only used by archived BT-SPAS-X scripts) |
| `.github/` workflows | |

Nothing outside SRJPEmodel calls an archived function. I checked runSRJPE, bt-spas-x, in-season, model-storage-app and SRJPEdashboard. The only calls are in comments in SRJPEdashboard and in `runSRJPE/scripts/original_code/`.

## How to archive

**Tag first, then move everything to a top-level `archive/` folder. Delete the folder after migration is done.**

1. **Tag `main` before touching anything.** For example `git tag legacy-pre-migration && git push origin legacy-pre-migration`. This is the permanent snapshot. It still installs (`pak::pak("SRJPE/SRJPEmodel@legacy-pre-migration")`) if anyone needs the old functions.
2. **`git mv` the archive column into `archive/`, keeping the original paths** (`archive/R/stock_recruit.R`, `archive/data-raw/adult_model/...`, etc.). `git mv` keeps history, which `git log --follow` can trace.
3. **Add `^archive$` to `.Rbuildignore`** so R CMD check and the package build skip it. (`^planning$` is already there.)
4. **Add a short `archive/README.md`**: one sentence on why it's here, the tag name, and that it gets deleted after migration.
5. **After the last submodel migrates (post Apr 2027), delete `archive/` in one commit.** The tag keeps it.

## Package fixes that go with the archive

Do these in the same PR or the package won't build cleanly:

- [ ] `devtools::document()` to regenerate `NAMESPACE` and `man/`. Delete orphaned `.Rd` files.
- [ ] Remove archived objects from `R/data.R`. Trim `data-raw/cache-models.R` and `cache-params.R` to BT-SPAS-X.
- [ ] `julian_week_to_date_lookup` is built from `data-raw/juvenile_abundance/archive/btspas_model_code/Jwk_Dates.txt`. Move that file to `data-raw/juvenile_abundance/` before archiving.
- [ ] Superseded vignettes (SR, survival, in-season, PLAD, and P2S marked as retired): add a one-line note at the top, and `knitr::opts_chunk$set(eval = FALSE)` in the setup chunk. The SR and in-season vignettes currently run code that calls functions being archived, so they'll fail to build without it.
- [ ] `_pkgdown.yml`: drop the reference sections for archived functions. Keep the articles, and add the superseded vignettes to the navbar. Also remove functions it lists that no longer exist (`prepare_mainstem_pCap_data`, `get_most_recent_model_results`, `generate_diagnostic_plot`, `pin_model_data`). Rebuild `docs/`, which is stale.
- [ ] `DESCRIPTION`: see the dependency items in the standards list below.
- [ ] Leave `storage.R`'s `.approved_input_model_names` alone. `stock_recruit`, `survival` and `inseason` are database model names, and other repos still store fits under them.
- [ ] `R CMD check`, then confirm `runSRJPE/scripts/01_run_btspasx.R` and `in-season/scripts/01_data_prep/prepare_model_data.R` still run.
- [ ] Bump the version (e.g. `0.1.0`) and start `NEWS.md` with the archive as its first entry.

## Bringing SRJPEmodel up to the R standards

In the order the infra standards suggest. Items marked *(changes behavior)* need your sign-off, because a caller could notice the difference. Everything else is mechanical.

**1. Secrets and `.gitignore`**
- [ ] Add `.RData` and `.httr-oauth` to `.gitignore`. `*.Rdata` only matches `.RData` on Macs, because git ignores case there; on Linux CI it wouldn't.
- [ ] Credentials are already read from environment variables (`AZ_CONTAINER_ACCESS_KEY`), and `.Renviron` and `config.yml` are ignored. Nothing to change.

**2. renv**: doesn't apply. SRJPEmodel is a package, so dependencies live in `DESCRIPTION`.

**3. Air, lintr and CI**
- [ ] Add `air.toml` and `.lintr`. The infra branch refers to an org `.lintr` template but doesn't include one yet.
- [ ] Run `air format .` as its own PR, with no other changes. This also fixes `if(` → `if (` everywhere.
- [ ] Add GitHub Actions for `R CMD check`, lint and the format check (`usethis::use_github_action()`). Tests don't fit models, so WinBUGS isn't needed in CI.

**4. Tests**
- [ ] Create `tests/testthat/`. Start with `test-juvenile-abundance.R`: one regression test of `prepare_pCap_inputs()` and `prepare_abundance_inputs()` on a small fixture passed through `input_catch_data` / `input_efficiency_data`. That way the test doesn't change with each SRJPEdata release.
- [ ] Add error tests for the argument checks (`expect_snapshot(error = TRUE)`).
- [ ] Add `Config/testthat/parallel: true` to `DESCRIPTION`.

**5. Everything else**

`DESCRIPTION` and namespace:
- [ ] Add every package `R/` uses to Imports (`ggplot2`, `DBI`, `pins`, `AzureStor`, `stringr`, `purrr`, `tibble`, `glue`, `tidyr`, `gridExtra`, `RPostgres`), with minimum versions, via `usethis::use_package()`. Remove `rstanarm` and `gtools` (unused).
- [ ] `SRJPEdata` isn't in `DESCRIPTION` at all, although `R/` calls it 30 times. It's GitHub-only, so adding it means `Remotes:`, which the standards say to avoid in shared packages. The alternative is publishing SRJPEdata to R-universe. **[decide]**
- [ ] `juvenile_abundance.R` calls `filter()`, `mutate()` etc. without `dplyr::`, and `NAMESPACE` has no imports. Add `@importFrom dplyr ...` once, so it no longer depends on the user having run `library(dplyr)`.
- [ ] `generate_lt_pCap_Us()` calls `rskew_normal()` without `brms::` in two places (`juvenile_abundance.R:909`, `:925`) (a bt-spas-x script loads brms to work around it). Prefix it.

Code in `juvenile_abundance.R` and `plot.R`:
- [ ] `fit_pCap_model()` sets `options(mc.cores = ...)` and never restores it. Use `on.exit()`.
- [ ] `T`/`F` → `TRUE`/`FALSE` (about 19 places).
- [ ] `sapply()` → `vapply()` in `plot.R`.
- [ ] Two `suppressWarnings(readr::parse_number(...))` calls in the extract functions: add a comment saying which warning is expected and why it's safe.
- [ ] *(changes behavior)* `run_abundance_workflow()` catches every error and returns only `error = TRUE`, so the reason is lost. Keep the error message in the returned row and report it with `cli::cli_warn()`.
- [ ] *(changes behavior)* `prepare_pCap_inputs()` warns and returns `NULL` on bad arguments. Switch to `cli::cli_abort()`.
- [ ] *(changes behavior)* Rename `generate_diagnostic_plot_juv()` → `plot_juv_diagnostics()`. No other repo calls it.
- [ ] **[decide]** Rename `R/juvenile_abundance.R` → `R/juvenile-abundance.R` (kebab-case file names).
- [ ] Remove `return()` from the last line of functions. Low priority; do it only when touching a function.
- [ ] Joins in `juvenile_abundance.R` don't declare `relationship =`. Add it when touching that code. If one then errors, that's a data issue to look into, not something to make go away.
- [ ] Rows dropped by `!is.na()` filters and `na.rm = TRUE`: review only, no change without a decision. `prepare_pCap_inputs()` already reports one of these.

Roxygen:
- [ ] Remove `@md` (`DESCRIPTION` already enables markdown).
- [ ] Use `@returns`, not `@return`, in `storage.R` and `canonical.R`.
- [ ] `R/data.R`: fill in the empty and TODO entries (e.g. `bt_spas_x_bayes_params`), with `@format` types and `@source`.

Repo files:
- [ ] Convert `README.md` to `README.Rmd`, covering install (`pak::pak("SRJPE/SRJPEmodel")`), test, lint/format, check and document.
- [ ] Add `NEWS.md`.
- [ ] Add a `CLAUDE.md` that points to the infra standards and `planning/style_standards_DRAFT.md` (once approved) and lists the check commands.
- [ ] Remove `cran-comments.md`; the package isn't going to CRAN.

## Things to consider for the migration (Jan-Apr 2027)

**Which code is the source.** For each submodel, decide which repo and commit is canonical (the submodel repo, or the version runSRJPE consumes), and put the commit SHA in the migration commit message. bt-spas-x is published and static; for the others, the submodel repo may keep changing after migration. Decide whether SRJPEmodel or the submodel repo is where changes happen from then on.

**Data has to land in SRJPEdata first.** Most submodel repos read local files (`stock-recruit/data/`, the survival export in `runSRJPE/data/survival/`). Package functions can't read local files, so each migration needs its inputs in an SRJPEdata release first.

**Match runSRJPE's hand-off formats.** runSRJPE consumes file hand-offs now (its README calls this Scenario A) and wants to run the code itself (Scenario B). If each `extract_*` function returns what the hand-off contains (SR: α, β, γ, sd_pro draws; survival: `survival_draws_<site>` + `survival_sites`; in-season: `timing.rds` + `nx_<issue_date>.csv`), runSRJPE's `compute_*` functions barely change.

**Stan interface.** SRJPEmodel uses `rstan`. The infra standards prefer `brms` or `cmdstanr`, and in-season already uses `cmdstanr` (installed from R-universe, which the standards allow). cmdstanr compiles executables, and IT security currently blocks those on at least one machine (see `in-season/CLAUDE.md`). Decide before migrating in-season. `store_model_fit()` only validates `stanfit` and `bugs` objects, so a cmdstanr fit would need a change there.

**Check it reproduces.** Before replacing a hand-off, check that the migrated function gives the same numbers as the submodel repo with a fixed seed. Then save a small version of that check as the submodel's regression test.

**Prep scripts in `archive/`.** Look at them at the break point. Pull in the useful parts as `prepare_*_inputs()` functions and leave the rest.

**PLAD.** PLAD has no fitting code in any repo yet; `read_plad_csv()` (archived with `integrated_model.R`) is the only function. Decide whether it gets functions here.

**Each submodel PR includes:**
- functions in `R/<submodel>.R`
- model code in `model_files/` and the cached data object
- one regression test plus argument-error tests
- the vignette rewritten, with the superseded note removed
- a `_pkgdown.yml` section
- a `NEWS.md` entry

Do one submodel per PR, in whatever order they reach their break point.

## Open questions

1. Kebab-case file names? Recommended yes. It changes one file today.
2. SRJPEdata as a dependency: `Remotes:` (simple, but against the standards), or publish SRJPEdata to R-universe?
3. Sign-off on the three behavior changes marked above (error reporting in `run_abundance_workflow()`, abort instead of `NULL` in `prepare_pCap_inputs()`, renaming the diagnostic plot function).
