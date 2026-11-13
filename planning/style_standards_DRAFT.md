# SRJPEmodel style standards (DRAFT)

Draft for review, 2026-11-12.

SRJPEmodel follows the FlowWest R standards in the `FlowWest/infra` repo (`standards/r-coding-standards.md` and `standards/r-data-analysis-data-science-standards.md`, currently a draft on the `add-r-standards` branch). They use **MUST/NEVER**, **SHOULD/AVOID** and **MAY/CONSIDER** the same way this file does.

This file covers only what's specific to SRJPEmodel: where the infra standards leave a choice open, and the conventions this package already has. Items marked **[decide]** still need your call.

## The main rule

Keep it simple. A new function should look like the ones next to it. Don't add helper layers, config objects, defensive wrappers or long comments unless something in the standards calls for them. In the infra repo's words: "Prefer boring, readable, validated, tested R over clever R."

## Files

- One file per submodel in `R/`, named for the submodel. Shared code goes in `storage.R`, `canonical.R`, `plot.R` and `data.R`.
- **[decide]** The infra standards name files in kebab-case (`juvenile-abundance.R`). SRJPEmodel uses snake_case. Only `juvenile_abundance.R` would change today, so I recommend switching to kebab-case during the cleanup.
- Tests mirror `R/` one file to one file: `R/juvenile-abundance.R` is tested by `tests/testthat/test-juvenile-abundance.R`.
- Stan/BUGS code goes in `model_files/`. `data-raw/cache-models.R` reads it into a data object with `usethis::use_data()` (e.g. `bt_spas_x_model_code`). Functions use the data object, not the file.
- Fixed MCMC settings and seeds go in a data object built in `data-raw/cache-params.R` (e.g. `bt_spas_x_bayes_params`).
- Scripts for exploring, testing or uploading go in `data-raw/<submodel>/`. They're never part of the package.

## Function names

Each submodel follows the same steps:

```
prepare_<x>_inputs()  ->  fit_<x>_model()  ->  extract_<x>_estimates()
plot_<x>_results()    plot_<x>_diagnostics()
```

- Plot functions start with `plot_`, never `generate_*_plot_*`.
- `run_*()` only for wrappers that chain the steps over many sites/years (e.g. `run_bt_spas_x_JPE_sites()`).
- Names are `snake_case` verbs. Never use dots in names. A leading dot (`.approved_model_names`) marks something internal.
- **[decide]** Strict snake_case would be `pcap`, not `pCap`. I recommend leaving existing exported names alone, since bt-spas-x, in-season and runSRJPE call them, and using lowercase in new names.

## Arguments and return values

- Use the same argument names as existing functions: `site`, `stream`, `run_year`, `con` for the database connection.
- `prepare_*` returns a named list that is passed straight to `fit_*`. Usually that's `inputs` (with `data`, `inits`, `parameters`) plus whatever `extract_*` needs later.
- `fit_*` returns the raw fit object.
- `extract_*` returns a tibble.
- Input data comes from `SRJPEdata::` by default. An optional argument can override it (e.g. `input_catch_data = NULL`), which also makes testing with small fixtures easy.
- Never read local files in package code. Never use `here::here()` in package code.
- New columns put units in the name (`flow_cfs`, `fork_length_mm`). Leave existing SRJPEdata column names as they are.

## Code

Air handles formatting, including `if (` spacing. Beyond that:

- Native pipe `|>`, never `%>%`. `<-` for assignment. `TRUE`/`FALSE`, never `T`/`F`.
- Prefix functions from other packages (`dplyr::filter()`), or use `@importFrom` for ones used constantly. Never use `@import` for a whole package, and never call `library()` inside a function.
- Check arguments at the top of the function with `if (...) cli::cli_abort()`, naming the bad argument and the allowed values. Don't warn and return `NULL`. Don't add checkmate; plain `if` checks are enough here.
- Use `cli::cli_inform()` / `cli::cli_alert_info()` for messages, never `print()` or `cat()`.
- Use `return()` only to exit early. The last line of the function is what it returns.
- Don't change global state. If a function has to set `options()` (e.g. `mc.cores`), restore it with `on.exit()`.
- Don't catch errors just to keep going. If a loop over sites has to continue past a failure, keep the error message and report it.
- Only suppress a warning if you know exactly what it is, with a comment saying why it's safe.
- Joins declare `relationship =` in new code.
- When rows are dropped (`!is.na()` filters, `na.rm = TRUE`), it's on purpose, and the function reports how many when it matters. Never replace `NA` with 0. A missing catch is not zero catch.
- Code that uses random numbers sets a seed. Fitting functions take the seed from the package's params data object.
- Comments say *why*, not *what*: `# keep knots away from the tails, which can have unsampled weeks`, not `# set knot positions`.

## Documentation (roxygen)

`prepare_stock_recruit_inputs()` (in `archive/` after the cleanup) is the reference for length. `juvenile_abundance.R` is the most you should write.

```r
#' Prepare inputs for the stock-recruit model
#'
#' @details One to three sentences on what the function does and any
#'   assumption that isn't obvious.
#' @param stream Character. The stream to fit the model for.
#' @param covariate Character. One of `"a"`, `"b"`, or `"null"`.
#' @returns A named list:
#' * **stock_recruit_table** Tibble of juvenile and adult abundance by brood year.
#' * **covariate_table** Tibble of covariates for the same years.
#' @family Prepare Model Inputs
#' @export
```

- Title line, `@details`, one `@param` per argument (type, units, allowed values), `@returns` (type, plus columns or list elements), `@family`, `@export`.
- No `@md`: `DESCRIPTION` already turns markdown on for the whole package.
- Use `@inheritParams` rather than copying a parameter's description.
- Document internal functions with `@noRd`, or not at all if they're short and obvious.
- **[decide]** Examples: the infra standards say examples should run under `R CMD check`. Most model functions need the database, WinBUGS or a long fit, so I recommend no `@examples` on those and putting usage in the vignette instead. Small functions that run quickly get a real example.
- Datasets in `R/data.R` get `@format` (each element or column, its type and units) and `@source`.
- Explain the model itself in the submodel's vignette, the manuscript or the submodel repo, not in roxygen. Link rather than repeat.
- File headers are one to four comment lines saying what the file is for. No lists of the functions in the file.
- Never edit `NAMESPACE` or `man/` by hand. Run `devtools::document()`.

## Tests

- testthat 3rd edition. Each migrated submodel gets **one small regression test**: run `prepare_<x>_inputs()` on a small fixture in `tests/testthat/fixtures/` and compare the result to a saved expected output. Tests don't run MCMC fitting.
- Each submodel also gets a test that bad arguments raise an error, using `expect_snapshot(error = TRUE)`.
- Set seeds with `withr::local_seed()` if a test uses randomness.

## Vignettes and pkgdown

- One vignette per submodel: what the model does, a diagram, how to run it with the functions. Link to the manuscript or submodel repo for full detail.
- A vignette for code that isn't in the package yet starts with a "Superseded" note and has `eval = FALSE` set for all chunks.
- Each submodel gets a `subtitle` section under "Model and Submodels Functions" in `_pkgdown.yml`, listing its functions in pipeline order.

## Commits and PRs

- Commit messages are lowercase and present tense, starting with a verb: `adds ...`, `fixes ...`, `updates ...`.
- Keep PRs small, with one submodel or one change each. Keep refactoring separate from changes that affect model results.
- The PR description says why. If results changed, it says how, with before and after values.
- Before a PR: `air format .`, `lintr::lint_package()`, `devtools::document()`, `devtools::check()`.

## Specific to Claude

The infra standards' rules for AI agents apply in full (Section 25 of the coding standards, Section 11 of the data-analysis standards). Above all:

- Change only what was asked. No renaming, reformatting or reorganizing nearby code.
- No new dependencies without approval.
- Never change thresholds, filters, priors, units or how `NA`s are handled unless asked. If something in the data looks wrong, point it out instead of fixing it.
- Never weaken a test, skip a check, or suppress a warning to make something pass.

Also for this package:

- No new helper functions, wrappers or abstractions unless the same code is repeated in several places.
- Roxygen and comments at the length above. No file headers longer than four lines, verification tables, or design-rationale essays in the package.
- When porting code from a submodel repo, keep the original variable names and math unless asked to change them. Put the source repo and commit SHA in the commit message, not in the code.
