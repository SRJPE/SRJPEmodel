# SRJPEmodel 0.1.0

* Archived the stock-recruit, survival, in-season, passage to spawner and integrated model code into `archive/`. The package before this change is tagged `legacy-pre-migration`. Passage to spawner is retired from the JPE.
* Renamed `generate_diagnostic_plot_juv()` to `plot_juv_diagnostics()`.
* `prepare_pCap_inputs()` now stops with an error on invalid arguments instead of warning and returning `NULL`.
* `run_abundance_workflow()` now reports why a site-year failed, as a warning and in a new `error_message` column.
* Added SRJPEdata to `Imports` (installed from GitHub via `Remotes`) and declared all other imported packages.
* Added regression tests for `prepare_pCap_inputs()`.
