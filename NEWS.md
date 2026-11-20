# SRJPEmodel 0.1.0

* Archived the stock-recruit, survival, in-season, passage to spawner and integrated model code into `archive/`. The package before this change is tagged `legacy-pre-migration`. Passage to spawner is retired from the JPE.
* Renamed `generate_diagnostic_plot_juv()` to `plot_juv_diagnostics()`.
* `prepare_pCap_inputs()` now stops with an error on invalid arguments instead of warning and returning `NULL`.
* Removed `run_abundance_workflow()` and `run_bt_spas_x_JPE_sites()`. BT-SPAS-X runs for all JPE sites are handled in runSRJPE.
* Added SRJPEdata to `Imports` (installed from GitHub via `Remotes`) and declared all other imported packages.
* Added regression tests for `prepare_pCap_inputs()`.
* Added `log_transform_flow` to `prepare_pCap_inputs()` and `prepare_abundance_inputs()` to fit BT-SPAS-X on log flow. Default `FALSE`.
