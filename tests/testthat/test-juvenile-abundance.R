# Fixtures are a small subset of SRJPEdata (ubc, deer creek, mill creek and
# tisdale, run years 2019-2023) so results don't change with SRJPEdata releases.
# Expected outputs were saved from SRJPEmodel at tag legacy-pre-migration.
catch <- readRDS(test_path("fixtures", "catch.rds"))
efficiency <- readRDS(test_path("fixtures", "efficiency.rds"))

test_that("prepare_pCap_inputs() matches saved output for all_sites", {
  result <- prepare_pCap_inputs(
    "all_sites",
    input_catch_data = catch,
    input_efficiency_data = efficiency
  )
  expect_equal(result, readRDS(test_path("fixtures", "expected_pcap_all_sites.rds")))
})

test_that("prepare_pCap_inputs() matches saved output for one_site", {
  result <- prepare_pCap_inputs(
    "one_site",
    skew = FALSE,
    site_selection = "tisdale",
    input_catch_data = catch,
    input_efficiency_data = efficiency
  )
  expect_equal(result, readRDS(test_path("fixtures", "expected_pcap_one_site.rds")))
})

test_that("prepare_pCap_inputs() errors on bad arguments", {
  expect_snapshot(prepare_pCap_inputs("all_sites", skew = TRUE), error = TRUE)
  expect_snapshot(prepare_pCap_inputs("one_site", skew = FALSE), error = TRUE)
  expect_snapshot(prepare_pCap_inputs("one_site", site_selection = "tisdale"), error = TRUE)
  expect_snapshot(prepare_pCap_inputs("all_sites", exclude_from_hyper = "not a site"), error = TRUE)
  expect_snapshot(prepare_pCap_inputs("all_sites", log_transform_flow = "yes"), error = TRUE)
})
