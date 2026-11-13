# prepare_pCap_inputs() errors on bad arguments

    Code
      prepare_pCap_inputs("all_sites", skew = TRUE)
    Condition
      Error in `prepare_pCap_inputs()`:
      ! `skew` must be NULL for the all_sites model type.

---

    Code
      prepare_pCap_inputs("one_site", skew = FALSE)
    Condition
      Error in `prepare_pCap_inputs()`:
      ! `site_selection` must be one site name for the one_site model type.

---

    Code
      prepare_pCap_inputs("one_site", site_selection = "tisdale")
    Condition
      Error in `prepare_pCap_inputs()`:
      ! `skew` must be TRUE or FALSE for the one_site model type.

---

    Code
      prepare_pCap_inputs("all_sites", exclude_from_hyper = "not a site")
    Condition
      Error in `prepare_pCap_inputs()`:
      ! `exclude_from_hyper` contains a site not in the efficiency data.

