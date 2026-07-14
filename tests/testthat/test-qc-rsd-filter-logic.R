test_that("QC RSD filter flag helper supports positive and legacy meanings", {
  expect_true(.qc_rsd_filter_enabled(list(remove_qc_rsd_filter = TRUE), Inf))
  expect_false(.qc_rsd_filter_enabled(list(remove_qc_rsd_filter = FALSE), 30))

  expect_true(.qc_rsd_filter_enabled(list(post_cor_filter = FALSE), Inf))
  expect_false(.qc_rsd_filter_enabled(list(post_cor_filter = TRUE), 30))

  expect_true(.qc_rsd_filter_enabled(list(), 30))
  expect_false(.qc_rsd_filter_enabled(list(), Inf))
})

test_that("report text mentions QC RSD filtering only when enabled", {
  enabled <- report_text_scatter_intro(
    list(
      rsd_cal = "met",
      remove_qc_rsd_filter = TRUE,
      rsd_cutoff = 30
    ),
    list()
  )
  disabled <- report_text_scatter_intro(
    list(
      rsd_cal = "met",
      remove_qc_rsd_filter = FALSE,
      rsd_cutoff = 30
    ),
    list()
  )

  enabled_text <- paste(as.character(enabled), collapse = " ")
  disabled_text <- paste(as.character(disabled), collapse = " ")

  expect_match(enabled_text, "filtered out of the post-corrected dataset")
  expect_false(grepl("filtered out of the post-corrected dataset", disabled_text))
})
