# test-loess.R

set.seed(1)  # stabilize any randomness in helpers/LOESS

# helpers
mk_df_single <- function() {
  data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("sample","QC","sample","sample","sample","QC"),
    order  = c(2,1,3,4,5,6),   # sorted -> QC ... QC
    M1 = c(110,100,120,130,140,150),
    M2 = c(55,  50, 52,  60,  58,  65),
    check.names = FALSE
  )
}

mk_df_batches <- function() {
  df <- data.frame(
    sample = paste0("b", rep(1:2, each = 6), "_s", 1:6),
    batch  = rep(c("A","B"), each = 6),
    class  = rep(c("QC","sample","sample","sample","sample","QC"), 2),
    order  = rep(1:6, 2),
    M1 = c(100,105,110,115,120,125, 200,195,190,185,180,175),
    M2 = c( 50, 48, 51, 49, 47, 46, 100,102,104,106,108,110),
    check.names = FALSE
  )
  df$class[df$batch == "B" & df$order %in% 2:5] <- "sample"
  df
}

test_that("loess_correction errors if first/last are not QC", {
  df <- mk_df_single()
  df$class[df$order == 1] <- "sample"
  expect_error(loess_correction(df, metab_cols = c("M1","M2")), "First and last samples must be QCs")
})

test_that("loess_correction sorts by order and returns same shape", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_single()
  out <- loess_correction(df, metab_cols = c("M1","M2"))
  expect_equal(out$order, sort(df$order))
  expect_setequal(names(out), names(df))
  expect_true(all(vapply(out[c("M1","M2")], is.numeric, TRUE)))
})

test_that("loess_correction centers QC near 1 and never negative", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_single()
  expect_silent({ out <- loess_correction(df, metab_cols = c("M1","M2")) })
  qc <- out[out$class == "QC", c("M1","M2")]
  expect_true(all(qc$M1 >= 0 & qc$M2 >= 0))
  expect_lt(abs(median(qc$M1) - 1), 0.35)
  expect_lt(abs(median(qc$M2) - 1), 0.35)
})

test_that("bw_loess_correction errors if a batch does not start/end with QC", {
  df <- mk_df_batches()
  df$class[df$batch == "A" & df$order == 6] <- "sample"
  expect_error(bw_loess_correction(df, metab_cols = c("M1","M2")), "must start and end with QC")
})

test_that("bw_loess_correction warns on too few QCs and corrects with degree=1", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_batches()
  p <- testthat::evaluate_promise(
    bw_loess_correction(df, metab_cols = c("M1","M2"), min_qc = 5)
  )
  # Assert we got the expected warnings, but don't leak them to testthat
  expect_true(any(grepl("Skipping batch 'A'", p$warnings)))
  expect_true(any(grepl("Skipping batch 'B'", p$warnings)))
  tmp <- p$result
  expect_setequal(names(tmp), names(df))
  expect_silent({
    out <- bw_loess_correction(df, metab_cols = c("M1","M2"), min_qc = 2, degree = 1)
  })
  qc_A <- subset(out, batch == "A" & class == "QC", select = c("M1","M2"))
  qc_B <- subset(out, batch == "B" & class == "QC", select = c("M1","M2"))
  expect_true(all(qc_A$M1 >= 0 & qc_A$M2 >= 0 & qc_B$M1 >= 0 & qc_B$M2 >= 0))
  expect_lt(abs(median(qc_A$M1) - 1), 0.35)
  expect_lt(abs(median(qc_A$M2) - 1), 0.35)
  expect_lt(abs(median(qc_B$M1) - 1), 0.35)
  expect_lt(abs(median(qc_B$M2) - 1), 0.35)
})

test_that("bw_loess_correction corrects per-batch; QC centered near 1", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_batches()
  out <- bw_loess_correction(df, metab_cols = c("M1","M2"), min_qc = 2, degree = 1)
  
  qc_A <- subset(out, batch == "A" & class == "QC", select = c("M1","M2"))
  qc_B <- subset(out, batch == "B" & class == "QC", select = c("M1","M2"))
  
  expect_lt(abs(median(qc_A$M1) - 1), 0.35)
  expect_lt(abs(median(qc_A$M2) - 1), 0.35)
  expect_lt(abs(median(qc_B$M1) - 1), 0.35)
  expect_lt(abs(median(qc_B$M2) - 1), 0.35)
})

test_that("bw_loess_correction outputs numeric metabolites with no NaN/Inf", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_batches()
  expect_silent({
    out <- bw_loess_correction(df, metab_cols = c("M1","M2"), min_qc = 2, degree = 1)
  })
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(vapply(out[c("M1","M2")], is.numeric, TRUE)))
})

test_that("loess_correction returns finite, non-negative, no-NA values", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_single()
  expect_silent({ out <- loess_correction(df, metab_cols = c("M1","M2")) })
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(is.na(mets)))
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(mets >= 0))
})

test_that("bw_loess_correction returns finite, non-negative, no-NA values", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_batches()
  expect_silent({
    out <- bw_loess_correction(df, metab_cols = c("M1","M2"), min_qc = 2, degree = 1)
  })
  mets <- as.matrix(out[c("M1","M2")])
  expect_false(any(is.na(mets)))
  expect_false(any(is.nan(mets)))
  expect_false(any(is.infinite(mets)))
  expect_true(all(mets >= 0))
})

test_that("cleanup imputes with smallest positive or 0 fallback (loess_correction)", {
  testthat::skip_if_not_installed("impute")
  df <- data.frame(
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC","sample","sample","sample","sample","QC"),
    order  = 1:6,
    M3     = c(0, 0, 0, 0, 0, 0),
    check.names = FALSE
  )
  out <- loess_correction(df, metab_cols = "M3")
  expect_true(all(out$M3 == 0))
})

test_that("cleanup imputes with smallest positive or 0 fallback (bw_loess_correction)", {
  testthat::skip_if_not_installed("impute")
  df <- data.frame(
    sample = paste0("b", rep(c("A","B"), each = 6), "_s", 1:6),
    batch  = rep(c("A","B"), each = 6),
    class  = rep(c("QC","sample","sample","sample","sample","QC"), 2),
    order  = rep(1:6, 2),
    M4     = c(1, 2, NA, 3, 4, 5,   0, 0, 0, 0, 0, 0),
    check.names = FALSE
  )
  out <- bw_loess_correction(df, metab_cols = "M4", min_qc = 2, degree = 1)
  expect_false(any(is.na(out$M4)))
  expect_true(all(out$M4 >= 0))
  mp <- suppressWarnings(min(out$M4[out$M4 > 0]))
  bvals <- out$M4[out$batch == "B"]
  if (is.finite(mp)) {
    expect_true(all(bvals %in% c(0, mp)))
  } else {
    expect_true(all(bvals == 0))
  }
  min_pos_A <- min(out$M4[out$batch == "A" & out$M4 > 0])
  expect_true(min_pos_A >= 0)
})
