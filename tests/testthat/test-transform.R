# test-transform.R

#  helpers 
mk_filtered_corrected <- function() {
  df <- data.frame(
    sample = paste0("s", 1:3),
    batch  = 1L,
    class  = c("QC","sample","sample"),
    order  = 1:3,
    A = c(1, 2, 0),
    B = c(1, NA, 0),
    C = c(2, 2, 0),
    ISTD_Met1 = c(100, 120, 140),
    ITSD_X    = c(50,  60,  70),
    check.names = FALSE
  )
  
  list(
    df_mv    = df,           # has NA in B row 2
    df_no_mv = df            # for unit tests it's fine if same; adjust if desired
  )
}

meta_cols  <- c("sample","batch","class","order")

test_that(".total_ratio_norm scales rows by non_missing_count / row_sum", {
  df <- mk_filtered_corrected()$df_mv
  # use only true metabolites (exclude ISTD here)
  mc <- c("A","B","C")
  
  out <- .total_ratio_norm(df, mc)
  
  # Row 1: A=1,B=1,C=2 -> sum=4, non-missing=3, ratio=3/4=0.75
  expect_equal(out$A[1], 1 * 0.75, tolerance = 1e-12)
  expect_equal(out$B[1], 1 * 0.75, tolerance = 1e-12)
  expect_equal(out$C[1], 2 * 0.75, tolerance = 1e-12)
  
  # Row 2: A=2,B=NA,C=2 -> sum=4, non-missing=2, ratio=0.5
  expect_equal(out$A[2], 2 * 0.5, tolerance = 1e-12)
  expect_true(is.na(out$B[2]))
  expect_equal(out$C[2], 2 * 0.5, tolerance = 1e-12)
  
  # Row 3: A=0,B=0,C=0 -> sum=0, non-missing=3, ratio=Inf
  # 0 * Inf -> NaN in R
  expect_true(is.nan(out$A[3]))
  expect_true(is.nan(out$B[3]))
  expect_true(is.nan(out$C[3]))
})

test_that(".istd_norm divides non-ISTD metabolites by row-wise ISTD mean", {
  df <- mk_filtered_corrected()$df_mv
  
  istd_cols <- c("ISTD_Met1", "ITSD_X")
  metab_cols <- c("A","B","C")
  
  # Row means of ISTDs:
  # r1: (100+50)/2 = 75
  # r2: (120+60)/2 = 90
  # r3: (140+70)/2 = 105
  out <- .istd_norm(df, metab_cols = metab_cols, istd_cols = istd_cols, na_action = "leave")
  
  expect_equal(out$A, df$A / c(75, 90, 105), tolerance = 1e-12)
  expect_equal(out$C, df$C / c(75, 90, 105), tolerance = 1e-12)
  
  # NA stays NA
  expect_true(is.na(out$B[2]))
  
  # ISTD cols not modified
  expect_equal(out$ISTD_Met1, df$ISTD_Met1)
  expect_equal(out$ITSD_X, df$ITSD_X)
})

test_that(".istd_norm errors if no ISTD columns are provided", {
  df <- mk_filtered_corrected()$df_mv
  expect_error(
    .istd_norm(df, metab_cols = c("A","B","C"), istd_cols = character(0)),
    "no ISTD/ITSD columns"
  )
})

test_that(".istd_norm na_action='leave' keeps original values when ISTD mean is NA", {
  df <- mk_filtered_corrected()$df_mv
  df$ISTD_Met1[2] <- NA
  df$ITSD_X[2]    <- NA
  
  out <- .istd_norm(
    df,
    metab_cols = c("A","B","C"),
    istd_cols  = c("ISTD_Met1","ITSD_X"),
    min_istd   = 1L,
    na_action  = "leave"
  )
  
  # row 2 should remain unchanged for A/B/C
  expect_equal(out$A[2], df$A[2])
  expect_true(is.na(out$B[2]))
  expect_equal(out$C[2], df$C[2])
})

test_that(".istd_norm na_action='na' sets normalized values to NA when ISTD mean is NA", {
  df <- mk_filtered_corrected()$df_mv
  df$ISTD_Met1[2] <- NA
  df$ITSD_X[2]    <- NA
  
  out <- .istd_norm(
    df,
    metab_cols = c("A","B","C"),
    istd_cols  = c("ISTD_Met1","ITSD_X"),
    min_istd   = 1L,
    na_action  = "na"
  )
  
  expect_true(is.na(out$A[2]))
  expect_true(is.na(out$B[2]))
  expect_true(is.na(out$C[2]))
})

test_that(".istd_norm na_action='error' throws when ISTD mean is NA", {
  df <- mk_filtered_corrected()$df_mv
  df$ISTD_Met1[2] <- NA
  df$ITSD_X[2]    <- NA
  
  expect_error(
    .istd_norm(
      df,
      metab_cols = c("A","B","C"),
      istd_cols  = c("ISTD_Met1","ITSD_X"),
      min_istd   = 1L,
      na_action  = "error"
    ),
    "Cannot compute ISTD mean"
  )
})


test_that("transform_data('none') returns unmodified dfs and sets withheld cols when eITSD_X=TRUE", {
  fc <- mk_filtered_corrected()
  
  res <- transform_data(
    filtered_corrected = fc,
    transform = "none",
    withheld_cols = character(0),
    ex_ISTD = TRUE
  )
  
  expect_true(startsWith(res$str, "After correction, no scaling"))
  
  expect_equal(res$df_mv, fc$df_mv)
  expect_equal(res$df_no_mv, fc$df_no_mv)
  
  expect_setequal(res$withheld_cols_mv, c("ISTD_Met1", "ITSD_X"))
  expect_setequal(res$withheld_cols_no_mv, c("ISTD_Met1", "ITSD_X"))
})

test_that("transform_data('none') does not add ISTDs to withheld cols when eITSD_X=FALSE", {
  fc <- mk_filtered_corrected()
  
  res <- transform_data(
    filtered_corrected = fc,
    transform = "none",
    withheld_cols = character(0),
    ex_ISTD = FALSE
  )
  
  expect_identical(res$withheld_cols_mv, character(0))
  expect_identical(res$withheld_cols_no_mv, character(0))
  expect_equal(res$df_mv, fc$df_mv)
})

test_that("transform_data('ISTD_norm') normalizes non-ISTD metabolites by ISTD mean and does not change ISTD cols", {
  fc <- mk_filtered_corrected()
  
  res <- transform_data(
    filtered_corrected = fc,
    transform = "ISTD_norm",
    withheld_cols = character(0),
    ex_ISTD = TRUE
  )
  
  # expected ISTD means
  istd_mean <- c(75, 90, 105)
  
  expect_equal(res$df_no_mv$A, fc$df_no_mv$A / istd_mean, tolerance = 1e-12)
  expect_equal(res$df_no_mv$C, fc$df_no_mv$C / istd_mean, tolerance = 1e-12)
  expect_true(is.na(res$df_no_mv$B[2]))
  expect_equal(res$df_no_mv$ISTD_Met1, fc$df_no_mv$ISTD_Met1)
  expect_equal(res$df_no_mv$ITSD_X, fc$df_no_mv$ITSD_X)
  expect_equal(res$df_no_mv[, meta_cols], fc$df_no_mv[, meta_cols])
  expect_true(grepl("normalized to the average", res$str))
})

test_that("transform_data('ISTD_norm') errors when no ISTD/ITSD columns exist", {
  fc <- mk_filtered_corrected()
  # remove ISTD cols
  fc$df_mv    <- fc$df_mv[, setdiff(names(fc$df_mv), c("ISTD_Met1","ITSD_X")), drop = FALSE]
  fc$df_no_mv <- fc$df_no_mv[, setdiff(names(fc$df_no_mv), c("ISTD_Met1","ITSD_X")), drop = FALSE]
  
  expect_error(
    transform_data(fc, transform = "ISTD_norm", withheld_cols = character(0), ex_ISTD = TRUE),
    "no ISTD/ITSD columns"
  )
})

test_that("transform_data('TRN') applies .total_ratio_norm to non-withheld columns and leaves withheld unchanged", {
  fc <- mk_filtered_corrected()
  
  # Withhold column C and exclude ISTDs from TRN
  res <- transform_data(
    filtered_corrected = fc,
    transform = "TRN",
    withheld_cols = "C",
    ex_ISTD = TRUE
  )
  
  exp_no_mv <- fc$df_no_mv
  exp_no_mv <- .total_ratio_norm(exp_no_mv, metab_cols = c("A","B"))
  
  expect_equal(res$df_no_mv$A, exp_no_mv$A, tolerance = 1e-12)
  expect_equal(res$df_no_mv$B, exp_no_mv$B, tolerance = 1e-12)
  
  # withheld unchanged
  expect_equal(res$df_mv$C, fc$df_mv$C)
  expect_equal(res$df_mv$ISTD_Met1, fc$df_mv$ISTD_Met1)
  expect_equal(res$df_mv$ITSD_X, fc$df_mv$ITSD_X)
  
  expect_setequal(res$withheld_cols_mv, c("C","ISTD_Met1","ITSD_X"))
  expect_true(grepl("ratiometrically normalized", res$str))
})

test_that("transform_data preserves all columns in df_mv/df_no_mv outputs", {
  fc <- mk_filtered_corrected()
  
  res <- transform_data(fc, transform = "ISTD_norm", withheld_cols = character(0), ex_ISTD = TRUE)
  
  expect_setequal(names(res$df_mv), names(fc$df_mv))
  expect_setequal(names(res$df_no_mv), names(fc$df_no_mv))
})

