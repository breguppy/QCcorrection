# test-transform.R

#  helpers 
mk_df <- function() {
  data.frame(
    sample = paste0("s", 1:3),
    batch  = 1L,
    class  = c("QC","sample","sample"),
    order  = 1:3,
    A = c(1, 2, 0),
    B = c(1, NA, 0),
    C = c(2, 2, 0),
    Met1_ISTD = c(100, 120, 140),
    X_ITSD    = c(50,  60,  70),
    check.names = FALSE
  )
}

meta_cols  <- c("sample","batch","class","order")
metab_cols <- c("A","B","C","Met1_ISTD","X_ITSD")

test_that(".total_ratio_norm scales rows by non_missing_count / row_sum", {
  df <- mk_df()
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

test_that("transform_data('none') returns unmodified values for metabolites, plus text", {
  df <- mk_df()
  # ex_ISTD=TRUE should exclude ISTD/ITSD from transformed set
  res <- transform_data(df, transform = "none", withheld_cols = character(), ex_ISTD = TRUE)
  
  expect_true(startsWith(res$str, "After correction, no scaling"))
  # ISTD columns are not in metab set of output frame
  expect_false("Met1_ISTD" %in% names(res$df))
  expect_false("X_ITSD"    %in% names(res$df))
  # original metabolite values preserved
  expect_equal(res$df$A, df$A)
  expect_equal(res$df$B, df$B)
  expect_equal(res$df$C, df$C)
  # withheld_cols augmented with ISTD/ITSD
  expect_setequal(res$withheld_cols, c("Met1_ISTD","X_ITSD"))
})

test_that("transform_data('log2') applies log2 to metabolites only and preserves metadata", {
  df <- mk_df()
  res <- transform_data(df, transform = "log2", withheld_cols = "pre_withheld", ex_ISTD = TRUE)
  
  # metadata present
  expect_setequal(names(res$df)[names(res$df) %in% meta_cols], meta_cols)
  # ISTD excluded
  expect_false("Met1_ISTD" %in% names(res$df))
  expect_false("X_ITSD"    %in% names(res$df))
  # log2 applied where finite
  idx <- which(!is.na(df$A))
  expect_equal(res$df$A[idx], log(df$A[idx], 2), tolerance = 1e-12)
  # NA stays NA
  expect_true(is.na(res$df$B[2]))
  # withheld includes prior + ISTD/ITSD
  expect_setequal(res$withheld_cols, c("pre_withheld","Met1_ISTD","X_ITSD"))
})

test_that("transform_data('TRN') uses .total_ratio_norm on non-ISTD metabolites", {
  df <- mk_df()
  # Compute expected by calling the helper on A,B,C only
  base <- df[, c(meta_cols, "A","B","C"), drop = FALSE]
  exp  <- .total_ratio_norm(base, c("A","B","C"))
  
  res <- transform_data(df, transform = "TRN", withheld_cols = character(), ex_ISTD = TRUE)
  
  # same columns as base (ISTD removed)
  expect_setequal(names(res$df), names(exp))
  # values equal to helper output
  expect_equal(res$df[, c("A","B","C")], exp[, c("A","B","C")])
  # description string present
  expect_true(grepl("ratiometrically normalized", res$str))
})

test_that("transform_data respects ex_ISTD=FALSE and keeps ISTD columns", {
  df <- mk_df()
  res <- transform_data(df, transform = "none", withheld_cols = character(), ex_ISTD = FALSE)
  expect_true(all(c("Met1_ISTD","X_ITSD") %in% names(res$df)))
  # withheld_cols unchanged
  expect_identical(res$withheld_cols, character(0))
})
