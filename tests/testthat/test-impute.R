# test-impute_missing.R
make_df <- function() {
  data.frame(
    rid    = seq_len(6),        
    sample = paste0("s", 1:6),
    batch  = 1L,
    class  = c("QC","QC","QC","sample","sample","sample"),
    order  = c(3,1,6,2,5,4),                      
    A = c(1, NA, 3,  NA,  5, NA),
    B = c(NA, 2, NA,  4, NA,  6),
    C = c(10,20,30, 40, 50, 60),
    stringsAsFactors = FALSE
  )
}

test_that("impute_missing counts NAs and preserves metadata and order", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "median", samImputeM = "mean"
  )
  
  # n_missv equals total NAs in selected metabolite columns
  expect_equal(out$n_missv, sum(is.na(df[metab_cols])))
  
  # output contains same columns, metadata unchanged
  expect_setequal(names(out$df), names(df))
  expect_true(all(out$df$class %in% c("QC","sample")))
  
  # rows are sorted by order
  expect_equal(out$df$order, sort(df$order))
})

test_that("median and mean strategies fill NAs appropriately", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "median", samImputeM = "mean"
  )
  imp <- out$df
  
  # map output rows back to input via rid
  map <- match(df$rid, imp$rid)
  
  # QC stats from QC subset only
  qc_idx <- df$class == "QC"
  qc_med_A <- median(df$A[qc_idx], na.rm = TRUE)
  qc_med_B <- median(df$B[qc_idx], na.rm = TRUE)
  
  # Sample stats from Sample subset only
  sam_idx <- !qc_idx
  sam_mean_A <- mean(df$A[sam_idx], na.rm = TRUE)
  sam_mean_B <- mean(df$B[sam_idx], na.rm = TRUE)
  
  # positions that were NA originally
  was_na_qc_A  <- which(qc_idx & is.na(df$A))
  was_na_qc_B  <- which(qc_idx & is.na(df$B))
  was_na_sam_A <- which(sam_idx & is.na(df$A))
  was_na_sam_B <- which(sam_idx & is.na(df$B))
  
  expect_true(all(imp$A[map[was_na_qc_A]]  == qc_med_A))
  expect_true(all(imp$B[map[was_na_qc_B]]  == qc_med_B))
  expect_true(all(imp$A[map[was_na_sam_A]] == sam_mean_A))
  expect_true(all(imp$B[map[was_na_sam_B]] == sam_mean_B))
  
  # no new NAs introduced in metabolite cols
  expect_false(any(is.na(imp[metab_cols])))
})

test_that("class_* strategies operate per class", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "class_median", samImputeM = "class_mean"
  )
  
  # strings returned
  expect_identical(out$qc_str,  "class-metabolite median")
  expect_identical(out$sam_str, "class-metabolite mean")
  
  # no NAs remain in metabolite columns
  expect_false(any(is.na(out$df[metab_cols])))
})

test_that("min and minHalf fill with subset minima", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "min", samImputeM = "minHalf"
  )
  
  imp <- out$df
  map <- match(df$rid, imp$rid)
  
  qc_idx  <- df$class == "QC"
  sam_idx <- !qc_idx
  
  # mins computed on each subset, as the function does
  qc_min  <- vapply(df[qc_idx,  metab_cols], min, numeric(1), na.rm = TRUE)
  sam_min <- vapply(df[sam_idx, metab_cols], min, numeric(1), na.rm = TRUE)
  
  # check representative NA positions
  expect_true(all(imp$A[map[which(qc_idx & is.na(df$A))]]  == qc_min["A"]))
  expect_true(all(imp$B[map[which(qc_idx & is.na(df$B))]]  == qc_min["B"]))
  expect_true(all(imp$A[map[which(sam_idx & is.na(df$A))]] == 0.5 * sam_min["A"]))
  expect_true(all(imp$B[map[which(sam_idx & is.na(df$B))]] == 0.5 * sam_min["B"]))
  
  expect_false(any(is.na(imp[metab_cols])))
})

test_that("zero strategy sets missing to 0", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "zero", samImputeM = "zero"
  )
  
  expect_true(all(out$df$A %in% c(0, na.omit(df$A))))
  expect_true(all(out$df$B %in% c(0, na.omit(df$B))))
})

test_that("nothing_to_impute leaves data unchanged and returns labels", {
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "nothing_to_impute", samImputeM = "nothing_to_impute"
  )
  
  # sort both by order and drop rowname attributes
  ord <- order(df$order)
  got <- out$df[order(out$df$order), names(df), drop = FALSE]
  exp <- df[ord, names(df), drop = FALSE]
  rownames(got) <- rownames(exp) <- NULL
  
  expect_equal(got, exp, ignore_attr = TRUE)
  expect_identical(out$qc_str,  "nothing to impute")
  expect_identical(out$sam_str, "nothing to impute")
})

test_that("KNN path runs when impute is installed", {
  testthat::skip_if_not_installed("impute")
  df <- make_df()
  metab_cols <- c("A","B","C")
  
  out <- impute_missing(
    df, metab_cols, qcImputeM = "KNN", samImputeM = "KNN"
  )
  
  # All metabolite NAs should be imputed
  expect_false(any(is.na(out$df[metab_cols])))
  
  # numeric matrix returned back to data.frame
  expect_true(all(vapply(out$df[metab_cols], is.numeric, TRUE)))
  
  # labels
  expect_identical(out$qc_str,  "KNN")
  expect_identical(out$sam_str, "KNN")
})
