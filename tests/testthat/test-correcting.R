# test-correct_data.R

# --- helpers ---------------------------------------------------------------
mk_df_single <- function(n = 11) {
  stopifnot(n %% 2 == 1)               # ensures last row is QC
  data.frame(
    sample = paste0("s", 1:n),
    batch  = "A",
    class  = ifelse((1:n) %% 2 == 1, "QC", "sample"),
    order  = 1:n,
    M1 = 100 + 2*(1:n),
    M2 =  50 + 1*(1:n),
    check.names = FALSE
  )
}

mk_df_batches_ok <- function(n = 11) {
  stopifnot(n %% 2 == 1)
  one <- function(tag) data.frame(
    sample = paste0(tag, "_s", 1:n),
    batch  = tag,
    class  = ifelse((1:n) %% 2 == 1, "QC", "sample"),
    order  = 1:n,
    M1 = 300 + 2*(1:n),
    M2 = 150 + 1.5*(1:n),
    check.names = FALSE
  )
  rbind(one("A"), one("B"))
}

meta_cols <- c("sample","batch","class","order")
met_cols  <- function(df) setdiff(names(df), meta_cols)

# --- LOESS -----------------------------------------------------------------
test_that("correct_data LOESS returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_single()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "LOESS")
  expect_type(out, "list")
  expect_true(all(c("df","str","parameters") %in% names(out)))
  expect_equal(names(out$df), names(df))
  expect_equal(out$str, "LOESS")
  expect_true(is.character(out$parameters))
  mets <- as.matrix(out$df[met_cols(df)])
  expect_false(any(is.na(mets) | is.nan(mets) | is.infinite(mets)))
  expect_true(all(vapply(out$df[met_cols(df)], is.numeric, TRUE)))
})

test_that("correct_data BW_LOESS returns clean, shaped output", {
  testthat::skip_if_not_installed("impute")
  df <- mk_df_batches_ok()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_LOESS")
  expect_equal(names(out$df), names(df))
  expect_equal(out$str, "Batchwise LOESS")
  mets <- as.matrix(out$df[met_cols(df)])
  expect_false(any(is.na(mets) | is.nan(mets) | is.infinite(mets)))
})

# --- RF --------------------------------------------------------------------
test_that("correct_data RF equals median across the three RF seeds", {
  testthat::skip_if_not_installed("randomForest")
  df <- mk_df_single()  # ≥5 QCs, first/last are QC
  seeds <- c(42, 31416, 272)
  # manual median across models
  dfs <- lapply(seeds, function(s) rf_correction(df, met_cols(df), ntree = 500, seed = s))
  # compute per-cell median for metabolite columns
  mdm <- dfs[[1]]
  mdm[met_cols(df)] <- Map(function(...) apply(cbind(...), 1, stats::median),
                           dfs[[1]][met_cols(df)], dfs[[2]][met_cols(df)], dfs[[3]][met_cols(df)])
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "RF")
  expect_equal(out$str, "Random Forest")
  expect_true(grepl("median value of the 3 models", out$parameters))
  expect_equal(out$df[met_cols(df)], mdm[met_cols(df)], tolerance = 1e-8)
  expect_equal(names(out$df), names(df))
})

test_that("correct_data BW_RF equals median across batch-wise RF seeds", {
  testthat::skip_if_not_installed("randomForest")
  df <- mk_df_batches_ok()
  seeds <- c(42, 31416, 272)
  
  dfs <- lapply(seeds, function(s) bw_rf_correction(df, met_cols(df), ntree = 500, seed = s))
  
  mdm <- dfs[[1]]
  mdm[met_cols(df)] <- Map(function(...) apply(cbind(...), 1, stats::median),
                           dfs[[1]][met_cols(df)], dfs[[2]][met_cols(df)], dfs[[3]][met_cols(df)])
  
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "BW_RF")
  
  # align row order before compare
  ord_out <- order(out$df$batch, out$df$order)
  ord_mdm <- order(mdm$batch, mdm$order)
  
  expect_equal(out$str, "Batchwise Random Forest")
  expect_true(grepl("median value of the 3 models", out$parameters))
  expect_equal(out$df[ord_out, met_cols(df)], mdm[ord_mdm, met_cols(df)], tolerance = 1e-8)
  expect_setequal(names(out$df), names(df))
})

# --- Generic properties ----------------------------------------------------
test_that("correct_data preserves metadata columns and types", {
  # use LOESS to avoid RF dependency here
  testthat::skip_if_not_installed("impute")
  df <- mk_df_single()
  out <- correct_data(df, metab_cols = met_cols(df), corMethod = "LOESS")
  expect_true(all(meta_cols %in% names(out$df)))
  expect_identical(out$df[meta_cols], df[meta_cols])
})
