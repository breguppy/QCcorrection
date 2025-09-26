# test-group_fold.R

make_df <- function() {
  data.frame(
    sample = paste0("s", 1:4),
    batch  = c(1,1,1,1),
    class  = c("QC","QC","sample","sample"),
    order  = 1:4,
    A = c(1,3,5,7),
    B = c(2,4,6,NA_real_),
    check.names = FALSE
  )
}

test_that("group_stats returns per-group data and summary stats", {
  df <- make_df()
  out <- group_stats(df)
  
  # structure
  expect_named(out, c("group_dfs","group_stats_dfs"))
  expect_setequal(names(out$group_dfs), c("QC","sample"))
  expect_setequal(names(out$group_stats_dfs), c("QC","sample"))
  
  # group_dfs: class renamed to Group; batch/order removed
  g_qc <- out$group_dfs$QC
  expect_true("Group" %in% names(g_qc))
  expect_false("class" %in% names(g_qc))
  expect_false("batch" %in% names(g_qc))
  expect_false("order" %in% names(g_qc))
  expect_setequal(names(g_qc), c("sample","Group","A","B"))
  
  # stats rows and columns
  s_qc <- out$group_stats_dfs$QC
  expect_identical(s_qc$` `, c("Mean","SE","CV"))
  expect_identical(colnames(s_qc), c(" ", "A", "B"))
  
  # QC expected values
  # A: values 1,3 -> mean=2, sd=sqrt(2), SE=1, CV=sd/mean
  expect_equal(s_qc$A[s_qc$` `=="Mean"], 2)
  expect_equal(s_qc$A[s_qc$` `=="SE"], 1, tolerance = 1e-12)
  expect_equal(s_qc$A[s_qc$` `=="CV"], sqrt(2)/2, tolerance = 1e-12)
  # B: values 2,4 -> mean=3, sd=sqrt(2), SE=1, CV=sd/mean
  expect_equal(s_qc$B[s_qc$` `=="Mean"], 3)
  expect_equal(s_qc$B[s_qc$` `=="SE"], 1, tolerance = 1e-12)
  expect_equal(s_qc$B[s_qc$` `=="CV"], sqrt(2)/3, tolerance = 1e-12)
  
  # sample group: A values 5,7 -> mean=6, SE=1, CV=sqrt(2)/6; B values 6,NA -> mean=6, SE=NA, CV=NA
  s_sam <- out$group_stats_dfs$sample
  expect_equal(s_sam$A[s_sam$` `=="Mean"], 6)
  expect_equal(s_sam$A[s_sam$` `=="SE"], 1, tolerance = 1e-12)
  expect_equal(s_sam$A[s_sam$` `=="CV"], sqrt(2)/6, tolerance = 1e-12)
  expect_equal(s_sam$B[s_sam$` `=="Mean"], 6)
  expect_true(is.na(s_sam$B[s_sam$` `=="SE"]))
  expect_true(is.na(s_sam$B[s_sam$` `=="CV"]))
})

test_that("fold_changes divides metabolite columns by control means and preserves metadata", {
  df <- make_df()
  # use QC means as control
  control_mean <- data.frame(
    A = mean(df$A[df$class=="QC"]),
    B = mean(df$B[df$class=="QC"]),
    check.names = FALSE
  )
  
  out <- fold_changes(df[df$class=="sample", ], control_mean)
  
  # metadata unchanged
  expect_setequal(colnames(out), colnames(df))
  expect_identical(out$sample, df$sample[df$class=="sample"])
  expect_identical(out$class,  df$class[df$class=="sample"])
  expect_identical(out$order,  df$order[df$class=="sample"])
  
  # A: 5/2=2.5, 7/2=3.5
  expect_equal(out$A, c(2.5, 3.5))
  # B: 6/3=2, NA stays NA
  expect_equal(out$B, c(2, NA_real_))
})

test_that("fold_changes errors if control_mean does not have exactly one row", {
  df <- make_df()
  bad0 <- df[0, c("A","B")]
  bad2 <- rbind(
    data.frame(A=1,B=1),
    data.frame(A=2,B=2)
  )
  expect_error(fold_changes(df, bad0), "exactly one row")
  expect_error(fold_changes(df, bad2), "exactly one row")
})

test_that("fold_changes handles zero control mean by producing Inf/NaN", {
  df <- make_df()
  control_zero <- data.frame(A = 0, B = 3)
  out <- fold_changes(df[df$class=="sample", ], control_zero)
  expect_true(all(is.infinite(out$A)))     # 5/0, 7/0 -> Inf
  expect_equal(out$B, c(2, NA_real_))      # 6/3=2, NA stays NA
})
