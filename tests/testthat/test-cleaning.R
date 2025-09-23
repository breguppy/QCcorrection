test_that("clean_data basic cleaning and outputs", {
  # toy input (unsorted; QC should start and end after sort + normalization)
  df <- data.frame(
    SampleID  = paste0("s", 1:6),
    BatchID   = c(1, 1, 1, 1, 1, 1),
    Type      = c(NA, "qc", "qc", "Qc", "sample", "QC"),
    Injection = c(3, 1, 6, 2, 5, 4),
    met1      = c("1", "a", "0", "5", "3", NA),
    met2      = c(2, "7", "foo", 0, NA, 1),
    met3      = c("0", "0", "0", "0", "0", "0"),
    note      = c("x", "y", "z", "w", "v", "u"),
    stringsAsFactors = FALSE
  )
  
  out <- clean_data(
    df,
    sample = "SampleID",
    batch  = "BatchID",
    class  = "Type",
    order  = "Injection",
    withheld_cols = c("note")
  )
  
  # structure
  expect_type(out, "list")
  expect_named(out, c("df", "replacement_counts", "withheld_cols"))
  expect_equal(out$withheld_cols, "note")
  
  # columns renamed, withheld removed, order applied
  expect_setequal(names(out$df),
                  c("sample", "batch", "class", "order", "met1", "met2", "met3"))
  expect_true(is.unsorted(df$Injection))
  expect_equal(out$df$order, sort(out$df$order))
  
  # class normalization to "QC"
  expect_true(all(out$df$class %in% c("QC", "sample")))
  expect_identical(out$df$class[1], "QC")
  expect_identical(out$df$class[nrow(out$df)], "QC")
  
  # numeric coercion + zero→NA
  expect_true(all(vapply(out$df[c("met1", "met2", "met3")], is.numeric, TRUE)))
  expect_true(any(is.na(out$df$met1)))
  expect_true(any(is.na(out$df$met2)))
  expect_true(all(is.na(out$df$met3)))
  
  # replacement counts
  rc <- out$replacement_counts
  # met1: one non-numeric ("a"), one zero ("0")
  expect_equal(rc$non_numeric_replaced[rc$metabolite == "met1"], 1)
  expect_equal(rc$zero_replaced[rc$metabolite == "met1"], 1)
  # met2: one non-numeric ("foo"), one zero
  expect_equal(rc$non_numeric_replaced[rc$metabolite == "met2"], 1)
  expect_equal(rc$zero_replaced[rc$metabolite == "met2"], 1)
  # met3: six zeros
  expect_equal(rc$non_numeric_replaced[rc$metabolite == "met3"], 0)
  expect_equal(rc$zero_replaced[rc$metabolite == "met3"], 6)
})

test_that("clean_data errors when first sample after sort is not QC", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3"),
    BatchID   = 1,
    Type      = c("sample", "QC", "QC"),
    Injection = c(1, 2, 3),
    met1      = c(1, 2, 3)
  )
  expect_error(
    clean_data(
      df,
      "SampleID",
      "BatchID",
      "Type",
      "Injection",
      withheld_cols = character()
    ),
    "begin with a QC sample"
  )
})

test_that("clean_data errors when last sample after sort is not QC", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3"),
    BatchID   = 1,
    Type      = c("QC", "QC", "sample"),
    Injection = c(1, 2, 3),
    met1      = c(1, 2, 3)
  )
  expect_error(
    clean_data(
      df,
      "SampleID",
      "BatchID",
      "Type",
      "Injection",
      withheld_cols = character()
    ),
    "end with a QC sample"
  )
})

test_that("withheld cols are not present in df but are tracked", {
  df <- data.frame(
    SampleID  = c("s1", "s2", "s3", "s4"),
    BatchID   = 1,
    Type      = c("QC", "QC", "sample", "QC"),
    Injection = c(1, 2, 3, 4),
    metA      = c(0, 1, 2, 3),
    keep_me   = c("a", "b", "c", "d")
  )
  out <- clean_data(df,
                    "SampleID",
                    "BatchID",
                    "Type",
                    "Injection",
                    withheld_cols = "keep_me")
  expect_false("keep_me" %in% names(out$df))
  expect_equal(out$withheld_cols, "keep_me")
})
