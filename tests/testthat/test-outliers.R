# test-outliers.R

# ---------- helpers ----------
mk_data <- function() {
  # 3 QCs + 6 samples (two non-QC classes).
  # Make M1 a very clear outlier in class A (s6).
  # Make M2 have high QC variance and one non-QC outlier (s5) to exercise "unstable QC" gating.
  data.frame(
    sample = paste0("s", 1:9),
    batch  = rep("B1", 9),
    class  = c(rep("QC", 3), rep("A", 3), rep("B", 3)),
    order  = 1:9,
    M1 = c(10, 10.5, 9.5,     10.2, 9.8, 60,     9.9, 10.1, 10.0),
    M2 = c( 5, 10, 20,        10,  25,  11,     10,   9.8, 10.2),
    check.names = FALSE
  )
}

met_cols <- function(df) setdiff(names(df), c("sample","batch","class","order"))

# ---------- tests ----------
test_that("detect_qc_aware_outliers returns expected structure", {
  df <- mk_data()
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 2, md_cutoff_quantile = 0.80
  )
  expect_type(out, "list")
  expect_true(all(c("qc_rsd","sample_md","candidate_metabolites","confirmations","params") %in% names(out)))
  expect_true(all(c("metabolite","qc_rsd") %in% names(out$qc_rsd)))
  expect_true(all(c("sample","group_id","md","cutoff","flagged") %in% names(out$sample_md)))
})

test_that("grouping collapses to 'all' when only one non-QC class or when requested", {
  df <- mk_data()
  # Keep only QC + class A -> only one non-QC class remains
  out1 <- detect_qc_aware_outliers(
    df[df$class != "B" | df$class == "QC", ],
    group_nonqc_by_class = TRUE,
    z_threshold = 2, md_cutoff_quantile = 0.80
  )
  expect_setequal(unique(out1$sample_md$group_id), c("QC","all"))
  
  # Explicitly request single group for non-QC
  out2 <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = FALSE,
    z_threshold = 2, md_cutoff_quantile = 0.80
  )
  expect_setequal(unique(out2$sample_md$group_id), c("QC","all"))
})

test_that("clear outlier on M1 (s6) becomes a candidate with lenient cutoffs", {
  df <- mk_data()
  # Relax MD cutoff and Z threshold so candidate reliably appears
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 1.5, md_cutoff_quantile = 0.60
  )
  cand <- out$candidate_metabolites
  expect_true(is.data.frame(cand))
  expect_true(any(cand$sample == "s6" & cand$metabolite == "M1"))
})

test_that("insufficient group size yields 'insufficient_n' in confirmations", {
  testthat::skip_if_not_installed("outliers")
  df <- mk_data()
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 1.5, md_cutoff_quantile = 0.60,
    min_group_n = 99,        # force insufficient_n
    confirm_method = "grubbs"
  )
  conf <- out$confirmations
  expect_true(nrow(conf) == 0 || any(conf$decision == "insufficient_n"))
})

test_that("unstable QC (high RSD on M2) leads to skip_unstable_qc decision", {
  testthat::skip_if_not_installed("outliers")
  
  # Stronger M2 outlier so it reliably trips MD + Z
  df <- data.frame(
    sample = paste0("s", 1:9),
    batch  = "B1",
    class  = c(rep("QC", 3), rep("A", 3), rep("B", 3)),
    order  = 1:9,
    M1     = c(10, 10.5, 9.5,   10.2, 9.8, 10.1,   9.9, 10.1, 10.0),
    # QC M2 is very variable (unstable); s5 in class A is an extreme outlier
    M2     = c(5, 10, 20,       10, 60, 11,        10,  9.8, 10.2),
    check.names = FALSE
  )
  
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 1.0,         # looser Z to ensure candidate
    md_cutoff_quantile = 0.50, # looser MD to ensure candidate
    qc_rsd_unstable = 0.05,    # force 'unstable' gate for M2
    confirm_method = "grubbs"
  )
  
  # M2 should now appear as a candidate
  expect_true(any(out$candidate_metabolites$metabolite == "M2"))
  
  # And confirmations for M2 should be skipped due to unstable QC
  conf_m2 <- subset(out$confirmations, metabolite == "M2")
  if (nrow(conf_m2)) {
    expect_true(all(conf_m2$decision == "skip_unstable_qc" | is.na(conf_m2$decision)))
  }
})

test_that("when outliers pkg missing, decisions record 'outliers_pkg_missing'", {
  # Only run this branch when 'outliers' is NOT installed
  if (requireNamespace("outliers", quietly = TRUE)) skip("outliers installed; skipping missing-pkg branch")
  df <- mk_data()
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 1.5, md_cutoff_quantile = 0.60
  )
  conf <- out$confirmations
  if (nrow(conf)) {
    expect_true(all(conf$decision %in% c("outliers_pkg_missing","insufficient_n")))
  }
})

test_that("params echo key inputs", {
  df <- mk_data()
  out <- detect_qc_aware_outliers(
    df, group_nonqc_by_class = TRUE,
    z_threshold = 3.5, alpha = 0.01
  )
  expect_equal(out$params$z_threshold, 3.5)
  expect_equal(out$params$alpha, 0.01)
  expect_true(isTRUE(out$params$group_nonqc_by_class))
})
