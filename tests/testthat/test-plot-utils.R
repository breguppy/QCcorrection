make_rsd_key_test_data <- function() {
  data.frame(
    Type = factor(
      c(rep("Samples", 4), rep("QC", 4)),
      levels = c("Samples", "QC")
    ),
    before = c(10, 20, 30, 40, 10, 20, 30, 40),
    after = c(12, 18, 30, 35, 8, 18, 25, 35),
    change = factor(
      c(
        "Increased", "Decreased", "No Change", "Decreased",
        "Decreased", "Decreased", "Decreased", "Decreased"
      ),
      levels = lab_levels
    )
  )
}

test_that("RSD key rows preserve per-facet percentages", {
  key_rows <- .rsd_key_rows(make_rsd_key_test_data(), "before", "after")

  testthat::expect_equal(nrow(key_rows), 6L)
  testthat::expect_equal(
    key_rows |>
      dplyr::filter(Type == "Samples") |>
      dplyr::pull(label),
    c("Increased 25%", "No change 25%", "Decreased 50%")
  )
  testthat::expect_equal(
    key_rows |>
      dplyr::filter(Type == "QC") |>
      dplyr::pull(label),
    c("Increased 0%", "No change 0%", "Decreased 100%")
  )
  testthat::expect_equal(as.character(key_rows$change), rep(lab_levels, 2L))
})

test_that("RSD key columns are separated and flush to the panel top", {
  d_all <- make_rsd_key_test_data()
  key_rows <- .rsd_key_rows(d_all, "before", "after")
  x_range <- range(d_all[["before"]])
  x_span <- diff(x_range)

  testthat::expect_equal(
    (key_rows$key_x[seq_len(3)] - x_range[[1L]]) / x_span,
    c(0.005, 0.41, 0.76)
  )
  testthat::expect_equal(
    (key_rows$label_x[seq_len(3)] - x_range[[1L]]) / x_span,
    c(0.025, 0.435, 0.785)
  )
  testthat::expect_equal(key_rows$label_hjust[seq_len(3)], c(0, 0, 0))
  testthat::expect_true(all(key_rows$band_y_max == key_rows$limit_y))

  plot <- mk_plot(d_all, "before", "after", facet_label_map(d_all), "Correction")
  y_scale <- plot$scales$get_scales("y")

  testthat::expect_true(inherits(y_scale, "ScaleContinuousPosition"))
  testthat::expect_equal(y_scale$expand, ggplot2::expansion(mult = c(0.05, 0)))
})
test_that("RSD comparison plot uses native point key markers", {
  d_all <- make_rsd_key_test_data()
  plot <- mk_plot(d_all, "before", "after", facet_label_map(d_all), "Correction")

  key_point_layers <- purrr::keep(
    plot$layers,
    function(layer) {
      inherits(layer$geom, "GeomPoint") &&
        isFALSE(layer$inherit.aes) &&
        identical(layer$aes_params$shape, 19) &&
        identical(layer$aes_params$size, 2)
    }
  )

  testthat::expect_length(key_point_layers, 1L)
  testthat::expect_equal(nrow(key_point_layers[[1]]$data), 6L)
  testthat::expect_equal(
    unname(plot$scales$get_scales("colour")$palette(3)),
    unname(color_values[lab_levels])
  )
})

test_that("RSD comparison key renders without Unicode marker warnings", {
  d_all <- make_rsd_key_test_data()
  plot <- mk_plot(d_all, "before", "after", facet_label_map(d_all), "Correction")
  warnings <- character()

  withCallingHandlers(
    ggplot2::ggplotGrob(plot),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  testthat::expect_false(any(grepl("U\\+|unicode|translate|conversion", warnings, ignore.case = TRUE)))
})
