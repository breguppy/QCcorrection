#' Plot entry points used by the app
#' @keywords internal
#' @noRd
# Helper for creating metabolite scatter plot
make_met_scatter <- function(d, met_col) {
  # choose the correct plotting function based on the correction method.
  cor_method <- d$corrected$str
  tryCatch({
    if (cor_method %in% c("Random Forest", "Batchwise Random Forest")) {
      met_scatter_rf(d$filtered$df, d$filtered_corrected$df, i = met_col)
    } else if (cor_method %in% c("LOESS", "Batchwise LOESS")) {
      met_scatter_loess(d$filtered$df, d$filtered_corrected$df, i = met_col)
    } else {
      ggplot2::ggplot() + ggplot2::labs(title = "No correction method selected.")
    }
  }, error = function(e) {
    shiny::showNotification(paste("Scatter failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "Scatter failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the RSD plot
make_rsd_plot <- function(p, d) {
  df_before <- d$filtered$df
  # Determine df_after based on rsd_compare selected by user.
  if (p$rsd_compare == "filtered_cor_data") {
    df_after <- d$filtered_corrected$df
    compared_to <- "Correction"
  } else {
    df_after <- d$transformed$df
    compared_to <- "Correction and Transformation"
  }
  
  # Need at least 1 metabolite column
  shiny::validate(
    shiny::need(ncol(df_before) > 4L, "No metabolites left before correction."),
    shiny::need(ncol(df_after)  > 4L, "No metabolites left after correction.")
  )
  
  tryCatch({
    if (identical(p$rsd_plot_type, "scatter")) {
      if (identical(p$rsd_cal, "met")) {
        plot_rsd_comparison(df_before, df_after, compared_to)
      } else {
        plot_rsd_comparison_class_met(df_before, df_after, compared_to)
      }
    } else {
      if (identical(p$rsd_cal, "met")) {
        plot_met_rsd_distributions(df_before, df_after, compared_to)
      } else {
        #distribution plot for class-met rsd comparison
      }
    }
  }, error = function(e) {
    shiny::showNotification(
      paste("RSD comparison failed:", e$message),
      type = "error",
      duration = 8
    )
    ggplot2::ggplot() + ggplot2::labs(title = "RSD comparison failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the PCA plot
make_pca_plot <- function(p, d) {
  # get after based on pca_compare selected by user.
  if (p$pca_compare == "filtered_cor_data") {
    df <- d$filtered_corrected$df
    after <- d$filtered_corrected
    compared_to <- "Correction"
  } else {
    df <- d$transformed$df
    after <- d$transformed
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(df), c("sample", "batch", "class", "order"))
  shiny::validate(
    shiny::need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    shiny::need(nrow(df) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- df[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  shiny::validate(need(
    any(keep),
    "All metabolite columns are constant/invalid after filtering."
  ))
  
  # before data cannot have any missing values.
  before <- d$imputed
  tryCatch({
    plot_pca(p, before, after, compared_to)
  }, error = function(e) {
    shiny::showNotification(paste("PCA failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "PCA failed \u2013 see notification")
  })
}

#' @keywords internal
#' @noRd
# Helper for creating the PCA plot
make_pca_loading_plot <- function(p, d) {
  # get after based on pca_compare selected by user.
  if (p$pca_compare == "filtered_cor_data") {
    df <- d$filtered_corrected$df
    after <- d$filtered_corrected
    compared_to <- "Correction"
  } else {
    df <- d$transformed$df
    after <- d$transformed
    compared_to <- "Correction and Transformation"
  }
  mets <- setdiff(names(df), c("sample", "batch", "class", "order"))
  shiny::validate(
    shiny::need(length(mets) >= 2, "Need at least 2 metabolite columns for PCA."),
    shiny::need(nrow(df) >= 3, "Need at least 3 samples for PCA.")
  )
  
  # Also ensure non-constant / non-NA columns
  X <- df[, mets, drop = FALSE]
  keep <- vapply(X, function(v) {
    v <- suppressWarnings(as.numeric(v))
    ok <- all(is.finite(v))
    nz <- (length(unique(v)) >= 2)
    ok && nz
  }, logical(1))
  shiny::validate(need(
    any(keep),
    "All metabolite columns are constant/invalid after filtering."
  ))
  
  # before data cannot have any missing values.
  before <- d$imputed
  tryCatch({
    plot_pca_loading(p, before, after, compared_to)
  }, error = function(e) {
    shiny::showNotification(paste("PCA failed:", e$message),
                            type = "error",
                            duration = 8)
    ggplot2::ggplot() + ggplot2::labs(title = "PCA Loading failed \u2013 see notification")
  })
}