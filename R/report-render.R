#' Writes HTML of correction report then writes PDF of correction report
#'
#' @keywords internal
#' @noRd
render_report <- function(p, d, out_dir,
                          template = system.file("app","report_templates","report.Rmd", package = "QCcorrection")) {
  .require_pkg("rmarkdown", "render reports")
  .require_pkg("pagedown",  "print to PDF")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  env <- new.env(parent = baseenv())
  
  # build plots
  met_candidates <- setdiff(names(d$filtered$df), c("sample","batch","class","order"))
  met1_plot <- make_met_scatter(d, met_candidates[1])
  met2_plot <- make_met_scatter(d, met_candidates[2])
  rsd_plot  <- make_rsd_plot(p, d)
  pca_plot  <- make_pca_plot(p, d)
  
  # plain strings only; no shiny::tagList here
  descriptions <- list(
    Withheld = if (isTRUE(p$withhold_cols) && length(d$cleaned$withheld_cols)) {
      paste0("Withheld columns from correction: ", paste(d$cleaned$withheld_cols, collapse = ", "))
    } else "No columns withheld from correction.",
    Imputation = paste0("QC: ", d$imputed$qc_str, "; Samples: ", d$imputed$sam_str,
                        if (isTRUE(p$remove_imputed)) ". Imputed values removed after correction." else ""),
    Correction = sprintf("Correction method: %s. %s", d$corrected$str, d$corrected$parameters),
    Transform  = d$transformed$str
  )
  
  params <- list(
    title = "QC Correction Report",
    notes = p$notes %||% "",
    plots = list(`Metabolite 1` = met1_plot, `Metabolite 2` = met2_plot,
                 `RSD Comparison` = rsd_plot, `PCA Comparison` = pca_plot),
    choices = list(
      raw_df = d$cleaned$df,
      replacement_counts = d$cleaned$replacement_counts,
      filtered = d$filtered,
      filtered_corrected = d$filtered_corrected
    ),
    descriptions = descriptions
  )
  
  html_out <- rmarkdown::render(
    input = template,
    output_format = "html_document",
    output_file   = file.path(out_dir, "correction_report.html"),
    params = params,
    envir  = env,
    quiet  = TRUE
  )
  
  chrome <- pagedown::find_chrome()
  if (is.null(chrome)) {
    warning("Chrome/Chromium not found. Returning HTML only.")
    return(normalizePath(html_out, winslash = "/"))
  }
  pdf_out <- file.path(out_dir, "correction_report.pdf")
  pagedown::chrome_print(input = html_out, output = pdf_out)
  normalizePath(pdf_out, winslash = "/")
}
