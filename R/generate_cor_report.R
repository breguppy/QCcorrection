# generate correction report

generate_cor_report <- function(input, rv, out_dir, template = "report.Rmd") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_pdf <- file.path(out_dir, "correction_report.pdf")
  env <- new.env(parent = globalenv())
  
  # create plots
  # get 2 columns for met scatter plots
  raw_cols <- setdiff(names(rv$filtered$df), c("sample","batch","class","order"))
  cor_cols <- setdiff(names(rv$filtered_corrected$df), c("sample","batch","class","order"))
  cols <- intersect(raw_cols, cor_cols)
  validate(need(length(cols) >= 1, "No overlapping metabolites between raw and corrected data."))
  met1 <- cols[1]
  met2 <- cols[ceiling(length(cols)/2)]
  
  met1_plot <- make_met_scatter(rv, met1)
  met2_plot <- make_met_scatter(rv, met2)
  rsd_plot <- make_rsd_plot(input, rv)
  pca_plot <- make_pca_plot(input, rv)
  
  params <- list(
    title   = "QC correction report",
    notes   = input$notes %||% "",
    plots = list(
      "Metabolite Scatter 1" = met1_plot,
      "Metabolite Scatter 2" = met2_plot,
      "RSD Comparison"       = rsd_plot,
      "PCA Comparison"       = pca_plot
    ),
    include = NULL
  )
  
  # Always render HTML first
  html_out <- rmarkdown::render(
    input         = template,
    output_format = "html_document",
    output_file   = file.path(out_dir, "correction_report.html"),
    params        = params,
    envir         = env
  )
  
  # Then print HTML to PDF with Chrome
  if (is.null(pagedown::find_chrome())) {
    stop("Chrome/Chromium not found. Please install Google Chrome or Microsoft Edge.")
  }
  pagedown::chrome_print(input = html_out, output = out_pdf)
  out_pdf
}