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
  
  choices <- list(
    rsd_cal         = input$rsd_cal,
    color_col       = input$color_col,
    fig_format      = input$fig_format,
    correction      = rv$corrected$str,
    cor_param       = rv$corrected$parameters,
    transformation  = rv$transformed$str,
    post_cor_filter = input$post_cor_filter,
    rsd_filter      = input$rsd_filter
  )
  
  # Get descriptions for plots
  descriptions <- list(
    "Correction Description" = sprintf(
      "Instrument drift in this data is corrected using %s. For each metabolite, this method %s This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
      choices$correction, choices$cor_param
    ),
    "Transformation Description" = sprintf(
      "%s", choices$transformation
    ),
    "Metabolite Scatter Plots" = sprintf(
      "These plots show a metabolites before and after signal drift correction."),
    "RSD Comparison" = sprintf(
      "To judge how well the correction method worked, we can visualize the change in variation with the relative standard deviation (RSD) comparison plots. \r\n RSD = (standard deviation / mean) * 100%%.\r\n For these figures RSD is calculated for each metabolite %s %s%s",
      if (choices$rsd_cal == "class_met") "grouping by sample class." else "",
      if (isTRUE(!choices$post_cor_filter)) "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above " else "",
      if (isTRUE(!choices$post_cor_filter)) sprintf("%s%%.", choices$rsd_filter) else ""
    ),
    "PCA Comparison" = sprintf(
      "PCA colored by %s. Correction method: %s.",
      choices$color_col, choices$correction
    )
  )
  
  params <- list(
    title   = "QC correction report",
    notes   = input$notes %||% "",
    plots = list(
      "Metabolite Scatter 1" = met1_plot,
      "Metabolite Scatter 2" = met2_plot,
      "RSD Comparison"       = rsd_plot,
      "PCA Comparison"       = pca_plot
    ),
    include = NULL,
    choices = choices,
    descriptions = descriptions
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