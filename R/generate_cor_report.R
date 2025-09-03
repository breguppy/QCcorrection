# generate correction report
source("R/processing_helpers.R")
source("R/plotting_helpers.R")

generate_cor_report <- function(p, d, out_dir, template = "report.Rmd") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_pdf <- file.path(out_dir, "correction_report.pdf")
  env <- new.env(parent = globalenv())
  
  # create plots
  # get 2 columns for met scatter plots
  if (p$rsd_cal == "met") {
    top2 <- metabolite_rsd(d$filtered$df)  %>%
    select(Metabolite, RSD_NonQC_before = RSD_NonQC) %>%
    inner_join(
      metabolite_rsd(d$filtered_corrected$df) %>%
        select(Metabolite, RSD_NonQC_after = RSD_NonQC),
      by = "Metabolite"
    ) %>%
    mutate(decrease = RSD_NonQC_before - RSD_NonQC_after) %>%
    filter(is.finite(decrease)) %>%
    arrange(desc(decrease)) %>%
    slice_head(n = 2) %>%
    pull(Metabolite)
    
    increased_qc <- metabolite_rsd(d$filtered$df) %>%
      select(Metabolite, RSD_QC_before = RSD_QC) %>%
      inner_join(
        metabolite_rsd(d$filtered_corrected$df) %>%
          select(Metabolite, RSD_QC_after = RSD_QC),
        by = "Metabolite"
      ) %>%
      filter(RSD_QC_after > RSD_QC_before) %>%
      arrange(desc(RSD_QC_after - RSD_QC_before)) %>%
      pull(Metabolite)
  } else {
    top2 <- class_metabolite_rsd(d$filtered$df) %>%
      filter(class != "QC") %>%
      select(Metabolite, RSD_before = RSD) %>%
      inner_join(
        class_metabolite_rsd(d$filtered_corrected$df) %>%
          filter(class != "QC") %>%
          select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) %>%
      mutate(decrease = RSD_before - RSD_after) %>%
      filter(is.finite(decrease)) %>%
      arrange(desc(decrease)) %>%
      distinct(Metabolite, .keep_all = TRUE) %>% 
      slice_head(n = 2) %>%
      pull(Metabolite)
    
    increased_qc <- class_metabolite_rsd(d$filtered$df) %>%
      filter(class == "QC") %>%
      select(Metabolite, RSD_before = RSD) %>%
      inner_join(
        class_metabolite_rsd(d$filtered_corrected$df) %>%
          filter(class == "QC") %>%
          select(Metabolite, RSD_after = RSD),
        by = "Metabolite"
      ) %>%
      filter(RSD_after > RSD_before) %>%
      arrange(desc(RSD_after - RSD_before)) %>%
      pull(Metabolite)
  }
  met1 <- top2[1]
  met2 <- top2[2]
  
  met1_plot <- make_met_scatter(d, met1)
  met2_plot <- make_met_scatter(d, met2)
  rsd_plot <- make_rsd_plot(p, d)
  pca_plot <- make_pca_plot(p, d)
  
  choices <- list(
    rsd_cal         = p$rsd_cal,
    color_col       = p$color_col,
    fig_format      = p$fig_format,
    correction      = d$corrected$str,
    cor_param       = d$corrected$parameters,
    transformation  = d$transformed$str,
    post_cor_filter = p$post_cor_filter,
    rsd_cutoff      = p$rsd_cutoff,
    rsd_compare     = p$rsd_compare,
    raw_df          = d$cleaned$df,
    replacement_counts = d$cleaned$replacement_counts,
    filtered        = d$filtered,
    Frule           = p$Frule,
    filtered_corrected = d$filtered_corrected
    
  )
  
  # Get descriptions for plots
  descriptions <- list(
    "Withheld Columns" = sprintf(
      "%s%s",
      tagList(
        tags$span(
          style = "font-weight:bold;", "The following columns are non-metabolite columns providing meta-information about the data:"
        ),
        tags$ul(
        lapply(c("sample = Identifies sample name", "batch = Identifies batch (large sample sets are separated into batches)", "class = Identifies sample type", "order = The order in which the sample was injected into the instrument."), function(name) {
          tags$li(name)
        })
      )
      ),
      if (isTRUE(p$withhold_cols) && !is.null(p$n_withhold)) {
        tagList(
          tags$span(
            style = "font-weight:bold;", "The following columns were withheld from correction:"
          ),
          tags$ul(
          lapply(d$cleaned$withheld_cols, function(name) {
            tags$li(name)
          })
          )
        )
      } else {""}
    ),
    "Imputation Description" = sprintf(
      "%s%s%s",
      if (d$imputed$qc_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.<br/>", d$imputed$qc_str)
      } else { "No missing QC values.<br/>"},
      if (d$imputed$sam_str != "nothing to impute") {
        sprintf("Missing QC values are imputed with %s.", d$imputed$sam_str)
      } else { "No missing sample values."},
      if (d$imputed$qc_str == "nothing to impute" && d$imputed$sam_str == "nothing to impute") {
        ""
      } else if (p$remove_imputed == TRUE) {
          "<br/>Imputed values are removed after correction."
      }
      else {""}
    ),
    "Correction Description" = sprintf(
      "Data was corrected using %s. For each metabolite, this method %s This model regresses peak areas in experimental samples, on an individual metabolite basis, against peak areas in pooled quality control samples.",
      choices$correction, choices$cor_param
    ),
    "Transformation Description" = sprintf(
      "%s <br/> %s", 
      choices$transformation,
      if (length(d$transformed$withheld_cols) > 0) {
        tagList(
          tags$span(
            style = "font-weight:bold;", "The following columns are withheld from the transformation:"
          ),
          tags$ul(
            lapply(d$transformed$withheld_cols, function(name) {
              tags$li(name)
          })
        )
        )
      } else {""}
    ),
    "Metabolite Scatter Plots" = sprintf(
      "These plots show a metabolites before and after signal drift correction before any transformation is applied. The two metabolites shown above have the largest decrease in sample variation. The change in variation was determined by calculating relative standard deviation (RSD) for each metabolite %s %s%s A full explanation of RSD is in the next section.",
      if (choices$rsd_cal == "class_met") "grouping by sample class." else "",
      if (isTRUE(!choices$post_cor_filter)) "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above " else "",
      if (isTRUE(!choices$post_cor_filter)) sprintf("%s%%.", choices$rsd_cutoff) else ""
      ),
    "RSD Comparison" = sprintf(
      "In these plots, the green indicates RSD decreased after %s, red indicates RSD increased after %s, and gray indicates no change in RSD. For these figures RSD is calculated for each metabolite %s %s%s <br/> %s ",
      if (choices$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation",
      if (choices$rsd_compare == "filtered_cor_data") "correction" else "correction and transformation",
      if (choices$rsd_cal == "class_met") "grouping by sample class." else "",
      if (isTRUE(!choices$post_cor_filter)) "Some metabolites may have been filtered out of the post-corrected dataset if the QC RSD is above " else "",
      if (isTRUE(!choices$post_cor_filter)) sprintf("%s%%.", choices$rsd_cutoff) else "",
      if (length(increased_qc) > 0) {
        sprintf(
          "<br/>The following metabolites increased QC RSD after correction: <br/> %s <br/> More investiagtion is needed to determine if these metabolites should be excluded from the data.",
          paste(increased_qc, collapse = ", ")
        )
      } else {
        ""
      }
    ),
    "PCA Comparison" = sprintf(
      "This PCA plot shows both the raw data and %s data colored by %s.",
      if (p$pca_compare == "filtered_cor_data") "corrected" else "corrected and transformed",
      choices$color_col
    )
  )
  
  params <- list(
    title   = "QC Correction Report",
    notes   = p$notes %||% "",
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