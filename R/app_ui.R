#' app_ui.R
#' @importFrom bslib card navset_tab nav_panel layout_sidebar sidebar bs_theme
#' @importFrom shiny tags icon fluidPage
#' @keywords internal
#' @noRd

app_ui <- function() {
  fluidPage(
    theme = bslib::bs_theme(preset = "cosmo"),
    shinyjs::useShinyjs(),
    titlePanel("QC Correction for Metabolomics Data"),
    bslib::navset_tab(
      id = "main_steps",
      mod_import_ui("import"),
      mod_correct_ui("correct"),
      mod_visualize_ui("viz"),
      mod_export_ui("export")
    )
  )
}
