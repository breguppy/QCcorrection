#app_ui.R

app_ui <- function() {
  fluidPage(
    theme = bs_theme(preset = "cosmo"),
    useShinyjs(),
    titlePanel("QC Correction for Metabolomics Data"),
    navset_tab(
      id = "main_steps",
      mod_import_ui("import"),
      #mod_correct_ui("correct"),
      #mod_visualize_ui("viz"),
      #mod_export_ui("export")
    )
  )
}
