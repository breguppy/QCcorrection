# app_server.R

app_server <- function(input, output, session) {
  session$onSessionEnded(function() { stopApp() })
  
  import <- mod_import_server("import")
  
  #correct <- mod_correct_server("correct",
  #                              filtered = import$filtered,
  #                              params   = import$params)
  
  #viz <- mod_visualize_server("viz",
  #                            filtered           = import$filtered,
  #                            filtered_corrected = correct$filtered_corrected,
  #                            transformed        = correct$transformed,
  #                            params             = import$params)
  
  #mod_export_server("export",
   #                 filtered_corrected = correct$filtered_corrected,
  #                  transformed        = correct$transformed,
  #                  params             = import$params)
}
