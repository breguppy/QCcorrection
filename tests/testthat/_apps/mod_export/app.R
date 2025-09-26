options(shiny.testmode = TRUE)

library(shiny)
if (!requireNamespace("QCcorrection", quietly = TRUE)) {
  if (!requireNamespace("pkgload", quietly = TRUE)) stop("pkgload needed")
  pkgload::load_all(path = "../../..", helpers = FALSE, quiet = TRUE)
}
library(QCcorrection)

# --- STUB EXPORTERS to avoid heavy logic ---
stub_exporters <- function() {
  ns <- asNamespace("QCcorrection")
  stab <- function(sym, fn) {
    if (bindingIsLocked(sym, ns)) unlockBinding(sym, ns)
    assign(sym, fn, envir = ns); lockBinding(sym, ns)
  }
  stab("export_xlsx", function(p, d) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "data"); wb
  })
  stab("export_stats_xlsx", function(p, d) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "stats"); wb
  })
  stab("export_outliers_xlsx", function(p, d) {
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, "outliers"); wb
  })
  stab("export_figures", function(p, d, out_dir = tempdir()) {
    dir <- file.path(out_dir, "figures")
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    writeLines("ok", file.path(dir, "dummy.txt"))
    list(fig_dir = dir)
  })
  stab("render_report", function(p, d, out_dir) {
    writeLines("ok", file.path(out_dir, "correction_report.html"))
    writeLines("%PDF-1.4", file.path(out_dir, "correction_report.pdf"))
    invisible(TRUE)
  })
}
stub_exporters()
# -------------------------------------------

df <- data.frame(
  sample=c("S1","S2","S3","S4"), batch="B1", class=c("QC","Sample","Sample","QC"),
  order=1:4, A=c(10,11,12,11), B=c(20,21,NA,21), check.names=FALSE
)
rv <- list(df=df)
data_stub <- reactiveVal(list(
  cleaned=rv, filtered=rv, imputed=rv, corrected=rv,
  filtered_corrected=rv, transformed=rv
))
params_stub <- reactiveVal(list(
  sample_col="sample", batch_col="batch", class_col="class", order_col="order",
  rsd_compare="QC vs Sample", rsd_cal="by_batch", pca_compare="before_after",
  color_col="class", fig_format="png"
))

ui <- fluidPage(
  tabsetPanel(id="main_steps",
              QCcorrection:::mod_export_ui("export")
  )
)

server <- function(input, output, session) {
  QCcorrection:::mod_export_server(
    "export",
    data   = function() data_stub(),
    params = function() params_stub()
  )
}

shinyApp(ui, server)
