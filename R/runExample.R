# Run the example shiny-app: shiny-examples/metaGEO

mypackage <- "CompGenomics-metaGEO"

#' @export
runExample <- function() {
  appDir <- system.file("shiny-examples", "metaGEO", package = mypackage)
  if (appDir == "") {
    stop(sprintf("Could not find example directory. Try re-installing `%s`.", mypackage), call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
