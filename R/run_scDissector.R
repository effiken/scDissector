#' @export
runApp <- function() {
  appDir <- system.file("scDissector", package = "scDissector")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `scDissector`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}
