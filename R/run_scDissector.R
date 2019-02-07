#' @export
run_scDissector <- function(preloaded_data=NULL,clustering_data_path=NULL) {
  appDir <- system.file("scDissector", package = "scDissector")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `scDissector`.", call. = FALSE)
  }
  
  .scDissector_preloaded_data<<-preloaded_data
  .scDissector_clustering_data_path<<-clustering_data_path
  shiny::runApp(appDir, display.mode = "normal")
}
