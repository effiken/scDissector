#' Runs scDissector shiny app
#'  
#' @param preloaded_data LDM object (optional)
#' @param clustering_data_path path to clustering data folder (optional)
#' 
#' @export
run_scDissector <- function(preloaded_data=NULL,clustering_data_path=NULL) {
  appDir <- system.file("scDissector", package = "scDissector")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `scDissector`.", call. = FALSE)
  }
  
  if (!is.null(preloaded_data)){
    .scDissector_preloaded_data<<-preloaded_data
  }
  if (!is.null(clustering_data_path)){
    .scDissector_clustering_data_path<<-clustering_data_path
  }
  shiny::runApp(appDir, display.mode = "normal")
  on.exit(rm(list= list(.scDissector_preloaded_data,.scDissector_clustering_data_path)))
}
