#' @name runSmartPhos
#' 
#' @title Launch the SmartPhos Shiny Application
#'
#' @description
#' \code{runSmartPhos} launches the SmartPhos \code{Shiny} application, which provides an interactive interface for analyzing phosphoproteomic data.
#'
#' @return The function does not return a value; it starts the Shiny application for \code{SmartPhos}.
#' @details
#' The \code{runSmartPhos} function locates the \code{Shiny} app directory within the \code{SmartPhos} package and launches the application. 
#' If the app directory cannot be found, the function will stop and prompt the user to re-install the \code{SmartPhos} package.
#' @export
#' @examples
#' # To run the SmartPhos Shiny application, simply call:
#' runSmartPhos()

runSmartPhos <- function() {
  # Locate the Shiny app directory within the SmartPhos package
  appDir <- system.file("shiny-app", package = "SmartPhos")
  if (appDir == "") {
    # If the directory is not found, stop the function and display an error message
    stop("Could not find example directory. Try re-installing `SmartPhos`.", call. = FALSE)
  }
  # Launch the Shiny application in normal display mode
  shiny::runApp(appDir, display.mode = "normal")
}