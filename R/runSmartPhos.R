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
#' # runSmartPhos()

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


#' @name makeSmartPhosDirectory
#'
#' @title Create SmartPhos Directory Structure
#'
#' @description
#' \code{makeSmartPhosDirectory} creates a directory for the SmartPhos shiny app,
#' and copies the necessary Shiny app files into the newly created directory.
#'
#' @param path A \code{character} string specifying the directory path where the SmartPhos folder should be created.
#'
#' @return None (invisible NULL). The function creates the necessary directories and copies files.
#'
#' @details
#' The function first creates the main directory at the specified path and a subdirectory named `"save"` for storing \code{MultiAssayExperiment} object.
#' It then locates the Shiny application files from the SmartPhos package and copies them into the new directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' makeSmartPhosDirectory("path/to/destination")
#' }
makeSmartPhosDirectory <- function(path) {
    # Create the directory folder and subdirectory save
    dir.create(path)
    dir.create(file.path(path, "save"))
    # Get location of current shiny app
    shinyPath <- system.file(package = "SmartPhos", "shiny-app")
    # Copy the contents to the newly created folder
    list_of_files <- list.files(shinyPath)
    file.copy(file.path(shinyPath, list_of_files), path, recursive = TRUE)
    message("Successful!!")
}
