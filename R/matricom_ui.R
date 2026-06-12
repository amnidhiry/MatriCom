#' matricom.app: launch matricom as an interactive, local app
#'
#' @return a local instance of the matricom shinyapp, without connection to the (online-only) OA data
#' @export
matricom.app <- function(){
  # where is webApp? Find and run from there
  app_dir <- system.file("webApp", package = "MatriCom")

  if (app_dir == "") {
    stop(
      "Could not find the app directory. Try re-installing `MatriCom`.",
      call. = FALSE
    )
  }

  shiny::runApp(appDir = app_dir, launch.browser = TRUE)
}

