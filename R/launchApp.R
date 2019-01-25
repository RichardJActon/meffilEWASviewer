#' launches the shinyApp
#'
#' @export ewas_viewer
#'
#' @return shiny application object
#'
#' @example \dontrun {ewas_viewer()}
#'
#' @import shiny
#'

# wrapper for shiny::shinyApp()
ewas_viewer <- function() {
	# shinyAppDir("./"
	#   #system.file("shinyApp", package="meffilEWASviewer")
	# )
	shinyApp(ui = ui, server = server)
}

# NB shinyAppDir bug means globals.r is not sourced to having to use single app.r file for now.
