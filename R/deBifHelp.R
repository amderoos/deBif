#' Opens the deBif manual
#'
#' \code{deBifHelp} opens the manual of the the deBif package in html format.
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' deBifHelp()
#' }
#'
#' @importFrom rstudioapi viewer
#' @export
deBifHelp <- function ()
{
  tempDir <- tempdir()
  file.copy(system.file("manual", package = "deBif"), tempDir, recursive=TRUE)
  htmlFile <- file.path(tempDir, "manual/index.html")
  rstudioapi::viewer(htmlFile, height="maximize")
}

