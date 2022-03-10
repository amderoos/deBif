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
#' @importFrom utils unzip
#' @export
deBifHelp <- function ()
{
  oldwd <- getwd()
  tempDir <- tempdir()
  unlink(paste0(tempDir, "/manual"), recursive = TRUE)
  dir.create(paste0(tempDir, "/manual"))
  setwd(paste0(tempDir, "/manual"))
  unzip(paste0(system.file("manual", package = "deBif"), "/deBif-manual.zip"))
  file.copy(paste0(system.file("doc", package = "deBif"), "/deBif-manual.pdf"), paste0(tempDir, "/manual"))
  # file.copy(system.file("manual", package = "deBif"), tempDir, recursive=TRUE)
  setwd(oldwd)
  htmlFile <- file.path(tempDir, "manual/index.html")
  rstudioapi::viewer(htmlFile, height="maximize")
}

