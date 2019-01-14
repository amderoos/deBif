#' Examples of phaseplane analysis of a system of ODEs
#'
#' \code{runExample}
#'
#'
#'   runExample(example)
#'
#'
#' Function runs one of the examples provided with the deBif package
#'
#' @param   example  (string, optional)
#' \preformatted{}
#'               Name of the example. If not provided a list of examples
#'               is returned
#'
#' @export
runExample <- function(example) {
  # locate all the shiny app examples that exist
  validExamples <- list.files(system.file("shiny-examples", package = "deBif"))

  validExamplesMsg <-
    paste0(
      "Valid examples are: '",
      paste(validExamples, collapse = "', '"),
      "'")
  if (nzchar(example) && (regexpr("\\.R", example) != (nchar(example)-1))) example <- paste0(example, ".R")

  # if an invalid example is given, throw an error
  if (missing(example) || !nzchar(example) ||
      !example %in% validExamples) {
    stop(
      'Please run `runExample()` with a valid example app as an argument.\n',
      validExamplesMsg,
      call. = FALSE)
  }

  # find and launch the app
  appDir <- system.file("shiny-examples", example, package = "deBif")
  shiny::runApp(appDir, display.mode = "normal")
}
