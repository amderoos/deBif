.onAttach <- function(libname, pkgname) {

  msg <- "\nWelcome to the deBif package for analysis of ordinary differential equations\n"

  options(stringsAsFactors=FALSE)
  packageStartupMessage(msg)
}
