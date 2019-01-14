index <- function(strings, names, error=TRUE) {   # Return indices of strings in names
  hit <- strings %in% names
  if (error & length(strings[!hit] > 0)) stop("Unknown: ", paste(strings[!hit], collapse=", "))
  m <- match(strings[hit], names)
  if (length(m) > 0) return(m)
  return(NULL)
}
