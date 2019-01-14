run <- function(tmin=0, tmax=100, tstep=1, state, parms, odes, dbopts, ...) {
  # run model and make a table
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% dbopts$args_run])
    if (length(unknown) > 0) warning(paste("Unknown argument(s):", unknown, sep=" "))
    dots_run <- dots[names(dots) %in% dbopts$args_run]
  } else dots_run <- NULL

  times <- seq(tmin, tmax, by=tstep)
  nsol <- as.data.frame(do.call('ode', c(list(times=times, func=odes, y=state, parms=parms), dots_run)))
  return(nsol)
}
