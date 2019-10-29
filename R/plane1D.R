plane1D <- function(xmin=0, xmax=1.1, ymin=0, ymax=1.1, log="", odes, state, parms, eps=0, grid=5, npixels=200, dbopts, ...) {
  # Make a phase plane for a single ODE with nullclines and/or phase portrait
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(dbopts$args_run, dbopts$args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):", unknown, sep=" "))
  }
  if (!is.null(dots)) dots_run <- dots[names(dots) %in% dbopts$args_run]
  else dots_run <- NULL
  x <- 1
  logx <- ifelse(grepl('x', log), TRUE, FALSE)
  logy <- ifelse(grepl('y', log), TRUE, FALSE)
  if (logx) xc <- 10^seq(log10(xmin), log10(xmax), length.out=npixels)
  else xc <- seq(xmin+eps, xmax, length.out=npixels)
  xvar <- names(state)[x];
  yvar <- paste0("d", xvar, "/dt")
  par(cex = dbopts$plotopts["cex"], mar = dbopts$plotmar)
  do.call('plot', c(list(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xvar, ylab=yvar, log=log,
                         cex.lab=dbopts$plotopts["cex.lab"], cex.axis=dbopts$plotopts["cex.axis"]),
                    dots[names(dots) %in% dbopts$args_plot]))
  lines(c(xmin,xmax), c(0,0), col="black", lwd=1, lty=2)
  legend("topright",legend=names(state)[1], col=dbopts$colors[1], lty=1, lwd=dbopts$plotopts["lwd"], cex=dbopts$plotopts["cex.legend"])

  xc <- seq(xmin, xmax, length.out = npixels)
  dxdt <- as.numeric(lapply(xc, function(i) {state[1] <- i; odes(0, state, parms)[[1]][[1]]}))

  lines(xc, dxdt, lwd=dbopts$plotopts["lwd"], col=dbopts$colors[1])
}

