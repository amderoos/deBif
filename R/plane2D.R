plane2D <- function(xmin=0, xmax=1.1, ymin=0, ymax=1.1, log="", odes, state, parms,
                    x=1, y=2, time=0, grid=5, eps=0, npixels=500, portrait=FALSE, vector=FALSE, dbopts, ...) {
  # Make a phase plane with nullclines and/or phase portrait
  dots <- list(...)
  if (!is.null(dots)) {
    unknown <- names(dots[!names(dots) %in% c(dbopts$args_run, dbopts$args_plot)])
    if (length(unknown)>0) warning(paste("Unknown argument(s):", unknown,sep=" "))
  }
  if (!is.null(dots)) dots_run <- dots[names(dots) %in% dbopts$args_run]
  else dots_run <- NULL
  if (!is.numeric(x)) x <- index(x, names(state))
  if (!is.numeric(y)) y <- index(y, names(state))
  ishows <- c(x, y)
  nvar <- length(state)
  lvec <- 50                         # length of vector
  logx <- ifelse(grepl('x', log), TRUE, FALSE)
  logy <- ifelse(grepl('y', log), TRUE, FALSE)
  if (logx) xc <- 10^seq(log10(xmin), log10(xmax), length.out=npixels)
  else xc <- seq(xmin+eps, xmax, length.out=npixels)
  if (logy) yc <- 10^seq(log10(ymin), log10(ymax), length.out=npixels)
  else yc <- seq(ymin+eps, ymax, length.out=npixels)
  xvar <- names(state)[x]; yvar <- names(state)[y]
  par(cex = as.numeric(dbopts$plotopts["cex"]), mar = dbopts$plotmar)
  do.call('plot', c(list(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), xlab=xvar, ylab=yvar, log=log,
                         cex.lab=as.numeric(dbopts$plotopts["cex.lab"]), cex.axis=as.numeric(dbopts$plotopts["cex.axis"])),
                    dots[names(dots) %in% dbopts$args_plot]))
  legend("topright", legend=names(state)[ishows], col=dbopts$colors[ishows], lty=1, lwd=dbopts$plotopts["lwd"], cex=as.numeric(dbopts$plotopts["cex.legend"]))

  vstate <- as.list(state)
  npixels2 <- npixels^2
  #vstate<-lapply(vstate,rep,vstate,npixels2);vstate[[x]]<-0;vstate[[y]]<-0
  for (j in seq(1, nvar)) if (j!=x & j!=y) vstate[[j]]<-rep(vstate[[j]],npixels2);
  FUN <- function(xc, yc, i){  # wrapper around model()
    vstate[[x]] <- xc; vstate[[y]] <- yc
    odes(time, vstate, parms)[[1]][seq((i-1)*npixels2+1, i*npixels2)]
  }
  for (i in ishows)
    contour(xc, yc, outer(xc, yc, FUN, i), levels=0, drawlabels=FALSE, add=TRUE, col=dbopts$colors[i], lwd=dbopts$plotopts["lwd"])

  if (portrait | vector) {
    if (logx) {dx <- (log10(xmax)-log10(xmin))/grid; vx <- 1+3.32*grid*dx/lvec}
    else {dx <- (xmax-xmin)/grid; vx = grid*dx/lvec}
    if (logy) {dy <- (log10(ymax)-log10(ymin))/grid; vy <- 1+3.32*grid*dy/lvec}
    else {dy <- (ymax-ymin)/grid; vy = grid*dy/lvec}

    for (i in seq(1,grid)) {
      if (logx) state[x] <- 10^((i-1)*dx + dx/2 + log10(xmin))
      else state[x] <- (i-1)*dx + dx/2 + xmin
      for (j in seq(1,grid,1)) {
        if (logy) state[y] <- 10^((j-1)*dy + dy/2 + log10(ymin))
        else state[y] <- (j-1)*dy + dy/2 + ymin
        if (portrait) {
          points(state[x], state[y], pch=as.numeric(dbopts$plotopts["pch"]))
          nsol <- do.call('run', c(list(state=state, parms=parms, odes=odes, dbopts=dbopts), dots_run))
          lines(cbind(nsol[x+1], nsol[y+1]), col="black")
        }
        if (vector) {
          dt <- sign(unlist(odes(time, state, parms)))
          if (logx) lines(c(state[x], state[x]*vx^dt[x]), c(state[y], state[y]))
          else lines(c(state[x], state[x]+vx*dt[x]), c(state[y], state[y]))
          if (logy) lines(c(state[x], state[x]), c(state[y], state[y]*vy^dt[y]))
          else lines(c(state[x], state[x]), c(state[y], state[y]+vy*dt[y]))
        }
      }
    }
  }
}
