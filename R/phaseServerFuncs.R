allequi <- function(curtab, odes, state, parms, plotopts, numopts) {
  xcol <- as.numeric(plotopts$xcol)
  ycol <- as.numeric(plotopts$ycol)
  xmin <- as.numeric(plotopts$xmin)
  xmax <- as.numeric(plotopts$xmax)
  ymin <- as.numeric(plotopts$ymin)
  ymax <- as.numeric(plotopts$ymax)
  logx <- (as.numeric(plotopts$logx) == 1)
  logy <- (as.numeric(plotopts$logy) == 1)
  grid <- numopts$ssgrid

  # find a steady state
  gridy <- ifelse((length(state) > 1), grid, 1)
  msg <- ""
  if (logx) dx <- (log10(xmax)-log10(xmin))/grid
  else dx <- (xmax-xmin)/grid
  if (logy) dy <- (log10(ymax)-log10(ymin))/gridy
  else dy <- (ymax-ymin)/gridy
  eqlst <- matrix(nrow = grid*gridy, ncol = length(state)); eqnr <- 0
  for (i in seq(1, grid)) {
    if (logx) state[xcol] <- 10^((i-1)*dx + dx/2 + log10(xmin))
    else state[xcol] <- (i-1)*dx + dx/2 + xmin
    for (j in seq(1,gridy,1)) {
      if (length(state) > 1) {
        if (logy) state[ycol] <- 10^((j-1)*dy + dy/2 + log10(ymin))
        else state[ycol] <- (j-1)*dy + dy/2 + ymin
      }
      q <- tryCatch(steady(y=state, func=odes, parms=parms))
      if (attr(q,"steady") && (!any(is.nan(q$y))) &&
          (q$y[xcol] >= xmin - 1e-8) && (q$y[xcol] <= xmax + 1e-8) &&
          ((length(state) == 1) || ((q$y[ycol] >= ymin - 1e-8) && (q$y[ycol] <= ymax + 1e-8)))) {
            equ <- q$y
        equ <- ifelse(abs(equ) < 1e-8, 0, equ)
        if (eqnr < 1) neweq <- TRUE
        else neweq <- (!any(sapply(1:eqnr, function(i) {any(c(all(abs(eqlst[i,] - equ) < 1e-4),all(abs(eqlst[i,] - equ) < 0.5e-4*(eqlst[i,] + equ))))})))
        if (neweq)
        {
          if (curtab >= 3) msg <- paste0(msg, paste(unlist(lapply(1:length(equ), function(i) {paste(names(equ[i]), "=", round(equ[i],6), sep=" ")})), collapse=', '), "\n")
          jac <- jacobian.full(y=equ,func=odes,parms=parms)
          eig <- eigen(jac)
          dom <- max(sort(Re(eig$values)))
          if (curtab >= 3) {
            if (dom < 0) msg <- paste0(msg, "Stable point\n")
            else msg <- paste0(msg, "Unstable point\n")
            msg <- paste0(msg, "Eigenvalues: ", paste(round(eig$values,5), collapse=', '), "\n")
          }
          if (curtab >= 3) msg <- paste0(msg, "\n")
          if (dom < 0) points(equ[xcol],ifelse((length(state) > 1),equ[ycol],0),pch=19,cex=2)
          else points(equ[xcol],ifelse((length(state) > 1),equ[ycol],0),pch=1,cex=2)
          eqnr <- eqnr + 1
          eqlst[eqnr,] <- equ
        }
      }
    }
  }
  return(msg)
}

nullclines <- function(odes, state, parms, plotopts, numopts) {
  xcol <- as.numeric(plotopts$xcol)
  ycol <- as.numeric(plotopts$ycol)
  xmin <- as.numeric(plotopts$xmin)
  xmax <- as.numeric(plotopts$xmax)
  ymin <- as.numeric(plotopts$ymin)
  ymax <- as.numeric(plotopts$ymax)
  logx <- (as.numeric(plotopts$logx) == 1)
  logy <- (as.numeric(plotopts$logy) == 1)
  eps  <- numopts$eps
  npixels <- numopts$npixels

  # Make a phase plane with nullclines and/or phase portrait
  nvar <- length(state)

  if (logx) {
    xc <- 10^seq(log10(xmin), log10(xmax), length.out=npixels)
  } else xc <- seq(xmin+eps, xmax, length.out=npixels)
  if (logy) {
    yc <- 10^seq(log10(ymin), log10(ymax), length.out=npixels)
  } else yc <- seq(ymin+eps, ymax, length.out=npixels)

  vstate <- as.list(state)
  npixels2 <- npixels^2

  xyc <- cbind(as.vector(matrix(xc, nrow = npixels, ncol = npixels, byrow = T)),
               as.vector(matrix(yc, nrow = npixels, ncol = npixels)))

  ztmp <- apply(xyc, 1, function(cst) {vstate[[xcol]] <- cst[1]; vstate[[ycol]] <- cst[2]; unlist(odes(0, vstate, parms))}, simplify = T)

  zlst <- list(xc = xc, yc = yc,
               dxc = matrix(ztmp[1,], npixels, npixels, byrow = T),
               dyc = matrix(ztmp[2,], npixels, npixels, byrow = T))

  return(zlst)
}
