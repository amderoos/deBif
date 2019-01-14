allequi <- function(odes, state, parms, time=0, x=1, xmin=0, xmax=1, y=2, ymin=0, ymax=1.1, log="", grid=10,
                    positive=FALSE, report=TRUE, eigenvalues=TRUE, jacobian=FALSE, vector=FALSE, plot=FALSE, ...) {
  # find a steady state
  if (!is.numeric(x)) x <- index(x, names(state))
  if ((length(state) > 1) & !is.numeric(y)) y <- index(y, names(state))
  logx <- ifelse(grepl('x', log), TRUE, FALSE)
  logy <- ifelse(grepl('y', log), TRUE, FALSE)
  gridy <- ifelse((length(state) > 1), grid, 1)
  if (logx) dx <- (log10(xmax)-log10(xmin))/grid
  else dx <- (xmax-xmin)/grid
  if (logy) dy <- (log10(ymax)-log10(ymin))/gridy
  else dy <- (ymax-ymin)/gridy
  eqlst <- matrix(nrow = grid*gridy, ncol = length(state)); eqnr <- 0
  for (i in seq(1, grid)) {
    if (logx) state[x] <- 10^((i-1)*dx + dx/2 + log10(xmin))
    else state[x] <- (i-1)*dx + dx/2 + xmin
    for (j in seq(1,gridy,1)) {
      if (length(state) > 1) {
        if (logy) state[y] <- 10^((j-1)*dy + dy/2 + log10(ymin))
        else state[y] <- (j-1)*dy + dy/2 + ymin
      }
      q <- steady(y=state, func=odes, parms=parms)
      if (attr(q,"steady")) {
        equ <- q$y
        equ <- ifelse(abs(equ) < 1e-8, 0, equ)
        if (eqnr < 1) neweq <- TRUE
        else neweq <- (!any(sapply(1:eqnr, function(i) {any(c(all(abs(eqlst[i,] - equ) < 1e-4),all(abs(eqlst[i,] - equ) < 0.5e-4*(eqlst[i,] + equ))))})))
        if (neweq)
        {
          if (report) print(equ)
          jac <- jacobian.full(y=equ,func=odes,parms=parms)
          eig <- eigen(jac)
          dom <- max(sort(Re(eig$values)))
          if (eigenvalues) {
            if (dom < 0) cat("Stable point, ")
            else cat("Unstable point, ")
            cat("eigenvalues: ",eig$values,"\n")
            if (vector) {cat("Eigenvectors:\n"); print(eig$vectors)}
            if (jacobian) {cat("Jacobian:\n"); print(jac)}
          }
          if (dom < 0) points(equ[x],ifelse((length(state) > 1),equ[y],0),pch=19,cex=2)
          else points(equ[x],ifelse((length(state) > 1),equ[y],0),pch=1,cex=2)
          eqnr <- eqnr + 1
          eqlst[eqnr,] <- equ
        }
      }
    }
  }
  return(NULL)
}
