bialt2AI <- function(A) {
  # Constructs the bialternate matrix product of 2A and I.
  #
  # See pages 485-487 in Kuznetsov, 1995; Elements of Applied Bifurcation
  # Analysis In particular, figure 10.8 and the equation at the top of page 487.
  n <- ncol(A)
  m <-n*(n-1)/2

  # p = 2, 3, 3, 4, 4, 4, ....
  p <- unlist(lapply((2:n), function(i) rep(i, (i-1))))
  # q = 1, 1, 2, 1, 2, 3, ....
  q <- unlist(lapply((2:n), function(i) (1:(i-1))))
  r <- p
  s <- q
  P <- matrix(p, m, m)
  Q <- matrix(q, m, m)
  R <- t(P)
  S <- t(Q)
  Aps <- matrix(c(A[p,]), length(p), n)[,s]
  Apr <- matrix(c(A[p,]), length(p), n)[,r]
  App <- matrix(c(A[p,]), length(p), n)[,p]
  Aqq <- matrix(c(A[q,]), length(q), n)[,q]
  Aqs <- matrix(c(A[q,]), length(q), n)[,s]
  Aqr <- matrix(c(A[q,]), length(q), n)[,r]

  twoAI <- (R == Q)*(-Aps) + ((R != P) & (S == Q))*Apr + ((R == P) & (S == Q))*(App + Aqq) + ((R == P) & (S != Q))*Aqs + (S == P)*(-Aqr)

  return(twoAI)
}

approxNullVec <- function(A) {
  # Find the eigenvector of the matrix A pertaining to the eigenvalue with
  # smallest absolute value
  eig <- eigen(A)
  minindx <- which.min(abs(Re(eig$values)))
  eigvec <- c(eig$vectors[,minindx])
  return (as.numeric(eigvec))
}

converty2y <- function(y, ymin1, ymax1, logy1, ymin2, ymax2, logy2) {
  # Convert values on the second y axis to corresponding values on the first y-axis
  if (logy2 == 0) {
    ytmp <- pmax((y - ymin2)/(ymax2 - ymin2), -0.5)
  } else {
    ytmp <- pmax(log(pmax(y, 1.0E-15)/ymin2)/(log(ymax2/ymin2)), -0.5)
  }

  if (logy1 == 0) {
    yval <- ymin1 + ytmp*(ymax1 - ymin1)
  } else {
    yval <- exp(log(ymin1) + ytmp*(log(ymax1/ymin1)))
  }
  return(yval)
}

