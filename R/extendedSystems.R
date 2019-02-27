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

TangentVecEQ <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts) {

  # Routine calculates the tangent vector to the curve determined by the equation
  #
  #                           F(x, p1) = 0
  #
  # This is a system of n equations, in which n is the length of the state variable
  # vector (given by statedim)

  #  When freeparsdim == 2 the Jacobian equals the following square (n+2)x(n+2)
  #  matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  # Otherwise, when freeparsdim == 1 the Jacobian equals the following square
  # (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # If freeparsdim == 2 delete the last row and 2nd column of the Jacobian, the
  # latter being the derivative w.r.t. second free parameters
  jac <- jac[(1:(statedim+1)), c(1, ((freeparsdim+1):(freeparsdim+statedim)))]

  tvdone <- FALSE
  # Append the current tangent vector as the last row to the jacobian to
  # preserve direction. See the matcont manual at
  # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
  # But only do this if the tangent vector passed as argument has the same
  # dimension as the equilibrium curve
  if (!is.null(tanvec) && (length(tanvec) == (statedim + 1))) {
    jac[nrow(jac),] <- tanvec
    if (rcond(jac) > nopts$atol) {
      # Solve for the tangent vector to the curve of the 1st free parameter and
      # the state variables.
      tvnew <- solve(jac, c(rep(0, (length(tanvec)-1)), 1))
      tvdone <- TRUE
    }
  }

  if (!tvdone) {
    # Set the last row with NA to 0
    jac[nrow(jac),] <- 0

    # Estimate the condition numbers
    condnrs <- unlist(lapply((1:ncol(jac)),
                             function(i){
                               jacc <- jac;
                               jacc[nrow(jac),i] <- 1;
                               if (abs(det(jacc)) > 1.0E-5) return(rcond(jacc))
                               else return(0.0)
                             }))
    maxind <- which.max(condnrs)
    jac[nrow(jac),maxind] <- 1

    # Solve for the tangent vector to the curve of the 1st free parameter and
    # the state variables.
    tvnew <- solve(jac, c(rep(0, (nrow(jac)-1)), 1))
  }
  tvnorm <- sqrt(sum(tvnew^2))
  tvnew <- tvnew/tvnorm
  names(tvnew) <- NULL
  jac[nrow(jac),] <- NA

  return (list(jacobian = jac, tanvec = tvnew))
}

ExtSystem <- function(t, state, parms, model, ystart = NULL, tanvec = NULL, condfun = NULL, statedim, freeparsdim, nopts = NULL) {
  # Evaluate the equations
  rhsval <- unlist(model(t, state, parms))
  names(rhsval) <- NULL

  # add possible conditions for freeparsdim > 1 continuation
  if (!is.null(condfun)) {
    for (i in (1:length(condfun))) {
      rhsval <- do.call(condfun[[i]], list(state, parms, model, tanvec, statedim, freeparsdim, nopts = nopts, rhsval))
    }
    names(rhsval) <- NULL
  }

  # add the dot product of the difference between the current state and
  # the initial guess with tangent vector for pseudo-arclength continuation
  if (!is.null(ystart) && !is.null(tanvec)) {
    rhsval <- c(unlist(rhsval), c((state - ystart) %*% tanvec))
  }

  return(list(rhsval))
}

BPcontinuation <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  # Continuation of a branching point is carried out using the defining system
  # eq. 41 on page 36 of the Matcont documentation (august 2011):
  #
  #  F(x, p) + b*v   = 0
  #  (F_x(x, p))^T v = 0
  #  v^T F_p(x, p)   = 0
  #  v^T v - 1       = 0
  #
  #  with initial conditions b = 0 and v the eigenvector of the matrix
  #  (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
  #  The unknowns are:
  #
  #  p:  the bifurcation parameter
  #  x:  the solution point
  #  b:  an additional value
  #  v:  the eigenvector (same dimension as x)
  #
  # This function is only called when freeparsdim == 2, in which case the
  # Jacobian equals the following square (n+2)x(n+2) matrix of partial
  # derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state[1:(freeparsdim + statedim)], func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  Fx <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # Extract F_p(x, p)
  Fp <- jac[(1:statedim), 1]

  # Extract the value of b and v from the state
  eigval = unlist(state[freeparsdim + statedim + 1]);
  eigvec = unlist(state[(freeparsdim + statedim + 2):length(state)]);

  #  F(x, p) + b*v   = 0
  rhsval <- unlist(rhsval)
  rhsval[1:length(eigvec)] <- rhsval[1:length(eigvec)] + eigval*eigvec

  # (F_x(x, p))^T v = 0
  rhsval <- c(rhsval, c(t(Fx) %*% eigvec))

  # v^T F_p(x, p)   = 0
  rhsval <- c(rhsval, c(eigvec %*% Fp))

  # v^T v - 1       = 0
  rhsval <- c(rhsval, (c(eigvec %*% eigvec) - 1))

  return(rhsval)
}

HPcontinuation <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  #  When freeparsdim == 2 the Jacobian equals the following square (n+2)x(n+2)
  #  matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  # Otherwise, when freeparsdim == 1 the Jacobian equals the following square
  # (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  A <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # And return the determinant of the bialternate product
  twoAI <-bialt2AI(A)

  return(c(unlist(rhsval), det(twoAI)))
}

LPcontinuation <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  # Continuation of a limitpoint is carried out using the defining system
  # (10.97) on page 515 of Kuznetsov (1996):
  #
  #  F(x, p)                 = 0
  #  F_x(x, p) q             = 0
  #  (F_x(x, p))^T p - eps*p = 0
  #  q^T q - 1               = 0
  #  p^T p - 1               = 0
  #
  #  with initial conditions b = 0 and v the eigenvector of the matrix
  #  (F_x(x, p))^T pertaining to the eigenvalue with the smallest norm.
  #  The unknowns are:
  #
  #  p:   the bifurcation parameter
  #  x:   the solution point
  #  eps: an additional value
  #  q:   the eigenvector of F_x pertaining to the 0 eigenvalue
  #  p:   the eigenvector of (F_x)^T pertaining to the 0 eigenvalue
  #
  # The advantage of this defining system is that detection of Bogdanov-Takens
  # points and cusp point are straightforward

  # This function is only called when freeparsdim == 2, in which case the
  # Jacobian equals the following square (n+2)x(n+2) matrix of partial
  # derivatives:
  # This function is only called when freeparsdim == 2, in which case the
  # Jacobian equals the following square (n+2)x(n+2) matrix of partial
  # derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  A <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # Extract the additional variables from the state
  eps = unlist(state[freeparsdim + statedim + 1]);
  q = unlist(state[(freeparsdim + statedim + 2):(freeparsdim + statedim + 1 + statedim)]);
  p = unlist(state[(freeparsdim + statedim + 2 + statedim):(freeparsdim + statedim + 1 + 2*statedim)]);

  # Add the additional values
  # A q = 0
  rhsval <- c(unlist(rhsval), c(A %*% q))

  # A^T p - eps*p = 0
  rhsval <- c(unlist(rhsval), c((t(A) %*% p) - eps*p))

  # <q, q) - 1 =
  rhsval <- c(unlist(rhsval), c((q %*% q) - 1))

  # <p, p) - 1 =
  rhsval <- c(unlist(rhsval), c((p %*% p) - 1))

  return(unlist(rhsval))
}

EQ_BPtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  res <- TangentVecEQ(state, parms, model, tanvec, statedim, freeparsdim, nopts)
  jac <- res$jacobian
  jac[nrow(jac),] <- res$tanvec

  return(c(unlist(rhsval), det(jac)))
}

EQ_HPtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  #  When freeparsdim == 2 the Jacobian equals the following square (n+2)x(n+2)
  #  matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFn/dp1 dFn/dp2 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #           |   NA      NA      NA   ...    NA  |
  #
  # Otherwise, when freeparsdim == 1 the Jacobian equals the following square
  # (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
  #           |   NA      NA   ...    NA  |
  #
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system
  A <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # And return the determinant of the bialternate product
  twoAI <-bialt2AI(A)

  return(c(unlist(rhsval), det(twoAI)))
}

EQ_LPtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  res <- TangentVecEQ(state, parms, model, tanvec, statedim, freeparsdim, nopts)

  return(c(unlist(rhsval), res$tanvec[1]))
  # return(c(unlist(rhsval), det(res$jac[(1:statedim),((freeparsdim+1):(freeparsdim+statedim))])))
}

HP_BTtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  res <- TangentVecEQ(state, parms, model, tanvec, statedim, freeparsdim, nopts)
  jac <- res$jac[(1:statedim),(2:(statedim+1))]

  return(c(unlist(rhsval), det(jac)))
}

LP_BTtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  # Extract the additional variables from the state
  q = unlist(state[(freeparsdim + statedim + 2):(freeparsdim + statedim + 1 + statedim)]);
  p = unlist(state[(freeparsdim + statedim + 2 + statedim):(freeparsdim + statedim + 1 + 2*statedim)]);

  # See pg 515 of Kuznetsov (1996), just below eq. (10.97)

  return(c(unlist(rhsval), as.numeric(p %*% q)))
}

LP_CPtest <- function(state, parms, model, tanvec, statedim, freeparsdim, nopts, rhsval) {

  # Extract the additional variables from the state
  q = unlist(state[(freeparsdim + statedim + 2):(freeparsdim + statedim + 1 + statedim)]);
  p = unlist(state[(freeparsdim + statedim + 2 + statedim):(freeparsdim + statedim + 1 + 2*statedim)]);

  # The quantity B(q, q) can be computed most easily using a directional derivative
  # see eq. (10.52) and its approximation some lines below (10.52) on page 490 of
  # Kuznetsov (1996)
  s <- state
  s[(freeparsdim+1):(freeparsdim + statedim)] <- s[(freeparsdim+1):(freeparsdim + statedim)] + sqrt(nopts$jacdif)*q
  rhsP <- unlist(model(0, s, parms))
  s <- state
  s[(freeparsdim+1):(freeparsdim + statedim)] <- s[(freeparsdim+1):(freeparsdim + statedim)] - sqrt(nopts$jacdif)*q
  rhsM <- unlist(model(0, s, parms))
  Bqq <- (rhsP + rhsM)/nopts$jacdif

  # See bottom of pg 514 of Kuznetsov (1996), just above eq. (10.97)

  return(c(unlist(rhsval), as.numeric(p %*% Bqq)))
}

initBP <- function(state, parms, model, tanvec, statedim, freeparsdim, starttype, nopts) {

  # Compute the Jacobian
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system.
  # Because the initial point is BP, freeparsdim is known to equal 2
  Fx <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # Initialize the eigenvalue estimate to 0
  eigval <- 0.0

  # Find the eigenvector of (F_x(x, p))^T pertaining to the eigenvalue with
  # smallest absolute value
  eigvec <- approxNullVec(t(Fx))

  return(list(y = c(state, eigval, eigvec), tanvec = tanvec))
}

initEQ <- function(state, parms, model, tanvec, statedim, freeparsdim, starttype, nopts) {
  if (starttype == "BP") {
    #  The Jacobian equals the following square (n+1)x(n+1) matrix of partial
    #  derivatives:
    #
    #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
    #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
    #           |   .       .    ...    .   |
    #      Df = |   .       .    ...    .   |
    #           |dFn/dp1 dFn/dx1 ... dFn/dxn|
    #           |   NA      NA   ...    NA  |
    #
    # Here is state the vector of state variables (length n)

    jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)
    # Append the current tangent vector as the last row to the jacobian to
    # obtain the matrix D(0) (see eq.(10.61)-(10.63) on pg. 499 of Kuznetsov
    # (1996))
    jac[nrow(jac),] <- tanvec
    D <- jac

    # Approximate the vector q1 with the stored tangent vector
    q1 <- as.numeric(c(tanvec))
    names(q1) <- names(state)

    # Solve the vector q2 from D q2 = 0
    eig <- eigen(D)
    minindx <- which.min(abs(Re(eig$values)))
    q2 <- as.numeric(c(eig$vectors[,minindx]))
    names(q2) <- names(state)

    # Extract the (n+1)xn matrix (J(0))^T (see eq.(10.61)-(10.63) on pg. 499 of
    # Kuznetsov (1996))
    JT <- t(jac[(1:(nrow(jac)-freeparsdim)), (1:ncol(jac))])

    # Solve the vector phi from (J(0))^T phi = 0. According to Kuznetsov (1996,
    # pg. 497, eq. (10.59)) there is a unique vector (up to scalar multiple) phi
    # of length n satisfying this linear equation. To handle the arbitrary
    # scaling vector, drop that row of the matrix which yields the least
    # singular restricted matrix, such that eigenvaues can be determined easily
    # Estimate the condition numbers
    condnrs <- unlist(lapply((1:nrow(JT)),
                             function(i){
                               incrows <- rep(TRUE, nrow(JT))
                               incrows[i] <- FALSE
                               jacc <- JT[incrows,]
                               if (abs(det(jacc)) > 1.0E-5) return(rcond(jacc))
                               else return(0.0)
                             }))
    maxind <- which.max(condnrs)
    incrows <- rep(TRUE, nrow(JT))
    incrows[maxind] <- FALSE

    # For the restricted matrix solve for the eigenvalues and select the
    # eigenvector belonging to the eigenvalue closest to 0
    eig <- eigen(JT[incrows,])
    minindx <- which.min(abs(Re(eig$values)))
    phi <- c(eig$vectors[,minindx])

    # Now compute B(q1,q2) and B(q2, q2) as explained. The function B(x,y) is
    # defined in eq. (10.56) on pg. 496 of Kuznetsov (1996). The vectors
    # B(q1,q2) and B(q2, q2) of length n (n equals the number of state
    # variables) are needed to compute the factors b12 and b22 (see below)

    # The quantity B(q2, q2) can be computed most easily using a directional derivative
    # see eq. (10.52) and its approximation some lines below (10.52) on page 490 of
    # Kuznetsov (1996)
    rhsP1 <- unlist(model(0, state + sqrt(nopts$jacdif)*(q1 + q2), parms))
    rhsM1 <- unlist(model(0, state - sqrt(nopts$jacdif)*(q1 + q2), parms))
    rhsP2 <- unlist(model(0, state + sqrt(nopts$jacdif)*(q1 - q2), parms))
    rhsM2 <- unlist(model(0, state - sqrt(nopts$jacdif)*(q1 - q2), parms))

    Bq1q2 <- (rhsP1 + rhsM1 - rhsP2 - rhsM2)/(4*nopts$jacdif)

    rhsP <- unlist(model(0, state + sqrt(nopts$jacdif)*q2, parms))
    rhsM <- unlist(model(0, state - sqrt(nopts$jacdif)*q2, parms))
    Bq2q2 <- (rhsP + rhsM)/nopts$jacdif

    # The inner product of phi with the vectors B(q1,q2) and B(q2, q2) of length
    # n (n equals the number of state variables) yields the factors b12 and b22
    b12 <- c(phi %*% Bq1q2)
    b22 <- c(phi %*% Bq2q2)
    v2 <- (-b22/(2*b12))*q1 + q2
    v2 <- v2/sqrt(sum(v2^2))

    # Determine the relative change in the components and the index of the
    # largest change Indices of zero elements of y. Ignore there relative change
    indx0s <- (1:length(state))[abs(state) < as.numeric(nopts$iszero)]
    dy <- abs(v2)/(pmax(abs(state), as.numeric(nopts$iszero)))
    dy[indx0s] <- 0
    # Index with maximum relative change
    dyind <- which.max(dy)
    dy <- as.numeric(nopts$stepsize)*v2
    dyscaled <- max(abs(as.numeric(nopts$stepsize)*state[dyind]), as.numeric(nopts$minstepsize))*(dy/abs(dy[dyind]))

    guess <- state + c(as.numeric(nopts$stepsize)*dyscaled)

    return(list(y = guess, tanvec = v2))
  } else return(NULL)
}

initLP <- function(state, parms, model, tanvec, statedim, freeparsdim, starttype, nopts) {

  # To continue an LP curve we use the extended system of eqs. (10.76) on pg. 504 in
  # Kuznetsov (1996), as the matrix of this system has full rank at a Bogdanov-Takens
  # point. As additional values we need then an approximation to the null vector of
  # the Jacobian q, which we also use as additional vector q0.

  # Compute the Jacobian
  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system.
  # Because the initial point is BP, freeparsdim is known to equal 2
  Fx <- jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))]

  # Initialize the eigenvalue estimate to 0
  eigval <- 0.0

  # Find the eigenvector of F_x(x, p) pertaining to the eigenvalue with
  # smallest absolute value
  q <- approxNullVec(Fx)

  # Find the eigenvector of (F_x(x, p))^T pertaining to the eigenvalue with
  # smallest absolute value
  p <- approxNullVec(t(Fx))

  return(list(y = c(state, eigval, as.numeric(q), as.numeric(p))))
}

analyseEQ <- function(state, parms, model, tanvec, statedim, freeparsdim, lastvals, nopts, session) {

  bpval <- EQ_BPtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(bpval) <- NULL
  hpval <- EQ_HPtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(hpval) <- NULL
  lpval <- EQ_LPtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(lpval) <- NULL

  testvals <- list()
  testvals$bpval <- unlist(bpval)
  testvals$hpval <- unlist(hpval)
  testvals$lpval <- unlist(lpval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    if (!is.null(lastvals$hpval) && ((lastvals$hpval)*hpval < -(nopts$atol*nopts$atol))) {
      biftype = "HP"
    }
    if (!is.null(lastvals$lpval) && ((lastvals$lpval)*lpval < -(nopts$atol*nopts$atol))) {
      biftype = "LP"
    }
    else if (!is.null(lastvals$bpval) && ((lastvals$bpval)*bpval < -(nopts$atol*nopts$atol))) {
      biftype = "BP"
    }
    if (!is.null(biftype)) {
      cfun <- get(paste0("EQ_", biftype, "test"), mode = "function")
      res <- tryCatch(stode(state, time = 0, func = ExtSystem, parms = parms, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                            maxiter = nopts$maxiter, verbose = FALSE, model = model, condfun = c(cfun),
                            tanvec = tanvec, statedim = statedim, freeparsdim = freeparsdim, nopts = nopts),
                      warning = function(e) {
                        msg <- gsub(".*:", "Warning in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      },
                      error = function(e) {
                        msg <- gsub(".*:", "Error in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      })
      if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
        y <- res$y[1:length(state)]
        names(y) <- names(state)

        # Compute the Jacobian w.r.t. to the free parameter and the state variables
        jac <- jacobian.full(y=y, func=model, parms=parms, pert = nopts$jacdif)

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:statedim), function(i){paste0("Eigenvalue", i)}))

        # Append the current tangent vector as the last row to the jacobian to
        # preserve direction. See the matcont manual at
        # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
        # Notice that jacobian.full returns a square matrix with NA values on the last row
        jac[nrow(jac),] <- tanvec
        if (rcond(jac) > nopts$atol) {
          tvnew <- solve(jac, c(rep(0, (length(tanvec)-1)), 1))
          tvnorm <- sqrt(sum(tvnew^2))
          tvnew <- tvnew/tvnorm
          names(tvnew) <- unlist(lapply((1:length(state)), function(i){paste0("d", names(state)[i])}))
        } else {
          tvnew <- tanvec
        }

        if (biftype == "BP") testvals$bpval <- NULL
        else if (biftype == "HP") testvals$hpval <- NULL
        else testvals$lpval <- NULL

        testvals[["y"]] <- y
        testvals[["tanvec"]] <- tvnew
        testvals[["eigval"]] <- eigval
        testvals[["biftype"]] <- biftype
      } else {
        if (biftype == "BP") msg <- "Locating branching point failed\n"
        else if (biftype == "HP") msg <- "Locating Hopf bifurcation point failed\n"
        else msg <- "Locating limit point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}

analyseHP <- function(state, parms, model, tanvec, statedim, freeparsdim, lastvals, nopts, session) {

  btval <- HP_BTtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(btval) <- NULL

  testvals <- list()
  testvals$btval <- unlist(btval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    if (!is.null(lastvals$btval) && ((lastvals$btval)*btval < -(nopts$atol*nopts$atol))) {
      biftype = "BT"
    }
    if (!is.null(biftype)) {
      cfun <- get(paste0("HP_", biftype, "test"), mode = "function")
      res <- tryCatch(stode(state, time = 0, func = ExtSystem, parms = parms, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                            maxiter = nopts$maxiter, verbose = FALSE, model = model, condfun = c(HPcontinuation, cfun),
                            tanvec = tanvec, statedim = statedim, freeparsdim = freeparsdim, nopts = nopts),
                      warning = function(e) {
                        msg <- gsub(".*:", "Warning in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      },
                      error = function(e) {
                        msg <- gsub(".*:", "Error in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      })
      if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
        y <- res$y[1:length(state)]
        names(y) <- names(state)

        # Compute the Jacobian w.r.t. to the free parameter and the state variables
        jac <- jacobian.full(y=y, func=model, parms=parms, pert = nopts$jacdif)

        # Discard all th rows and columns that pertain to additional variables
        jac <- jac[(1:(freeparsdim+statedim)), (1:(freeparsdim+statedim))]

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:statedim), function(i){paste0("Eigenvalue", i)}))

        # Append the current tangent vector as the last row to the jacobian to
        # preserve direction. See the matcont manual at
        # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
        # Notice that jacobian.full returns a square matrix with NA values on the last row
        jac[nrow(jac),] <- tanvec[(1:(freeparsdim+statedim))]
        if (rcond(jac) > nopts$atol) {
          tvnew <- solve(jac, c(rep(0, ((freeparsdim+statedim)-1)), 1))
          tvnorm <- sqrt(sum(tvnew^2))
          tvnew <- tvnew/tvnorm
          names(tvnew) <- unlist(lapply((1:length(state)), function(i){paste0("d", names(state)[i])}))
        } else {
          tvnew <- tanvec[(1:(freeparsdim+statedim))]
        }

        if (biftype == "BT") testvals$btval <- NULL

        testvals[["y"]] <- y
        testvals[["tanvec"]] <- tvnew
        testvals[["eigval"]] <- eigval
        testvals[["biftype"]] <- biftype
      } else {
        if (biftype == "BT") msg <- "Locating Bogdanov-Takens point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}

analyseLP <- function(state, parms, model, tanvec, statedim, freeparsdim, lastvals, nopts, session) {

  btval <- LP_BTtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(btval) <- NULL
  cpval <- LP_CPtest(state, parms, model, tanvec, statedim, freeparsdim, nopts, NULL)
  names(cpval) <- NULL

  testvals <- list()
  testvals$btval <- unlist(btval)
  testvals$cpval <- unlist(cpval)
  biftype <- NULL

  if (!is.null(lastvals)) {
    if (!is.null(lastvals$cpval) && ((lastvals$cpval)*cpval < -(nopts$atol*nopts$atol))) {
      biftype = "CP"
    }
    if (!is.null(lastvals$btval) && ((lastvals$btval)*btval < -(nopts$atol*nopts$atol))) {
      biftype = "BT"
    }
    if (!is.null(biftype)) {
      cfun <- get(paste0("LP_", biftype, "test"), mode = "function")
      res <- tryCatch(stode(state, time = 0, func = ExtSystem, parms = parms, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                            maxiter = nopts$maxiter, verbose = FALSE, model = model, condfun = c(LPcontinuation, cfun),
                            tanvec = tanvec, statedim = statedim, freeparsdim = freeparsdim, nopts = nopts),
                      warning = function(e) {
                        msg <- gsub(".*:", "Warning in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      },
                      error = function(e) {
                        msg <- gsub(".*:", "Error in rootSolve:", e)
                        if (!is.null(session)) updateConsoleLog(session, msg)
                        else cat(msg)
                        return(NULL)
                      })
      if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
        y <- res$y[1:length(state)]
        names(y) <- names(state)

        # Compute the Jacobian w.r.t. to the free parameter and the state variables
        jac <- jacobian.full(y=y, func=model, parms=parms, pert = nopts$jacdif)

        # Discard all th rows and columns that pertain to additional variables
        jac <- jac[(1:(freeparsdim+statedim)), (1:(freeparsdim+statedim))]

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:statedim), ((freeparsdim+1):(freeparsdim+statedim))])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:statedim), function(i){paste0("Eigenvalue", i)}))

        # Append the current tangent vector as the last row to the jacobian to
        # preserve direction. See the matcont manual at
        # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
        # Notice that jacobian.full returns a square matrix with NA values on the last row
        jac[nrow(jac),] <- tanvec[(1:(freeparsdim+statedim))]
        if (rcond(jac) > nopts$atol) {
          tvnew <- solve(jac, c(rep(0, ((freeparsdim+statedim)-1)), 1))
          tvnorm <- sqrt(sum(tvnew^2))
          tvnew <- tvnew/tvnorm
          names(tvnew) <- unlist(lapply((1:length(state)), function(i){paste0("d", names(state)[i])}))
        } else {
          tvnew <- tanvec[(1:(freeparsdim+statedim))]
        }

        if (biftype == "BT") testvals$btval <- NULL
        else testvals$cpval <- NULL

        testvals[["y"]] <- y
        testvals[["tanvec"]] <- tvnew
        testvals[["eigval"]] <- eigval
        testvals[["biftype"]] <- biftype
      } else {
        if (biftype == "BT") msg <- "Locating Bogdanov-Takens point failed\n"
        else msg <- "Locating cusp point failed\n"

        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
      }
    }
  }
  return(testvals)
}
