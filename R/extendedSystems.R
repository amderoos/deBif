extsys <- function(t, state, parms, model, ystart = NULL, tanvec = NULL, condfun = NULL, codim = 1, nopts = NULL) {
  # Evaluate the equations
  rhsval <- unlist(model(t, state, parms))
  names(rhsval) <- NULL

  # add possible conditions for codim > 1 continuation
  if (!is.null(condfun)) {
    rhsval <- condfun(state, parms, model, codim, nopts = nopts, rhsval)
    names(rhsval) <- NULL
  }

  # add the dot product of the difference between the current state and
  # the initial guess with tangent vector for pseudo-arclength continuation
  if (!is.null(ystart) && !is.null(tanvec)) {
    rhsval <- c(unlist(rhsval), c((state - ystart) %*% tanvec))
  }

  return(list(rhsval))
}

BPcondition <- function(state, parms, model, codim, nopts, rhsval) {

  #  When codim == 2 the Jacobian equals the following square
  #  (n+2)x(n+2) matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFm/dp1 dFm/dp2 dFm/dx1 ... dFm/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #
  # Here is the vector of state variables (lenght n) and
  # m = (n+1)
  #
  # Otherwise, when codim == 1 the Jacobian equals the following
  # square (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFm/dp1 dFm/dx1 ... dFm/dxn|
  #
  # In which m = n+1 (i.e. equal to the number of state variables plus 1).

  rhsdim <- length(unlist(rhsval))
  jac <- jacobian.full(y=state[1:(codim + rhsdim)], func=model, parms=parms, pert = nopts$jacdif)

  # If codim == 2 the Jacobian has two rows of NA at the end, because
  # the state contains 2 free parameters and hence the dimension of the vector
  # y=state is 2 larger than the right-hand side of the system of ODEs. Delete
  # the last row and 2nd column of the Jacobian, the latter being the derivative
  # w.r.t. second free parameters

  # Extract the restricted Jacobian of the system
  Fx <- jac[(1:rhsdim), ((codim+1):(codim+rhsdim))]

  # Extract F_p(x, p)
  Fp <- jac[(1:rhsdim), 1]

  # The extended system to solve for is see the Matcont documentation (Branch
  # point locator, page 36, eq. 41)
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

  # The value of ncol(jac) equals the number of free parameters plus the
  # the number of state variables
  eigval = unlist(state[codim + rhsdim + 1]);
  eigvec = unlist(state[(codim + rhsdim + 2):length(state)]);

  #  F(x, p) + b*v   = 0
  rhsval <- unlist(rhsval) + eigval*eigvec

  # (F_x(x, p))^T v = 0
  rhsval <- c(rhsval, c(t(Fx) %*% eigvec))

  # v^T F_p(x, p)   = 0
  rhsval <- c(rhsval, c(eigvec %*% Fp))

  # v^T v - 1       = 0
  rhsval <- c(rhsval, (c(eigvec %*% eigvec) - 1))

  return(rhsval)
}

HPcondition <- function(state, parms, model, codim, nopts, rhsval) {

  #  When codim == 2 the Jacobian equals the following square
  #  (n+2)x(n+2) matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFm/dp1 dFm/dp2 dFm/dx1 ... dFm/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #
  # Here is the vector of state variables (length n) and
  # m = (n+1)
  #
  # Otherwise, when codim == 1 the Jacobian equals the following
  # square (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFm/dp1 dFm/dx1 ... dFm/dxn|
  #
  # In which m = n+1 (i.e. equal to the number of state variables plus 1).

  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # If codim == 2 the Jacobian has two rows of NA at the end, because the state
  # contains 2 free parameters and hence the dimension of the vector y=state is
  # 2 larger than the right-hand side of the system of ODEs.
  # Extract the restricted Jacobian of the system
  A <- jac[(1:(nrow(jac)-codim)), ((codim+1):ncol(jac))]

  # Now construct the bialternate matrix product of 2A and I.
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

  # And return the determinant of the bialternate product

  return(c(unlist(rhsval), det(twoAI)))
}

LPcondition <- function(state, parms, model, codim, nopts, rhsval) {

  #  When codim == 2 the Jacobian equals the following square
  #  (n+2)x(n+2) matrix of partial derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFm/dp1 dFm/dp2 dFm/dx1 ... dFm/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #
  # Here is the vector of state variables (length n) and
  # m = (n+1)
  #
  # Otherwise, when codim == 1 the Jacobian equals the following
  # square (n+1)x(n+1) matrix
  #
  #           |dF1/dp1 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dx1 ... dF2/dxn|
  #           |   .       .    ...    .   |
  #      Df = |   .       .    ...    .   |
  #           |   .       .    ...    .   |
  #           |dFm/dp1 dFm/dx1 ... dFm/dxn|
  #
  # In which m = n+1 (i.e. equal to the number of state variables plus 1).

  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # If codim == 2 the Jacobian has two rows of NA at the end, because
  # the state contains 2 free parameters and hence the dimension of the vector
  # y=state is 2 larger than the right-hand side of the system of ODEs. Delete
  # the last row and 2nd column of the Jacobian, the latter being the derivative
  # w.r.t. second free parameters
  if (codim == 2) jac <- jac[(1:(nrow(jac)-1)), c(1, (3:ncol(jac)))]

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

  # Solve for the tangent vector to the curve of the 1st free parameter and
  # the state variables.
  jac[nrow(jac),maxind] <- 1
  tvlp <- solve(jac, c(rep(0, (nrow(jac)-1)), 1))
  names(tvlp) <- NULL

  return(c(unlist(rhsval), tvlp[1]))
}

initBP <- function(state, parms, model, tanvec, codim, starttype, nopts) {
  #  The Jacobian equals the following square (n+2)x(n+2) matrix of partial
  #  derivatives:
  #
  #           |dF1/dp1 dF1/dp2 dF1/dx1 ... dF1/dxn|
  #           |dF2/dp1 dF2/dp2 dF2/dx1 ... dF2/dxn|
  #           |   .       .       .    ...    .   |
  #      Df = |   .       .       .    ...    .   |
  #           |dFm/dp1 dFm/dp2 dFm/dx1 ... dFm/dxn|
  #           |   NA      NA      NA   ...    NA  |
  #
  # Here is state the vector of state variables (length n) and m = (n+1)

  jac <- jacobian.full(y=state, func=model, parms=parms, pert = nopts$jacdif)

  # Extract the restricted Jacobian of the system.
  # Because the initial point is BP, codim is known to equal 2
  Fx <- jac[(1:(nrow(jac)-codim)), ((codim+1):ncol(jac))]

  # Find the eigenvector of (F_x(x, p))^T pertaining to the eigenvalue with
  # smallest absolute value
  eig <- eigen(t(Fx))
  minindx <- which.min(abs(Re(eig$values)))
  eigvec <- c(eig$vectors[,minindx])
  # Initialize the eigenvalue estimate to 0
  eigval <- 0.0

  return(list(y = c(state, eigval, eigvec), tanvec = tanvec))
}

initEQ <- function(state, parms, model, tanvec, codim, starttype, nopts) {
  if(starttype == "BP") {
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
    q1 <- tanvec

    # Solve the vector q2 from D q2 = 0
    eig <- eigen(D)
    minindx <- which.min(abs(Re(eig$values)))
    q2 <- as.numeric(c(eig$vectors[,minindx]))

    # Extract the (n+1)xn matrix (J(0))^T (see eq.(10.61)-(10.63) on pg. 499 of
    # Kuznetsov (1996))
    JT <- t(jac[(1:(nrow(jac)-codim)), (1:ncol(jac))])

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
    Bq1q2 <- unlist(lapply((1:(length(state)-codim)),
                           function(i) {
                             fun <- function(x) {res <- unlist(model(0, x, parms)); return(res[i])}
                             bb <- hessian(fun, state, pert = nopts$jacdif)
                             q1m <- matrix(as.numeric(q1), nrow(bb), ncol(bb), byrow = FALSE)
                             q2m <- matrix(as.numeric(q2), nrow(bb), ncol(bb), byrow = TRUE)
                             return(sum(sum((bb * q1m) * q2m)))
                           }))
    Bq2q2 <- unlist(lapply((1:(length(state)-codim)),
                           function(i) {
                             fun <- function(x) {res <- unlist(model(0, x, parms)); return(res[i])}
                             bb <- hessian(fun, state, pert = nopts$jacdif)
                             q1m <- matrix(as.numeric(q2), nrow(bb), ncol(bb), byrow = FALSE)
                             q2m <- matrix(as.numeric(q2), nrow(bb), ncol(bb), byrow = TRUE)
                             return(sum(sum((bb * q1m) * q2m)))
                           }))

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

testEQ <- function(state, parms, model, tanvec, extjac, lastvals, nopts, output) {
  if (exists("deBifdebug", envir = .GlobalEnv)) debug <- get("deBifdebug", envir = .GlobalEnv)
  else debug <- FALSE

  rhsval <- unlist(model(0, state, parms))
  rhsdim <- length(rhsval)
  names(rhsval) <- NULL

  hpval <- HPcondition(state, parms, model, codim = 1, nopts, NULL)
  lpval <- tanvec[1]
  names(lpval) <- NULL

  # See pg. 499, eq. (10.62) in Kuznetsov (1996). The variable extjac is the
  # extended Jacobian with the tangent to the solution curve added as additional
  # row
  bpval <- det(extjac)

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
      cfun <- get(paste0(biftype, "condition"), mode = "function")
      if (biftype == "BP") {
        initres <- initBP(state, parms, model, tanvec, 1, "EQ", nopts)
        if (!is.null(initres)) guess <- initres$y
      } else guess <- state

      res <- tryCatch(stode(guess, time = 0, func = extsys, parms = parms, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                            maxiter = nopts$maxiter, verbose = FALSE, model = model, condfun = cfun, codim = 1, nopts = nopts),
                      warning = function(e) {
                        msg <- gsub(".*:", "Warning in rootSolve:", e)
                        if (debug) cat(msg)
                        else if (!is.null(output)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
                        return(NULL)
                      },
                      error = function(e) {
                        msg <- gsub(".*:", "Error in rootSolve:", e)
                        if (debug) cat(msg)
                        else if (!is.null(output)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
                        return(NULL)
                      })
      if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
        y <- res$y[1:length(state)]
        names(y) <- names(state)

        # Compute the Jacobian w.r.t. to the free parameter and the state variables
        jac <- jacobian.full(y=y, func=model, parms=parms, pert = nopts$jacdif)

        # Compute the eigenvalues of the restricted Jacobian
        eig <- eigen(jac[(1:rhsdim), (2:ncol(jac))])

        # Sort them on decreasing real part
        eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
        names(eigval) <-  unlist(lapply((1:rhsdim), function(i){paste0("Eigenvalue", i)}))

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

        if (debug) cat(msg)
        else if (!is.null(output)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
      }
    }
  }
  return(testvals)
}
