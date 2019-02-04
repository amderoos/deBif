extsysEQ <- function(t, state, parms, model, ystart, tanvec) {
  # Evaluate the equations
  fval <- model(0, state, parms)
  print(c(unlist(fval), c((state - ystart) %*% tanvec)))

  # and add the dot product of the difference between the current state and
  # the initial guess with tangent vector for pseudo-arclength continuation
  return(list(c(unlist(fval), c((state - ystart) %*% tanvec))))
}

computeCurve <- function(model, state, parms, freepars, popts, nopts, curvetype, tanvec = NULL, output = NULL, updateProgress = NULL) {
  if (exists("deBifdebug", envir = .GlobalEnv)) debug <- get("deBifdebug", envir = .GlobalEnv)
  else debug <- FALSE

  if (exists("deBifverbose", envir = .GlobalEnv)) verbose <- get("deBifverbose", envir = .GlobalEnv)
  else verbose <- FALSE

  if (!is.numeric(freepars)) x <- index(freepars, names(parms))
  codim <- length(freepars)

  y <- c(parms[freepars], state)
  fixedpars <- parms[!(1:length(parms)) %in% freepars]
  extsys <- get(paste0("extsys", curvetype), mode = "function")
  if (is.null(tanvec))  tanvec <- c(1.0, rep(0, (length(y)-1)))

  guess <- y
  pntnr <- 1
  eignames <- unlist(lapply((1:length(state)), function(i){paste0("Eigenvalue", i)}))
  tvnames <- unlist(lapply((1:length(y)), function(i){paste0("d", names(y)[i])}))

  allsols <- NULL
  alltvs <- NULL
  alleigs <- NULL

  if (debug) {
    assign("state", state, env = .GlobalEnv)
    assign("parms", parms, env = .GlobalEnv)
    assign("freepars", freepars, env = .GlobalEnv)
    assign("popts", popts, env = .GlobalEnv)
    assign("nopts", nopts, env = .GlobalEnv)
    assign("curvetype", curvetype, env = .GlobalEnv)
    assign("tanvec", tanvec, env = .GlobalEnv)
  }
  stepscalefac <- 1
  while (pntnr < as.numeric(nopts$maxpoints)) {
    if (verbose) {
      res <- rootSolve::stode(guess, time = 0, func = extsys, parms = fixedpars, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol, maxiter = nopts$maxiter, verbose = TRUE,
                              model = model, ystart = guess, tanvec = tanvec)
    } else {
      {sink("/dev/null");
        res <- rootSolve::stode(guess, time = 0, func = extsys, parms = fixedpars, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol, maxiter = nopts$maxiter,
                                model = model, ystart = guess, tanvec = tanvec);
        sink()}
    }
    if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
      y <- res$y

      # Compute the Jacobian w.r.t. to free parameters (the first codim columns) and the state variables
      jac <- jacobian.full(y=y, func=model, parms=fixedpars)

      # Compute the eigenvalues of the restricted Jacobian (exclude the first codim columns)
      eig <- eigen(jac[(1:length(state)),((codim+1):ncol(jac))])
      # Sort them on decreasing real part
      eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]

      # Report the solution point and the eigenvalues
      # msg <- paste0("\nSolution ", pntnr, " found:\n\n", paste(unlist(lapply(1:length(y), function(i) {paste0(names(y[i]), "=", round(y[i], 4))})), collapse=', '), "\n\n")
      # msg <- paste0(msg, "Eigenvalues:\n\n", paste(unlist(lapply(1:length(eigval), function(i) {paste0(round(eigval[i], 4))})), collapse=' '), "\n\n")
      msg <- paste0("\nSolution ", pntnr, " found:\n\n", paste(unlist(lapply(1:length(y), function(i) {paste0(names(y[i]), "=", sprintf("%12.5E", y[i]))})), collapse=', '), "\n\n")
      msg <- paste0(msg, "Eigenvalues:\n\n", paste(unlist(lapply(1:length(eigval), function(i) {paste0(sprintf("%12.5E", eigval[i]))})), collapse=' '), "\n\n")
      if (debug) cat(msg)
      else if (!is.null(updateProgress)) updateProgress(detail = msg)

      names(eigval) <- eignames

      # Append the current tangent vector as the last row to the jacobian to
      # preserve direction. See the matcont manual at
      # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
      # Notice that jacobian.full returns a square matrix with NA values on the last row
      jac[nrow(jac),] <- tanvec
      tvnew <- solve(jac, c(rep(0, (length(tanvec)-1)), 1))
      tvnorm <- sqrt(sum(tvnew^2))
      tvnew <- tvnew/tvnorm
      names(tvnew) <- tvnames

      # Store the results
      allsols <- rbind(allsols, c(y))
      alltvs <- rbind(alltvs, c(tvnew))
      alleigs <- rbind(alleigs, c(eigval))

      # Determine the relative change in the components and the index of the largest change
      indx0s <- (1:length(y))[abs(y) < as.numeric(nopts$iszero)]   # Indices of zero elements of y. Ignore there relative change
      dy <- abs(tvnew)/(pmax(abs(y), as.numeric(nopts$iszero)))
      dy[indx0s] <- 0
      dyind <- which.max(dy)                      # Index with maximum relative change

      dy <- as.numeric(nopts$stepsize)*tvnew

      dyscaled <- max(abs(as.numeric(nopts$stepsize)*y[dyind]), as.numeric(nopts$minstepsize))*(dy/abs(dy[dyind]))

      guess <- y + dyscaled/stepscalefac
      yold <- y
      tanvec <- tvnew

      # Stop the curve if outside the visible plotting region
      if (pntnr > 10) {
        if (pntnr > as.numeric(nopts$maxpoints) ||
            (y[1] < (as.numeric(popts$xmin) - as.numeric(nopts$iszero))) ||
            (y[1] > (as.numeric(popts$xmax) + as.numeric(nopts$iszero)))) {
          msg <- "Minimum or maximum of x-axis domain reached: Computation halted"
          if (debug) cat(msg)
          else if (!is.null(output)) output[["console"]] <- renderText({msg})
          break
        }

        if ((curvetype == "EQ") && (popts$ycol == 1)) {
          if (all(y[(2:(length(state)+1))] < (as.numeric(popts$ymin) - as.numeric(nopts$iszero)))) {
            msg <- "Minimum of y-axis domain reached for all y-axis variables: Computation halted"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- renderText({msg})
            break
          }
          if (all(y[(2:(length(state)+1))] > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Maximum of y-axis domain reached for all y-axis variables: Computation halted"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- renderText({msg})
            break
          }
        } else {
          if ((y[popts$ycol] < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
              (y[popts$ycol] > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Minimum or maximum of y-axis domain reached for 1st y-axis variable: Computation halted"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- renderText({msg})
            break
          }
          if ((popts$y2col > 1) &&
              ((y[popts$y2col] < (as.numeric(popts$y2min) - as.numeric(nopts$iszero))) ||
               (y[popts$y2col] > (as.numeric(popts$y2max) + as.numeric(nopts$iszero))))) {
            msg <- "Maximum or maximum of y-axis domain reached for 2nd y-axis variable: Computation halted"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- renderText({msg})
            break
          }
        }
      }
    } else {                                      # Solution not found
      if (pntnr == 1) {
        msg <- "No convergence at initial point"
        cat(msg)
        return(NULL)
      }
      if (stepscalefac > 4096) {
        msg <- "Step size too small"
        cat(msg)
        return(NULL)
      }
      stepscalefac <- stepscalefac/2
      guess <- y + dyscaled/stepscalefac
    }
    pntnr <- pntnr + 1

    if (pntnr >= as.numeric(nopts$maxpoints)) {
      msg <- "Maximum number of points along the curve reached: Computation halted"
      if (debug) cat(msg)
      else if (!is.null(output)) output[["console"]] <- renderText({msg})
    }

    # Take a breath
    Sys.sleep(as.numeric(nopts$computedelay))
  }
  return(list("points" = allsols, "eig.vals" = alleigs, "tangent" = alltvs))
}
