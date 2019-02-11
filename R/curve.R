computeCurve <- function(model, state, parms, freepars, popts, nopts, starttype, curvetype, tanvec = NULL, report2console = TRUE, session = NULL, output = NULL) {

  if (exists("deBifdebug", envir = .GlobalEnv)) debug <- get("deBifdebug", envir = .GlobalEnv)
  else debug <- FALSE

  if (exists("deBifverbose", envir = .GlobalEnv)) verbose <- get("deBifverbose", envir = .GlobalEnv)
  else verbose <- FALSE

  if (!(curvetype %in% c("EQ", "BP", "HP", "LP"))) {
    msg <- paste0("Computation aborted:\nContinuation for curve type ", curvetype, " not implemented\n")
    if (debug) cat(msg)
    else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
    return(NULL)
  }

  if (!is.numeric(freepars)) x <- index(freepars, names(parms))
  codim <- length(freepars)
  rhsdim <- length(state)

  y <- c(parms[freepars], state)
  fixedpars <- parms[!(1:length(parms)) %in% freepars]

  varnames <- names(y)
  eignames <- unlist(lapply((1:rhsdim), function(i){paste0("Eigenvalue", i)}))
  tvnames <- unlist(lapply((1:length(y)), function(i){paste0("d", varnames[i])}))

  if (curvetype == "EQ") {
    condfun <- NULL
  } else {
    if (exists(paste0(curvetype, "condition"), mode = "function")) condfun <- get(paste0(curvetype, "condition"), mode = "function")
    else {
      msg <- paste0("Computation aborted:\nAdditional condition function for curve type ", curvetype, " not found\n")
      if (debug) cat(msg)
      else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
      return(NULL)
    }
  }

  # Test functions will only be implemented for EQ curves
  if (exists(paste0("test", curvetype), mode = "function")) {
    testfun <- get(paste0("test", curvetype), mode = "function")
    testvals <- NULL
  }
  else testfun <- NULL

  # If present run the initializer
  if (exists(paste0("init", curvetype), mode = "function")) {
    initfun <- get(paste0("init", curvetype), mode = "function")
    names(y) <- varnames
    initres <- initfun(y, fixedpars, model, tanvec, codim, starttype, nopts)
    if (!is.null(initres)) {
      y <- c(as.numeric(initres$y))
      names(y) <- varnames
      tanvec <- c(as.numeric(initres$tanvec))
    }
  }

  if (is.null(tanvec) || (length(tanvec) != length(y))) tanvec <- c(1.0, rep(0, (length(y)-1)))
  else if (abs(tanvec[1]) > nopts$iszero) tanvec <- sign(tanvec[1])*tanvec

  guess <- y
  pntnr <- 1

  allsols <- NULL
  alltvs <- NULL
  alleigs <- NULL
  specialsols <- NULL
  specialtvs <- NULL
  specialeigs <- NULL
  specialtags <- NULL

  stepscalefac <- 1
  while (pntnr < as.numeric(nopts$maxpoints)) {
    if (verbose) {
      res <- tryCatch(stode(guess, time = 0, func = extsys, parms = fixedpars, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                            maxiter = nopts$maxiter, verbose = TRUE, model = model, ystart = guess, tanvec = tanvec, condfun = condfun, codim = codim, nopts = nopts),
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
    } else {
      res <- tryCatch(stode(guess, time = 0, func = extsys, parms = fixedpars, rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol, maxiter = nopts$maxiter,
                            model = model, ystart = guess, tanvec = tanvec, condfun = condfun, codim = codim, nopts = nopts),
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
                      });
    }
    if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {                    # Solution found
      y <- res$y

      # Compute the Jacobian w.r.t. to free parameters (the first codim columns) and the state variables

      ##### This has to be resolved. Do we call the routine of the model equations or the extended system?
      # I think the latter, because we want to compute the derivatives of the additional conditions w.r.t.
      # all free parameters and variables to compute the tangent to the curve
      jac <- jacobian.full(y=y, func=extsys, parms=fixedpars, pert = nopts$jacdif, model = model, condfun = condfun, codim = codim, nopts = nopts)

      # Compute the eigenvalues of the restricted Jacobian (exclude the first codim columns
      # and takink only as many rows as there are state variables)
      eig <- eigen(jac[(1:rhsdim),((codim+1):(codim+rhsdim))])
      # Sort them on decreasing real part
      eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
      names(eigval) <- eignames

      # Append the current tangent vector as the last row to the jacobian to
      # preserve direction. See the matcont manual at
      # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
      # Notice that jacobian.full returns a square matrix with NA values on the last row
      jac[nrow(jac),] <- tanvec
      if (rcond(jac) > nopts$atol) {
        tvnew <- solve(jac, c(rep(0, (length(tanvec)-1)), 1))
        tvnorm <- sqrt(sum(tvnew^2))
        tvnew <- tvnew/tvnorm
        names(tvnew) <- tvnames
      } else {
        tvnew <- tanvec
      }

      ############## Execute the test functions and store the results
      if (!is.null(testfun)) {
        testvals <- testfun(state = y, parms = fixedpars, model = model, tanvec = tvnew, extjac = jac, lastvals = testvals, nopts, output = output)
        if (!is.null(testvals) && ("y" %in% names(testvals) > 0) && (length(testvals$y) > 0)) {

          allsols <- rbind(allsols, c(testvals$y[1:length(y)]))
          alltvs <- rbind(alltvs, c(testvals$tanvec[1:length(tanvec)]))
          alleigs <- rbind(alleigs, c(testvals$eigval[1:length(eigval)]))

          specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
          specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
          specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))

          dscp <- paste(unlist(lapply(1:length(c(testvals$y)), function(i) {paste0(names(y[i]), "=", round(c(testvals$y)[i], 3))})), collapse=', ')
          specialtags <- rbind(specialtags, c("Type" = testvals$biftype, "Description" = paste0(testvals$biftype, ": ", dscp)))

          if (testvals$biftype == "BP") msg <- paste("Branching point found:\n", dscp, "\n", sep=" ")
          else if (testvals$biftype == "HP") msg <- paste("Hopf bifurcation point found:\n", dscp, "\n", sep=" ")
          else msg <- paste("Limit point found:\n", dscp, "\n", sep=" ")
          msg <- paste0(msg, "Eigenvalues:\n",
                        paste(unlist(lapply(1:length(eigval),
                                            function(i) {paste0(ifelse(is.complex(testvals$eigval[i]), sprintf("%12.5E + %12.5Ei", Re(testvals$eigval[i]), Im(testvals$eigval[i])), sprintf("%12.5E", testvals$eigval[i])))})), collapse=' '), "\n")

          if (debug) cat(msg)
          else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)

          testvals$y <- NULL
          testvals$tanvec <- NULL
          testvals$eigval <- NULL
          testvals$biftype <- NULL

          pntnr <- pntnr + 1
        }
      }

      ############## Report the solution point and the eigenvalues
      # msg <- paste0("\nSolution ", pntnr, " found:\n\n", paste(unlist(lapply(1:length(y), function(i) {paste0(names(y[i]), "=", round(y[i], 4))})), collapse=', '), "\n\n")
      # msg <- paste0(msg, "Eigenvalues:\n\n", paste(unlist(lapply(1:length(eigval), function(i) {paste0(round(eigval[i], 4))})), collapse=' '), "\n\n")
      msg <- paste0("Solution ", pntnr, " found:\n", paste(unlist(lapply(1:length(y), function(i) {paste0(names(y[i]), "=", sprintf("%12.5E", y[i]))})), collapse=', '), "\n")
      msg <- paste0(msg, "Eigenvalues:\n",
                    paste(unlist(lapply(1:length(eigval),
                                        function(i) {paste0(ifelse(is.complex(eigval[i]), sprintf("%12.5E + %12.5Ei", Re(eigval[i]), Im(eigval[i])), sprintf("%12.5E", eigval[i])))})), collapse=' '), "\n")
      if (debug) cat(msg)
      else if (!is.null(output)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))

      if (report2console) {
        yy <- c(y, eigval)
        if (pntnr == 1) names(yy) <- c(names(y), eignames)
        else names(yy) <- NULL
        print(c(yy), print.gap = 2, digits = 6, right = TRUE)
      }

      ############## Store the results
      allsols <- rbind(allsols, c(y))
      alltvs <- rbind(alltvs, c(tvnew))
      alleigs <- rbind(alleigs, c(eigval))

      ############## Take a new step along the curve
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

      ############## Stop the curve if outside the visible plotting region
      pntnr <- pntnr + 1

      if (pntnr >= as.numeric(nopts$maxpoints)) {
        msg <- "Computation halted:\nMaximum number of points along the curve reached\n"
        if (debug) cat(msg)
        else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
      }

      if (pntnr > 10) {
        if (pntnr > as.numeric(nopts$maxpoints) ||
            (y[1] < (as.numeric(popts$xmin) - as.numeric(nopts$iszero))) ||
            (y[1] > (as.numeric(popts$xmax) + as.numeric(nopts$iszero)))) {
          msg <- "Computation halted:\nMinimum or maximum of x-axis domain reached\n"
          if (debug) cat(msg)
          else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
          break
        }

        if ((curvetype == "EQ") && (popts$ycol == 1)) {
          if (all(y[(2:(rhsdim+1))] < (as.numeric(popts$ymin) - as.numeric(nopts$iszero)))) {
            msg <- "Computation halted:\nMinimum of y-axis domain reached for all y-axis variables\n"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
            break
          }
          if (all(y[(2:(rhsdim+1))] > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Computation halted:\nMaximum of y-axis domain reached for all y-axis variables\n"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
            break
          }
        } else {
          if ((y[popts$ycol] < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
              (y[popts$ycol] > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
            break
          }
          if ((popts$y2col > 1) &&
              ((y[popts$y2col] < (as.numeric(popts$y2min) - as.numeric(nopts$iszero))) ||
               (y[popts$y2col] > (as.numeric(popts$y2max) + as.numeric(nopts$iszero))))) {
            msg <- "Computation halted:\nMaximum or maximum of y-axis domain reached for 2nd y-axis variable\n"
            if (debug) cat(msg)
            else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
            break
          }
        }

        if (any(y[(1:(rhsdim+1))] < - as.numeric(nopts$iszero))) {
          msg <- "Computation halted:\nOne of the variables has become negative\n"
          if (debug) cat(msg)
          else if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
          break
        }
      }
    } else {                                      # Solution not found
      if (pntnr == 1) {
        msg <- "No convergence at initial point\n"
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        return(NULL)
      }
      if (stepscalefac > 4096) {
        msg <- "Step size too small\n"
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        return(NULL)
      }
      stepscalefac <- stepscalefac*2
      guess <- yold + dyscaled/stepscalefac
    }

    # Take a breath
    Sys.sleep(as.numeric(nopts$computedelay))
  }
  return(list("points" = allsols, "eigvals" = alleigs, "tangent" = alltvs,
              "special.points" = specialsols, "special.eigvals" = specialeigs, "special.tangent" = specialtvs, "special.tags" = specialtags))
}
