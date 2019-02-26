initCurveContinuation <- function(session, model, initstate, initparms, tanvec, curtabname, clist,
                                  curvetype, inittype, popts, nopts, reportlevel) {

  if (exists("deBifverbose", envir = .GlobalEnv)) verbose <- get("deBifverbose", envir = .GlobalEnv)
  else verbose <- FALSE

  if (!(curvetype %in% c("EQ", "BP", "HP", "LP"))) {
    msg <- paste0("Computation aborted:\nContinuation for curve type ", curvetype, " not implemented\n")
    if (!is.null(session)) updateConsoleLog(session, msg)
    else cat(msg)
    return(NULL)
  }

  if (curtabname == 'BifurcationCurves') {
    freepars = c(as.numeric(popts[["xcol"]]))
  } else {
    freepars = c(as.numeric(popts[["xcol"]]), as.numeric(popts[["ycol"]]))
  }

  freeparsdim <- length(freepars)
  statedim <- length(initstate)

  y <- c(initparms[freepars], initstate)
  fixedpars <- initparms[!(1:length(initparms)) %in% freepars]

  varnames <- names(y)
  eignames <- unlist(lapply((1:statedim), function(i){paste0("Eigenvalue", i)}))
  tvnames <- unlist(lapply((1:length(y)), function(i){paste0("d", varnames[i])}))

  if (curvetype == "EQ") {
    condfun <- NULL
  } else {
    if (exists(paste0(curvetype, "continuation"), mode = "function"))
      condfun <- c(get(paste0(curvetype, "continuation"), mode = "function"))
    else {
      msg <- paste0("Computation aborted:\nAdditional condition function for curve type ", curvetype, " not found\n")
      if (!is.null(session)) updateConsoleLog(session, msg)
      else cat(msg)
      return(NULL)
    }
  }

  # Test functions will only be implemented for EQ curves
  if (exists(paste0("analyse", curvetype), mode = "function")) {
    analysefun <- get(paste0("analyse", curvetype), mode = "function")
    testvals <- NULL
  }
  else analysefun <- NULL

  # If present run the initializer
  if (exists(paste0("init", curvetype), mode = "function")) {
    initfun <- get(paste0("init", curvetype), mode = "function")
    names(y) <- varnames
    initres <- initfun(y, fixedpars, model, tanvec, statedim, freeparsdim, inittype, nopts)
    if (!is.null(initres)) {
      if (!is.null(initres$y)) {
        y <- c(as.numeric(initres$y))
        names(y) <- varnames
      }
      if (!is.null(initres$tanvec)) tanvec <- c(as.numeric(initres$tanvec))
    }
  }

  if (is.null(tanvec) || (length(tanvec) != length(y)))
    tanvec <- c(1.0, rep(0, (length(y)-1)))
  else if (abs(tanvec[1]) > nopts$iszero) tanvec <- sign(tanvec[1])*tanvec

  computeSpecs <- list()
  computeSpecs$guess <- y
  computeSpecs$yold <- y
  computeSpecs$tanvec <- tanvec
  computeSpecs$testvals <- NULL
  computeSpecs$pntnr <- 1
  computeSpecs$stepscalefac <- 1

  computeSpecs$tabname <- curtabname
  computeSpecs$model <- model
  computeSpecs$fixedpars <- fixedpars
  computeSpecs$condfun <- condfun
  computeSpecs$analysefun <- analysefun
  computeSpecs$statedim <- statedim
  computeSpecs$freeparsdim <- freeparsdim
  computeSpecs$curvetype <- curvetype
  computeSpecs$eignames <- eignames
  computeSpecs$tvnames <- tvnames
  computeSpecs$reportlevel <- reportlevel
  computeSpecs$newcurvenr <- length((clist[[curtabname]]))+1

  nsol <- tryCatch(nextCurvePoints(1, computeSpecs, popts, nopts, session = session),
                   warning = function(e) {
                     msg <- gsub(".*:", "Warning in nextCurvePoints:", e)
                     if (!is.null(session)) updateConsoleLog(session, msg)
                     else cat(msg)
                     return(NULL)
                   },
                   error = function(e) {
                     msg <- gsub(".*:", "Error in nextCurvePoints:", e)
                     if (!is.null(session)) updateConsoleLog(session, msg)
                     else cat(msg)
                     return(NULL)
                   })

  if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {
    startPnt <- c("Type" = curvetype,
                  "Description" = paste(unlist(lapply(1:length(nsol$points[1,]),
                                                      function(i) {paste0(names(nsol$points[1,i]), "=",
                                                                          round(nsol$points[1, i], 3))})),
                                        collapse=', '))
    lbl <- paste0(curvetype, sprintf("%02d", (clist$TotalCurves + 1)),": ", startPnt["Description"])
    startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])

    newcurve <- list(label = lbl, type = curvetype, initstate = initstate, parameters = initparms, bifpars = freepars,
                     points = nsol$points, eigvals = nsol$eigvals, tangent = nsol$tangent,
                     special.points = nsol$points, special.eigvals = nsol$eigvals,
                     special.tangent = nsol$tangent, special.tags = rbind(NULL, c(startPnt)))

    clist[[computeSpecs$tabname]][[computeSpecs$newcurvenr]] <- newcurve
    clist$TotalCurves <- (clist$TotalCurves + 1)
    return(clist)
  }

  return(NULL)
}

nextCurvePoints <- function(maxpoints, computeSpecs, popts, nopts, session = NULL) {

  if (exists("deBifverbose", envir = .GlobalEnv) && (get("deBifverbose", envir = .GlobalEnv)))
    verbose <- TRUE
  else verbose <- FALSE

  cmpSpecs <- computeSpecs
  guess <- cmpSpecs$guess
  yold <- cmpSpecs$yold
  tanvec <- cmpSpecs$tanvec
  testvals <- cmpSpecs$testvals
  pntnr <- cmpSpecs$pntnr
  stepscalefac <- cmpSpecs$stepscalefac

  model <- cmpSpecs$model
  fixedpars <- cmpSpecs$fixedpars
  condfun <- cmpSpecs$condfun
  statedim <- cmpSpecs$statedim
  freeparsdim <- cmpSpecs$freeparsdim
  pointdim <- freeparsdim + statedim
  eignames <- cmpSpecs$eignames
  tvnames <- cmpSpecs$tvnames
  analysefun <- cmpSpecs$analysefun
  reportlevel <- cmpSpecs$reportlevel
  curvetype <- cmpSpecs$curvetype

  allsols <- NULL
  alltvs <- NULL
  alleigs <- NULL
  specialsols <- NULL
  specialtvs <- NULL
  specialeigs <- NULL
  specialtags <- NULL

  while ((pntnr - cmpSpecs$pntnr) < as.numeric(maxpoints)) {
    res <- tryCatch(stode(guess, time = 0, func = ExtSystem, parms = fixedpars,
                          rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                          maxiter = nopts$maxiter, verbose = verbose, model = model,
                          ystart = guess, tanvec = tanvec, condfun = condfun,
                          statedim = statedim, freeparsdim = freeparsdim, nopts = nopts),
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
    # Solution found
    if (!is.null(res) && !is.null(attr(res, "steady")) && attr(res, "steady")) {
      y <- res$y

      # Compute the Jacobian w.r.t. to free parameters (the first freeparsdim
      # columns) and the state variables
      jac <- jacobian.full(y=y, func=ExtSystem, parms=fixedpars, pert = nopts$jacdif, model = model,
                           condfun = condfun, statedim = statedim, freeparsdim = freeparsdim, nopts = nopts)

      # Compute the eigenvalues of the restricted Jacobian (exclude the first freeparsdim columns
      # and takink only as many rows as there are state variables)
      eig <- eigen(jac[(1:statedim),((freeparsdim+1):pointdim)])
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
      if (!is.null(analysefun)) {
        testvals <- analysefun(state = y, parms = fixedpars, model = model, tanvec = tvnew,
                               statedim = statedim, freeparsdim = freeparsdim,
                               lastvals = testvals, nopts, session = session)
        if (!is.null(testvals) && ("y" %in% names(testvals) > 0) && (length(testvals$y) > 0)) {

          allsols <- rbind(allsols, c(testvals$y[1:pointdim]))
          alltvs <- rbind(alltvs, c(testvals$tanvec[1:pointdim]))
          alleigs <- rbind(alleigs, c(testvals$eigval[1:length(eigval)]))

          specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
          specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
          specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))

          dscp <- paste(unlist(lapply(1:pointdim, function(i) {paste0(names(y[i]), "=", round(c(testvals$y)[i], 3))})),
                        collapse=', ')
          specialtags <- rbind(specialtags, c("Type" = testvals$biftype, "Description" = paste0(testvals$biftype, ": ", dscp)))

          if (testvals$biftype == "BP") msg <- paste("Solution ", pntnr, ": Branching point found:\n", dscp, "\n", sep=" ")
          else if (testvals$biftype == "HP") msg <- paste("Solution ", pntnr, ": Hopf bifurcation point found:\n", dscp, "\n", sep=" ")
          else if (testvals$biftype == "LP") msg <- paste("Solution ", pntnr, ": Limit point found:\n", dscp, "\n", sep=" ")
          else if (testvals$biftype == "BT") msg <- paste("Solution ", pntnr, ": Bogdanov-Takens point found:\n", dscp, "\n", sep=" ")
          else msg <- paste("Cusp point found:\n", dscp, "\n", sep=" ")
          msg <- paste0(msg, "Eigenvalues:\n",
                        paste(unlist(lapply(1:length(eigval),
                                            function(i) {paste0(ifelse(is.complex(testvals$eigval[i]),
                                                                       sprintf("%12.5E + %12.5Ei", Re(testvals$eigval[i]),
                                                                               Im(testvals$eigval[i])),
                                                                       sprintf("%12.5E", testvals$eigval[i])))})),
                              collapse=' '), "\n")

          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)

          if (reportlevel == 2) {
            msg <- paste0(msg, "Test values:\n")
            if ("bpval" %in% names(testvals)) msg <- paste0(msg, "BP:", sprintf("%12.5E", testvals$bpval), "; ")
            if ("hpval" %in% names(testvals)) msg <- paste0(msg, "HP:", sprintf("%12.5E", testvals$hpval), "; ")
            if ("lpval" %in% names(testvals)) msg <- paste0(msg, "LP:", sprintf("%12.5E", testvals$lpval), "; ")
            if ("btval" %in% names(testvals)) msg <- paste0(msg, "BT:", sprintf("%12.5E", testvals$btval), "; ")
            if ("cpval" %in% names(testvals)) msg <- paste0(msg, "CP:", sprintf("%12.5E", testvals$cpval), "; ")
            msg <- paste0(msg, "\n")
            cat(msg)
          }
          else if (reportlevel == 1) {
            yy <- c(as.numeric(y), eigval)
            if (pntnr == 1) names(yy) <- c(names(y), eignames)
            else names(yy) <- NULL
            print(c(yy), print.gap = 2, digits = 6, right = TRUE)
          }

          testvals$y <- NULL
          testvals$tanvec <- NULL
          testvals$eigval <- NULL
          testvals$biftype <- NULL

          pntnr <- pntnr + 1
        }
      }

      ############## Report the solution point and the eigenvalues
      msg <- paste0("Solution ", pntnr, " found:\n",
                    paste(unlist(lapply(1:pointdim,
                                        function(i) {paste0(names(y[i]), "=", sprintf("%12.5E", y[i]))})), collapse=', '), "\n")
      msg <- paste0(msg, "Eigenvalues:\n",
                    paste(unlist(lapply(1:length(eigval),
                                        function(i) {paste0(ifelse(is.complex(eigval[i]),
                                                                   sprintf("%12.5E + %12.5Ei", Re(eigval[i]), Im(eigval[i])),
                                                                   sprintf("%12.5E", eigval[i])))})), collapse=' '), "\n")
      if (!is.null(session)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
      else cat(msg)

      if (reportlevel == 2) {
        specvar <- c("bpval", "hpval", "lpval", "btval", "cpval")
        speclbl <- c("BP", "HP", "LP", "BT", "CP")
        if (any(specvar %in% names(testvals)))
          msg <- paste0(msg, "Test values:\n",
                        paste(unlist(lapply((1:length(specvar)),
                                            function(i) {if (specvar[i] %in% names(testvals))
                                              paste0(speclbl[i], ": ", sprintf("%12.5E", testvals[[specvar[i]]]))})),
                              collapse='; '), "\n")
        cat(msg)
      }
      else if (reportlevel == 1) {
        yy <- c(as.numeric(y), eigval)
        if (pntnr == 1) names(yy) <- c(names(y), eignames)
        else names(yy) <- NULL
        print(c(yy), print.gap = 2, digits = 6, right = TRUE)
      }

      ############## Store the results
      allsols <- rbind(allsols, (c(y)[1:pointdim]))
      alltvs <- rbind(alltvs, (c(tvnew)[1:pointdim]))
      alleigs <- rbind(alleigs, c(eigval))

      ############## Take a new step along the curve
      # Determine the relative change in the components and the index of the largest change
      indx0s <- (1:length(y))[abs(y) < as.numeric(nopts$iszero)]   # Indices of zero elements of y. Ignore there relative change
      dy <- abs(tvnew)/(pmax(abs(y), as.numeric(nopts$iszero)))
      dy[indx0s] <- 0
      dyind <- which.max(dy)                      # Index with maximum relative change

      dy <- as.numeric(nopts$stepsize)*tvnew

      dyscaled <- max(abs(as.numeric(nopts$stepsize)*y[dyind]), as.numeric(nopts$minstepsize))*(dy/abs(dy[dyind]))

      stepscalefac <- max(1, stepscalefac/2)
      guess <- y + dyscaled/stepscalefac
      yold <- y
      tanvec <- tvnew

      ############## Stop the curve if outside the visible plotting region
      pntnr <- pntnr + 1

      if (pntnr > as.numeric(nopts$maxpoints)) {
        msg <- "Computation halted:\nMaximum number of points along the curve reached\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cmpSpecs <- NULL
        break
      }

      if (pntnr > 10) {
        if (pntnr > as.numeric(nopts$maxpoints) ||
            (as.numeric(y[1]) < (as.numeric(popts$xmin) - as.numeric(nopts$iszero))) ||
            (as.numeric(y[1]) > (as.numeric(popts$xmax) + as.numeric(nopts$iszero)))) {
          msg <- "Computation halted:\nMinimum or maximum of x-axis domain reached\n"
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
          cmpSpecs <- NULL
          break
        }
        if (curvetype == "EQ") {
          if (popts$ycol == 1) {
            if (all(as.numeric(y[(2:(statedim+1))]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMinimum of y-axis domain reached for all y-axis variables\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cmpSpecs <- NULL
              break
            }
            if (all(as.numeric(y[(2:(statedim+1))]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMaximum of y-axis domain reached for all y-axis variables\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cmpSpecs <- NULL
              break
            }
          } else {
            if ((as.numeric(y[popts$ycol]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
                (as.numeric(y[popts$ycol]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cmpSpecs <- NULL
              break
            }
            if ((popts$y2col > 1) &&
                ((as.numeric(y[popts$y2col]) < (as.numeric(popts$y2min) - as.numeric(nopts$iszero))) ||
                 (as.numeric(y[popts$y2col]) > (as.numeric(popts$y2max) + as.numeric(nopts$iszero))))) {
              msg <- "Computation halted:\nMaximum or maximum of y-axis domain reached for 2nd y-axis variable\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cmpSpecs <- NULL
              break
            }
          }
        } else {
          if (curvetype == "HP") {
            if (!any(is.complex(eigval))) {
              msg <- "Computation halted:\nEigenvalues have become real valued\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cmpSpecs <- NULL
              break
            }
          }
          if ((as.numeric(y[2]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
              (as.numeric(y[2]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
            if (!is.null(session)) updateConsoleLog(session, msg)
            else cat(msg)
            cmpSpecs <- NULL
            break
          }
        }

        if (any(as.numeric(y[(1:(freeparsdim + statedim))]) < 0 - as.numeric(nopts$iszero))) {
          msg <- "Computation halted:\nOne of the variables has become negative\n"
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
          cmpSpecs <- NULL
          break
        }
      }
    } else {                                      # Solution not found
      if (pntnr == 1) {
        msg <- "No convergence at initial point\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        return(NULL)
      }
      if (stepscalefac > 4096) {
        msg <- "Step size too small\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cmpSpecs <- NULL
        break
      }
      stepscalefac <- stepscalefac*2
      guess <- yold + dyscaled/stepscalefac
    }
  }
  if (is.null(cmpSpecs)) {
    if (!is.null(allsols)) {
      endPnt <- c("Type" = curvetype,
                  "Description" = paste(unlist(lapply(1:length(allsols[1,]),
                                                      function(i) {paste0(names(allsols[1,i]), "=",
                                                                          round(allsols[nrow(allsols), i], 3))})),
                                        collapse=', '))
      updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
      endPnt["Description"] <- paste0(sprintf("%04d: ", pntnr-1), endPnt["Description"])

      specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
      specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
      specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))
      specialtags <- rbind(specialtags, c(endPnt))
    }
  } else {
    cmpSpecs$guess <- guess
    cmpSpecs$yold <- yold
    cmpSpecs$tanvec <- tanvec
    cmpSpecs$testvals <- testvals
    cmpSpecs$pntnr <- pntnr
    cmpSpecs$stepscalefac <- stepscalefac
  }
  session$userData$computeSpecs <- cmpSpecs

  if (!is.null(allsols)) {
    return(list("points" = allsols, "eigvals" = alleigs, "tangent" = alltvs,
                "special.points" = specialsols, "special.eigvals" = specialeigs,
                "special.tangent" = specialtvs, "special.tags" = specialtags))
  } else {
    return(NULL)
  }
}
