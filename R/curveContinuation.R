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
      condfun <- list(get(paste0(curvetype, "continuation"), mode = "function"))
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

  curveData <- list()
  curveData$guess <- y
  curveData$yold <- y
  curveData$tanvec <- tanvec
  curveData$testvals <- NULL
  curveData$pntnr <- 1
  curveData$stepscalefac <- 1

  curveData$tabname <- curtabname
  curveData$model <- model
  curveData$fixedpars <- fixedpars
  curveData$condfun <- condfun
  curveData$extSys <- switch(curvetype, "BP" = ExtSystemEQ, "EQ" = ExtSystemEQ, "HP" = ExtSystemEQ, "LP" = ExtSystemEQ,
                             "LC" = ExtSystemLC)
  curveData$analysefun <- analysefun
  curveData$statedim <- statedim
  curveData$freeparsdim <- freeparsdim
  curveData$pointdim <- freeparsdim + statedim
  curveData$inittype <- inittype
  curveData$curvetype <- curvetype
  curveData$eignames <- eignames
  curveData$tvnames <- tvnames
  curveData$reportlevel <- reportlevel
  curveData$newcurvenr <- length((clist[[curtabname]]))+1

  # If present run the initializer
  if (exists(paste0("init", curvetype), mode = "function")) {
    initfun <- get(paste0("init", curvetype), mode = "function")
    names(y) <- varnames
    initres <- initfun(y, fixedpars, curveData, nopts)
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

  curveData$guess <- y
  curveData$yold <- y
  curveData$tanvec <- tanvec

  nsol <- tryCatch(nextCurvePoints(1, curveData, popts, nopts, session = session),
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

    clist[[curveData$tabname]][[curveData$newcurvenr]] <- newcurve
    clist$TotalCurves <- (clist$TotalCurves + 1)
    return(clist)
  }

  return(NULL)
}

nextCurvePoints <- function(maxpoints, curveData, popts, nopts, session = NULL) {

  if (exists("deBifverbose", envir = .GlobalEnv) && (get("deBifverbose", envir = .GlobalEnv)))
    verbose <- TRUE
  else verbose <- FALSE

  cData <- curveData
  curvetype <- cData$curvetype
  pntnr <- cData$pntnr

  allsols <- NULL
  alltvs <- NULL
  alleigs <- NULL
  specialsols <- NULL
  specialtvs <- NULL
  specialeigs <- NULL
  specialtags <- NULL

  while ((pntnr - cData$pntnr) < as.numeric(maxpoints)) {
    res <- tryCatch(stode(cData$guess, time = 0, func = cData$extSys, parms = cData$fixedpars,
                          rtol = nopts$rtol, atol = nopts$atol, ctol = nopts$ctol,
                          maxiter = nopts$maxiter, verbose = verbose, curveData = cData, nopts = nopts),
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

      # Compute the Jacobian w.r.t. to free parameters (the first cData$freeparsdim
      # columns) and the state variables
      jac <- jacobian.full(y=y, func=cData$extSys, parms=cData$fixedpars, pert = nopts$jacdif, curveData = cData, nopts = nopts)

      # Compute the eigenvalues of the restricted Jacobian (exclude the first cData$freeparsdim columns
      # and takink only as many rows as there are state variables)
      eig <- eigen(jac[(1:cData$statedim),((cData$freeparsdim+1):cData$pointdim)])
      # Sort them on decreasing real part
      eigval <- eig$values[order(Re(eig$values), decreasing = TRUE)]
      names(eigval) <- cData$eignames

      # Append the current tangent vector as the last row to the jacobian to
      # preserve direction. See the matcont manual at
      # http://www.matcont.ugent.be/manual.pdf, page 10 & 11
      # Notice that jacobian.full returns a square matrix with NA values on the last row
      jac[nrow(jac),] <- cData$tanvec
      if (rcond(jac) > nopts$atol) {
        tvnew <- solve(jac, c(rep(0, (length(cData$tanvec)-1)), 1))
        tvnorm <- sqrt(sum(tvnew^2))
        tvnew <- tvnew/tvnorm
        names(tvnew) <- cData$tvnames
        cData$tanvec <- tvnew
      }

      ############## Execute the test functions and store the results
      if (!is.null(cData$analysefun)) {
        testvals <- cData$analysefun(state = y, parms = cData$fixedpars, cData, nopts, session = session)
        if (!is.null(testvals) && ("y" %in% names(testvals) > 0) && (length(testvals$y) > 0)) {

          allsols <- rbind(allsols, c(testvals$y[1:cData$pointdim]))
          alltvs <- rbind(alltvs, c(testvals$tanvec[1:cData$pointdim]))
          alleigs <- rbind(alleigs, c(testvals$eigval[1:length(eigval)]))

          specialsols <- rbind(specialsols, c(allsols[nrow(allsols),]))
          specialtvs <- rbind(specialtvs, c(alltvs[nrow(alltvs),]))
          specialeigs <- rbind(specialeigs, c(alleigs[nrow(alleigs),]))

          dscp <- paste(unlist(lapply(1:cData$pointdim, function(i) {paste0(names(y[i]), "=", round(c(testvals$y)[i], 3))})),
                        collapse=', ')
          specialtags <- rbind(specialtags, c("Type" = testvals$biftype, "Description" = paste0(testvals$biftype, ": ", dscp)))

          msg <- switch(testvals$biftype,
                        BP = "Branching", HP = "Hopf bifurcation", LP = "Limit", BT = "Bogdanov-Takens", CP = "Cusp")
          msg <- paste0("Solution ", pntnr, ": ", msg, " point found:\n", dscp, "\n")
          msg <- paste0(msg, "Eigenvalues:\n",
                        paste(unlist(lapply(1:length(eigval),
                                            function(i) {paste0(ifelse(is.complex(testvals$eigval[i]),
                                                                       sprintf("%12.5E + %12.5Ei", Re(testvals$eigval[i]),
                                                                               Im(testvals$eigval[i])),
                                                                       sprintf("%12.5E", testvals$eigval[i])))})),
                              collapse=' '), "\n")

          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)

          if (cData$reportlevel == 2) {
            msg <- paste0(msg, "Test values:\n")
            if ("bpval" %in% names(testvals)) msg <- paste0(msg, "BP:", sprintf("%12.5E", testvals$bpval), "; ")
            if ("hpval" %in% names(testvals)) msg <- paste0(msg, "HP:", sprintf("%12.5E", testvals$hpval), "; ")
            if ("lpval" %in% names(testvals)) msg <- paste0(msg, "LP:", sprintf("%12.5E", testvals$lpval), "; ")
            if ("btval" %in% names(testvals)) msg <- paste0(msg, "BT:", sprintf("%12.5E", testvals$btval), "; ")
            if ("cpval" %in% names(testvals)) msg <- paste0(msg, "CP:", sprintf("%12.5E", testvals$cpval), "; ")
            msg <- paste0(msg, "\n")
            cat(msg)
          }
          else if (cData$reportlevel == 1) {
            yy <- c(as.numeric(y), eigval)
            if (pntnr == 1) names(yy) <- c(names(y), cData$eignames)
            else names(yy) <- NULL
            print(c(yy), print.gap = 2, digits = 6, right = TRUE)
          }

          testvals$y <- NULL
          testvals$tanvec <- NULL
          testvals$eigval <- NULL
          testvals$biftype <- NULL

          pntnr <- pntnr + 1
        }
        cData$testvals <- testvals
      }

      ############## Report the solution point and the eigenvalues
      msg <- paste0("Solution ", pntnr, " found:\n",
                    paste(unlist(lapply(1:cData$pointdim,
                                        function(i) {paste0(names(y[i]), "=", sprintf("%12.5E", y[i]))})), collapse=', '), "\n")
      msg <- paste0(msg, "Eigenvalues:\n",
                    paste(unlist(lapply(1:length(eigval),
                                        function(i) {paste0(ifelse(is.complex(eigval[i]),
                                                                   sprintf("%12.5E + %12.5Ei", Re(eigval[i]), Im(eigval[i])),
                                                                   sprintf("%12.5E", eigval[i])))})), collapse=' '), "\n")
      if (!is.null(session)) shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
      else cat(msg)

      if (cData$reportlevel == 2) {
        specvar <- c("bpval", "hpval", "lpval", "btval", "cpval")
        speclbl <- c("BP", "HP", "LP", "BT", "CP")
        if (any(specvar %in% names(cData$testvals)))
          msg <- paste0(msg, "Test values:\n",
                        paste(unlist(lapply((1:length(specvar)),
                                            function(i) {if (specvar[i] %in% names(cData$testvals))
                                              paste0(speclbl[i], ": ", sprintf("%12.5E", cData$testvals[[specvar[i]]]))})),
                              collapse='; '), "\n")
        cat(msg)
      }
      else if (cData$reportlevel == 1) {
        yy <- c(as.numeric(y), eigval)
        if (pntnr == 1) names(yy) <- c(names(y), cData$eignames)
        else names(yy) <- NULL
        print(c(yy), print.gap = 2, digits = 6, right = TRUE)
      }

      ############## Store the results
      allsols <- rbind(allsols, (c(y)[1:cData$pointdim]))
      alltvs <- rbind(alltvs, (c(tvnew)[1:cData$pointdim]))
      alleigs <- rbind(alleigs, c(eigval))

      ############## Take a new step along the curve
      # Determine the relative change in the components and the index of the largest change
      indx0s <- (1:length(y))[abs(y) < as.numeric(nopts$iszero)]   # Indices of zero elements of y. Ignore there relative change
      dy <- abs(tvnew)/(pmax(abs(y), as.numeric(nopts$iszero)))
      dy[indx0s] <- 0
      dyind <- which.max(dy)                      # Index with maximum relative change

      dy <- as.numeric(nopts$stepsize)*tvnew

      dyscaled <- max(abs(as.numeric(nopts$stepsize)*y[dyind]), as.numeric(nopts$minstepsize))*(dy/abs(dy[dyind]))

      cData$stepscalefac <- max(1, cData$stepscalefac/2)
      cData$guess <- y + dyscaled/cData$stepscalefac
      cData$yold <- y
      tanvec <- tvnew

      ############## Stop the curve if outside the visible plotting region
      pntnr <- pntnr + 1

      if (pntnr > as.numeric(nopts$maxpoints)) {
        msg <- "Computation halted:\nMaximum number of points along the curve reached\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cData <- NULL
        break
      }

      if (pntnr > 10) {
        if (pntnr > as.numeric(nopts$maxpoints) ||
            (as.numeric(y[1]) < (as.numeric(popts$xmin) - as.numeric(nopts$iszero))) ||
            (as.numeric(y[1]) > (as.numeric(popts$xmax) + as.numeric(nopts$iszero)))) {
          msg <- "Computation halted:\nMinimum or maximum of x-axis domain reached\n"
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
          cData <- NULL
          break
        }
        if (curvetype == "EQ") {
          if (popts$ycol == 1) {
            if (all(as.numeric(y[(2:(cData$statedim+1))]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMinimum of y-axis domain reached for all y-axis variables\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cData <- NULL
              break
            }
            if (all(as.numeric(y[(2:(cData$statedim+1))]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMaximum of y-axis domain reached for all y-axis variables\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cData <- NULL
              break
            }
          } else {
            if ((as.numeric(y[popts$ycol]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
                (as.numeric(y[popts$ycol]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
              msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cData <- NULL
              break
            }
            if ((popts$y2col > 1) &&
                ((as.numeric(y[popts$y2col]) < (as.numeric(popts$y2min) - as.numeric(nopts$iszero))) ||
                 (as.numeric(y[popts$y2col]) > (as.numeric(popts$y2max) + as.numeric(nopts$iszero))))) {
              msg <- "Computation halted:\nMaximum or maximum of y-axis domain reached for 2nd y-axis variable\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cData <- NULL
              break
            }
          }
        } else {
          if (curvetype == "HP") {
            if (!any(is.complex(eigval))) {
              msg <- "Computation halted:\nEigenvalues have become real valued\n"
              if (!is.null(session)) updateConsoleLog(session, msg)
              else cat(msg)
              cData <- NULL
              break
            }
          }
          if ((as.numeric(y[2]) < (as.numeric(popts$ymin) - as.numeric(nopts$iszero))) ||
              (as.numeric(y[2]) > (as.numeric(popts$ymax) + as.numeric(nopts$iszero)))) {
            msg <- "Computation halted:\nMinimum or maximum of y-axis domain reached for 1st y-axis variable\n"
            if (!is.null(session)) updateConsoleLog(session, msg)
            else cat(msg)
            cData <- NULL
            break
          }
        }

        if (any(as.numeric(y[(1:(cData$freeparsdim + cData$statedim))]) < 0 - as.numeric(nopts$iszero))) {
          msg <- "Computation halted:\nOne of the variables has become negative\n"
          if (!is.null(session)) updateConsoleLog(session, msg)
          else cat(msg)
          cData <- NULL
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
      if (cData$stepscalefac > 4096) {
        msg <- "Step size too small\n"
        if (!is.null(session)) updateConsoleLog(session, msg)
        else cat(msg)
        cData <- NULL
        break
      }
      cData$stepscalefac <- cData$stepscalefac*2
      cData$guess <- cData$yold + dyscaled/cData$stepscalefac
    }
  }
  if (is.null(cData)) {
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
    cData$pntnr <- pntnr
  }
  session$userData$curveData <- cData

  if (!is.null(allsols)) {
    return(list("points" = allsols, "eigvals" = alleigs, "tangent" = alltvs,
                "special.points" = specialsols, "special.eigvals" = specialeigs,
                "special.tangent" = specialtvs, "special.tags" = specialtags))
  } else {
    return(NULL)
  }
}
