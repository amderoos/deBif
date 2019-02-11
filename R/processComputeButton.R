processComputeButton <- function(input, output, session, model, state, parms, bifenv, direction) {
  if (exists("deBifdebug", envir = .GlobalEnv)) debug <- get("deBifdebug", envir = .GlobalEnv)
  else debug <- FALSE

  curtab <- as.numeric(input$plottab)
  popts <- get("plotopts", envir=bifenv)
  nopts <- get("numopts", envir=bifenv)
  clist <- get("curveList", envir=bifenv)
  report2console <- get("report2console", envir=bifenv)
  pointid <- as.numeric(input$selectpoint)

  nopts$stepsize <- direction*abs(nopts$stepsize)
  curvescomputed <- as.numeric(clist[["TotalCurves"]])

  for (i in names(state)) state[i] <- input[[paste0(i, "_", curtab)]]
  for (i in names(parms)) parms[i] <- input[[paste0(i, "_", curtab)]]

  if (pointid > 0) {
    ind1 <- round(pointid/1000000)
    ind2 <- round((pointid-ind1*1000000)/1000)
    ind3 <- round(pointid-ind1*1000000-ind2*1000)
    ii <- ifelse((ind1 == 3), 2, 1)   # 2 parameter bifurcation points have 2 columns before the state, otherwise 1 only

    initstate <- as.numeric(clist[[ind1]][[ind2]]$special.points[ind3, (ii + (1:length(state)))])
    initparms <- as.numeric(clist[[ind1]][[ind2]]$parameters)
    inittype <- clist[[ind1]][[ind2]]$special.tags[ind3, "Type"]
    if (inittype != "TS")
      initparms[as.numeric(clist[[ind1]][[ind2]]$bifpars)] <- as.numeric(clist[[ind1]][[ind2]]$special.points[ind3, (1:ii)])
    tanvec <- clist[[ind1]][[ind2]]$tangent[ind3,]

    names(initstate) <- names(state)
    names(initparms) <- names(parms)
  } else {
    initstate <- state
    initparms <- parms
    inittype <- "US"
    tanvec <- NULL
  }

  newcurvenr <- length((clist[[curtab]]))+1

  if (curtab == 1) {
    nsol <- run(tmax=nopts$tmax, tstep=nopts$tstep, odes=model, state=initstate, parms=initparms, dbopts=nopts, method=nopts$odemethod)

    startPnt <- c("Type" = "TS",
                  "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[1, (1+i)], 3))})), collapse=', '))
    endPnt <- c("Type" = "TS",
                "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 3))})), collapse=', '))

    output[["console"]] <- updateConsoleText(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
    curvescomputed <- curvescomputed + 1

    lbl <- paste0("TS", sprintf("%02d", curvescomputed),": ", startPnt["Description"])
    startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])
    endPnt["Description"] <- paste0(sprintf("%04d: ", nrow(nsol)), endPnt["Description"])

    newcurve <- list(label = lbl, type = "TS", initstate = initstate, parameters = initparms, points = nsol,
                     special.points = rbind(c(nsol[1,]), c(nsol[nrow(nsol),])), special.tags = rbind(startPnt, endPnt))

    clist$Orbits[[newcurvenr]] <- newcurve
    clist$TotalCurves <- curvescomputed
  } else {
    if (curtab == 2) {
      curvetype = "EQ"
      freepars = c(as.numeric(popts[[curtab]][["xcol"]]))
    } else {
      curvetype = input$curvetype
      freepars = c(as.numeric(popts[[curtab]][["xcol"]]), as.numeric(popts[[curtab]][["ycol"]]))
    }

    nsol <- tryCatch(computeCurve(model, initstate, initparms, freepars, popts[[curtab]], nopts, inittype, curvetype,
                                  tanvec = tanvec, report2console, session = session, output = output),
                    warning = function(e) {
                      msg <- gsub(".*:", "Warning in computeCurve:", e)
                      if (debug) cat(msg)
                      else shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
                      return(NULL)
                    },
                    error = function(e) {
                      msg <- gsub(".*:", "Error in computeCurve:", e)
                      if (debug) cat(msg)
                      else shinyjs::html(id = "progress", html = HTML(gsub("\n", "<br>", msg)))
                      return(NULL)
                    })

    if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {
      curvescomputed <- curvescomputed + 1

      startPnt <- c("Type" = curvetype,
              "Description" = paste(unlist(lapply(1:length(nsol$points[1,]), function(i) {paste0(names(nsol$points[1,i]), "=", round(nsol$points[1, i], 3))})), collapse=', '))
      endPnt <- c("Type" = curvetype,
              "Description" = paste(unlist(lapply(1:length(nsol$points[1,]), function(i) {paste0(names(nsol$points[1,i]), "=", round(nsol$points[nrow(nsol$points), i], 3))})), collapse=', '))

      output[["console"]] <- updateConsoleText(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))

      lbl <- paste0(curvetype, sprintf("%02d", curvescomputed),": ", startPnt["Description"])
      startPnt["Description"] <- paste0(sprintf("%04d: ", 1), startPnt["Description"])
      endPnt["Description"] <- paste0(sprintf("%04d: ", nrow(nsol$points)), endPnt["Description"])

      if ("special.points" %in% names(nsol)) {
        newcurve <- list(label = lbl, type = curvetype, initstate = initstate, parameters = initparms, bifpars = freepars, points = nsol$points,
                         eigvals = nsol$eigvals, tangent = nsol$tangent,
                         special.points = rbind(c(nsol$points[1,]), nsol$special.points, c(nsol$points[nrow(nsol$points),])),
                         special.tags = rbind(startPnt, nsol$special.tags, endPnt),
                         special.eigvals = rbind(nsol$eigvals[1,], nsol$special.eigvals, nsol$eigvals[nrow(nsol$points),]),
                         special.tangent = rbind(nsol$tangent[1,], nsol$special.tangent, nsol$tangent[nrow(nsol$points),]))
      } else {
        newcurve <- list(label = lbl, type = curvetype, initstate = initstate, parameters = initparms, bifpars = freepars, points = nsol$points,
                         eigvals = nsol$eigvals, tangent = nsol$tangent, special.points = rbind(c(nsol$points[1,]), c(nsol$points[nrow(nsol$points),])),
                         special.tags = rbind(startPnt, endPnt))
      }

      clist[[curtab]][[newcurvenr]] <- newcurve
      clist$TotalCurves <- curvescomputed
    }
  }

  lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
  ids <- c((0:length(clist[[curtab]])), -1)
  names(ids) <- c("None", lbls[,1], "All")[1:length(ids)]
  updateSelectInput(session, "deletecurve", choices=ids, selected=0)
  updateSelectInput(session, "savecurve", choices=ids, selected=0)

  rm("curveList", envir=bifenv)
  assign("curveList", clist, envir=bifenv)

  if (exists("CurveList", envir = .GlobalEnv)) {
    rm("CurveList", envir = .GlobalEnv)
  }
  assign("CurveList", clist, envir = .GlobalEnv)

  output[[paste0("plot", curtab)]] <- renderPlot({
    if (curtab == 1) biforbitplot(clist[[1]], popts[[1]])
    else if (curtab == 2) bif1parplot(clist[[2]], popts[[2]])
    else bif2parplot(clist[[3]], popts[[3]])
  },
  height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
  width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
}
