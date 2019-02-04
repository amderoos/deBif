processComputeButton <- function(input, output, session, model, state, parms, bifenv, updateProgress) {
  curtab <- as.numeric(input$plottab)
  popts <- get("plotopts", envir=bifenv)
  nopts <- get("numopts", envir=bifenv)
  clist <- get("curveList", envir=bifenv)
  pointid <- as.numeric(input$selectpoint)

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
    inittype <- clist[[ind1]][[ind2]]$special.points[[ind3, "Type"]]
    tanvec <- clist[[ind1]][[ind2]]$tangent[ind3,]
    }
    names(initstate) <- names(state)
    names(initparms) <- names(parms)
  } else {
    initstate <- state
    initparms <- parms
    tanvec <- c(1.0, rep(0, length(state)))
  }

  newcurvenr <- length((clist[[curtab]]))+1

  msg <- ""
  output[["console"]] <- renderText({msg})

  if (curtab == 1) {
    nsol <- run(tmax=nopts$tmax, tstep=nopts$tstep, odes=model, state=initstate, parms=initparms, dbopts=nopts, method=nopts$odemethod)

    OS <- c(nsol[1,], "Type" = "OS", "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[1, (1+i)], 3))})), collapse=', '))
    OE <- c(nsol[nrow(nsol),], "Type" = "OE", "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 3))})), collapse=', '))

    output[["console"]] <- renderText({paste("Ended in", OE["Description"], "\n", sep=" ")})
    curvescomputed <- curvescomputed + 1

    lbl <- paste0("TS", sprintf("%02d", curvescomputed),": ", OS["Description"])
    OS["Description"] <- paste0("OS: ", OS["Description"])
    OE["Description"] <- paste0("OE: ", OE["Description"])

    newcurve <- list(label = lbl, type = "TS", initstate = initstate, parameters = initparms, points = nsol, special.points = rbind(OS, OE))

    clist$Orbits[[newcurvenr]] <- newcurve
    clist$TotalCurves <- curvescomputed
  } else if (curtab == 2) {
    if (inittype == "BP") {
      # Something has to be done here to jump to the other branch
      # tanvec <- clist[[ind1]][[ind2]]$tangent[ind3,]
    }

    nsol <- computeCurve(model, initstate, initparms, as.numeric(popts[[curtab]][["xcol"]]), popts[[curtab]], nopts, "EQ", tanvec = tanvec, output = output, updateProgress = updateProgress)
    if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {
      curvescomputed <- curvescomputed + 1

      ES <- c(nsol$points[1,], "Type" = "EQ",
              "Description" = paste(unlist(lapply(1:length(nsol$points[1,]), function(i) {paste0(names(nsol$points[1,i]), "=", round(nsol$points[1, i], 3))})), collapse=', '))
      EE <- c(nsol$points[nrow(nsol$points),], "Type" = "EQ",
              "Description" = paste(unlist(lapply(1:length(nsol$points[1,]), function(i) {paste0(names(nsol$points[1,i]), "=", round(nsol$points[nrow(nsol$points), i], 3))})), collapse=', '))

      lbl <- paste0("EQ", sprintf("%02d", curvescomputed),": ", ES["Description"])
      ES["Description"] <- paste0("ES: ", ES["Description"])
      EE["Description"] <- paste0("EE: ", EE["Description"])

      newcurve <- list(label = lbl, type = "EQ", initstate = initstate, parameters = initparms, points = nsol$points, eig.vals = nsol$eig.vals, tangent = nsol$tangent, special.points = rbind(ES, EE))

      clist$BifurcationCurves[[newcurvenr]] <- newcurve
      clist$TotalCurves <- curvescomputed
    } else {
      output[["console"]] <- renderText({"Computation of curve failed"})
    }
  }

  lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
  ids <- (0:length(clist[[curtab]]))
  names(ids) <- c("None", lbls[,1])[1:length(ids)]
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
