processTabSwitch <- function(input, output, session, state, parms, bifenv) {
  curtab <- as.numeric(input$plottab)

  # Update the save and delete curve menu
  clist <- get("curveList", envir=bifenv)
  lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
  ids <- c((0:length(clist[[curtab]])), -1)
  names(ids) <- c("None", lbls[,1], "All")[1:length(ids)]
  updateSelectInput(session, "deletecurve", choices=ids, selected=0)
  updateSelectInput(session, "savecurve", choices=ids, selected=0)

  # Update the plot options
  popts <- get("plotopts", envir=bifenv)
  if (curtab == 1){
    updateSelectInput(session, "xcol",  label=h4('Variable(s) on X-axis'),
                      choices=c("Time" = 1, setNames((2:(length(state)+1)), names(state))), selected=popts[[curtab]]$xcol)
    updateSelectInput(session, "ycol",  label=h4('Variable(s) on Y-axis'),
                      choices=c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=popts[[curtab]]$ycol)
    updateSelectInput(session, "y2col", label=h4("Variable on 2nd Y-axis"),
                      choices=c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=popts[[curtab]]$y2col)
  } else if (curtab == 2){
    updateSelectInput(session, "xcol",  label=h4('Bifurcation parameter'),
                      choices=c(setNames((1:(length(parms))), names(parms))), selected=popts[[curtab]]$xcol)
    updateSelectInput(session, "ycol",  label=h4('Variable(s) on Y-axis'),
                      choices=c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=popts[[curtab]]$ycol)
    updateSelectInput(session, "y2col", label=h4("Variable on 2nd Y-axis"),
                      choices=c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=popts[[curtab]]$y2col)
  } else {
    updateSelectInput(session, "xcol",  label=h4('1st bifurcation parameter'),
                      choices=c(setNames((1:(length(parms))), names(parms))), selected=popts[[curtab]]$xcol)
    updateSelectInput(session, "ycol",  label=h4('2nd bifurcation parameter'),
                      choices=c(setNames((1:(length(parms))), names(parms))), selected=popts[[curtab]]$ycol)
  }
  updateSelectInput(session,  "logx", selected=popts[[curtab]]$logx)
  updateNumericInput(session, "xmin", value=popts[[curtab]]$xmin)
  updateNumericInput(session, "xmax", value=popts[[curtab]]$xmax)
  popts[[curtab]]$xlab <- ifelse(curtab == 1, (c("Time", names(state)))[popts[[curtab]]$xcol], (names(parms))[popts[[curtab]]$xcol])

  updateSelectInput(session,  "logy", selected=popts[[curtab]]$logy)
  updateNumericInput(session, "ymin", value=popts[[curtab]]$ymin)
  updateNumericInput(session, "ymax", value=popts[[curtab]]$ymax)
  popts[[curtab]]$ylab <- ifelse(curtab < 3, (c("State variables", names(state)))[popts[[curtab]]$ycol],
                                 (names(parms))[popts[[curtab]]$ycol])

  updateSelectInput(session,  "logy2", selected=popts[[curtab]]$logy2)
  updateNumericInput(session, "y2min", value=popts[[curtab]]$y2min)
  updateNumericInput(session, "y2max", value=popts[[curtab]]$y2max)
  popts[[curtab]]$y2lab <- (c("None", names(state)))[popts[[curtab]]$y2col]

  output[[paste0("plot", curtab)]] <- renderPlot({
    if (curtab == 1) biforbitplot(clist[[1]], output, session, popts[[1]])
    else if (curtab == 2) bif1parplot(clist[[2]], output, session, popts[[2]])
    else bif2parplot(clist[[3]], output, session, popts[[3]])
  },
  height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
  width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
}

processDeleteCurve <- function(input, output, session, bifenv) {
  curtab <- as.numeric(input$plottab)
  clist <- get("curveList", envir=bifenv)
  popts <- get("plotopts", envir=bifenv)
  deletenr <- as.numeric(input$deletecurve)
  totalcurves <- as.numeric(length((clist[[curtab]])))

  if ((totalcurves > 0) && ((deletenr > 0) || (deletenr == -1)) && (deletenr < (totalcurves + 1))) {

    if (deletenr == -1) clist[[curtab]] <- list()
    else clist[[curtab]][[deletenr]] <- NULL

    # Update the save and delete curve menu
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
      if (curtab == 1) biforbitplot(clist[[1]], output, session, popts[[1]])
      else if (curtab == 2) bif1parplot(clist[[2]], output, session, popts[[2]])
      else bif2parplot(clist[[3]], output, session, popts[[3]])
    },
    height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
    width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
  }
}

processSaveCurve <- function(input, output, session, bifenv) {
  curtab <- as.numeric(input$plottab)
  clist <- get("curveList", envir=bifenv)
  savenr <- as.numeric(input$savecurve)
  totalcurves <- as.numeric(length((clist[[curtab]])))

  if ((totalcurves > 0) && ((savenr > 0) || (savenr == -1)) && (savenr < (totalcurves + 1))) {
    varname <- make.names(input$curvename, unique = TRUE)
    if (exists(varname, envir = .GlobalEnv)) {
      rm(varname, envir = .GlobalEnv)
    }
    if (savenr == -1) assign(varname, clist[[curtab]], envir = .GlobalEnv)
    else assign(varname, clist[[curtab]][[savenr]], envir = .GlobalEnv)
  }
}

processLoadCurve <- function(input, output, session, bifenv, replace = FALSE) {
  curtab <- as.numeric(input$plottab)
  clist <- get("curveList", envir=bifenv)
  popts <- get("plotopts", envir=bifenv)
  snames <- names(get("state", envir=bifenv))
  pnames <- names(get("parms", envir=bifenv))

  varname <- input$loadcurve
  if (exists(varname, envir = .GlobalEnv)) {
    inlist <- get(varname, envir = .GlobalEnv)
  } else return(NULL)

  # If inlist is a list containing elements 'Orbits', 'BifurcationCurves',
  # 'BifurcationBounds', and 'TotalCurves' assume that it is a correct list
  # generated by an earlier call to bifurcation()
  if (is.list(inlist) &&
      all(c("Orbits", "BifurcationCurves", "BifurcationBounds", "TotalCurves") %in% names(inlist))) {
    clist <- inlist
  } else {
    # If inlist is a list containing elements that each contain elements
    # 'points' and 'special.points' with columns that have appropriate names
    # assume that it is a correct list generated by an earlier call to
    # bifurcation()
    if (curtab == 1) {
      if (is.list(inlist) && !is.null(inlist[[1]]) &&
          all(unlist(lapply(inlist, function(x) {((c("points", "special.points") %in% names(x)) &&
                                               all((colnames(x$points))[1] == "Time") &&
                                               all((colnames(x$points))[2:ncol(x$points)] %in% snames))})))) {

        if (replace) clist[[curtab]] <- inlist
        else clist[[curtab]] <- c(clist[[curtab]], inlist)
      } else {
        msg <- paste0("Curves not loaded: variable '", varname, "' does not have the correct structure for this plot\n")
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        inlist <- NULL
      }
    } else if (curtab == 2) {
      if (is.list(inlist) && !is.null(inlist[[1]]) &&
          all(unlist(lapply(inlist, function(x) {((c("points", "eigvals", "tangent", "special.points", "special.tags",
                                                  "special.eigvals", "special.tangent") %in% names(x)) &&
                                               all((colnames(x$points))[1] %in% pnames) &&
                                               all((colnames(x$points))[2:ncol(x$points)] %in% snames))})))) {
        if (replace) clist[[curtab]] <- inlist
        else clist[[curtab]] <- c(clist[[curtab]], inlist)
      } else {
        msg <- paste0("Curves not loaded: variable '", varname, "' does not have the correct structure for this plot\n")
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        inlist <- NULL
      }
    } else if (curtab == 3) {
      if (is.list(inlist) && !is.null(inlist[[1]]) &&
          all(unlist(lapply(inlist, function(x) {((c("points", "eigvals", "tangent", "special.points", "special.tags",
                                                  "special.eigvals", "special.tangent") %in% names(x)) &&
                                               all((colnames(x$points))[1:2] %in% pnames) &&
                                               all((colnames(x$points))[3:ncol(x$points)] %in% snames))})))) {
        if (replace) clist[[curtab]] <- inlist
        else clist[[curtab]] <- c(clist[[curtab]], inlist)
      } else {
        msg <- paste0("Curves not loaded: variable '", varname, "' does not have the correct structure for this plot\n")
        if (!is.null(output)) output[["console"]] <- updateConsoleText(session, msg)
        else cat(msg)
        inlist <- NULL
      }
    } else inlist <- NULL
  }

  if (!is.null(inlist)) {
    # Update the save and delete curve menu
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

    # If called via the load curve button, replot the current plot
    output[[paste0("plot", curtab)]] <- renderPlot({
      if (curtab == 1) biforbitplot(clist[[1]], output, session, popts[[1]])
      else if (curtab == 2) bif1parplot(clist[[2]], output, session, popts[[2]])
      else bif2parplot(clist[[3]], output, session, popts[[3]])
    },
    height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
    width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
  }
}

updateSelectedPoint <- function(input, session, state, parms, bifenv) {
  curtab <- as.numeric(input$plottab)
  clist <- get("curveList", envir=bifenv)

  pointid <- as.numeric(input$selectpoint)
  if (pointid > 0) {
    ind1 <- round(pointid/1000000)
    ind2 <- round((pointid-ind1*1000000)/1000)
    ind3 <- round(pointid-ind1*1000000-ind2*1000)
    ii <- ifelse((ind1 == 3), 2, 1)   # 2 parameter bifurcation points have 2 columns before the state, otherwise 1 only
    lapply(1:length(state), function(i){updateNumericInput(session, paste0(names(state[i]), "_", curtab), value=as.numeric(clist[[ind1]][[ind2]]$special.points[[ind3,(i+ii)]]))})
    lapply(1:length(parms), function(i){updateNumericInput(session, paste0(names(parms[i]), "_", curtab), value=as.numeric(clist[[ind1]][[ind2]]$parameters[[i]]))})
    if (ind1 > 1) {
      cnames <- colnames(clist[[ind1]][[ind2]]$special.points)
      updateNumericInput(session, paste0(names(parms[cnames[1]]), "_", curtab), value=as.numeric(clist[[ind1]][[ind2]]$special.points[[ind3,1]]))
      if (ind1 == 3) {
        updateNumericInput(session, paste0(names(parms[cnames[2]]), "_", curtab), value=as.numeric(clist[[ind1]][[ind2]]$special.points[[ind3,2]]))
      }
    }
  }
}

updateSpecialPointsList <- function(session, selected, bifenv) {
  clist <- get("curveList", envir=bifenv)

  splist <- list()
  listlbls <- NULL
  for (i in (1:3)) {
    if (length(clist[[i]]) > 0) {
      listlbls <- c(listlbls, unlist(lapply((1:length(clist[[i]])), function(j){return(clist[[i]][[j]]$label)})))
      splist <- c(splist, lapply((1:length(clist[[i]])),
                                 function(j){
                                   lbls <- unlist(clist[[i]][[j]]$special.tags[,"Description"], use.names=FALSE);
                                   ids <- ((i*1000000)+j*1000)+(1:length(lbls)); names(ids) <- lbls;
                                   return(ids)}))
      }
  }
  if (length(splist) > 0) {
    names(splist) <- listlbls
    updateSelectInput(session, "selectpoint", choices=c(list("User specified" = 0), splist), selected=selected)
  } else {
    updateSelectInput(session, "selectpoint", choices=c(list("User specified" = 0)), selected=0)
  }
}

processOptionsApply <- function(input, output, session, state, parms, bifenv) {
  curtab <- as.numeric(input$plottab)
  popts <- get("plotopts", envir=bifenv)
  nopts <- get("numopts", envir=bifenv)
  clist <- get("curveList", envir=bifenv)

  popts[[curtab]]$xcol <- as.numeric(input[["xcol"]])
  popts[[curtab]]$logx <- as.numeric(input[["logx"]])
  popts[[curtab]]$xmin <- ifelse(input[["logx"]] == 1,
                                 max(as.numeric(input[["xmin"]]), 1.0E-10), as.numeric(input[["xmin"]]))
  popts[[curtab]]$xmax <- as.numeric(input[["xmax"]])
  popts[[curtab]]$xlab <- ifelse(curtab == 1, (c("Time", names(state)))[popts[[curtab]]$xcol], (names(parms))[popts[[curtab]]$xcol])
  popts[[curtab]]$ycol <- as.numeric(input[["ycol"]])
  popts[[curtab]]$logy <- as.numeric(input[["logy"]])
  popts[[curtab]]$ymin <- ifelse(input[["logy"]] == 1,
                                 max(as.numeric(input[["ymin"]]), 1.0E-10), as.numeric(input[["ymin"]]))
  popts[[curtab]]$ymax <- as.numeric(input[["ymax"]])
  popts[[curtab]]$ylab <- ifelse(curtab < 3, (c("State variables", names(state)))[popts[[curtab]]$ycol],
                                 (names(parms))[popts[[curtab]]$ycol])

  if (popts[[curtab]]$xmax < 1.0001*popts[[curtab]]$xmin) {
    cat("Maximum of x-axis not significantly different from its minimum\n")
    popts[[curtab]]$xmax <- 1.0001*popts[[curtab]]$xmin
  }
  if (popts[[curtab]]$ymax < 1.0001*popts[[curtab]]$ymin) {
    cat("Maximum of y-axis not significantly different from its minimum\n")
    popts[[curtab]]$ymax <- 1.0001*popts[[curtab]]$ymin
  }

  if ((curtab < 3) && (popts[[curtab]]$ycol > 1)) {
    popts[[curtab]]$y2col <- as.numeric(input[["y2col"]])
    popts[[curtab]]$logy2 <- as.numeric(input[["logy2"]])
    popts[[curtab]]$y2min <- ifelse(input[["logy2"]] == 1,
                                    max(as.numeric(input[["y2min"]]), 1.0E-10), as.numeric(input[["y2min"]]))
    popts[[curtab]]$y2max <- as.numeric(input[["y2max"]])
    popts[[curtab]]$y2lab <- (c("None", names(state)))[popts[[curtab]]$y2col]

    if (popts[[curtab]]$y2max < 1.0001*popts[[curtab]]$y2min) {
      cat("Maximum of y-axis not significantly different from its minimum\n")
      popts[[curtab]]$y2max <- 1.0001*popts[[curtab]]$y2min
    }
  }

  text2numeric <- function(oldval, newval){
    newval2 <- gsub("[^0-9.E+-]*", "", newval)
    if (!is.na(suppressWarnings(as.numeric(newval2)))) return (as.numeric(newval2))
    else return(as.numeric(oldval))
  }

  if (curtab == 1) {
    nopts$tmax <- as.numeric(input[["tmax"]])
    nopts$tstep <- as.numeric(input[["tstep"]])
    nopts$odemethod <- input[["method"]]
  }
  else {
    nopts$rtol <- max(text2numeric(nopts$rtol, input[["rtol"]]), 1.0E-10)
    updateTextInput(session, "rtol", value=sprintf("%.1E", nopts$rtol))
    nopts$atol <- max(text2numeric(nopts$atol, input[["atol"]]), 1.0E-10)
    updateTextInput(session, "atol", value=sprintf("%.1E", nopts$atol))
    nopts$iszero <- max(text2numeric(nopts$iszero, input[["iszero"]]), 1.0E-10)
    updateTextInput(session, "iszero", value=sprintf("%.1E", nopts$iszero))
    nopts$jacdif <- max(text2numeric(nopts$jacdif, input[["jacdif"]]), 1.0E-10)
    updateTextInput(session, "jacdif", value=sprintf("%.1E", nopts$jacdif))

    nopts$stepsize <- as.numeric(input[["stepsize"]])
    nopts$minstepsize <- max(as.numeric(input[["minstepsize"]]), 1.0E-10)
    nopts$maxiter <- max(as.numeric(input[["maxiter"]]), 1)
    nopts$maxpoints <- max(as.numeric(input[["maxpoints"]]), 1)
    nopts$computedelay <- max(as.numeric(input[["computedelay"]]), 0.0)
  }

  rm("numopts", envir=bifenv)
  assign("numopts", nopts, envir=bifenv)

  rm("plotopts", envir=bifenv)
  assign("plotopts", popts, envir=bifenv)

  output[[paste0("plot", curtab)]] <- renderPlot({
    if (curtab == 1) biforbitplot(clist[[1]], output, session, popts[[1]])
    else if (curtab == 2) bif1parplot(clist[[2]], output, session, popts[[2]])
    else bif2parplot(clist[[3]], output, session, popts[[3]])
  },
  height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
  width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
}

updateConsoleText <- function(session, addtext) {
  if ("alltext" %in% names(session$userData)) session$userData$alltext <- paste(session$userData$alltext, addtext, sep = "")
  else session$userData$alltext <- addtext
  return(renderText({session$userData$alltext}))
}
