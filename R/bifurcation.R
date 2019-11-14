#' Phaseplane analysis of a system of ODEs
#'
#' \code{bifurcation}
#'
#'
#'   bifurcation(model, state, parms, resume = TRUE, ...)
#'
#'
#' @param   model  (function, required)
#' \preformatted{}
#'               An R-function that computes the values of the derivatives
#'               in the ODE system (the model definition) at time t.
#'               The model must be defined as: model <- function(t, state, parms),
#'               where t is the current time point in the integration, state is
#'               the current value of the variables in the ODE #' system and
#'               parms is a vector or list of parameters.
#'               The return value of func should be a list, whose first and single
#'               element is a vector containing the derivatives of y with respect
#'               to time. The derivatives must be specified in the same order as
#'               the state variables state. The vector state and parms should both
#'               have name attributes for all their elements
#'
#' @param   state  (numeric vector, required)
#' \preformatted{}
#'               The initial (state) values for the ODE system. This vector should
#'               have name attributes for all its elements
#'
#' @param   parms  (numeric vector, required)
#' \preformatted{}
#'               The values of the parameters in the ODE system. This vector should
#'               have name attributes for all its elements
#'
#' @param   resume  (boolean, optional)
#' \preformatted{}
#'               If TRUE the program will try to load the curves computed during
#'               the last session from global variable 'deBifCurves' and try to
#'               restore the numerical and plot settings by importing them from
#'               the global variable 'deBifSettings'.
#'               The program saves the curves computed during a session and the
#'               numerical and plot settings of this last session in these global
#'               variables 'deBifCurves' and 'deBifSettings'.
#'
#' @param   ...  (optional arguments)
#' \preformatted{}
#'               Additional arguments that can be included at the command line to tweak
#'               graphical default values used by the application.
#'               Valid arguments are:
#' \preformatted{}
#'               \code{lwd}:         Line width (default 3)
#' \preformatted{}
#'               \code{cex}:         Base font size (default 1.2)
#' \preformatted{}
#'               \code{tcl.len}:     Length of axes ticks (default 0.03)
#' \preformatted{}
#'               \code{bifsym}:      Symbol used to mark a bifurcation point
#'                                     in an equilibrium curve (default: 8)
#' \preformatted{}
#'               \code{biflblpos}:   Location of label of a bifurcation point. Values
#'                                   of 1, 2, 3 and 4, respectively, indicate positions
#'                                   below, to the left of, above and to the right of
#'                                   the symbol marking the bifurcation point (default: 3)
#' \preformatted{}
#'               \code{unstablelty}: Line style of curve section representing unstable
#'                                   equilibrium points (default: 3 (refers to dotted lines))
#'
#' @return None.
#'
#' @examples
#' \dontrun{
#' # The initial state of the system has to be specified as a named vector of state values.
#' state <- c(R=1, N=0.01)
#'
#' # Parameters has to be specified as a named vector of parameters.
#' parms <- c(r=1, K=1, a=1, c=1, delta=0.5)
#'
#' # The model has to be specified as a function that returns
#' # the derivatives as a list.
#' model <- function(t, state, parms) {
#'   with(as.list(c(state,parms)), {
#'
#'     dR <- r*R*(1 - R/K) - a*R*N
#'     dN <- c*a*R*N - delta*N
#'
#'     return(list(c(dR, dN)))
#'   })
#' }
#'
#' bifurcation(model, state, parms)
#' }
#' @importFrom graphics contour legend lines par plot points text title axis mtext persp axTicks segments
#' @importFrom grDevices trans3d
#' @import deSolve rootSolve shiny shinydashboard shinydashboardPlus
#' @importFrom shinyjs useShinyjs click removeClass html
#' @importFrom stats setNames
#' @importFrom grDevices dev.off png
#' @export
bifurcation <- function(model, state, parms, resume = TRUE, ...) {

  if (interactive()) {
    # Get the names of the state variables and the parameters
    statenames <- names(state)
    parmsnames <- names(parms)

    # Initialize numerical options
    initnopts <- list(odemethod = "lsoda", tmax = 1000, tstep = 0.1,
                      args_run = unique(names(c(formals(deSolve::ode), formals(deSolve::lsoda)))),
                      methods_run = as.character(formals(deSolve::ode)$method),
                      rtol = 1e-7, atol = 1e-9, ctol = 1e-8, jacdif = 1.0E-6, maxiter = 100,
                      maxpoints = 500, iszero = 1.0E-5, stepsize = 0.01, minstepsize = 1.0E-5, replotfreq = 10,
                      ninterval = 10, glorder = 4, lcampl = 1.0E-6
    )

    # Initialize options for plotting etc.
    initpopts <- vector(mode = "list", 3)
    for (i in 1:length(initpopts)) {
      initpopts[[i]] <- list(xcol = 1, xmin = 0, xmax = 1, logx = 0, xlab = "",
                             ycol = 1, ymin = 0, ymax = 1, logy = 0, ylab = "",
                             y2col = 1, y2min = 0, y2max = 1, logy2 = 0, y2lab = "None", plot3d = 0,
                             lwd = 3, bifsym = 8, unstablelty = 3, tcl.len = 0.03, theta = -35,
                             cex = 1.2, cex.lab = 1.25, cex.axis = 1, cex.legend = 1, cex.sym = 1, biflblpos = 3,
                             colors = c("red","blue","darkgreen","darkorange","darkmagenta",
                                        "gold","darkorchid","aquamarine","deeppink","gray",seq(2,991)))
    }
    initpopts[[1]]$xmax <- initnopts$tmax
    initpopts[[3]]$ycol <- 2

    # Read options from the environment
    if (resume && exists("deBifSettings", envir = .GlobalEnv)) {
      inlist    <- get("deBifSettings", envir = .GlobalEnv)
      initnopts <- checkNumSettings(initnopts, inlist)
      initpopts <- checkPlotSettings(initpopts, inlist, state, parms)
    }

    # Read options from the command line
    adjustableopts <- c("lwd", "cex", "tcl.len", "bifsym", "biflblpos", "unstablelty")
    dots <- list(...)
    if (!is.null(dots)) {
      useropts <- dots[names(dots) %in% adjustableopts]
      if (!is.null(useropts)) {
        for (j in 1:length(useropts)) {
          for (i in 1:length(initpopts)) initpopts[[i]][names(useropts)[j]] <- useropts[j]
        }
      }
    }

    initpopts[[1]]$xlab  <- (c("Time", statenames))[initpopts[[1]]$xcol]
    initpopts[[1]]$ylab  <- ifelse(initpopts[[1]]$ycol  == 1, "State variables", statenames[initpopts[[1]]$ycol -1])
    initpopts[[1]]$y2lab <- ifelse(initpopts[[1]]$y2col == 1, "None",            statenames[initpopts[[1]]$y2col-1])
    initpopts[[2]]$xlab  <- parmsnames[initpopts[[2]]$xcol]
    initpopts[[2]]$ylab  <- ifelse(initpopts[[2]]$ycol  == 1, "State variables", statenames[initpopts[[2]]$ycol -1])
    initpopts[[2]]$y2lab <- ifelse(initpopts[[2]]$y2col == 1, "None",            statenames[initpopts[[2]]$y2col-1])
    initpopts[[3]]$xlab  <- parmsnames[initpopts[[3]]$xcol]
    initpopts[[3]]$ylab  <- parmsnames[initpopts[[3]]$ycol]

    # Read the curves from the environment
    initCurves <- list(Orbits = list(), BifurcationCurves = list(), BifurcationBounds = list(), TotalCurves = 0)
    if (resume && exists("deBifCurves", envir = .GlobalEnv)) {
      inlist     <- get("deBifCurves", envir = .GlobalEnv)
      initCurves <- checkInputCurves(NULL, inlist, statenames, parmsnames)
    }

    ui <- buildUI(state, parms, initpopts, initnopts)

    ############################################# BEGIN SERVER FUNCTION ####################################################
    server <- function(input, output, session) {

      # Hide the button to collapse the left sidebar
      shinyjs::runjs("document.getElementsByClassName('sidebar-toggle')[0].style.visibility = 'hidden';")

      # Create the variable curveList as a reactive value, such that the plots will be updated when curveList changes
      curveList <- reactiveValues()
      curveList$Orbits <- initCurves$Orbits
      curveList$BifurcationCurves <- initCurves$BifurcationCurves
      curveList$BifurcationBounds <- initCurves$BifurcationBounds
      curveList$TotalCurves <- initCurves$TotalCurves
      curveListNames <- c('Orbits', 'BifurcationCurves', 'BifurcationBounds', 'TotalCurves')

      updateSpecialPointsList(session, initCurves, 0)

      # Create the variable plotopts as a reactive value, such that the plots will be updated when plotopts changes
      plotopts <- reactiveValues()
      plotopts$Orbits <- initpopts[[1]]
      plotopts$BifurcationCurves <- initpopts[[2]]
      plotopts$BifurcationBounds <- initpopts[[3]]

      # Create the variable numopts as a reactive value, even though computations will only start after a button push
      numopts <- reactiveValues()
      lapply((1:length(initnopts)), function(i) {numopts[[(names(initnopts)[[i]])]] <- initnopts[[i]]})

      # Create the variable consoleLog as a reactive value, such that the console will be updated when text is added to it
      consoleLog <- reactiveVal()

      busyComputing <- reactiveVal(0)
      updatePlot <- reactiveVal(0)
      curveDirection <- reactiveVal(0)
      changeCurveMenu <- reactiveVal(0)

      # React to changes in curveList or plotopts, by replotting the current plot
      observe({
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        if (as.numeric(updatePlot()) == 0) return(NULL)

        output[[paste0("plot", curtab)]] <- renderPlot({
          if (curtab == 1) biforbitplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else if (curtab == 2) bif1parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else bif2parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
        },
        height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
        width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
        updatePlot(0)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      output$saveplot <- downloadHandler(
        filename = function() {
          tempfile(pattern = "Rplot", tmpdir = '', fileext = ".png")
        },
        content = function(file) {
          curtab <- as.numeric(isolate(input$plottab))
          curtabname <- curveListNames[curtab]
          png(file,
              height = 0.75*session$clientData[[paste0("output_plot", curtab, "_width")]],
              width = 0.99*session$clientData[[paste0("output_plot", curtab, "_width")]])
          if (curtab == 1) biforbitplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else if (curtab == 2) bif1parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          else bif2parplot(session, curveList[[curtabname]], plotopts[[curtabname]])
          dev.off()
        },
        contentType = "image/png")

      # React to changes in curveList by updating the special point selection menu
      observe({
        if (as.numeric(isolate(busyComputing())) == 1) return(NULL)

        if (as.numeric(changeCurveMenu()) == 0) return(NULL)
        else if (as.numeric(changeCurveMenu()) == 1)
          updateSpecialPointsList(session, reactiveValuesToList(curveList), as.numeric(isolate(input$selectpoint)))
        else
          updateSpecialPointsList(session, reactiveValuesToList(curveList), 0)

        # Updating the save and delete curve menu
        curtabname <- curveListNames[as.numeric(isolate(input$plottab))]
        updateCurveMenu(session, curveList[[curtabname]])
        changeCurveMenu(0)
      })

      # React to changes in consoleLog, by rendering the new text
      observe({
        shinyjs::html(id = "console", html = HTML(gsub("\n", "<br>", consoleLog())))
        session$sendCustomMessage(type = "scrollCallback", 1)
      })

      observe({
        pointid <- as.numeric(input$selectpoint)
        if (pointid > 0) {
          isolate({
            curtab <- as.numeric(isolate(input$plottab))
            curtabname <- curveListNames[curtab]
            clist <- reactiveValuesToList(curveList)
            ind1 <- round(pointid/1000000)
            ind2 <- round((pointid-ind1*1000000)/1000)
            ind3 <- round(pointid-ind1*1000000-ind2*1000)
            cln1 <- curveListNames[ind1]
            inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
            if ((curtabname == 'BifurcationBounds') && (inittype %in% c("BP", "HP", "LP"))) {
              updateSelectInput(session, "curvetype3", selected = inittype)
            } else if (curtabname == 'BifurcationCurves') {
              if (inittype == "HP") updateSelectInput(session, "curvetype2", selected = "LC")
              else updateSelectInput(session, "curvetype2", selected = "EQ")
            }
          })
        }
      })

      observe({
        if (as.numeric(curveDirection()) == 0) return(NULL)

        isolate({
          # Close the rightSidebar
          shinyjs::removeClass(selector = "body.skin-blue.sidebar-mini", class = "control-sidebar-open")
          # Collapse the State variables and Parameters stacks
          shinyjs::removeClass(selector = "li.treeview", class = "active")
          shinyjs::hide(selector = "ul.menu-open");
          updateSelectInput(session, "deletecurve", selected = 0)

          curtab <- as.numeric(input$plottab)
          curtabname <- curveListNames[curtab]
          pointid <- as.numeric(input$selectpoint)
          if (curtabname == 'BifurcationCurves') curvetype <- input$curvetype2
          else if (curtabname == 'BifurcationBounds') curvetype <- input$curvetype3
          clist <- reactiveValuesToList(curveList)
          popts <- reactiveValuesToList(plotopts)
          numopts$stepsize <- as.numeric(curveDirection())*abs(numopts$stepsize)

          # Get the starting point
          if (pointid > 0) {
            ind1 <- round(pointid/1000000)
            ind2 <- round((pointid-ind1*1000000)/1000)
            ind3 <- round(pointid-ind1*1000000-ind2*1000)
            cln1 <- curveListNames[ind1]
            ii <- ifelse((ind1 == 3), 2, 1) # 2 parameter bifurcation points have 2 columns before state, otherwise only 1

            initstate <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (ii + (1:length(state)))])
            initparms <- as.numeric(clist[[cln1]][[ind2]]$parameters)
            inittype <- clist[[cln1]][[ind2]]$special.tags[ind3, "Type"]
            if (inittype != "TS")
              initparms[as.numeric(clist[[cln1]][[ind2]]$bifpars)] <- as.numeric(clist[[cln1]][[ind2]]$special.points[ind3, (1:ii)])
            inittanvec <- clist[[cln1]][[ind2]]$tangent[ind3,]
          } else {
            initstate <- state
            initparms <- parms
            for (i in statenames) initstate[i] <- input[[paste0(i, "_", curtab)]]
            for (i in parmsnames) initparms[i] <- input[[paste0(i, "_", curtab)]]
            inittype <- "US"
            inittanvec <- NULL
          }
          names(initstate) <- names(state)
          names(initparms) <- names(parms)

          newlist <- initCurveContinuation(session, model, initstate, initparms, inittanvec, curtabname, clist,
                                           curvetype, inittype, popts[[curtabname]], numopts, as.numeric(input$reportlevel))
          if (!is.null(newlist)) {
            lapply((1:length(newlist)), function(i) {curveList[[(curveListNames[[i]])]] <- newlist[[(curveListNames[[i]])]]})
            updatePlot(1)
            busyComputing(1)
            updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
            shinyjs::show("pausebtn")
            shinyjs::show("stopbtn")
          } else {
            curveDirection(0)
          }

          # Update the console log
          consoleLog(session$userData$alltext)
        })
      })

      observe({
        if ((as.numeric(busyComputing()) != 1) || is.null(session$userData$curveData)) return(NULL)

        isolate({
          curveData <- session$userData$curveData
          nsol <- tryCatch(nextCurvePoints(isolate(round(as.numeric(numopts$replotfreq))), session$userData$curveData,
                                           plotopts[[curveData$tabname]], numopts, session = session),
                           warning = function(e) {
                             msg <- gsub(".*:", "Warning in nextCurvePoints:", e)
                             if (!is.null(session)) updateConsoleLog(session, msg)
                             else cat(msg)
                             session$userData$curveData <- NULL
                             return(NULL)
                           },
                           error = function(e) {
                             msg <- gsub(".*:", "Error in nextCurvePoints:", e)
                             if (!is.null(session)) updateConsoleLog(session, msg)
                             else cat(msg)
                             session$userData$curveData <- NULL
                             return(NULL)
                           })

          if (!is.null(nsol) && (length(nsol) > 0) && !is.null(nsol$points)) {

            newcurve <- curveList[[curveData$tabname]][[curveData$newcurvenr]]
            newcurve$points <- rbind(newcurve$points, nsol$points)
            newcurve$eigvals <- rbind(newcurve$eigvals, nsol$eigvals)
            newcurve$tangent <- rbind(newcurve$tangent, nsol$tangent)
            newcurve$special.points <- rbind(newcurve$special.points, nsol$special.points)
            newcurve$special.tags <- rbind(newcurve$special.tags, nsol$special.tags)
            newcurve$special.eigvals <- rbind(newcurve$special.eigvals, nsol$special.eigvals)
            newcurve$special.tangent <- rbind(newcurve$special.tangent, nsol$special.tangent)

            curveList[[curveData$tabname]][[curveData$newcurvenr]] <- newcurve
          }
        })

        # Update the console log
        consoleLog(session$userData$alltext)

        # Invalidate this for later if computation has not ended
        if (!is.null(session$userData$curveData)) {
          updatePlot(1)
          invalidateLater(10, session)
        } else {
          changeCurveMenu(1)
          updatePlot(1)
          curveDirection(0)
          busyComputing(0)
          shinyjs::hide("pausebtn")
          shinyjs::hide("stopbtn")
        }
      })

      # observeEvent is non-reactive, it only reacts to invalidation of the specified event
      observeEvent(input$pausebtn, {
        busycomp <- as.numeric(isolate(busyComputing()))
        if (busycomp == 0) return(NULL)
        else if (busycomp == 1) {
          updateActionButton(session, "pausebtn", label = "Continue", icon = icon("forward"))
          busyComputing(-1)
          changeCurveMenu(1)
          updatePlot(1)
        }
        else {
          updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
          busyComputing(1)
        }
      })

      observeEvent(input$stopbtn, {
        if (as.numeric(isolate(busyComputing())) == 0) return(NULL)
        updateConsoleLog(session, "Computation interrupted by the user\n")

        curveData <- session$userData$curveData

        # Add the final points as special point
        newcurve <- curveList[[curveData$tabname]][[curveData$newcurvenr]]
        if (!is.null(newcurve) && !is.null(newcurve$points)) {
          if (curveData$curvetype == "LC") {
            statedim <- length(newcurve$initstate)
            freeparsdim <- length(newcurve$bifpars)

            vals <- lapply((1:statedim), function(i) {
              indxrange <- statedim*(1:(numopts$ninterval*numopts$glorder))
              yname <- names(newcurve$points[1, freeparsdim+i])
              y <- newcurve$points[nrow(newcurve$points), freeparsdim+i+indxrange]
              return(paste0("Min.", yname, "=", round(min(y), 3), ", Max.", yname, "=", round(max(y), 3)))
            })
            endPnt <- c("Type" = curveData$curvetype,
                        "Description" = paste0(names(newcurve$points[1, 1]), "=",
                                               round(newcurve$points[nrow(newcurve$points), 1], 3), ", ",
                                               names(newcurve$points[1, ncol(newcurve$points)]), "=",
                                               round(newcurve$points[nrow(newcurve$points), ncol(newcurve$points)], 3), ", ",
                                               paste(unlist(vals), collapse = ', ')))
          } else {
            endPnt <- c("Type" = curveData$curvetype,
                        "Description" = paste(unlist(lapply(1:length(newcurve$points[1,]),
                                                            function(i) {
                                                              paste0(names(newcurve$points[1,i]), "=",
                                                                     round(newcurve$points[nrow(newcurve$points), i], 3))
                                                            })),
                                              collapse=', '))
          }
          updateConsoleLog(session, paste("Ended in", endPnt["Description"], "\n", sep=" "))
          endPnt["Description"] <- paste0(sprintf("%04d: ", (curveData$pntnr-1)), endPnt["Description"])

          newcurve$special.points <- rbind(newcurve$special.points, newcurve$points[nrow(newcurve$points),])
          newcurve$special.eigvals <- rbind(newcurve$special.eigvals, newcurve$special.eigvals[nrow(newcurve$special.eigvals),])
          newcurve$special.tangent <- rbind(newcurve$special.tangent, newcurve$special.tangent[nrow(newcurve$special.tangent),])
          newcurve$special.tags <- rbind(newcurve$special.tags, c(endPnt))

          curveList[[curveData$tabname]][[curveData$newcurvenr]] <- newcurve
        }

        session$userData$curveData <- NULL
        changeCurveMenu(1)
        updatePlot(1)
        curveDirection(0)
        busyComputing(0)
        updateActionButton(session, "pausebtn", label = "Pause", icon = icon("pause-circle"))
        shinyjs::hide("pausebtn")
        shinyjs::hide("stopbtn")

        # Update the console log
        consoleLog(session$userData$alltext)
      }, ignoreInit = TRUE)

      observeEvent(input$plottab, {
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        clist <- reactiveValuesToList(curveList)
        popts <- reactiveValuesToList(plotopts)
        plotopts[[curtabname]] <- updatePlotOptionEntries(session, curtab, popts[[curtabname]], statenames, parmsnames)
        updatePlot(1)
        updateCurveMenu(session, curveList[[curtabname]])
      })

      observeEvent(input$selectpoint, {
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        clist <- reactiveValuesToList(curveList)
        updateSelectedPoint(session, curtab, clist, as.numeric(input$selectpoint), statenames, parmsnames)
      })

      observeEvent(c(input$plotoptsapply, input$numoptsapply), {
        # Close the rightSidebar
        # shinyjs::removeClass(selector = "body.skin-blue.sidebar-mini", class = "control-sidebar-open")

        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        popts <- reactiveValuesToList(plotopts)
        plotopts[[curtabname]] <- processPlotOptionsApply(session, input, curtab, popts[[curtabname]], statenames, parmsnames)
        processNumOptionsApply(session, input, curtab, numopts)
        updatePlot(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      observeEvent(input$computebtn, {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        # Close the rightSidebar
        shinyjs::removeClass(selector = "body.skin-blue.sidebar-mini", class = "control-sidebar-open")
        # Collapse the State variables and Parameters stacks
        shinyjs::removeClass(selector = "li.treeview", class = "active")
        shinyjs::hide(selector = "ul.menu-open");
        updateSelectInput(session, "deletecurve", selected = 0)

        clist <- reactiveValuesToList(curveList)

        initstate <- state
        initparms <- parms
        for (i in statenames) initstate[i] <- input[[paste0(i, "_1")]]
        for (i in parmsnames) initparms[i] <- input[[paste0(i, "_1")]]

        newlist <- processComputeButton(session, model, initstate, initparms, clist, as.numeric(input$selectpoint), numopts)
        if (!is.null(newlist)) lapply((1:length(newlist)),
                                      function(i) {curveList[[(curveListNames[[i]])]] <- newlist[[(curveListNames[[i]])]]})
        changeCurveMenu(1)

        # Update the console log
        consoleLog(session$userData$alltext)
      }, ignoreInit = TRUE)

      observeEvent(c(input$computefwrd2,input$computefwrd3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curveDirection(1)
      }, ignoreInit = TRUE)

      observeEvent(c(input$computebwrd2,input$computebwrd3), {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curveDirection(-1)
      }, ignoreInit = TRUE)

      observeEvent(input$deletebtn, {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        curveList[[curtabname]] <- processDeleteCurve(session, curtab, curveList[[curtabname]], as.numeric(input$deletecurve))
        changeCurveMenu(-1)
      })

      observeEvent(input$savebtn, {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        curtabname <- curveListNames[curtab]
        processSaveCurve(curtab, curveList[[curtabname]], as.numeric(input$savecurve), make.names(input$curvename, unique = TRUE))
      })

      observeEvent(input$appendbtn, {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, input$loadcurve, statenames, parmsnames, replace = FALSE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}

        updateCurveMenu(session, curveList[[curveListNames[curtab]]])
        updateSpecialPointsList(session, reactiveValuesToList(curveList), 0)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      observeEvent(input$replacebtn, {
        if (as.numeric(isolate(busyComputing())) != 0) return(NULL)
        curtab <- as.numeric(input$plottab)
        clist <- reactiveValuesToList(curveList)
        newlist <- processLoadCurve(session, clist, input$loadcurve, statenames, parmsnames, replace = TRUE)
        if (!is.null(newlist)) for (x in curveListNames) {curveList[[x]] <- newlist[[x]]}

        updateCurveMenu(session, curveList[[curveListNames[curtab]]])
        updateSpecialPointsList(session, reactiveValuesToList(curveList), 0)

        # Update the console log
        consoleLog(session$userData$alltext)
      })

      onStop(fun = function() {
        isolate({
          cat("Saving curves and programs settings")
          # Save the current curve list
          if (exists("deBifCurves", envir = .GlobalEnv)) {
            rm("deBifCurves", envir = .GlobalEnv)
          }
          # assign("deBifCurves", list(Orbits = curveList$Orbits, BifurcationCurves = curveList$BifurcationCurves,
          #                            BifurcationBounds = curveList$BifurcationBounds, TotalCurves = curveList$TotalCurves), envir = .GlobalEnv)
          # global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
          (function(key, val, pos) assign(key,val, envir=as.environment(pos)))("deBifCurves", list(Orbits = curveList$Orbits, BifurcationCurves = curveList$BifurcationCurves,
                                                                                                   BifurcationBounds = curveList$BifurcationBounds, TotalCurves = curveList$TotalCurves), 1L)
          # Save the plot and numerical settings in the global environment
          if (exists("deBifSettings", envir = .GlobalEnv)) {
            rm("deBifSettings", envir = .GlobalEnv)
          }
          # assign("deBifSettings", list(plotopts = reactiveValuesToList(plotopts), numopts = reactiveValuesToList(numopts)), envir = .GlobalEnv)
          # global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
          (function(key, val, pos) assign(key,val, envir=as.environment(pos)))("deBifSettings", list(plotopts = reactiveValuesToList(plotopts), numopts = reactiveValuesToList(numopts)), 1L)
        })
      })
    }
    ############################################## END SERVER FUNCTION #####################################################

    # output[["console"]] <- renderText({
    #   cnames <- names(cdata)
    #
    #   allvalues <- lapply(cnames, function(name) {
    #     paste(name, cdata[[name]], sep = " = ")
    #   })
    #   paste(allvalues, collapse = "\n")
    # })

    shinyApp(ui = ui, server = server)
  }
}

# global env set hack (function(key, val, pos) assign(key,val, envir=as.environment(pos)))(myKey, myVal, 1L) `
