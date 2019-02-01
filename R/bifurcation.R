#' Phaseplane analysis of a system of ODEs
#'
#' \code{bifurcation}
#'
#'
#'   bifurcation(model, state, parms)
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
#' @return None.
#'
#' @examples
#' \dontrun{
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
#' # The initial state of the system has to be specified as a named vector of state values.
#' state <- c(R=1, N=0.01)
#'
#' # Parameters has to be specified as a named vector of parameters.
#' parms <- c(r=1, K=1, a=1, c=1, delta=0.5)
#'
#' phaseplane(model, state, parms)
#' }
#' @importFrom graphics contour legend lines par plot points
#' @import deSolve rootSolve shiny shinydashboard shinydashboardPlus
#' @importFrom shinyjs useShinyjs click removeClass
#' @export
bifurcation <- function(model, state, parms) {

  bifenv <- environment()

  curveList <- list()
  curveList[["Orbits"]] <- list()
  curveList[["BifurcationCurves"]] <- list()
  curveList[["BifurcationBounds"]] <- list()
  curveList[["TotalCurves"]] <- 0

  # Options for plotting etc.
  plotopts <- vector(mode = "list", 3)
  for (i in 1:length(plotopts)) {
    plotopts[[i]] <- list(xcol = 1, xmin = 0, xmax = 1, logx = 0, xlab = "",
                          ycol = 1, ymin = 0, ymax = 1, logy = 0, ylab = "",
                          y2col = 1, y2min = 0, y2max = 1, logy2 = 0, y2lab = "",
                          lwd = 3, pch = 20, sizeLegend = 1,
                          font.main = 2, font.sub = 1,
                          cex.main = 2, cex.lab = 1.4, cex.axis = 1.2,
                          colors = c("red","blue","darkgreen","darkorange","darkmagenta","gold","darkorchid","aquamarine","deeppink","gray",seq(2,991)),
                          plotmar = c(5,5,4,4))
  }

  numopts <- list(odemethod = "lsoda", tmax = 100, tstep = 0.1,
                  args_run = unique(names(c(formals(deSolve::ode), formals(deSolve::lsoda)))),
                  methods_run = as.character(formals(deSolve::ode)$method)
  )

  choices <- c("Time series" = 1, "1 parameter bifurcation" = 2, "2 parameter bifurcation" = 3)
  myTabs <- lapply(1:length(choices), function(i) {tabPanel(title = names(choices)[i], plotOutput(outputId=paste0("plot", choices[i]), height = "100%"), value = choices[i])})

  ui <- dashboardPagePlus(
    shinyjs::useShinyjs(),
    ########## Header of window
    header = dashboardHeaderPlus(
      title = tagList(
        span(class = "logo-lg", "Bifurcation analysis"),
        icon("compass"), tags$style(".fa-compass {color:#E87722}")), titleWidth = 220,
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "gears"
    ),
    ########## Left side-bar
    sidebar = dashboardSidebar(
      width = 220,
      # conditionPanels inside sidebarMenu do not really work well, so the entire sidebarMenu should be wrapped inside a conditionaPanel
      conditionalPanel(
        condition = "input.plottab == 1",
        sidebarMenu(
          h4("Initial values", align = "center"),
          selectizeInput('selectpoint', "", c("User specified" = 0), selected=0),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_1"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_1"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          ),
          id = "leftsbmenu"
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 2",
        sidebarMenu(
          h4("Initial values", align = "center"),
          selectizeInput('selectpoint', "", c("User specified" = 0), selected=0),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_2"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_2"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          ),
          id = "leftsbmenu"
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 3",
        sidebarMenu(
          h4("Initial values", align = "center"),
          selectizeInput('selectpoint', "", c("User specified" = 0), selected=0),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),
                   function(i){numericInput(inputId=paste0(names(state[i]), "_3"), label=h6(names(state[i])),
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_3"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          ),
          id = "leftsbmenu"
        )
      ),
      sidebarMenu(
        br(),
        actionButton("computebtn", "Compute", icon("refresh")),
        br(),
        h4("Curve management", align = "center"),
        menuItem(h4("Delete curve"),
                 selectInput('deletecurve', 'Select curve', c("None" = 0), selected=0),
                 actionButton("deletebtn", "Delete curve"),
                 tabName = "deletetab"
        ),
        menuItem(h4("Save curve"),
                 selectInput('savecurve', 'Select curve', c("None" = 0), selected=0),
                 textInput("curvename", "Give a valid R variable name"),
                 actionButton("savebtn", "Save curve"),
                 tabName = "savetab"
        )
      )
    ),
    ########## Main panels
    body = dashboardBody(
      do.call(tabsetPanel, c(myTabs, id = "plottab")),
      box(width = NULL, height = NULL, verbatimTextOutput("console"))
    ),
    ########## Right side-bar
    rightsidebar = rightSidebar(
      background = "dark",
      rightSidebarTabContent(
        id = "plotopttab",
        title = "Plot options",
        icon = "chart-line",
        active = TRUE,
        selectInput('xcol', h4('Variable(s) on X-axis'),
                    c("Time" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$xcol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="xmin", label="Minimum", value=plotopts[[1]]$xmin),
                    numericInput(inputId="xmax", label="Maximum", value=plotopts[[1]]$xmax)),
        selectInput('logx', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
        selectInput('ycol', h4('Variable(s) on Y-axis'),
                    c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$ycol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="ymin", label="Minimum", value=plotopts[[1]]$ymin),
                    numericInput(inputId="ymax", label="Maximum", value=plotopts[[1]]$ymax)),
        selectInput('logy', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
        conditionalPanel(
          condition = "input.plottab != 3 && input.ycol > 1",
          selectInput('y2col', h4("Variable on 2nd Y-axis"),
                      c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$y2col),
          conditionalPanel(
            condition = "input.y2col > 1",
            splitLayout(cellWidths = c("50%", "50%"),
                        numericInput(inputId="y2min", label="Minimum", value=plotopts[[1]]$y2min),
                        numericInput(inputId="y2max", label="Maximum", value=plotopts[[1]]$y2max)),
            selectInput('logy2', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy2)
          )
        ),
        br(),
        actionButton("plotoptsapply", "Apply", icon("refresh"))
      ),
      # sidebarMenu() and menuItem() do not work in the right sidebar. Only conditionpanels can be used
      rightSidebarTabContent(
        id = "numopttab",
        title = "Numerical options",
        icon = "dashboard",
        h4("Time integration"),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="tmax", label="Maximum time", value=numopts$tmax),
                    numericInput(inputId="tstep", label="Time step", value=numopts$tstep,step=0.1)),
        selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4"),selected=numopts$odemethod),
        br(),
        actionButton("numoptsapply", "Apply", icon("refresh"))
      )
    )
  )

  # server = function(input, output) { }
  server <- function(input, output, session) {

    serverenv <- environment()

    observeEvent(input$plottab, {
      curtab <- as.numeric(input$plottab)

      # Update the delete and save menus
      clist <- get("curveList", envir=bifenv)
      lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
      ids <- (0:length(clist[[curtab]]))
      names(ids) <- c("None", lbls[,1])[1:length(ids)]
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
      popts[[curtab]]$y2lab <- ifelse(curtab < 2, (names(state))[popts[[curtab]]$y2col-1], (names(parms))[popts[[curtab]]$y2col])

      serverenv$output[[paste0("plot", curtab)]] <- renderPlot({
        if (curtab == 1) biforbitplot(clist[[1]], popts[[1]])
        else if (curtab == 2) bif1parplot(clist[[2]], popts[[2]])
        else bif2parplot(clist[[3]], popts[[3]])
      },
      height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
      width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
    })

    observeEvent(input$selectpoint, {
      curtab <- as.numeric(input$plottab)
      curveid <- as.numeric(input$selectpoint)
      if (curveid > 0) {
        ind1 <- round(curveid/1000000)
        ind2 <- round((curveid-ind1*1000000)/1000)
        ind3 <- round(curveid-ind1*1000000-ind2*1000)
        ii <- ifelse((ind1 == 3), 2, 1)   # 2 parameter bifurcation points have 2 columns before the state, otherwise 1 only
        lapply(1:length(state), function(i){updateNumericInput(session, paste0(names(state[i]), "_", curtab), value=as.numeric(curveList[[ind1]][[ind2]]$points[[ind3,(i+ii)]]))})
        lapply(1:length(parms), function(i){updateNumericInput(session, paste0(names(parms[i]), "_", curtab), value=as.numeric(curveList[[ind1]][[ind2]]$parameters[[i]]))})
        if (ind1 > 1) {
          cnames <- colnames(CurveList[[ind1]][[ind2]]$points)
          updateNumericInput(session, paste0(names(parms[cnames[1]]), "_", curtab), value=as.numeric(curveList[[ind1]][[ind2]]$points[[ind3,1]]))
          if (ind1 == 3) {
            updateNumericInput(session, paste0(names(parms[cnames[2]]), "_", curtab), value=as.numeric(curveList[[ind1]][[ind2]]$points[[ind3,2]]))
          }
        }
      }
    })

    observeEvent(input$computebtn, {
      shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")
      updateSelectInput(session, "deletecurve", selected = 0)

      curtab <- as.numeric(input$plottab)
      popts <- get("plotopts", envir=bifenv)
      clist <- get("curveList", envir=bifenv)
      curvescomputed <- as.numeric(clist[["TotalCurves"]])

      for (i in names(state)) state[i] <- input[[paste0(i, "_", curtab)]]
      for (i in names(parms)) parms[i] <- input[[paste0(i, "_", curtab)]]

      newcurvenr <- length((clist[[curtab]]))+1

      msg <- ""
      serverenv$output[["console"]] <- renderText({msg})

      if (curtab == 1) {
        nsol <- run(tmax=numopts$tmax, tstep=numopts$tstep, odes=model, state=state, parms=parms, dbopts=numopts, method=numopts$odemethod)

        OS <- c(nsol[1,], "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[1, (1+i)], 2))})), collapse=', '))
        OE <- c(nsol[nrow(nsol),], "Description" = paste(unlist(lapply(1:length(state), function(i) {paste0(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 2))})), collapse=', '))

        serverenv$output[["console"]] <- renderText({paste("Ended in", OE["Description"], "\n", sep=" ")})
        curvescomputed <- curvescomputed + 1

        lbl <- paste0("TS", sprintf("%02d", curvescomputed),": ", OS["Description"])
        OS["Description"] <- paste0("OS: ", OS["Description"])
        OE["Description"] <- paste0("OE: ", OE["Description"])

        newcurve <- list(label = lbl, initstate = state, parameters = parms, curve = nsol, points = rbind(OS, OE))
        clist$Orbits[[newcurvenr]] <- newcurve
        clist$TotalCurves <- curvescomputed

        lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
        ids <- (0:length(clist[[curtab]]))
        names(ids) <- c("None", lbls[,1])[1:length(ids)]
        updateSelectInput(session, "deletecurve", choices=ids, selected=0)
        updateSelectInput(session, "savecurve", choices=ids, selected=0)

        # Update the special point selection menu
        splist <- list()
        for (i in (1:3)) {
          if (length(clist[[i]]) > 0) {
            listlbls <- NULL
            splist <- c(splist, lapply((1:length(clist[[i]])),
                                       function(j){
                                         listlbls <<- c(listlbls, clist[[i]][[j]]$label)
                                         lbls <- unlist(clist[[i]][[j]]$points[,"Description"], use.names=FALSE);
                                         ids <- ((i*1000000)+j*1000)+(1:length(lbls)); names(ids) <- lbls;
                                         return(ids)}))}
          names(splist) <- listlbls
        }
        updateSelectInput(session, "selectpoint", choices=c(list("User specified" = 0), splist), selected=0)
      }

      serverenv$output[[paste0("plot", curtab)]] <- renderPlot({
        if (curtab == 1) biforbitplot(clist[[1]], popts[[1]])
        else if (curtab == 2) bif1parplot(clist[[2]], popts[[2]])
        else bif2parplot(clist[[3]], popts[[3]])
      },
      height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
      width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})

      rm("curveList", envir=bifenv)
      assign("curveList", clist, envir=bifenv)

      if (exists("CurveList", envir = .GlobalEnv)) {
        rm("CurveList", envir = .GlobalEnv)
      }
      assign("CurveList", clist, envir = .GlobalEnv)
    })

    observeEvent(input$deletebtn, {
      updateSelectInput(session, "deletecurve", selected = 0)

      curtab <- as.numeric(input$plottab)
      clist <- get("curveList", envir=bifenv)
      popts <- get("plotopts", envir=bifenv)
      deletenr <- as.numeric(input$deletecurve)
      totalcurves <- as.numeric(length((clist[[curtab]])))

      if ((totalcurves > 0) && (deletenr > 0) && (deletenr < (totalcurves + 1))) {

        clist[[curtab]][[deletenr]] <- NULL

        lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
        ids <- (0:length(clist[[curtab]]))
        names(ids) <- c("None", lbls[,1])[1:length(ids)]
        updateSelectInput(session, "deletecurve", choices=ids, selected=0)
        updateSelectInput(session, "savecurve", choices=ids, selected=0)

        # Update the special point selection menu
        splist <- list()
        for (i in (1:3)) {
          if (length(clist[[i]]) > 0) {
            listlbls <- NULL
            splist <- c(splist, lapply((1:length(clist[[i]])),
                                       function(j){
                                         listlbls <<- c(listlbls, clist[[i]][[j]]$label)
                                         lbls <- unlist(clist[[i]][[j]]$points[,"Description"], use.names=FALSE);
                                         ids <- ((i*1000000)+j*1000)+(1:length(lbls)); names(ids) <- lbls;
                                         return(ids)}))}
          names(splist) <- listlbls
        }
        updateSelectInput(session, "selectpoint", choices=c(list("User specified" = 0), splist), selected=0)

        serverenv$output[[paste0("plot", curtab)]] <- renderPlot({
          if (curtab == 1) biforbitplot(clist[[1]], popts[[1]])
          else if (curtab == 2) bif1parplot(clist[[2]], popts[[2]])
          else bif2parplot(clist[[3]], popts[[3]])
        },
        height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
        width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})

        rm("curveList", envir=bifenv)
        assign("curveList", clist, envir=bifenv)

        if (exists("CurveList", envir = .GlobalEnv)) {
          rm("CurveList", envir = .GlobalEnv)
        }
        assign("CurveList", clist, envir = .GlobalEnv)
      }
    })

    observeEvent(input$savebtn, {
      curtab <- as.numeric(input$plottab)
      clist <- get("curveList", envir=bifenv)
      savenr <- as.numeric(input$savecurve)
      totalcurves <- as.numeric(length((clist[[curtab]])))

      if ((totalcurves > 0) && (savenr > 0) && (savenr < (totalcurves + 1))) {
        varname <- make.names(input$curvename, unique = TRUE)
        if (exists(varname, envir = .GlobalEnv)) {
          rm(varname, envir = .GlobalEnv)
        }
        assign(varname, clist[[curtab]][[savenr]], envir = .GlobalEnv)
      }
    })

    observeEvent(input$plotoptsapply, {
      shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")

      isolate({
        curtab <- as.numeric(input$plottab)
        clist <- get("curveList", envir=bifenv)
        popts <- get("plotopts", envir=bifenv)

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
          popts[[curtab]]$y2lab <- ifelse(curtab < 2, (names(state))[popts[[curtab]]$y2col-1], (names(parms))[popts[[curtab]]$y2col])

          if (popts[[curtab]]$y2max < 1.0001*popts[[curtab]]$y2min) {
            cat("Maximum of y-axis not significantly different from its minimum\n")
            popts[[curtab]]$y2max <- 1.0001*popts[[curtab]]$y2min
          }
        }

        rm("plotopts", envir=bifenv)
        assign("plotopts", popts, envir=bifenv)
      })

      serverenv$output[[paste0("plot", curtab)]] <- renderPlot({
        if (curtab == 1) biforbitplot(clist[[1]], popts[[1]])
        else if (curtab == 2) bif1parplot(clist[[2]], popts[[2]])
        else bif2parplot(clist[[3]], popts[[3]])
      },
      height = function() {0.75*session$clientData[[paste0("output_plot", curtab, "_width")]]},
      width = function() {0.99*session$clientData[[paste0("output_plot", curtab, "_width")]]})
    })
  }

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

