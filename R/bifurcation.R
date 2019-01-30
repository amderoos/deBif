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
      conditionalPanel(
        condition = "input.plottab == 1",
        sidebarMenu(
          h4("Initial values", align = "center"),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),function(i){numericInput(inputId=paste0(names(state[i]), "_1"), label=h6(names(state[i])), value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),function(i){numericInput(inputId=paste0(names(parms[i]), "_1"), label=h6(names(parms[i])),value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          )
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 2",
        sidebarMenu(
          h4("Initial values", align = "center"),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),function(i){numericInput(inputId=paste0(names(state[i]), "_2"), label=h6(names(state[i])), value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),function(i){numericInput(inputId=paste0(names(parms[i]), "_2"), label=h6(names(parms[i])),value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          )
        )
      ),
      conditionalPanel(
        condition = "input.plottab == 3",
        sidebarMenu(
          h4("Initial values", align = "center"),
          menuItem(
            h4("State variables"),
            lapply(1:length(state),function(i){numericInput(inputId=paste0(names(state[i]), "_3"), label=h6(names(state[i])), value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")})
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),function(i){numericInput(inputId=paste0(names(parms[i]), "_3"), label=h6(names(parms[i])),value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")})
          )
        )
      ),
      br(),
      sidebarMenu(
        actionButton("computebtn", "Compute", icon("refresh")),
        br(),
        h4("Curve management", align = "center"),
        menuItem(
          h4("Delete curve"),
          sidebarMenuOutput("curvemenu")
        ),
        id = "leftsbmenu"
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
        conditionalPanel(
          condition = "input.plottab == 1",
          selectInput('xcol1', h4('Variable(s) on X-axis'), c("Time" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$xcol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="xmin1", label="Minimum", value=plotopts[[1]]$xmin),
                      numericInput(inputId="xmax1", label="Maximum", value=plotopts[[1]]$xmax)),
          selectInput('logx1', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
          selectInput('ycol1', h4('Variable(s) on Y-axis'), c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$ycol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="ymin1", label="Minimum", value=plotopts[[1]]$ymin),
                      numericInput(inputId="ymax1", label="Maximum", value=plotopts[[1]]$ymax)),
          selectInput('logy1', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
          conditionalPanel(
            condition = "input.ycol1 > 1",
            selectInput('y2col1', h4("Variable on 2nd Y-axis"), c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$y2col),
            conditionalPanel(
              condition = "input.y2col1 > 1",
              splitLayout(cellWidths = c("50%", "50%"),
                          numericInput(inputId="y2min1", label="Minimum", value=plotopts[[1]]$y2min),
                          numericInput(inputId="y2max1", label="Maximum", value=plotopts[[1]]$y2max)),
              selectInput('logy21', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy2)
            )
          )
        ),
        conditionalPanel(
          condition = "input.plottab == 2",
          selectInput('xcol2', h4('Bifurcation parameter'), c(setNames((1:(length(parms))), names(parms))), selected=plotopts[[2]]$xcol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="xmin2", label="Minimum", value=plotopts[[1]]$xmin),
                      numericInput(inputId="xmax2", label="Maximum", value=plotopts[[1]]$xmax)),
          selectInput('logx2', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
          selectInput('ycol2', h4('Variable(s) on Y-axis'), c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[2]]$ycol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="ymin2", label="Minimum", value=plotopts[[1]]$ymin),
                      numericInput(inputId="ymax2", label="Maximum", value=plotopts[[1]]$ymax)),
          selectInput('logy2', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
          conditionalPanel(
            condition = "input.ycol2 > 1",
            selectInput('y2col2', h4("Variable on 2nd Y-axis"), c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[2]]$y2col),
            conditionalPanel(
              condition = "input.y2col2 > 1",
              splitLayout(cellWidths = c("50%", "50%"),
                          numericInput(inputId="y2min2", label="Minimum", value=plotopts[[2]]$y2min),
                          numericInput(inputId="y2max2", label="Maximum", value=plotopts[[2]]$y2max)),
              selectInput('logy22', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[2]]$logy2)
            )
          )
        ),
        conditionalPanel(
          condition = "input.plottab == 3",
          selectInput('xcol3', h4('1st bifurcation parameter'), c(setNames((1:(length(parms))), names(parms))), selected=plotopts[[3]]$xcol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="xmin3", label="Minimum", value=plotopts[[1]]$xmin),
                      numericInput(inputId="xmax3", label="Maximum", value=plotopts[[1]]$xmax)),
          selectInput('logx3', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
          selectInput('ycol3', h4('2nd bifurcation parameter'), c(setNames((1:(length(parms))), names(parms))), selected=plotopts[[3]]$ycol),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="ymin3", label="Minimum", value=plotopts[[1]]$ymin),
                      numericInput(inputId="ymax3", label="Maximum", value=plotopts[[1]]$ymax)),
          selectInput('logy3', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy)
        ),
        br(),
        actionButton("plotoptsapply", "Apply", icon("refresh"))
      ),
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

    observeEvent(input$computebtn, {
      shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")
      updateSelectInput(session, "deletecurve", selected = 0)

      curtab <- as.numeric(input$plottab)
      popts <- get("plotopts", envir=bifenv)[[1]]
      clist <- get("curveList", envir=bifenv)

      for (i in names(state)) state[i] <- input[[paste0(i, "_", curtab)]]
      for (i in names(parms)) parms[i] <- input[[paste0(i, "_", curtab)]]

      newcurvenr <- length((clist[[curtab]]))+1

      msg <- ""
      serverenv$output[["console"]] <- renderText({msg})

      if (curtab == 1) {
        nsol <- run(tmax=numopts$tmax, tstep=numopts$tstep, odes=model, state=state, parms=parms, dbopts=numopts, method=numopts$odemethod)
        msg <- paste(unlist(lapply(1:length(state), function(i) {paste(names(state[i]), "=",round(nsol[nrow(nsol), (1+i)],5), sep=" ")})), collapse=', ')
        msg <- paste("Ended in ", msg, "\n", sep=" ")
        serverenv$output[["console"]] <- renderText({msg})

        lbl <- paste0("TS: ", paste(unlist(lapply(1:length(state), function(i) {paste(names(state[i]), "=", round(state[i],5), sep="")})), collapse=', '))
        newcurve <- list(label = lbl, initstate = state, parameters = parms, curve = nsol)
        clist$Orbits[[newcurvenr]] <- newcurve
      }

      rm("curveList", envir=bifenv)
      assign("curveList", clist, envir=bifenv)

      if (exists("CurveList", envir = .GlobalEnv)) {
        rm("CurveList", envir = .GlobalEnv)
      }
      assign("CurveList", clist, envir = .GlobalEnv)

      shinyjs::click("plotoptsapply")
    })

    observeEvent(input$deletebtn, {
      updateSelectInput(session, "deletecurve", selected = 0)

      curtab <- as.numeric(input$plottab)
      clist <- get("curveList", envir=bifenv)
      deletenr <- as.numeric(input$deletecurve)
      totalcurves <- as.numeric(length((clist[[curtab]])))

      if ((totalcurves > 0) && (deletenr > 0) && (deletenr < (totalcurves + 1))) {

        clist[[curtab]][[deletenr]] <- NULL

        rm("curveList", envir=bifenv)
        assign("curveList", clist, envir=bifenv)

        if (exists("CurveList", envir = .GlobalEnv)) {
          rm("CurveList", envir = .GlobalEnv)
        }
        assign("CurveList", clist, envir = .GlobalEnv)

        shinyjs::click("plotoptsapply")
      }
    })

    # render the curve menu
    output$curvemenu <- renderMenu({
      input$plotoptsapply
      curtab <- as.numeric(input$plottab)
      clist <- get("curveList", envir=bifenv)

      if (length((clist[[curtab]])) > 0) {
        lbls <- do.call("rbind", lapply(clist[[curtab]], "[[", "label"))
        ids <- (1:length(clist[[curtab]]))
        names(ids) <- lbls[,1]
        sidebarMenu(
          selectInput('deletecurve', 'Select curve', c("None" = 0, ids), selected=0),
          actionButton("deletebtn", "Delete curve", icon("trash-alt"))
        )
      } else {
        sidebarMenu(
          selectInput('deletecurve', 'Select curve', c("None" = 0), selected=0),
          actionButton("deletebtn", "Delete curve", icon("trash-alt"))
        )
      }
    })

    lapply(1:3, function(nr){
      serverenv$output[[paste0("plot", nr)]] <- renderPlot({
        curtab <- as.numeric(input$plottab)
        input$plotoptsapply
        input$curveClicked

        observeEvent(input$plotoptsapply, {
          shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")
        })

        isolate({
          popts <- get("plotopts", envir=bifenv)

          popts[[nr]]$xcol <- as.numeric(input[[paste0("xcol", nr)]])
          popts[[nr]]$logx <- as.numeric(input[[paste0("logx", nr)]])
          popts[[nr]]$xmin <- ifelse(input[[paste0("logx", nr)]] == 1, max(as.numeric(input[[paste0("xmin", nr)]]), 1.0E-10), as.numeric(input[[paste0("xmin", nr)]]))
          popts[[nr]]$xmax <- as.numeric(input[[paste0("xmax", nr)]])
          popts[[nr]]$xlab <- ifelse(nr == 1, (c("Time", names(state)))[popts[[nr]]$xcol], (names(parms))[popts[[nr]]$xcol])
          popts[[nr]]$ycol <- as.numeric(input[[paste0("ycol", nr)]])
          popts[[nr]]$logy <- as.numeric(input[[paste0("logy", nr)]])
          popts[[nr]]$ymin <- ifelse(input[[paste0("logy", nr)]] == 1, max(as.numeric(input[[paste0("ymin", nr)]]), 1.0E-10), as.numeric(input[[paste0("ymin", nr)]]))
          popts[[nr]]$ymax <- as.numeric(input[[paste0("ymax", nr)]])
          popts[[nr]]$ylab <- ifelse(nr < 3, (c("State variables", names(state)))[popts[[nr]]$ycol], (names(parms))[popts[[nr]]$ycol])

          if (popts[[nr]]$xmax < 1.0001*popts[[nr]]$xmin) {
            cat("Maximum of x-axis not significantly different from its minimum\n")
            popts[[nr]]$xmax <- 1.0001*popts[[nr]]$xmin
          }
          if (popts[[nr]]$ymax < 1.0001*popts[[nr]]$ymin) {
            cat("Maximum of y-axis not significantly different from its minimum\n")
            popts[[nr]]$ymax <- 1.0001*popts[[nr]]$ymin
          }

          if ((nr < 3) && (popts[[nr]]$ycol > 1)) {
            popts[[nr]]$y2col <- as.numeric(input[[paste0("y2col", nr)]])
            popts[[nr]]$logy2 <- as.numeric(input[[paste0("logy2", nr)]])
            popts[[nr]]$y2min <- ifelse(input[[paste0("logy2", nr)]] == 1, max(as.numeric(input[[paste0("y2min", nr)]]), 1.0E-10), as.numeric(input[[paste0("y2min", nr)]]))
            popts[[nr]]$y2max <- as.numeric(input[[paste0("y2max", nr)]])
            popts[[nr]]$y2lab <- ifelse(nr < 2, (names(state))[popts[[nr]]$y2col-1], (names(parms))[popts[[nr]]$y2col])

            if (popts[[nr]]$y2max < 1.0001*popts[[nr]]$y2min) {
              cat("Maximum of y-axis not significantly different from its minimum\n")
              popts[[nr]]$y2max <- 1.0001*popts[[nr]]$y2min
            }
          }
        })

        if (curtab == 1) biforbitplot(curveList[[1]], popts[[1]])
        else if (curtab == 2) bif1parplot(curveList[[2]], popts[[2]])
        else bif2parplot(curveList[[3]], popts[[3]])

        rm("plotopts", envir=bifenv)
        assign("plotopts", popts, envir=bifenv)
      }, height = function() {0.75*session$clientData[[paste0("output_plot", nr, "_width")]]}, width = function() {0.99*session$clientData[[paste0("output_plot", nr, "_width")]]})
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

