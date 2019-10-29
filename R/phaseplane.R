#' Phaseplane analysis of a system of ODEs
#'
#' \code{phaseplane}
#'
#'
#'   phaseplane(model, state, parms, ...)
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
phaseplane <- function(model, state, parms, ...) {
  debifOpts <- list(plotopts = c(lwd = 2, pch = 20, cex = 1.2, cex.lab = 1.25, cex.axis = 1, cex.legend = 1),
                    colors = c("red","blue","darkgreen","darkorange","darkmagenta","gold","darkorchid","aquamarine","deeppink","gray",seq(2,991)),
                    #colors = seq(2,101),  # Use standard R colors
                    args_plot = names(formals(graphics::plot.default)),
                    args_run = unique(names(c(formals(run),formals(deSolve::ode),formals(deSolve::lsoda)))),
                    methods_run = as.character(formals(deSolve::ode)$method),
                    plotdefaults = graphics::par(no.readonly = TRUE),
                    plotmar = c(5,5,4,4))

  adjustableopts <- c("lwd", "pch", "cex")
  dots <- list(...)
  if (!is.null(dots)) {
    useropts <- dots[names(dots) %in% adjustableopts]

    if (!is.null(useropts)) {
      for (j in 1:length(useropts)) {
        debifOpts$plotopts[names(useropts)[j]] <- useropts[j]
      }
    }
  }

  if (length(state) == 1)
    choices <- c("Time plot"=0, "Nullclines"=1, "Steady state"=2)
  else choices <- c("Time plot"=0, "Nullclines"=1, "Steady state"=2, "Vector field"=3, "Trajectory"=4, "Portrait"=5)
  myTabs <- lapply(1:length(choices), function(i) {tabPanel(title = names(choices)[i], plotOutput(outputId=paste0("plot", choices[i]), height = "100%"), value = choices[i])})

  ui <- dashboardPagePlus(
    shinyjs::useShinyjs(),
    header = dashboardHeaderPlus(
      title = tagList(
        span(class = "logo-lg", "Phaseplane analysis"),
        icon("compass"), tags$style(".fa-compass {color:#E87722}")), titleWidth = 220,
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "gears"
    ),
    sidebar = dashboardSidebar(
      width = 220,
      sidebarMenu(
        menuItem(
          h4("Initial state"),
          lapply(1:length(state), function(i){numericInput(inputId=names(state[i]),label=h6(names(state[i])),value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
          tabName = "initstate"
        ),
        menuItem(
          h4("Parameters"),
          lapply(1:length(parms), function(i){numericInput(inputId=names(parms[i]),label=h6(names(parms[i])),value=as.numeric(parms[i]),step=0.1*as.numeric(parms[i]),width="95%")}),
          tabName = "parvals"
        ),
        menuItem(
          actionButton("lapply", "Apply", icon("refresh")),
          tabName = "lapplubtn"),
        id = "leftsbmenu"
      )
    ),
    body = dashboardBody(
      do.call(tabsetPanel, c(myTabs, id = "plottab")),
      box(width = NULL, height = NULL, verbatimTextOutput("console"))
    ),
    rightsidebar = rightSidebar(
      background = "dark",
      rightSidebarTabContent(
        id = "plotopttab",
        title = "Plot options",
        icon = "chart-line",
        active = TRUE,
        h4("X-axis"),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="xmin", label="Minimum", value=0),
                    numericInput(inputId="xmax", label="Maximum", value=100)),
        selectInput('xscale', 'Scale type', c('Linear', "Logarithmic")),
        h4("Y-axis"),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="ymin", label="Minimum", value=0),
                    numericInput(inputId="ymax", label="Maximum", value=1)),
        selectInput('yscale', 'Scale type', c('Linear', "Logarithmic")),
        br(),
        actionButton("plotoptsapply", "Apply", icon("refresh"))
      ),
      rightSidebarTabContent(
        id = "numopttab",
        title = "Numerical options",
        icon = "dashboard",
        h4("Time integration"),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="tmax", label="Maximum time", value=100),
                    numericInput(inputId="tstep", label="Time step", value=0.1,step=0.1)),
        selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4")),
        br(),
        actionButton("numoptsapply", "Apply", icon("refresh"))
      )
    )
  )

  # server = function(input, output) { }
  server <- function(input, output, session) {
    cdata <- session$clientData

    lapply(0:5, function(nr){
      output[[paste0("plot", nr)]] <- renderPlot({
        plottype <- input$plottab
        input$lapply
        input$plotoptsapply
        input$numoptsapply

        observeEvent(input$plotoptsapply, {
          # Close the rightSidebar
          shinyjs::removeClass(selector = "body.skin-blue.sidebar-mini", class = "control-sidebar-open")
        })
        observeEvent(input$numoptsapply, {
          # Close the rightSidebar
          shinyjs::removeClass(selector = "body.skin-blue.sidebar-mini", class = "control-sidebar-open")
        })

        isolate({
          for (i in names(state)) state[i] <- input[[i]]
          for (i in names(parms)) parms[i] <- input[[i]]
          tmax <- input$tmax
          tstep <- input$tstep
          xmin <- input$xmin
          xmax <- input$xmax
          ymin <- input$ymin
          ymax <- input$ymax
          logxy <- ""
          if (input$xscale == "Logarithmic") {
            logxy <- paste0(logxy, "x")
            xmin <- max(xmin, 1.0E-10)
          }
          if (input$yscale == "Logarithmic") {
            logxy <- paste0(logxy, "y")
            ymin <- max(ymin, 1.0E-10)
          }
          if (xmax < 1.0001*xmin) {
            cat("Maximum of x-axis not significantly different from its minimum\n")
            xmax <- 1.0001*xmin
          }
          if (ymax < 1.0001*ymin) {
            cat("Maximum of y-axis not significantly different from its minimum\n")
            ymax <- 1.0001*ymin
          }

          msg <- ""
          output[["console"]] <- renderText({msg})

          if (plottype == 0) {
            nsol <- run(tmax=tmax, tstep=tstep, odes=model, state=state, parms=parms, dbopts=debifOpts, method=input$method)
            par(cex = debifOpts$plotopts["cex"], mar = debifOpts$plotmar)
            plot(NULL, type='n', xlim=c(xmin,xmax), ylim=c(ymin,ymax), log=logxy, xlab="Time", ylab="State variables",
                 cex.lab=debifOpts$plotopts["cex.lab"], cex.axis=debifOpts$plotopts["cex.axis"])
            lapply(2:ncol(nsol), function(i) {
              lines(nsol[,1], nsol[,i], col=debifOpts$colors[min(i-1, length(debifOpts$colors))], lwd=debifOpts$plotopts["lwd"], pch=(i-1))
              })
            legend("topright", legend=names(state), col=debifOpts$colors[1:(ncol(nsol)-1)], lty=1, lwd=debifOpts$plotopts["lwd"], cex=debifOpts$plotopts["cex.legend"])

            msg <- paste(unlist(lapply(1:length(state), function(i) {paste(names(state[i]), "=",round(nsol[nrow(nsol), (1+i)],5), sep=" ")})), collapse=', ')
            msg <- paste("Ended in ", msg, "\n", sep=" ")
          }
          else if (plottype == 1) {
            if (length(state) == 1) plane1D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, odes=model, state=state, parms=parms, eps=-.001, dbopts=debifOpts)
            else plane2D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, odes=model, state=state, parms=parms, eps=-.001, dbopts=debifOpts)
          }
          else if (plottype == 2) {
            if (length(state) == 1) {
              plane1D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, odes=model, state=state, parms=parms, eps=-.001, dbopts=debifOpts)
              msg <- allequi(state=state, parms=parms, odes=model, x=1, xmin=xmin, xmax=xmax, log=logxy, grid=10, plot=TRUE)
            }
            else {
              plane2D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, odes=model, state=state, parms=parms, eps=-.001, dbopts=debifOpts)
              msg <- allequi(state=state, parms=parms, odes=model, x=1, xmin=xmin, xmax=xmax, y=2, ymin=ymin, ymax=ymax, log=logxy, grid=10, plot=TRUE)
            }
          }
          else if (plottype == 3) {
            plane2D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, tmax=tmax, tstep=tstep, odes=model, state=state, parms=parms, eps=-.001,
                    vector=TRUE, grid=8, dbopts=debifOpts)
            allequi(state=state, parms=parms, odes=model, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, grid=10, report=FALSE, eigenvalues=FALSE, plot=TRUE)
          }
          else if (plottype == 4) {
            plane2D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, odes=model, state=state, parms=parms, eps=-.001, dbopts=debifOpts)
            allequi(state=state, parms=parms, odes=model, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, grid=10, report=FALSE, eigenvalues=FALSE, plot=TRUE)
            nsol <- run(tmax=tmax, tstep=tstep, odes=model, state=state, parms=parms, dbopts=debifOpts, method=input$method)

            points(nsol[1, 2], nsol[1, 3], pch=debifOpts$plotopts["pch"])
            lines(nsol[, 2], nsol[, 3], lwd=debifOpts$plotopts["lwd"], col="black")

            msg <- paste(unlist(lapply(1:length(state), function(i) {paste(names(state[i]), "=", round(nsol[nrow(nsol), (1+i)], 5), sep=" ")})), collapse=', ')
            msg <- paste("Ended in ", msg, "\n", sep=" ")
          }
          else if (plottype == 5) {
            plane2D(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, tmax=tmax, tstep=tstep, odes=model, state=state, parms=parms, eps=-.001,
                    portrait=TRUE, dbopts=debifOpts, method=input$method)
            allequi(state=state, parms=parms, odes=model, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, log=logxy, grid=10, report=FALSE, eigenvalues=FALSE, plot=TRUE)
          }
          output[["console"]] <- renderText({msg})
        })
      }, height = function() {0.75*session$clientData[[paste0("output_plot", nr, "_width")]]}, width = function() {0.99*session$clientData[[paste0("output_plot", nr, "_width")]]})

    })

    # output[["console"]] <- renderText({
    #   cnames <- names(cdata)
    #
    #   allvalues <- lapply(cnames, function(name) {
    #     paste(name, cdata[[name]], sep = " = ")
    #   })
    #   paste(allvalues, collapse = "\n")
    # })
  }


  shinyApp(ui = ui, server = server)
}

