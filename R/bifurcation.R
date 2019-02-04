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
                          lwd = 3, stablesym = 20, unstablesym = 1, unstablelty = "dotted", sizeLegend = 1,
                          font.main = 2, font.sub = 1,
                          cex.main = 2, cex.lab = 1.4, cex.axis = 1.2,
                          colors = c("red","blue","darkgreen","darkorange","darkmagenta","gold","darkorchid","aquamarine","deeppink","gray",seq(2,991)),
                          plotmar = c(5,5,4,4))
  }

  numopts <- list(odemethod = "lsoda", tmax = 100, tstep = 0.1,
                  args_run = unique(names(c(formals(deSolve::ode), formals(deSolve::lsoda)))),
                  methods_run = as.character(formals(deSolve::ode)$method),
                  rtol = 1e-7, atol = 1e-9, ctol = 1e-8, maxiter = 100,
                  maxpoints = 500, iszero = 1.0E-5, stepsize = 0.01, minstepsize = 1.0E-3, computedelay = 0.01
  )

  ui <- buildUI(state, parms, plotopts, numopts)

  # server = function(input, output) { }
  server <- function(input, output, session) {

    serverenv <- environment()

    observeEvent(input$plottab, {
      processTabSwitch(input, serverenv$output, session, state, parms, bifenv)
    })

    observeEvent(input$selectpoint, {
      updateSelectedPoint(input, session, state, parms, bifenv)
    })

    observeEvent(input$computebtn, {
      shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")
      updateSelectInput(session, "deletecurve", selected = 0)

      # Create a Progress object
      progress <- shiny::Progress$new()
      progress$set(message = "", value = 0)
      # Close the progress when this reactive exits (even if there's an error)
      on.exit(progress$close())

      updateProgress <- function(value = NULL, detail = NULL) {
        if (is.null(value)) {
          value <- progress$getValue()
          value <- value + (progress$getMax() - value) / numopts$maxpoints
        }
        progress$set(value = value, detail = detail)
      }

      processComputeButton(input, serverenv$output, session, model, state, parms, bifenv, updateProgress)

      # Update the special point selection menu
      updateSpecialPointsList(session, bifenv)
    })

    observeEvent(input$deletebtn, {
      processDeleteCurve(input, serverenv$output, session, bifenv)
      updateSelectInput(session, "deletecurve", selected = 0)
    })

    observeEvent(input$savebtn, {
      processSaveCurve(input, serverenv$output, session, bifenv)
    })

    observeEvent(c(input$plotoptsapply, input$numoptsapply), {
      shinyjs::removeClass(selector = "aside.control-sidebar", class = "control-sidebar-open")

      processOptionsApply(input, serverenv$output, session, state, parms, bifenv)
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

