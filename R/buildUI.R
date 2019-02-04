buildUI <- function(state, parms, plotopts, numopts) {
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
      sidebarMenu(
        h4("Initial values", align = "center"),
        selectizeInput('selectpoint', "", c("User specified" = 0), selected=0)
      ),
      conditionalPanel(
        condition = "input.plottab == 1",
        sidebarMenu(
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
      box(width = NULL, height = NULL, verbatimTextOutput("console")),
      tags$head(
        tags$style(
          HTML(".shiny-notification {
               font-size: 16px;
               height: 100px;
               width: calc(100% - 250px);
               position:fixed;
               top: calc(100% - 120px);;
               left: calc(230px);;
}
"
          )
        )
      )
    ),
    ########## Right side-bar
    rightsidebar = rightSidebar(
      background = "dark",
      rightSidebarTabContent(
        id = "plotopttab",
        title = h3("Plot options"),
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
        title = h3("Numerical options"),
        icon = "dashboard",
        conditionalPanel(
          condition = "input.plottab == 1",
          h4("Time integration"),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="tmax", label="Maximum time", value=numopts$tmax),
                      numericInput(inputId="tstep", label="Time step", value=numopts$tstep,step=0.1)),
          selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4"),selected=numopts$odemethod)),
        conditionalPanel(
          condition = "input.plottab > 1",
          h4("Curve continuation"),
          h5("Tolerances"),
          numericInput(inputId="rtol", label="Relative", value=numopts$rtol),
          numericInput(inputId="atol", label="Absolute", value=numopts$atol),
          numericInput(inputId="iszero", label="Zero identity", value=numopts$iszero),
          h5("Step size"),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="stepsize", label="Target", value=numopts$stepsize),
                      numericInput(inputId="minstepsize", label="Minimum", value=numopts$minstepsize)),
          numericInput(inputId="maxiter", label="Maximum number of iterations", value=numopts$maxiter),
          numericInput(inputId="maxpoints", label="Maximum number of points", value=numopts$maxpoints),
          numericInput(inputId="computedelay", label="Computational delay", min=0.0, max=1.0, value=numopts$computedelay, step=0.01)),

        actionButton("numoptsapply", "Apply", icon("refresh"))
      )
    )
  )
  return(ui)
}
