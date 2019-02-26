buildUI <- function(state, parms, plotopts, numopts) {
  choices <- c("Time series" = 1, "1 parameter bifurcation" = 2, "2 parameter bifurcation" = 3)
  myTabs <- lapply(1:length(choices), function(i) {tabPanel(title = names(choices)[i], plotOutput(outputId=paste0("plot", choices[i]), height = "100%"), value = choices[i])})

  ui <- dashboardPagePlus(
    shinyjs::useShinyjs(),
    ########## Header of window
    header = dashboardHeaderPlus(
      # Set height of dashboardHeader
      title = tagList(
        span(class = "logo-lg", "Bifurcation analysis"),
        icon("compass"), tags$style(".fa-compass {color:#E87722}")), titleWidth = 220,
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "gears"
    ),
    ########## Left side-bar
    sidebar = dashboardSidebar(
      width = 220,
      tags$style(type='text/css', "#computebtn { font-size: 13px;}"),
      tags$style(type='text/css', "#computefwrd2 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computefwrd3 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computebwrd2 { font-size: 13px;}"),
      tags$style(type='text/css', "#computebwrd3 { font-size: 13px;}"),
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
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tabName = "state1tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_1"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tabName = "pars1tab"
          ),
          br(),
          actionButton("computebtn", "Compute", icon("forward")),
          br(),
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
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tabName = "state2tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_2"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tabName = "pars2tab"
          ),
          br(),
          splitLayout(cellWidths = c("50%", "50%"),
                      actionButton("computebwrd2", "Compute", icon("backward")),
                      actionButton("computefwrd2", "Compute", icon("forward"))),
          br(),
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
                                            value=as.numeric(state[i]),step=0.1*as.numeric(state[i]),width="95%")}),
            tabName = "state3tab"
          ),
          menuItem(
            h4("Parameters"),
            lapply(1:length(parms),
                   function(i){numericInput(inputId=paste0(names(parms[i]), "_3"), label=h6(names(parms[i])),
                                            value=as.numeric(parms[i]), step=0.1*as.numeric(parms[i]),width="95%")}),
            tabName = "pars3tab"
          ),
          br(),
          splitLayout(cellWidths = c("50%", "50%"),
                      actionButton("computebwrd3", "Compute", icon("backward")),
                      actionButton("computefwrd3", "Compute", icon("forward"))),
          br(),
          id = "leftsbmenu"
        ),
        radioButtons("curvetype", "Curve type to compute", choices = c("BP", "HP", "LP"), inline = TRUE)
      ),
      sidebarMenu(
        h4("Curve management", align = "center"),
        menuItem(h4("Load curves"),
                 textInput("loadcurve", "Give a valid R variable name"),
                 splitLayout(cellWidths = c("50%", "50%"),
                             actionButton("appendbtn", "Append"),
                             actionButton("replacebtn", "Replace")),
                 tabName = "loadtab"
        ),
        menuItem(h4("Save curves"),
                 selectInput('savecurve', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                 textInput("curvename", "Give a valid R variable name"),
                 actionButton("savebtn", "Save curve"),
                 tabName = "savetab"
        ),
        menuItem(h4("Delete curves"),
                 selectInput('deletecurve', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                 actionButton("deletebtn", "Delete curve"),
                 tabName = "deletetab"
        ),
        radioButtons("reportlevel", "Console report level", choiceNames = c("None", "Normal", "Full"), choiceValues = (0:2), inline = TRUE),
        div(style="line-height: 18px !important", br()),
        conditionalPanel(
          condition = "input.plottab == 2 | input.plottab == 3",
          shinyjs::hidden(actionButton("pausebtn", "Pause", icon("pause-circle"), style = "color: white;
                       background-color: green;
                       font-size: 18px;
                       position: relative;
                       left: 18px;
                       height: 45px;
                       width: 150px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 6px;
                       border-width: 2px")),
          div(style="line-height: 8px !important", br()),
          shinyjs::hidden(actionButton("stopbtn", "Stop", icon("stop-circle"), style = "color: white;
                       background-color: red;
                       font-size: 18px;
                       position: relative;
                       left: 18px;
                       height: 45px;
                       width: 150px;
                       text-align:center;
                       text-indent: -2px;
                       border-radius: 6px;
                       border-width: 2px"))
        )
      )
    ),
    ########## Main panels
    body = dashboardBody(
      do.call(tabsetPanel, c(myTabs, id = "plottab")),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#console { font-size: 11px; width: calc(100%); left: calc(242px); height: 149px; overflow: auto; }"
      ))),
      tags$script(
        '
        Shiny.addCustomMessageHandler("scrollCallback",
        function(color) {
        var objDiv = document.getElementById("console");
        objDiv.scrollTop = objDiv.scrollHeight;
        }
        );'
      ),
      tags$head(tags$style(HTML('.box {margin-bottom: 0px; margin-top: 0px;}'))),
      box(width = NULL, height = "170px", verbatimTextOutput("console", placeholder = TRUE)),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#progress { font-size: 10px; width: calc(100%); left: calc(242px); height: 72px; overflow: auto;}"
      ))),
      tags$head(tags$style(HTML('.box {margin-bottom: 0px; margin-top: 0px;}'))),
      box(width = NULL, height = "90px", verbatimTextOutput("progress", placeholder = TRUE))
    ),
    ########## Right side-bar
    rightsidebar = rightSidebar(
      background = "dark",
      # I copied here the function rightSidebarTabContent and changed it to get
      # rid of the additional blank heading line
      shiny::tags$div(
        tags$head(
          tags$style(
            HTML(
              "
              .form-group {
              margin-bottom: 0 !important;
              }
              "
            )
          )
        ),
        class = "tab-pane active",
        id = "control-sidebar-plotopttab-tab",
        icon = "chart-line",
        h3("Plot options"),
        selectInput('xcol', ('Variable(s) on X-axis'),
                    c("Time" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$xcol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="xmin", label="Minimum", value=plotopts[[1]]$xmin),
                    numericInput(inputId="xmax", label="Maximum", value=plotopts[[1]]$xmax)),
        selectInput('logx', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
        div(style="line-height: 6px !important", br()),
        selectInput('ycol', ('Variable(s) on Y-axis'),
                    c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$ycol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="ymin", label="Minimum", value=plotopts[[1]]$ymin),
                    numericInput(inputId="ymax", label="Maximum", value=plotopts[[1]]$ymax)),
        selectInput('logy', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
        div(style="line-height: 6px !important", br()),
        conditionalPanel(
          condition = "input.plottab != 3 && input.ycol > 1",
          selectInput('y2col', ("Variable on 2nd Y-axis"),
                      c("None" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$y2col),
          conditionalPanel(
            condition = "input.y2col > 1",
            splitLayout(cellWidths = c("50%", "50%"),
                        numericInput(inputId="y2min", label="Minimum", value=plotopts[[1]]$y2min),
                        numericInput(inputId="y2max", label="Maximum", value=plotopts[[1]]$y2max)),
            selectInput('logy2', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy2),
            splitLayout(cellWidths = c("50%", "50%"),
                        strong("Plot type"),
                        radioButtons("plot3d", NULL, choices = c("2D" = 0, "3D" = 1), selected = plotopts[[1]]$plot3d,
                                     inline = TRUE)),
            sliderInput(inputId="theta", label="Viewing angle", min=-90, max=90, value=plotopts[[1]]$theta, step=1,
                        ticks = FALSE, round=TRUE)
          )
        ),
        div(style="line-height: 12px !important", br()),
        actionButton("plotoptsapply", "Apply", icon("refresh"))
      ),
      # sidebarMenu() and menuItem() do not work in the right sidebar. Only conditionpanels can be used
      shiny::tags$div(
        class = "tab-pane",
        id = "control-sidebar-numopttab-tab",
        icon = "dashboard",
        h3("Numerical options"),
        conditionalPanel(
          condition = "input.plottab == 1",
          h4("Time integration"),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="tmax", label="Maximum time", value=numopts$tmax),
                      numericInput(inputId="tstep", label="Time step", value=numopts$tstep,step=0.1)),
          selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4"),selected=numopts$odemethod)
          ),
        conditionalPanel(
          condition = "input.plottab > 1",
          h4("Curve continuation"),
          h4("Tolerances"),
          splitLayout(cellWidths = c("50%", "50%"),
                      textInput(inputId="rtol", label="Relative", value=sprintf("%.1E", numopts$rtol)),
                      textInput(inputId="atol", label="Absolute", value=sprintf("%.1E", numopts$atol))),
          textInput(inputId="iszero", label="Zero identity", value=sprintf("%.1E", numopts$iszero)),
          textInput(inputId="jacdif", label="Jacobian pertubation", value=sprintf("%.1E", numopts$jacdif)),
          div(style="line-height: 6px !important", br()),
          h4("Step size"),
          splitLayout(cellWidths = c("50%", "50%"),
                      numericInput(inputId="stepsize", label="Target", value=numopts$stepsize),
                      numericInput(inputId="minstepsize", label="Minimum", value=numopts$minstepsize)),
          h4("Iterations"),
          numericInput(inputId="maxiter", label="Maximum number of iterations", value=numopts$maxiter),
          numericInput(inputId="maxpoints", label="Maximum number of points", value=numopts$maxpoints),
          numericInput(inputId="replotfreq", label="Points between plot updates", min=1, max=10000,
                       value=numopts$replotfreq, step=1)
          ),
        div(style="line-height: 12px !important", br()),
        actionButton("numoptsapply", "Apply", icon("refresh"))
      )
    )
  )
  return(ui)
}
