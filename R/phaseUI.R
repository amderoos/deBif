phaseUI <- function(state, parms, plotopts, numopts) {
  if (length(state) == 1)
    choices <- c("Time plot"=1, "Nullclines"=2, "Steady states"=3)
  else choices <- c("Time plot"=1, "Nullclines"=2, "Steady states"=3, "Vector field"=4, "Trajectories"=5, "Portrait"=6)
  myTabs <- lapply(1:length(choices), function(i) {tabPanel(title = names(choices)[i], plotOutput(outputId=paste0("plot", choices[i]), height = "100%"), value = choices[i])})

  ui <- dashboardPagePlus(
    shinyjs::useShinyjs(),
    ########## Header of window
    header = dashboardHeaderPlus(
      # Set height of dashboardHeader
      title = tagList(
        span(class = "logo-lg", "Phaseplane analysis"),
        icon("compass"), tags$style(".fa-compass {color:#E87722}")),
      left_menu = tagList(span(class = "help-button", icon("question-circle"),
                               tags$style(".fa-question-circle {font-size: 24px; color:#66CC66; left: 10px; top: 13px; position: absolute;}"))),
      tags$li(class = "dropdown", actionButton("showODEs", "Show ODEs", class = "show-odes"),
              tags$style(".show-odes {font-size: 13px;
                                      border-width:2px;
                                      width: 85px; height: 30px;
                                      text-indent: -8px;
                                      left: 50px; top: 12px; position: absolute;}")),
      titleWidth = 220,
      enable_rightsidebar = TRUE,
      rightSidebarIcon = "gears"
    ),
    ########## Left side-bar
    sidebar = dashboardSidebar(
      width = 220,
      tags$style(type='text/css', "#computefwrd1 { font-size: 13px; margin-left: 2px; }"),
      tags$style(type='text/css', "#computebwrd1 { font-size: 13px;}"),
      h4("Initial values", align = "center"),
      selectizeInput('selectpoint1', "", c("User specified" = 0), selected=0),
      sidebarMenu(
        id = "varsparscurvesmenu",
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
        actionButton("lapply", "Apply", icon("refresh")),
        br(),
        conditionalPanel(condition = "input.plottab == 1 || input.plottab == 5",
                         splitLayout(cellWidths = c("52%", "48%"),
                                     actionButton("computebwrd1", "Compute", icon("backward")),
                                     actionButton("computefwrd1", "Compute", icon("forward")))),
        br(),
        h4("Curve management", align = "center"),
        menuItem(h4("Load curves"),
                 textInput("loadcurve1", "Give a valid R variable name"),
                 splitLayout(cellWidths = c("50%", "50%"),
                             actionButton("appendbtn1", "Append"),
                             actionButton("replacebtn1", "Replace")),
                 tabName = "loadtab1"
        ),
        menuItem(h4("Save curves"),
                 selectInput('savecurve1', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                 textInput("curvename1", "Give a valid R variable name"),
                 actionButton("savebtn1", "Save curve"),
                 tabName = "savetab1"
        ),
        menuItem(h4("Delete curves"),
                 selectInput('deletecurve1', 'Select curve', c("None" = 0, "All" = -1), selected=0),
                 actionButton("deletebtn1", "Delete curve"),
                 tabName = "deletetab1"
        )
      )),
    body = dashboardBody(
      tags$script(
        "
        Shiny.addCustomMessageHandler('scrollCallback',
        function(color) {
        var objDiv = document.getElementById('console');
        objDiv.scrollTop = objDiv.scrollHeight;
        }
        );
        // Bind function to respond to the help button
        $('i.fa.fa-question-circle').on('click',function(){
          Shiny.onInputChange('helpClicked', Math.random());
        });
        // Bind function to the toggle sidebar button
        $('i.fa.fa-gears').on('click',function(){
          $(window).trigger('resize'); // Trigger resize event
          $(window).trigger('resize'); // Trigger resize event
        });
        var dimension = [0, 0];
        $(document).on('shiny:connected', function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange('dimension', dimension);
        });
        $(window).resize(function(e) {
            dimension[0] = window.innerWidth;
            dimension[1] = window.innerHeight;
            Shiny.onInputChange('dimension', dimension);
        });
        "
      ),
      do.call(tabsetPanel, c(myTabs, id = "plottab")),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#saveplot { left: calc(100% - 125px); top: 108px; position: absolute;}"
      ))),
      downloadButton("saveplot", label = "Save png"),
      shiny::tags$head(shiny::tags$style(shiny::HTML(
        "#console { font-size: 11px; width: calc(100%); left: calc(242px); height: 149px; overflow: auto; }"
      ))),
      tags$head(tags$style(HTML('.box {margin-bottom: 0px; margin-top: 0px;}'))),
      box(width = NULL, height = "170px", verbatimTextOutput("console", placeholder = TRUE)),
    ),
    rightsidebar = rightSidebar(
      background = "dark",
      # I copied here the function rightSidebarTabContent and changed it to get
      # rid of the additional blank heading line
      shiny::tags$div(
        tags$head(
          tags$style(
            HTML(
              '
              .form-group {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              .shiny-split-layout {
              margin-top: 0 !important;
              margin-bottom: 2px !important;
              }
              #xmin{height: 30px}
              #xmax{height: 30px}
              #ymin{height: 30px}
              #ymax{height: 30px}
              #y2min{height: 30px}
              #y2max{height: 30px}
              #tmax{height: 30px}
              #tstep{height: 30px}
              '
            )
          )
        ),
        class = "tab-pane active",
        id = "control-sidebar-plotopttab-tab",
        icon = "chart-line",
        h3("Plot options"),
        selectInput('xcol', 'Variable(s) on X-axis', c("Time" = 1), selected=plotopts[[1]]$xcol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="xmin", label="Minimum", value=plotopts[[1]]$xmin),
                    numericInput(inputId="xmax", label="Maximum", value=plotopts[[1]]$xmax)),
        selectInput('logx', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logx),
        div(style="line-height: 6px !important", br()),
        selectInput('ycol', 'Variable(s) on Y-axis',
                    c("All" = 1, setNames((2:(length(state)+1)), names(state))), selected=plotopts[[1]]$ycol),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="ymin", label="Minimum", value=plotopts[[1]]$ymin),
                    numericInput(inputId="ymax", label="Maximum", value=plotopts[[1]]$ymax)),
        selectInput('logy', 'Scale type', c("Linear" = 0, "Logarithmic" = 1), selected=plotopts[[1]]$logy),
        div(style="line-height: 6px !important", br()),
        conditionalPanel(
          condition = "input.plottab == 1 && input.ycol > 1",
          selectInput('y2col', 'Variable on 2nd Y-axis',
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
            conditionalPanel(condition = "input.plot3d == 1",
                             sliderInput(inputId="theta", label="Viewing angle",
                                         min=-90, max=90, value=plotopts[[1]]$theta, step=1,
                                         ticks = FALSE, round=TRUE))
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
        h4("Time integration"),
        splitLayout(cellWidths = c("50%", "50%"),
                    numericInput(inputId="tmax", label="Maximum time", value=numopts$tmax),
                    numericInput(inputId="tstep", label="Time step", value=numopts$tstep, min=0, step=0.1)),
        selectInput('method', 'Integrator', c("lsoda", "ode23", "ode45", "rk4"),selected=numopts$odemethod),
        div(style="line-height: 12px !important", br()),
        actionButton("numoptsapply", "Apply", icon("refresh")
        )
      )
    )
  )

  return(ui)
}
