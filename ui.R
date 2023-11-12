library(shinycssloaders)

ui <- fluidPage(
  
  titlePanel("Bayesian optimal group sequential design with time to event endpoint"),
  "More details and instructions to use this app can be found ",
  tags$a(href="https://github.com/Stat-Li/BayesianOptimalDesign", 
         "here."),
  fluidRow(
    column(3,
           wellPanel(
             h4("Design Type:"),
             radioButtons("design_arm", label = ("Single-arm or two-arm RCT"),
                          choices = c("Single-arm",
                                      "Two-arm RCT"),
                          selected = "Single-arm"),
             radioButtons("design_stop", label = ("Stopping"),
                          choices = c("Futility stopping only",
                                         "Futility and superiority stopping"),
                          selected = "Futility stopping only"),
             numericInput("delta0",
                          label = ("Noninferiority Margin"),
                          value = 1.0, min = 1, max = 5, step = 0.1),
             textInput("interim_timing", label = ("Information time"),
                       value = "0.5"),
             helpText('Note: Seperate timing by space.', style = "font-size:10px"),
             conditionalPanel(
               condition = "input.design_arm == 'Two-arm RCT'",
               numericInput("ss_ratio",
                            label = ("Experimental arm and control arm sample size ratio"),
                            value = 1, min = 0.1, max = 10, step = 0.1)
             )
           ),
           wellPanel(
             h4("Simulation Setting:"),
             numericInput("delta1",
                          label = ("Desired hazard ratio"),
                          value = 0.6, min = 0, max = 1, step = 0.1),
             numericInput("type1",
                          label = ("Type I Error"),
                          value = 0.1, min = 0, max = 1, step = 0.1),
             numericInput("type2",
                          label = ("Type II Error"),
                          value = 0.2, min = 0, max = 1, step = 0.1),
             numericInput("weibull_shape",
                          label = ("Weibull shape parameter"),
                          value = 0.5, min = 0, max = 3, step = 0.1),
             conditionalPanel(
               condition = "input.design_arm == 'Single-arm'",
               textInput("gamma_a", label = ("Gamma prior shape parameter"),
                         value = "2 5 10")
             ),
             conditionalPanel(
               condition = "input.design_arm == 'Two-arm RCT'",
               textInput("sigma2_0", label = ("Normal prior variance parameter"),
                         value = "0.5 0.2 0.1")
             ),
             # textInput("gamma_a", label = ("Gamma prior shape parameter"),
             #           value = "2 5 10"),
             helpText('Note: Multiple values allowed. Seperate each shape by space.', style = "font-size:10px"),
             numericInput("ta",
                          label = ("Accrual duration (years)"),
                          value = 4, min = 0, max = 10, step = 1),
             numericInput("tf",
                          label = ("Follow-up duration (years)"),
                          value = 2, min = 0, max = 10, step = 1),
             # numericInput("surv_3year", 
             #              label = ("3-year survival probability for control"), 
             #              value = 0.53, min = 0, max = 1, step = 0.01),
             h5('x-year survival probability for control'),
             div(style="display: inline-block;vertical-align:top; width: 100px;",
                 numericInput("x_year", label = "Year", value = 3, min = 0, max = 10, step = 0.1)),
             div(style="display: inline-block;vertical-align:top; width: 100px;",
                 numericInput("surv_3year", label = ("Survival prob."), value = 0.53, min = 0, max = 1, step = 0.01)),
             numericInput("search_gamma",
                          label = withMathJax("Tuning parameter \\(\\gamma\\) upper bound"),
                          value = 0.5, min = 0, max = 10, step = 0.1),
             numericInput("n_sims", 
                          label = ("Number of simulation trials"), 
                          value = 2000, min = 0, max = 100000, step = 100),
             actionButton("go", "Submit", icon("redo"))
           )
    ),
    column(9,
           # wellPanel(
           #   h4("Frequentist operating characteristics:"),
           #   # dataTableOutput('out_table')
           #   shinycssloaders::withSpinner(
           #     dataTableOutput('out_table')
           #   ),
           #   span('PRN: probability of rejecting null. PET: probability of early stopping.
           #        ES: expected sample size. EE expected number of events.')
           # )
           mainPanel(
             
             # Output: Tabset w/ plot, summary, and table ----
             tabsetPanel(type = "tabs",
                         tabPanel("Operating characteristics", shinycssloaders::withSpinner(
                           dataTableOutput('out_table_ops')
                         ),
                         span('PRN: probability of rejecting null. PET: probability of early stopping.
                  ES: expected sample size. EE expected number of events.')),
                         tabPanel("Stopping cutoff", dataTableOutput('out_table_stop'), 
                         span('Optimal stopping boundary that minimize the expected sample size 
                                       under null hypothesis.')),
                         
             )
             
           )
    )
  )
)