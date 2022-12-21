source('MPATH.R')
library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinydashboard)

options(shiny.maxRequestSize=300*1024^2)

ui <- dashboardPage(
  dashboardHeader(title = 'MPATH'),
  dashboardSidebar(

    sidebarMenu(
      useShinyjs(),
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Tables", tabName = "tables", icon = icon("table"),

                menuSubItem('Expression data', tabName = 'expression_data', icon = icon('table')),
                menuSubItem('Log(2FC) data', tabName = 'log2fc_pval_dataframe', icon = icon('table'))

              ),
      menuItem("Plots", tabName = "plots", icon = icon("chart-line"),

               menuSubItem('Volcano plots', tabName = 'volcano_plots', icon = icon('chart-line')),
               menuSubItem('Gene regulation', tabName = 'regulation_barplot', icon = icon('chart-line')),
               menuSubItem('PCA', tabName = 'pca_plots', icon = icon('chart-line')),
               menuSubItem('Heatmaps', tabName = 'heatmaps', icon = icon('chart-line'))

      )
    )

  ),
  dashboardBody(

    useShinyjs(),

    tabItems(
      tabItem(tabName = "dashboard",
              fluidRow(

                box(
                  title = 'Get started', status = 'primary', solidHeader = T,
                  fileInput("expression_file", "Input file with expression data (.csv or .tsv)",
                            multiple = F,
                            accept = c(".tsv",
                                       ".csv")),

                  p('The first column should be labeled "Gene" and contain Ensembl identifiers;
         all other columns should contain expression values and be labeled with the
         sample and replicate IDs separated by an underscore (e.g. S1_R1, S1_R2...).',
                    style = 'margin-top:-2%'),

                  div(style="float: left; position: relative; left: 40%; margin: 0px",actionButton("load_data","Load data"))
                )),

              fluidRow(
                div(id = 'data_loaded',
                  box(title = 'Parameters', status = 'warning', solidHeader = T,
                    uiOutput('sampleSelection') %>% withSpinner(type = 8, size = 0.5),
                    textInput('log2fc_threshold', 'Log(2FC) threshold', value = '1.5'),
                    textInput('pval_threshold', 'P-value threshold', value = '0.01'),

                    div(style="float: left; position: relative; left: 40%; margin: 0px",actionButton("run_pipeline","Run pipeline"))
                  )
                ) %>% shinyjs::hidden()
              ),

              fluidRow(
                div(id = 'pipeline_run',
                    box(title = 'Run status', status = 'warning', solidHeader = T,
                          h2('Pipeline run complete') %>% withSpinner(type = 8, size = 0.5),
                        )
                ) %>% shinyjs::hidden()
              )
      ),

      tabItem(tabName = "expression_data",
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                  tableOutput("expression_data")))
      ),

      tabItem(tabName = "log2fc_pval_dataframe",
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   tableOutput("log2fc_pval_dataframe")))
      ),

      tabItem(tabName = "volcano_plots",
              uiOutput("volcano_plots")
      ),

      tabItem(tabName = "regulation_barplot",
              plotOutput("regulation_barplot")

      ),

      tabItem(tabName = "pca_plots",
              plotOutput('elbow_plot'),
              plotOutput('pca_plot'),
              plotOutput('loadings_plot')
      ),

      tabItem(tabName = "heatmaps",
              plotlyOutput("hm1")
      )


    )

  )
)


server <- function(input, output, session) {

  inst <- reactiveVal()

  observeEvent(input$load_data, {
    inst(MPATH_Pipeline(expression_file = input$expression_file$datapath))

    shinyjs::toggle("data_loaded")
  })

  output$volcano_plots <- renderUI({
    plot_output_list <- lapply(unique(inst()$expression_data$Sample[inst()$expression_data$Sample != input$benchmark_sample]), function(i) {
      plotOutput(i)
    })

    do.call(tagList, plot_output_list)
  })

  observeEvent(input$run_pipeline, {
    inst()$Volcano(benchmark_sample = input$benchmark_sample)

    for (i in 1:length(unique(inst()$expression_data$Sample))) {
        local({
          samp = unique(inst()$expression_data$Sample)[i]
          output[[samp]] <- renderPlot({
            inst()$volcano_plots[samp]
          })
        })
    }

    #inst()$GeneRegulation(as.integer(input$log2fc_threshold), as.integer(input$pval_threshold))
    inst()$GeneRegulation()
    inst()$PCA(clusters = 3)
    inst()$Pathways()

    shinyjs::toggle("pipeline_run")
  })

  output$sampleSelection <- renderUI({
    selectInput("benchmark_sample", "Benchmark sample", choices = as.character(unique(inst()$expression_data$Sample)))
  })

  output$expression_data <- renderTable(head(inst()$expression_data, 50))
  output$log2fc_pval_dataframe <- renderTable(head(inst()$log2fc_pval_dataframe, 50))
  output$regulation_barplot <- renderPlot(inst()$regulation_barplot$plot)
  output$pca_plot <- renderPlot(inst()$pca$pca_plot)
  output$elbow_plot <- renderPlot(inst()$pca$elbow_plot)
  output$loadings_plot <- renderPlot(inst()$pca$loadings_plot)
  output$hm1 <- renderPlotly(inst()$pathways$plots$heatmaps$`Lipid metabolism`)

}

shinyApp(ui = ui, server = server)












# ui <- dashboardPage(
#   dashboardHeader(title = 'MPATH'),
#   dashboardSidebar(
#
#       fileInput("expression_file", "Input file with expression data (.csv or .tsv)",
#         multiple = F,
#         accept = c(".tsv",
#                    ".csv")),
#
#       p('The first column should be labeled "Gene" and contain Ensembl identifiers;
#         all other columns should contain expression values and be labeled with the
#         sample and replicate IDs separated by an underscore (e.g. S1_R1, S1_R2...).',
#         style = 'margin-top:-5%;'),
#
#       actionButton("run_pipeline","Run pipeline")
#
#
#     ),
#
#     dashboardBody(
#       div(id = 'display-panels',
#         h2('Display plots'),
#
#           # selectInput("display1", "Display:",
#           #             c("Expression data" = "expression_data",
#           #               "Volcano plot" = "volcano_plots",
#           #               "Log(2FC) data" = "log2fc_pval_dataframe")),
#
#          (div(style='height:50%;overflow-x: scroll;overflow-y: scroll;',
#              tableOutput("Table1") %>% withSpinner(type = 8, size = 0.5)))
#         #tableOutput("Table1")# %>% withSpinner(type = 8, size = 0.5)
#       ) %>% shinyjs::hidden()
#     )
# )
#
#
#
#




