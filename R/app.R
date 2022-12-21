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
                menuSubItem('Log(2FC) data', tabName = 'log2fc_pval_dataframe', icon = icon('table')),
                menuSubItem('Pathway frequencies', tabName = 'pathway_frequencies', icon = icon('table'))

              ),
      menuItem("Plots", tabName = "plots", icon = icon("chart-line"),

               menuSubItem('Volcano plots', tabName = 'volcano_plots', icon = icon('chart-line')),
               menuSubItem('Gene regulation', tabName = 'regulation_barplot', icon = icon('chart-line')),
               menuSubItem('PCA', tabName = 'pca_plots', icon = icon('chart-line')),
               menuSubItem('Heatmaps', tabName = 'heatmaps', icon = icon('chart-line'))

      ),
      menuItem("Documentation", tabName = "documentation", icon = icon("book"))
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
                ),

                div(id = 'pipeline_run',
                    box(title = 'Run status', status = 'success', solidHeader = T, height = '50%',
                        h2('Pipeline run complete')# %>% withSpinner(type = 8, size = 0.5)
                    )
                ) %>% shinyjs::hidden()

                ),

              fluidRow(
                div(id = 'data_loaded',
                  box(title = 'Parameters', status = 'warning', solidHeader = T,
                    uiOutput('sampleSelection') %>% withSpinner(type = 8, size = 0.5),
                    textInput('log2fc_threshold', 'Log(2FC) threshold', value = '1.5'),
                    textInput('pval_threshold', 'P-value threshold', value = '0.01'),

                    div(style="float: left; position: relative; left: 40%; margin: 0px",actionButton("run_pipeline","Run pipeline"))
                  )
                ) %>% shinyjs::hidden()
              )
      ),

      tabItem(tabName = "documentation",
              h2('MPATH documentation'),
              br(),
              h3('Install as an R library'),
              br(),
              h4('The code can be installed using devtools:'),
              br(),
              p('install.packages("devtools")', style = 'text-indent: 30px;font-size: 15px'),
              p('library(devtools)', style = 'text-indent: 30px;font-size: 15px'),
              p('install_github("AudenCote/MPATH")', style = 'text-indent: 30px;font-size: 15px'),
      ),

      tabItem(tabName = "expression_data",
              h3('Reformatted expression data'),
              downloadButton("download_expression_data", "Download as .csv"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                  tableOutput("expression_data")))
      ),

      tabItem(tabName = "pathway_frequencies",
              h3('Up/downregulated genes by pathway'),
              uiOutput('pathway_tool_selection'),
              downloadButton("download_pathway", "Download as .csv"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                  tableOutput("pathway_frequency")))
      ),

      tabItem(tabName = "log2fc_pval_dataframe",
              h3('Log(2FC) and p-values relative to the benchmark sample'),
              downloadButton("download_log2fc_pval_dataframe", "Download as .csv"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   tableOutput("log2fc_pval_dataframe")))
      ),

      tabItem(tabName = "volcano_plots",
              uiOutput('volcano_plot_selection'),
              plotOutput("volcano_plot")
      ),

      tabItem(tabName = "regulation_barplot",
              plotOutput("regulation_barplot")

      ),

      tabItem(tabName = "pca_plots",
              plotOutput('sil_plot'),
              plotOutput('pca_plot'),
              plotOutput('loadings_plot')
      ),

      tabItem(tabName = "heatmaps",
              uiOutput('heatmap_pathway_selection'),
              plotlyOutput("heatmap")
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

  output$sampleSelection <- renderUI({
    selectInput("benchmark_sample", "Benchmark sample", choices = as.character(unique(inst()$expression_data$Sample)))
  })

  observeEvent(input$run_pipeline, {
    inst()$Volcano(benchmark_sample = input$benchmark_sample)
    inst()$GeneRegulation(as.numeric(input$log2fc_threshold), as.numeric(input$pval_threshold))
    inst()$PCA()
    inst()$Pathways()

    shinyjs::toggle("pipeline_run")
  })

  output$volcano_plot_selection <- renderUI({
    selectInput("volcano_plot_sample", "Sample to display",
                choices = as.character(unique(inst()$expression_data$Sample)[unique(inst()$expression_data$Sample) != input$benchmark_sample]))
  })

  output$volcano_plot <- renderPlot({
    inst()$volcano_plots[input$volcano_plot_sample]
  })

  output$heatmap_pathway_selection <- renderUI({
    selectInput("heatmap_pathway", "Pathway to display",
                choices = as.character(unique(inst()$top_pathway_genes_l2fcp$MitoPathway)))
  })

  output$heatmap <- renderPlotly({
    inst()$pathways$plots$heatmaps[[input$heatmap_pathway]]
  })

  output$pathway_tool_selection <- renderUI({
    selectInput("pathway_tool", "Pathway data type to display",
                choices = c('MitoCarta frequencies', 'Panther frequencies', 'Mitocarta Fisher results'),
                selected = 'MitoCarta')
  })

  output$expression_data <- renderTable(head(inst()$expression_data, 50))
  output$log2fc_pval_dataframe <- renderTable(head(inst()$log2fc_pval_dataframe, 50))
  output$pathway_frequency <- renderTable({
      inst()$pathways[input$pathway_tool]
  })

  output$regulation_barplot <- renderPlot(inst()$regulation_barplot$plot)
  output$pca_plot <- renderPlot(inst()$pca$pca_plot)
  output$sil_plot <- renderPlot(inst()$pca$sil_plot)
  output$loadings_plot <- renderPlot(inst()$pca$loadings_plot)

  output$download_expression_data <- downloadHandler(
    filename = 'MPATH_ExpressionData.csv',
    content = function(file) {
      write.csv(inst()$expression_data, file, row.names = FALSE)
    }
  )

  output$download_log2fc_pval_dataframe <- downloadHandler(
    filename = 'MPATH_Log2FC_PVals.csv',
    content = function(file) {
      write.csv(inst()$log2fc_pval_dataframe, file, row.names = FALSE)
    }
  )

  output$download_log2fc_pval_dataframe <- downloadHandler(
    filename = function(){
      paste(paste('MPATH_', str_replace(input$pathway_tool, ' ', '_')), '.csv', sep = '')
    },
    content = function(file) {
      write.csv(inst()$pathways[input$pathway_tool], file, row.names = FALSE)
    }
  )

}

shinyApp(ui = ui, server = server)

