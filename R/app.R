source('MPATH.R')
require(shiny)
require(shinycssloaders)
require(shinyjs)
require(shinydashboard)
require(DT)
require(writexl)
require(openxlsx)


options(shiny.maxRequestSize=300*1024^2)
options(spinner.color="#0275D8", spinner.color.background="#ffffff", spinner.size=2)

ui <- dashboardPage(
  dashboardHeader(title = 'MPATH'),
  dashboardSidebar(

    sidebarMenu(
      useShinyjs(),
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Tables", tabName = "tables", icon = icon("table"),

                menuSubItem('Expression data', tabName = 'expression_data', icon = icon('table')),
                menuSubItem('Log(2FC) data', tabName = 'log2fc_pval_dataframe', icon = icon('table')),
                menuSubItem('Genes of interest', tabName = 'goi_log2fc', icon = icon('table')),
                menuSubItem('Pathway frequencies', tabName = 'pathway_frequencies', icon = icon('table'))

              ),
      menuItem("Plots", tabName = "plots", icon = icon("chart-line"),

               menuSubItem('Volcano plots', tabName = 'volcano_plots', icon = icon('chart-line')),
               menuSubItem('Gene regulation', tabName = 'regulation_barplot', icon = icon('chart-line')),
               menuSubItem('PCA', tabName = 'pca_plots', icon = icon('chart-line')),
               menuSubItem('Pathway heatmaps', tabName = 'heatmaps', icon = icon('chart-line')),
               menuSubItem('Pathway regulation', tabName = 'path_reg_bars', icon = icon('chart-line'))

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
                        h2('Pipeline run complete'),
                        downloadButton("download_all", "Download all data")
                    )
                ) %>% shinyjs::hidden()

                ),

              fluidRow(
                div(id = 'data_loaded',
                  box(title = 'Parameters', status = 'warning', solidHeader = T,
                    uiOutput('sampleSelection') %>% withSpinner(type = 8, size = 0.5),
                    textInput('log2fc_threshold', 'Log(2FC) threshold', value = '1.5'),
                    textInput('pval_threshold', 'P-value threshold', value = '0.01'),

                    fileInput("gene_list_file", "Optional: input file with a list of genes (.csv or .tsv)",
                              multiple = F,
                              accept = c(".tsv",
                                         ".csv")),

                    div(style="float: left; position: relative; left: 40%; margin: 0px", actionButton("run_pipeline","Run pipeline"))
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
              downloadButton("download_expression_data", "Download as .xlsx"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   DT::dataTableOutput("expression_data")))
      ),

      tabItem(tabName = "pathway_frequencies",
              h3('Up/downregulated genes by pathway'),
              uiOutput('pathway_tool_selection'),
              downloadButton("download_pathway", "Download as .xlsx"),
              br(),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   DT::dataTableOutput("pathway_frequency")))
      ),

      tabItem(tabName = "log2fc_pval_dataframe",
              h3('Log(2FC) and p-values relative to the benchmark sample'),
              downloadButton("download_log2fc_pval_dataframe", "Download as .xlsx"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   DT::dataTableOutput("log2fc_pval_dataframe")))
      ),

      tabItem(tabName = "goi_log2fc",
              h3('Log(2FC) and p-values for genes of interest only'),
              downloadButton("download_goi_log2fc", "Download as .xlsx"),
              (div(style='overflow-x: scroll;overflow-y: scroll;',
                   DT::dataTableOutput("goi_log2fc_pval_dataframe")))
      ),

      tabItem(tabName = "volcano_plots",
              uiOutput('volcano_plot_selection'),
              plotlyOutput("volcano_plot")
      ),

      tabItem(tabName = "regulation_barplot",
              plotOutput("regulation_barplot")

      ),

      tabItem(tabName = "pca_plots",
              plotOutput('sil_plot'),
              uiOutput('pca_slider'),
              plotOutput('pca_plot'),
              plotOutput('loadings_plot')
      ),

      tabItem(tabName = "heatmaps",
              uiOutput('heatmap_pathway_selection'),
              plotlyOutput("heatmap")
      ),

      tabItem(tabName = "path_reg_bars",
              uiOutput('n_path_freq_slider'),
              uiOutput('path_sample_consider'),
              h3('Number of upregulated genes by pathway'),
              plotOutput("path_reg_bars_up"),
              h3('Number of downregulated genes by pathway'),
              plotOutput("path_reg_bars_down")
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
    inst()$filter_genes(input$gene_list_file$datapath)
    inst()$Volcano(benchmark_sample = input$benchmark_sample, as.numeric(input$log2fc_threshold), as.numeric(input$pval_threshold))
    inst()$GeneRegulation(as.numeric(input$log2fc_threshold), as.numeric(input$pval_threshold))
    inst()$silhouette()
    inst()$Pathways(method = 'MitoCarta')

    shinyjs::toggle("pipeline_run")
  })

  output$volcano_plot_selection <- renderUI({
    selectInput("volcano_plot_sample", "Sample to display",
                choices = as.character(unique(inst()$expression_data$Sample)[unique(inst()$expression_data$Sample) != input$benchmark_sample]))
  })

  output$volcano_plot <- renderPlotly({
    inst()$volcano_plots[[input$volcano_plot_sample]]
  })

  output$pca_slider <- renderUI({
    sliderInput("slider_clusters", label = h3("Change number of clusters"), min = 1,
              max = 10, value = as.integer(inst()$silhouette_clusters))
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
                selected = 'MitoCarta frequencies')
  })

  output$expression_data <- DT::renderDataTable(head(inst()$expression_data, 50))
  output$log2fc_pval_dataframe <- DT::renderDataTable({
    df <- head(inst()$log2fc_pval_dataframe, 50)
    df$P <- format(df$P, scientific = T)
    df
    })
  output$goi_log2fc_pval_dataframe <- DT::renderDataTable({
    df <- head(inst()$goi_log2fc_pval_dataframe, 50)
    df$P <- as.character(format(df$P, scientific = T))
    df
    })
  output$pathway_frequency <- DT::renderDataTable({
      data.frame(inst()$pathways[input$pathway_tool])
  })

  output$regulation_barplot <- renderPlot(inst()$regulation_barplot$plot)

  output$n_path_freq_slider <- renderUI({
    sliderInput("n_path_freq_slider", label = h6("Number of pathways to view"), min = 1,
                max = length(unique(inst()$path_freqs$MitoPathway)), value = 5)
  })

  output$path_sample_consider <- renderUI({
    selectInput("path_sample_consider", label = h6("Sample(s) to compare to reference"), choices = append(c('all'), as.character(unique(inst()$expression_data$Sample))[as.character(unique(inst()$expression_data$Sample)) != input$benchmark_sample]))
  })

  output$path_reg_bars_up <- renderPlot({
    inst()$pathway_freq_plots(input$path_sample_consider, input$n_path_freq_slider)
    inst()$pathways$plots$mitocarta.up
    }
  )

  output$path_reg_bars_down <- renderPlot({
    inst()$pathway_freq_plots(input$path_sample_consider, input$n_path_freq_slider)
    inst()$pathways$plots$mitocarta.down
  }
  )

  output$pca_plot <- renderPlot({
    inst()$PCA(as.integer(input$slider_clusters))
    inst()$pca$pca_plot
    }
  )

  output$loadings_plot <- renderPlot({
    inst()$PCA(as.integer(input$slider_clusters))
    inst()$pca$loadings_plot
  }
  )

  output$sil_plot <- renderPlot(inst()$pca$sil_plot)

  output$download_expression_data <- downloadHandler(
    filename = 'MPATH_ExpressionData.xlsx',
    content = function(file) {
      writexl::write_xlsx(inst()$expression_data, file)
    }
  )

  output$download_log2fc_pval_dataframe <- downloadHandler(
    filename = 'MPATH_Log2FC_PVals.xlsx',
    content = function(file) {
      writexl::write_xlsx(inst()$log2fc_pval_dataframe, file)
    }
  )

  output$download_pathway <- downloadHandler(
    filename = function(){
      paste(paste('MPATH_', str_replace(input$pathway_tool, ' ', '_')), '.xlsx', sep = '')
    },
    content = function(file) {
      writexl::write_xlsx(inst()$pathways[input$pathway_tool], file)
    }
  )

  output$download_goi_log2fc <- downloadHandler(
    filename = 'MPATH_GOI_Log2FC_PVals.xlsx',
    content = function(file) {
      writexl::write_xlsx(inst()$goi_log2fc_pval_dataframe, file)
    }
  )

  output$download_all <- downloadHandler(
    filename = 'MPATH_AllData.xlsx',
    content = function(file) {

      mitocarta_towrite <- data.frame(inst()$pathways['MitoCarta frequencies'])
      if(length(names(mitocarta_towrite)) == 4){
        names(mitocarta_towrite) <- c('MitoPathway', 'Sample', 'Sig.Up', 'Sig.Down')
      }

      panther_towrite <- data.frame(inst()$pathways['Panther frequencies'])
      if(length(names(panther_towrite)) == 4){
        names(panther_towrite) <- c('Pathway', 'Sample', 'Sig.Up', 'Sig.Down')
      }

      fisher_towrite <- data.frame(inst()$pathways['Mitocarta Fisher results'])
      if(length(names(fisher_towrite)) == 3){
        names(fisher_towrite) <- c('MitoPathway', 'P', 'Reg.Direction')
      }

      openxlsx::write.xlsx(list('Expression data' = inst()$expression_data,
                            'All genes Log2(FC) and P-values' = inst()$log2fc_pval_dataframe,
                            'Genes of interest' = inst()$goi_log2fc_pval_dataframe,
                            'MitoCarta frequencies' = mitocarta_towrite,
                            'Panther frequencies' = panther_towrite,
                            'MitoCarta Fisher results' = fisher_towrite), file)
    }
  )


}

shinyApp(ui = ui, server = server)

