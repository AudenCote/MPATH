library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinydashboard)

setwd('/Users/acl/Desktop/Khrapko/MPATH/Dev/MPATH/R')

source('MPATH.R')


inst <- MPATH_Pipeline(expression_file = '../../OvarianStage.tsv')
#inst$filter_genes(gene_list_file = '../../SampleGeneList.csv')
inst$Volcano(benchmark_sample = 'S1')
inst$GeneRegulation()
inst$silhouette()
inst$PCA(3)
inst$Pathways(method = 'MitoCarta')
#inst$pathway_freq_plots('S3', 5)
