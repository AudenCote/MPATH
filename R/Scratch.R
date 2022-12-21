library(shiny)
library(shinycssloaders)
library(shinyjs)
library(shinydashboard)

setwd('/Users/acl/Desktop/Khrapko/MPATH/Dev/MPATH/R')

source('MPATH.R')


inst <- MPATH_Pipeline(expression_file = '../../WholeCellRNASeq.tsv')
inst$Volcano(benchmark_sample = 'S1')
inst$GeneRegulation()
inst$PCA(clusters = 3)
#inst$Pathways(method = 'MitoCarta')

# sample_order = c('S3', 'S5', 'S7', 'S9', 'S12', 'S15', 'S18', 'S20')
# inst$test_hmdb <- inst$test_hmdb[, sample_order]
# hm <- heatmaply(inst$test_hmdb, Colv = NA, col = colorRampPalette(brewer.pal(3, "Blues"))(25))
