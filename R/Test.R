source('R/MPATH.R')

mpath <- MPATH_Pipeline(expression_file = '/Users/acl/Desktop/Khrapko/MPATH/Dev/ExampleExpressionData.tsv')

mpath$Volcano(benchmark_sample = '1')
mpath$GeneRegulation(log2fc_threshold = 1.5, pval_threshold = 0.01)
mpath$PCA()
mpath$Pathways()


