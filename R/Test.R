source('R/MPATH.R')

mpath <- MPATH_Pipeline(expression_file = '/Users/acl/Desktop/Khrapko/MPATH/Dev/ReformattedExpressionData.tsv')

#MAYBE SEPARATE DATA INITIALIZATION FROM THE VOLCANO FUNCTION TO ENABLE MODULARIZATION???
mpath$Volcano()
mpath$GeneRegulation(log2fc_threshold = 1.5, pval_threshold = 0.01)
mpath$PCA()

mpath$Pathways(method = 'MitoCarta')
