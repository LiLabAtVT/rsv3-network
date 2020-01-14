# This script uses the Random Forest network inference method to predict interactions between TF genes and gene modules.

'''
Random Forest network inference: 
Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods. PloS one. 2010;5(9):e12776.
GENIE3 package:
Huynh-Thu VA, Irrthum A, Wehenkel L, Geurts P. Inferring regulatory networks from expression data using tree-based methods. PloS one. 2010;5(9):e12776.
'''

setwd("/path/to/analysis/directory/")

# Import and format data - 5 gene clusters and 131 TFs
exprdat <- read.csv("cluster_TF_expression.csv", row.names = 1)
rows <- row.names(exprdat)
exprdat <- as.matrix(sapply(exprdat, as.numeric))
rownames(exprdat) <- paste(rows)
# Import DEG transcription factors and format
tfs <- read.csv("TF_expression.csv")
tfs <- tfs$X
tfs <- as.character(tfs)

# Run GENIE3
library(GENIE3)
set.seed(123) 
net <- GENIE3(exprdat, regulators = tfs, treeMethod = "RF", nTrees = 1000)

# Get list all interactions
net0 <- getLinkList(net)

# Get TF-gene cluster interactions
net1 <- net0[grep("module*", net0$targetGene), ]
net2 <- net1[order(-net1$weight),]  # sort by highest weight
write.csv(net2, file = "randomforest_network.csv")

