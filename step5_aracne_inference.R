# This script uses the ARACNE (mutual information) network inference method to predict interactions between TF genes and gene modules.

'''
ARACNE network inference: 
Margolin AA, Nemenman I, Basso K, Wiggins C, Stolovitzky G, Dalla Favera R, et al. ARACNE: an algorithm for the reconstruction of gene regulatory networks in a mammalian cellular context. BMC Bioinformatics. 2006;7(1):S7.
minet package:
Meyer PE, Lafitte F, Bontempi G. minet: AR/Bioconductor package for inferring large transcriptional networks using mutual information. BMC bioinformatics. 2008;9(1):461.
'''

setwd("/path/to/analysis/directory/")

# Import and format data - 5 gene clusters and 131 TFs
exprdat <- read.csv("cluster_TF_expression.csv", row.names = 1)
samples <- colnames(exprdat)
features <- rownames(exprdat)
exprdat <- t(exprdat)

# Run minet with ARACNE inference
library(minet)
set.seed(123)
mim <- build.mim(exprdat, estimator = "pearson")
net <- aracne(mim, eps = .1)
net_df <- data.frame(row=rownames(net)[row(net)], col=colnames(net)[col(net)], corr=c(net))   # Convert correlation matrix to list
net_df<-netp_df[!(net_df$corr==0),]    # Logic index to remove zero correlation
net_df0 <- net_df[grep("module*", net_df$col), ]
net_df0 <- net_df0[grep("Glyma*", net_df0$row), ]
net_df0 <- net_df0[order(-net_df0$corr),]  # sort by highest weight
write.csv(net_df0, file = "aracne_network.csv")