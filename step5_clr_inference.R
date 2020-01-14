# This script uses the context likelihood of relatedness (CLR) network inference method to predict interactions between TF genes and gene modules.

'''
CLR network inference: 
Faith JJ, Hayete B, Thaden JT, Mogno I, Wierzbowski J, Cottarel G, et al. Large-scale mapping and validation of Escherichia coli transcriptional regulation from a compendium of expression profiles. PLoS biology. 2007;5(1):e8.
minet package:
Meyer PE, Lafitte F, Bontempi G. minet: AR/Bioconductor package for inferring large transcriptional networks using mutual information. BMC bioinformatics. 2008;9(1):461.
'''

setwd("/path/to/analysis/directory/")

# Import and format data - 5 gene clusters and 131 TFs
expdat <- read.csv("cluster_TF_expression.csv", row.names = 1)
expdat <- as.matrix(expdat)

# Run minet with CLR inference
library(minet)
set.seed(345)
mim <- build.mim(t(expdat), estimator = "pearson")
net <- minet::clr(mim)
net_df <- data.frame(row=rownames(net)[row(net)], col=colnames(net)[col(net)], corr=c(net))   # Convert correlation matrix to list
net_df<-net_df[!(net_df$corr==0),]    # Logic index to remove zero correlation
net_df0 <- net_df[grep("module*", net_df$col), ]
net_df0 <- net_df0[grep("Glyma*", net_df0$row), ]
net_df0 <- net_df0[order(-net_df0$corr),]  # sort by highest weight
write.csv(net_df0, file = "clr_network.csv")