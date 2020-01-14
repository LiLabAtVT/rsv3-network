# This script uses the least angle regression (LARS) network inference method to predict interactions between TF genes and gene modules.

'''
LARS network inference: 
Haury A-C, Mordelet F, Vera-Licona P, Vert J-P. TIGRESS: trustful inference of gene regulation using stability selection. BMC systems biology. 2012;6(1):145
Tigress package:
Haury A-C, Mordelet F, Vera-Licona P, Vert J-P. TIGRESS: trustful inference of gene regulation using stability selection. BMC systems biology. 2012;6(1):145.
'''

setwd("/path/to/analysis/directory/")

# Import and format data - 5 gene clusters and 131 TFs
exprdat <- read.csv("cluster_TF_expression.csv", row.names = 1)
expdat <- t(as.matrix(expdat))

# Run TIGRESS
library(tigress)
set.seed(100)
nstepsLARS = 4
net <- tigress(expdat, tflist=colnames(expdat[,6:136]), targetlist=colnames(expdat[,1:5]),
                    nstepsLARS = nstepsLARS, normalizeexp = TRUE, verb = TRUE)
net_df <- net[[4]]
net_df <- data.frame(regulator=rownames(net_df)[row(net_df)], target=colnames(net_df)[col(net_df)], probability=c(net_df))   # Convert correlation matrix to list
net_df <- net_df[!(net_df$probability==0),]   # Logic index to remove 0 probability
write.csv(net_df, file = "lars_network.csv")
