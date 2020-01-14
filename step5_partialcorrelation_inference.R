# This script uses the partial correlation network inference method to predict interactions between TF genes and gene modules.

'''
Partial correlation network inference: 
Schäfer J, Strimmer K. An empirical Bayes approach to inferring large-scale gene association networks. Bioinformatics. 2004;21(6):754-64.
GeneNet package:
Schäfer J, Opgen-Rhein R, Strimmer K. Reverse engineering genetic networks using the GeneNet package. The Newsletter of the R Project Volume 6/5, December 2006. 2006;6(9):50.
'''

setwd("/path/to/analysis/directory/")

# Import and format data - 5 gene clusters and 131 TFs
exprdat <- read.csv("cluster_TF_expression.csv", row.names = 1)
IDs <- row.names(exprdat)
samples <- colnames(exprdat)
exprdat <- t(exprdat)
row.names(exprdat) <- samples
colnames(exprdat) <- IDs

# Run GeneNet
library("GeneNet")
set.seed(300)

# Create longitudinal object
explong <- as.longitudinal(exprdat, repeats = c(2,2,2,2,2), time = c(0, 2, 4, 6, 8))
summary(explong)

# Estimate partial correlations
pc = ggm.estimate.pcor(explong, method="dynamic")           # dynamic, with shrinkage

# Find edges
pc.edges = network.test.edges(pc, direct = TRUE, fdr = TRUE)
pc.net = extract.network(pc.edges, method.ggm = "qval", cutoff.ggm = 1)

# Add node labels to node numbers network
node.labels <- data.frame(colnames(explong)) 
features <- data.frame("node"=rownames(node.labels), node.labels)  # geneTF/genecluster number key
pc.net0 <- merge(x = pc.net, y = features, by.x = c("node1"), by.y = c("node"))  # match key to node1 #s in output
names(pc.net0)[names(pc.net0)=="colnames.explong."] <- "node1_ID"                # change column name of node 1 IDs
pc.net1 <- merge(x = pc.net0, y = features, by.x = c("node2"), by.y = c("node"))  # match key to node2 #s in output
names(pc.net1)[names(pc.net1)=="colnames.explong."] <- "node2_ID"                 # change column name of node 2 IDs

# Remove interactions that are undirected
pc.net2 <- pc.net1[!grepl("undirected", pc.net1$directions),]

# Filter for interactions where target node is a gene cluster
pc.net3 <- pc.net2[grep("module*", pc.net2$node1_ID),]

# Filter for interactions where TF gene is regulator and target is cluster
direction1 <- pc.net2[grepl("1to2", pc.net2$directions) 
                      & grepl("Glyma.*", pc.net2$node1_ID) 
                      & grepl("module*", pc.net2$node2_ID), ]
direction2 <- pc.net2[grepl("2to1", pc.net2$directions) 
                      & grepl("Glyma.*", pc.net2$node2_ID)
                      & grepl("cluster*", pc.net2$node1_ID), ]
pc.net4 <- rbind(direction1, direction2)
pc.net4 <- pc.net4[order(-pc.net4$pval),]  # sort by highest weight
write.csv(pc.net4, file = "partialcorrelation_network")




