# This script is for performing differential expression analyses (with DESeq2) between L29 and Williams82 at each time point using read count data from featurecounts.

'''
DESeq2 from:
Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology. 2014;15(12):550.
featurecounts from:
Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics. 2013;30(7):923-30.
'''


setwd("/path/to/analysis/directory/")

# Setting up count matrix - rows as genes and columns as samples
countdata <- read.csv("featurecounts_output.csv", header=T, row.names=1)
countdata <- as.data.frame(countdata)

# Remove Wm82 string at end of gene names
rownames(countdata)<-gsub(".Wm82.a2.v1", "",rownames(countdata))

# Sample information - genotypes and stages
sampleTable <- read.csv("deseq2_sampletable.csv", row.names=1, header=T, stringsAsFactors = TRUE)  # 'resistant' corresponds to L29, and 'susceptible' corresponds to Williams82.

# Set reference levels
library("magrittr")
sampleTable$condition %<>% relevel("susceptible")
sampleTable$time %<>% relevel("0h")

#-------------------- Data Quality Check with PCA thru DESeq2 -------------------#
library(DESeq2)
coldata <- data.frame(sampleTable)
ddsTC <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design= ~ condition + time + condition:time)
ddsTC

# Prefiltering Dataset - removing rows (genes) with < 1 count (little gene info)
nrow(ddsTC)    
ddsTC <- ddsTC[rowSums(counts(ddsTC)) > 1, ]
nrow(ddsTC)

# Variance Stabilizing Transformation (VST normalize) - use this expression data for downstream analysis
vsd <- vst(ddsTC)
write.csv(assay(vsd), file = "deseq2_VSTnormalized_expressiondata.csv")

# PCA for Outlier Detection
library(ggplot2)
data <- plotPCA(vsd, intgroup=c("condition", "time"), returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
pdf("Rsv3_PCA_OutlierDetection.pdf")
ggplot(data, aes(PC1, PC2, color=condition, shape=time)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) #+
#  geom_text(aes_string(label = "name"), color = "black") # this line labels the sample 
dev.off()



#~~~~~~~~~~~~~ Get DESeq2 differential expression results ~~~~~~~~~~~~#
# Run the DESeq pipeline
ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ condition + time)

## Get differential expression results
resTC <- results(ddsTC)
resTC$symbol <- mcols(ddsTC)$symbol 



#~~~~~~~~~~~~~~~~~~~~~~ Genotype-Stage Contrasts ~~~~~~~~~~~~~~~~~~~~~#
# Note, I have created a separate directory for saving my DESeq2 contrast data.
# R = L29 (resistant) and S = Williams82 (susceptible)

# R vs S at 0 hours
RS0<-results(ddsTC, name="condition_resistant_vs_susceptible", test="Wald")
write.csv(RS0, file="/deseq2/ContrastResults/RS0.csv")

# R vs S at 2 hours
RS2<-results(ddsTC, contrast=list(c("condition_resistant_vs_susceptible","conditionresistant.time2h")), test="Wald")
write.csv(RS2, file="RS2.csv")

# R vs S at 4 hours
RS4<-results(ddsTC, contrast=list(c("condition_resistant_vs_susceptible","conditionresistant.time4h")), test="Wald")
write.csv(RS4, file="RS4.csv")

# R vs S at 6 hours
RS6<-results(ddsTC, contrast=list(c("condition_resistant_vs_susceptible","conditionresistant.time6h")), test="Wald")
write.csv(RS6, file="RS6.csv")

# R vs S at 8 hours
RS8<-results(ddsTC, contrast=list(c("condition_resistant_vs_susceptible","conditionresistant.time8h")), test="Wald")
write.csv(RS8, file="RS8.csv")
