# This script extracts significantly differentially expressed genes (DEGs) between L29 and Williams82 at each time point from contrasts made with DESeq2. 
# It then compiles DEGs from contrasts while excluding those found at 0 hpi (control) and then merges DEGs to their VST normalized expressions (done in deseq2.R).
# Finally, script identifies transcription factor (TF) genes and then merges them to their VST expressions.

setwd("/path/to/deseq2/contrast/results/directory/")

all_files = list.files(pattern="*.csv")   # read in DESeq2 contrast files. Mine are in directory of their own.

# Extract significantly different DEGs.
pco=.05  # FDR adjusted p-value of .05
bmco=10  # base mean of 10
fcuco=1  # log2 fold change of 1
fcdco=-1 # log2 fold change of -1

# This loops through each contrast file to filter for DEGs defined by the parameters above. Those with FDR adjusted p-value less than 0.05, base mean above 10, and log2 fold change less than or greater than 1.
for (file in all_files){
  print(file)
  filename <- tools::file_path_sans_ext(file)
  data <-read.csv(file)
  sigGenes <- subset(data, padj<=pco & baseMean>=bmco & 
                       (log2FoldChange >= fcuco | log2FoldChange <= fcdco))
  # exporting results
  summary_sigGenes <- capture.output(summary(sigGenes))
  # Summary files
  cat(summary_sigGenes, file= paste0("path/to/deseq2/DEGs/summary/directory/",filename,"_DEGs_summary.txt"), sep = "\n")
  # Full files
  write.csv(as.data.frame(sigGenes), file = paste0("/path/to/deseq2/DEGs/directory/",filename,"_DEGs.csv"))
}


#----------------------------------------------------------------------------------------#
# At the point, compile DEGs, exclude those at 0 hpi, and merge to VST expression data. 

setwd("/path/to/deseq2/DEGs/directory/")

# Read in all DEG files and name them by file name
all_DEG_files <- list.files(pattern="*_DEGs.csv")  # read in all DEG files for each contrast.
for (file in all_DEG_files){
  filename <- tools::file_path_sans_ext(file)
  data <-read.csv(file)
  assign(paste(filename), data)
  print(filename)
  print(nrow(data))
}

# Get DEG gene IDs only
RS0_DEGs0 <- as.vector(RS0_DEGs$X)    # 213 DEGs between L29 and Williams82 at 0 hpi. These genes are control; remove them from compiled DEG list
RS2_DEGs0 <- as.vector(RS2_DEGs$X)
RS4_DEGs0 <- as.vector(RS4_DEGs$X)
RS6_DEGs0 <- as.vector(RS6_DEGs$X)
RS8_DEGs0 <- as.vector(RS8_DEGs$X)

# Compile DEGs from 2 to 8 hpi
degs2to8 <- data.frame(c(RS2_DEGs0, RS4_DEGs0, RS6_DEGs0, RS8_DEGs0))
nrow(degs2to8)  #1261 DEGs including duplicate instances
degs2to8 <- degs2to8[!duplicated(degs2to8),]
length(degs2to8) # 1164 unique DEGs 

# Give genes identifier to determine whether they came from 2to8 hpi set or 0 hpi set
times2to8 <- rep("2to8", each=1164)
degs2to8 <- data.frame("X" = degs2to8, "time" = times2to8)
time0 <- rep("0hr", each = 213)
degs0hr <- data.frame("X" = RS0_DEGs0, "time" = time0)

# Combine 0hr and 2to8 hpi DEGs
degs0to8 <- rbind(degs2to8, degs0hr)
nrow(degs0to8)  # 1377 DEGs from 0-8 hpi including duplicate instances

# Remove gene entirely if it's duplicated. That is, if gene appears in 2to8 hpi and 0 hpi, remove it entirely since 0 hpi is control.
degs_0to8 <- degs0to8[!(duplicated(degs0to8$X) | duplicated(degs0to8$X, fromLast = TRUE)), ]
nrow(degs_0to8)   # 1305 DEGs from 2to8hr only and 0hr only

# Remove DEGs that are DEGs at 0 hpi only by finding row index of 0hr DEGs. 
which(degs_0to8$time == "0hr")
degs0to8[1128,]
degs_0to8[1129,]
# Getting DEGs specific to 2 to 8 hpi only
DEGs_excluding0hr <- degs_0to8[1:1128,]         # 1128 DEGs specific to 2 to 8 hpi; USE THESE


### Merge DEGs with VST normalized expression data from DESeq2.
expr <- read.csv("path/to/deseq2_VSTnormalized_expressiondata.csv")
merged <- merge(x = DEGs_excluding0hr, y = expr, by.x = c("X"), by.y = c("X"))


### Average Replicates
# Preparing to average replicates
mergedt <- t(merged)
colnames(mergedt) <- mergedt[1, ] # the first row (gene names) will be the header
mergedt <- mergedt[c(-1, -2), ]   # remove first (gene id) and second (time identifier) rows
# Convert dataframes to numeric matrices
mergedt <- as.matrix(mergedt)
mergedt <- apply(mergedt, 2, as.numeric) 
# Average expression of replicates
n=2    # two biological replicates
merged <- t(aggregate(mergedt,list(rep(1:(nrow(mergedt)%/%n+1),each=n,len=nrow(mergedt))),mean)[-1])
colnames(merged) <- c("L29.0hpi", "L29.2hpi", "L29.4hpi", "L29.6hpi", "L29.8hpi",
                      "W82.0hpi", "W82.2hpi", "W82.4hpi", "W82.6hpi", "W82.8hpi")
write.csv(merged, file = "/path/to/analysis/directory/DEGs_normalized_data.csv")


#----------------------------------------------------------------------------------------#
# Lastly, identify TF genes and merge them to their VST expressions.
# Note, you will need to download a list of soybean's TF genes. This can be obtained from the Plant Transcription Factor Database (PlantTFDB).
'''
PlantTFDB from:
Jin J, Tian F, Yang D-C, Meng Y-Q, Kong L, Luo J, et al. PlantTFDB 4.0: toward a central hub for transcription factors and regulatory interactions in plants. Nucleic Acids Research. 2017;45(D1):D1040–D5.
'''

tf <- read.table("/soybean/TF/list/from/PlantTFDB/download", header = T)

# Merge DEGs with Transcription Factor annotations
expr_tf <- merge(x = expr, y = tf, by.x = c("X"), by.y = c("Gene_ID"))
expr_tf <- subset(expr_tf, select = -c(TF_ID, Family))

# Remove duplicates (some genes have multiple TF isoforms)
expr_tf0 <- expr_tf[!duplicated(expr_tf$X), ]  # 131 TF genes
write.csv(expr_tf0, file = "/path/to/analysis/directory/TF_expression.csv")


