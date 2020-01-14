# This script combines all TF-module interactions predicted by each network inference method and converts them to TF-DEG interactions
# Then it filters for interactions that were detected by four out of the five inference methods. No interactions were found by all five methods.
# Finally, at the end of this script are notes on computational methods for validating the network predicted interactions.

setwd("/path/to/analysis/directory/")

# Import network files from each inference method
anet <- read.csv("aracne_network.csv")
cnet <- read.csv("clr_network.csv")
lnet <- read.csv("lars_network.csv")
pnet <- read.csv("partialcorrelation_network.csv")
rnet <- read.csv("randomforest_network.csv")

# Merge regulator and target columns so have unique key for counting duplicates
anet$regulatortarget <- paste(anet$row, anet$col)
cnet$regulatortarget <- paste(cnet$row, cnet$col)
lnet$regulatortarget <- paste(lnet$regulator, lnet$target)
pnet$regulatortarget <- paste(pnet$node2_ID, pnet$node1_ID)
rnet$regulatortarget <- paste(rnet$regulatoryGene, rnet$targetGene)

# Combine all interaction columns
mergnet <- c(anet$regulatortarget, cnet$regulatortarget, lnet$regulatortarget, pnet$regulatortarget, rnet$regulatortarget)
mergnet <- data.frame("regulatortarget" = mergnet)

# Count number of duplicates (interactions found in more than one method)
intcounts <- table(mergnet)               # No interactions predicted by all 5 methods
intcounts <- data.frame(intcounts)
intcounts <- intcounts[order(intcounts$Freq, decreasing = T),]
unqint <- mergnet[!duplicated(mergnet),]    # 654 unique TF-module interactions

# Split and rename columns
library(stringr)
unqint <- data.frame(str_split_fixed(int_4methods$mergnet, " ", 2))  # For all 654 unique interactions, split interaction column into TF and module
colnames(unqint) <- c("regulator", "target")


### Match all TF genes to all DEGs belonging in their target modules
# Subset by module assignment
modlist <- list()
for (i in c(1:5)){
  module = paste("module",i, sep = "")
  print(module)
  moddat <- unqint[which(unqint$target == module),]
  modlist[[module]] <- moddat
}

# Subset gene module assignment by module number
degs <- read.csv("DEG_modules.csv")
degs <- data.frame("gene" = degs$X, "module" = degs$module)

moddeg <- list()
for (i in c(1:5)){
  print(i)
  moddat <- degs[degs$module == i,]
  moddeg[[i]] <- moddat
}

# Loop matches TFs to DEGs belonging in target module 
allmodules <- list()
for (i in c(1:5)){         
  tfdf = modlist[[i]]    
  moddf = moddeg[[i]]
  module = paste("module",i, sep = "")
  singlemodule <- list()
  for (row in 1:nrow(tfdf)){  
    reg = tfdf$regulator         
    reg = levels(droplevels(reg))
    tf = reg[row]
    target = degdf$gene
    tfrep = rep(tf, length(target))
    reg_tar = data.frame("TF" = tfrep, "target" = target)
    singlemodule[[row]] <- reg_tar
  }
  singlemodule <- do.call(rbind, singlemodule)
  allmodules[[module]] <- singlemodule
}

write.csv(allmodules, file = "TF-DEG_interactions_all.csv") # All TF-DEG interactions without filtering


#---------------------------------------------------------------------------------#
### Filter for interactions detected by 4 out of 5 inference methods. None were found by all five methods.

# Get interactions predicted by 4 methods; 
int_4methods <- intcounts[grep("4", intcounts$Freq),]   # 56 TF-module interactions predicted by 4 out of 5 methods 
tfmod_4methods <- data.frame(str_split_fixed(int_4methods$mergnet, " ", 2))  # split columns
colnames(tfmod_4methods) <- c("regulator", "target")

### Match all TF genes to all DEGs belonging in their target modules
# Subset by module assignment
modlist <- list()
for (i in c(1:5)){
  module = paste("module",i, sep = "")
  print(module)
  moddat <- tfmod_4methods[which(tfmod_4methods$target == module),]
  modlist[[module]] <- moddat
}

# Subset gene module assignment by module number
moddeg <- list()
for (i in c(1:5)){
  print(i)
  moddat <- degs[degs$module == i,]
  moddeg[[i]] <- moddat
}

# Loop matches TFs to DEGs belonging in target module 
allmodules <- list()
for (i in c(1:5)){         
  tfdf = modlist[[i]]    
  moddf = moddeg[[i]]
  module = paste("module",i, sep = "")
  singlemodule <- list()
  for (row in 1:nrow(tfdf)){  
    reg = tfdf$regulator         
    reg = levels(droplevels(reg))
    tf = reg[row]
    target = degdf$gene
    tfrep = rep(tf, length(target))
    reg_tar = data.frame("TF" = tfrep, "target" = target)
    singlemodule[[row]] <- reg_tar
  }
  singlemodule <- do.call(rbind, singlemodule)
  allmodules[[module]] <- singlemodule
}

write.csv(allmodules, file = "TF-DEG_interactions_4methods.csv") # TF-DEG interactions predicted by 4 out of 5 inference methods

'''
Now that you have high confidence predictions of TF-DEG interactions (those found by 4 out of 5 interactions), 
you can be validate them by comparison to published Arabidopsis interactions found by DNA affinity purification sequencing 
(DAP-seq)*. This will require you to idenitfy Arabidopsis homologs for soybean genes. We did this with BLAST using an 
E-value of 1e-5. The predicted interactions can be further validated by motif sequence analysis using Meme Suite**. 
You will need to download soybean promoter region sequences for all of the target genes. This can be done on Soybase***, 
where we defined the promoter regions as the 1000 base pairs flanking the 5 prime end of a gene. Keep these promoter sequences 
with their genes respective module assignments. Submit each set of module promoter sequences to the motif discovery tool, 
MEME****, to identify enriched motif sequences in a module. Here motif width was set from 6 to 12 and a maximum of 6 motifs 
were allowed to be found. The enriched motif sequences discovered can then me submitted directly to the TomTom tool*****, 
where they can be compared to motif sequences recognized by specific TFs. For TomTom, select the Arabidopsis DAP-seq motif 
database from OMalley 2016. The results from this will allow you to identify which motif sequences are putatively recognized 
by which types of TFs, which can then be compared back to the TFs predicted to regulate genes with those motif sequences.

Thus, with the regulatory interactions detected by 4 out of five inference methods and validated by both published Arabidopsis 
interactions and motif sequence analysis, a set of high confidence putative regulatory interactions can be obtained to aid
in hypothesis-driven studies.

*O’Malley RC, Huang S-sC, Song L, Lewsey MG, Bartlett A, Nery JR, et al. Cistrome and epicistrome features shape the regulatory DNA landscape. Cell. 2016;165(5):1280-92
**Bailey TL, Boden M, Buske FA, Frith M, Grant CE, Clementi L, et al. MEME SUITE: tools for motif discovery and searching. Nucleic acids research. 2009;37(suppl_2):W202-W8.
***Grant D, Nelson RT, Cannon SB, Shoemaker RC. SoyBase, the USDA-ARS soybean genetics and genomics database. Nucleic acids research. 2009;38(suppl_1):D843-D6
****Bailey TL, Elkan C. Fitting a mixture model by expectation maximization to discover motifs in bipolymers.  Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology; Menlo Park, CA: AAAl Press; 1994. p. 28-36.
*****Gupta S, Stamatoyannopoulos JA, Bailey TL, Noble WS. Quantifying similarity between motifs. Genome biology. 2007;8(2):R24.
'''
