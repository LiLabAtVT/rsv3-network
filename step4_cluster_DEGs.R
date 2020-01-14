# This script clusters 1128 DEGs'variance stabilized transformed expressions using Gaussian finite mixture models by mclust.
# Then it calculates each module's mean expression and merges this with the TF gene expression data (TF_expression.csv created in extract_DEGs.R).

'''
mclust from:
Scrucca L, Fop M, Murphy TB, Raftery AE. mclust 5: clustering, classification and density estimation using Gaussian finite mixture models. The R Journal. 2016;8(1):289- 317.
'''

setwd("/path/to/analysis/directory/")

exprdat <- read.csv("DEGs_normalized_data.csv", row.names = 1, header = TRUE)
genes <- row.names(exprdat)

# Run Mclust
library(mclust)
set.seed(100)       # get reproducible results
myMclust <- Mclust(exprdat)  # get parameters for most optimal model

# Bayesian Inforamtion Criteria to determine best model
BIC <- mclustBIC(exprdat)
ICL <- mclustICL(exprdat)

# Summary gives top models based on BIC and ICL and best models of each
myMclust
summary(myMclust)

BICsummary <- summary(myMclust$BIC, data=exprdat)
BICsummary

summary(BIC)
summary(ICL)

# Add column of cluster groups called "module" and reorder columns
exprdat$module <- myMclust$classification
mods <- exprdat[,c("module","L29.0hpi","L29.2hpi","L29.4hpi","L29.6hpi", "L29.8hpi", "W82.0hpi", "W82.2hpi","W82.4hpi", "W82.6hpi","W82.8hpi")]   # reorder columns                          
write.csv(mods, file="DEG_modules.csv", row.names=TRUE, quote=FALSE)


### Loop calculates module medians for each module at each stage and then puts them together in dataframe. 
mods <- mods[order(mods$module),]  # put modules in order
umods <- unique(mods$module)   # unique modules
module_means=list()
for (column in names(mods[,3:ncol(mods)])){
  print(column)
  coldat=mods[,column]
  colmod_means = list()
  for (module in umods){
    print(module)
    mod_mean=mean(coldat[mods$module==module])
    colmod_means[[length(colmod_means)+1]] = mod_mean
  }
  colmod_means=as.data.frame(t(colmod_means))
  module_means[[column]] <- colmod_means
}
module_means <- do.call(rbind, module_means)
module_means <- t(module_means)
module_means <- apply(module_means, 2, as.numeric)  

# Add module assignment so know which row is which module:
row.names(module_means) <- paste("module", umods, sep = "") # make row names module number


### Create dataframe with module expression means and TF gene expressions.
tf <- read.csv("TF_expression.csv")
module_tf <- rbind(module_means, tf)
write.csv(module_tf, file = "cluster_TF_expression.csv")

