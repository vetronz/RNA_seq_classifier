library(Biobase)
library(GEOquery)
library(ggplot2)
library(tidyr)
# library(rsample)
library(caret)
# library(h2o)
library(dplyr)
library(tidyverse)  # data manipulation and visualization
library(modelr)     # provides easy pipeline modeling functions
library(broom)

## GET DATA FROM GEO
# gset <- getGEO('GSE72809')
# head(gset)
# saveRDS(gset, file = 'gset_GSE72809')

## LOAD DATA FROM rds object
setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
getwd()
gset <-readRDS(file = "gset_GSE72809")
head(gset)
# exprs(gset_t[[1]])
n <- ncol(gset[[1]])


gset_df <- exprs(gset[[1]])
gset_df[1:5,1:5]

colnames(phenoData(gset[[1]]))
phenoData(gset[[1]])$'category:ch1'
class(phenoData(gset[[1]])$'category:ch1')



## ADD DEF BACT AND DEF VIRAL TO LISTS
bct.index <- c()
vrl.index <- c()
for (i in 1:n){
  if (phenoData(gset[[1]])$'category:ch1'[i] == 'Definite Bacterial'){
    bct.index[[i]] <- i
  }
  if (phenoData(gset[[1]])$'category:ch1'[i] == 'Definite Viral'){
    vrl.index[[i]] <- i
  }
}

vrl.index <- vrl.index[!is.na(vrl.index)]


# bct_df <- gset_df[, bct]
# vrl_df <- gset_df[, vrl]
bct_df <- gset.sig[, bct]
vrl_df <- gset.sig[, vrl]
dim(bct_df)
dim(vrl_df)


## DIFFERENTIAL EXPRESSION
bct.mean = apply(bct_df, 1, mean)
vrl.mean = apply(vrl_df, 1, mean)

# bct.mean[1:5]
# vrl.mean[1:5]

# Just get the maximum of all the means
limit = max(bct.mean, vrl.mean)

# Scatter plot
plot(bct.mean ~ vrl.mean, xlab = "WT", ylab = "KO",
     main = "GSE72809 - Scatter", xlim = c(0, limit), ylim = c(0, limit))
# Diagonal line
abline(0, 1, col = "red")


# Compute fold-change (biological significance)
# Difference between the means of the conditions
fold = bct.mean - vrl.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : nrow(gset_df)) { # For each gene : 
  x = bct_df[i,] # vector of gene expression for gene i in bct
  y = vrl_df[i,] # vector of gene expression for gene i in vrl
  
  # Compute t-test between the two conditions
  t = t.test(x, y)
  
  # Put the current p-value in the pvalues list
  pvalue[i] = t$p.value
  # Put the current t-statistic in the tstats list
  tstat[i] = t$statistic
}

# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot(fold, -log10(pvalue), main = "GSE72809 - Volcano")

fold_cutoff = 2
pvalue_cutoff = 0.01
abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)

# Screen for the genes that satisfy the filtering criteria

# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(gset_df[filter_by_fold, ])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(gset_df[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = gset_df[filter_combined,]
dim(filtered)
filtered[1:5, 1:5]

# Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "GSE72089 - Volcano #2")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")

sig.genes <- rownames(filtered)
gset.sig <- as.data.frame(t(gset_binary[sig.genes,]))
gset.sig[1:5, 1:5]
dim(gset.sig)