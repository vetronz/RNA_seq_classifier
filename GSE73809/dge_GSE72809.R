# Read the data into R
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

vrl.index <- vrl[!is.na(vrl)]


# bct_df <- gset_df[, bct]
# vrl_df <- gset_df[, vrl]
bct_df <- gset.sig[, bct]
vrl_df <- gset.sig[, vrl]
dim(bct_df)
dim(vrl_df)

## SUBSET THE GSET TO GET BACT VIRAL DF
index_binary <- c(bct, vrl)
gset_binary <- gset_df[,index_binary]
gset_binary_t <- as.data.frame(t(gset_binary))

gset_binary[1:5,1:5]
gset_binary_t[1:5,1:5]

dim(gset_binary)
dim(gset_binary_t)

## CONSTRUCTS GROUND TRUTH LIST 
truth <- character(0)
index_binary <- list(bct, vrl)
for (i in 1:length(index_binary)){
  for (j in 1:length(index_binary[[i]])){
    print(paste(i, j))
    truth = append(truth,i)
  }
}
truth
gset_binary_t$truth <- truth
gset_binary_t[, 47320:47324]

## SET TRUTH TO FACTOR
# sapply(gset_binary_t[, 47320:47324], class)
gset_binary_t$truth <- as.factor(gset_binary_t$truth)
sapply(gset_binary_t[, 47320:47324], class)


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

## CONSTRUCTS GROUND TRUTH LIST 
truth <- character(0)
index_binary <- list(bct, vrl)
for (i in 1:length(index_binary)){
  for (j in 1:length(index_binary[[i]])){
    print(paste(i, j))
    truth = append(truth,i)
  }
}

gset.sig$truth <- truth
gset.sig[,34:36]

## SET TRUTH TO FACTOR
# sapply(gset_binary_t[, 47320:47324], class)
gset.sig$truth <- as.factor(gset.sig$truth)
sapply(gset.sig[, 30:36], class)



## TRAIN TEST SPLIT
set.seed(123)
index <- sample(1:nrow(gset.sig), round(nrow(gset.sig) * 0.7))
index
train_1 <- gset.sig[index, ]
test_1  <- gset.sig[-index, ]

dim(train_1)
dim(test_1)


bac_prop <- summary(gset.sig$truth)[[1]]/nrow(gset.sig)
train_bac_prop <- summary(train_1$truth)[[1]]/dim(train_1)[1]
test_bac_prop <- summary(test_1$truth)[[1]]/dim(test_1)[1]

print(paste('bacterial proportion:', bac_prop))
print(paste('bacterial train proportion:', train_bac_prop))
print(paste('bacterial test proportion:', test_bac_prop))



## NORMALIZATION
# identify only the predictor variables
features <- setdiff(names(train_1), 'truth')

# pre-process estimation based on training features
pre_process <- preProcess(
  x      = train_1[, features],
  method = c("center", "scale")    
)

# apply to both training & test
train_x <- predict(pre_process, train_1[, features])
test_x  <- predict(pre_process, test_1[, features])

train_x[1:5, 47320:47323]
test_x[1:5, 1:5]
# remove(train_x, test_x)

## LOAD DATA FROM rds object
# saveRDS(train_x, file = 'train_x_GSE72809')
# saveRDS(test_x, file = 'test_x_GSE72809')
train_x <- readRDS(file = 'train_x_GSE72809')
test_x <- readRDS(file = 'test_x_GSE72809')

train_x$truth <- truth

model1 <- glm(truth ~ ILMN_1657871, family = "binomial", data = train_1)
summary(model1)

pred <- predict(model1, test_1, type='response')

pred
pred.class <- ifelse(pred > 0.5, 1, 0)
pred.class

test_1


