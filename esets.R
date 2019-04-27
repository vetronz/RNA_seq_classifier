# library(Biobase)
# library(GEOquery)
library(dplyr)
library(ggfortify)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()

rm(list=setdiff(ls(), "x"))
load('esets.RData')


dim(e.set) # megaExp data
dim(status) # megaExp labels

dim(e.set.i) # iris data
dim(status.iris) # iris labels

# attributes(status)
# status$Diagnosis
# status$category

# transpose
e.set.t <- t(e.set)
e.set.i.t <- t(e.set.i)

dim(e.set.t)
dim(e.set.i.t)

# rm(targets, targets.iris, e.set, e.set.i)
# e.set.t[1:5,(ncol(e.set.t)-2):ncol(e.set.t)]

# scale
e.set.s <- scale(e.set.t)
e.set.i.s <- scale(e.set.i.t)

# e.set.i.s <- as.data.frame(scale(e.set.i.t))
# e.set.s <- as.data.frame(scale(e.set.t))
# e.set.i.s <- as.data.frame(scale(e.set.i.t))

# apply(e.set.s, 2, var)
# apply(e.set.s, 2, mean)


# labels
label.i <- rownames(e.set.i.t)
label.s.i <- status.iris$My_code
common <- match(label.i, label.s.i)
# common

e.set.l <- as.character(status$most_general)
e.set.i.l <- as.character(status.iris$most_general[common]) # pass common index to extract common iris labels

e.set.df <- as.data.frame(e.set.s)
e.set.i.df <- as.data.frame(e.set.i.s)

e.set.df$label <- e.set.l
e.set.df$exp <- 0

e.set.i.df$label <- e.set.i.l
e.set.i.df$exp <- 1

e.set.df[1:5, (ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:5, (ncol(e.set.i.df)-4):ncol(e.set.i.df)]


# combine datasets
common.genes <-intersect(colnames(e.set.df), colnames(e.set.i.df))
length(common.genes) # 42004

dim(e.set.df[,common.genes])
dim(e.set.i.df[,common.genes])

e.set.df[1:10,(ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:10,(ncol(e.set.i.df)-4):ncol(e.set.i.df)]

e.set.c <- rbind(e.set.df[,common.genes], e.set.i.df[,common.genes])

dim(e.set.c)
e.set.c[1:10,(ncol(e.set.c)-5):ncol(e.set.c)]
e.set.c[(nrow(e.set.c)-10):nrow(e.set.c),(ncol(e.set.c)-5):ncol(e.set.c)]

e.set.r <- e.set.c[,-((ncol(e.set.c)-1):ncol(e.set.c))]
dim(e.set.r)

# View(status)
e.set.c$label == 'bacterial'
e.set.c$label == 'viral'

e.set.c$label <- ifelse(e.set.c$label == 'bacterial', 1, 0)


# PCA
gset.pca <- prcomp(e.set.c, scale = FALSE)
summary(gset.pca)
plot(gset.pca, type = 'l')
autoplot(gset.pca, data = e.set.c, colour = 'label')
autoplot(gset.pca, data = e.set.c, colour = 'exp')



## DGE
bct <- e.set.c[e.set.c$label == 'bacterial',]
vrl <- e.set.c[e.set.c$label == 'viral',]

# strip the labels for raw data
bct.r <- bct[-((ncol(bct)-1):ncol(bct))]
vrl.r <- vrl[-((ncol(vrl)-1):ncol(vrl))]
dim(bct.r)
dim(vrl.r)

# 
# dim(bct)
# bct[1:10,(ncol(bct)-5):ncol(bct)]
# bct[(nrow(bct)-10):nrow(bct),(ncol(bct)-5):ncol(bct)]
# 
# vrl[1:10,(ncol(vrl)-5):ncol(vrl)]
# vrl[(nrow(vrl)-10):nrow(vrl),(ncol(vrl)-5):ncol(vrl)]

bct.mean <- apply(bct.r, 2, mean)
vrl.mean <- apply(vrl.r, 2, mean)

head(bct.mean)
head(vrl.mean)

# Just get the maximum of all the means
limit = max(bct.mean, vrl.mean)

# Scatter plot
plot(vrl.mean ~ bct.mean, xlab = "bct", ylab = "vrl",
     main = "Bct vs Vrl - Scatter", xlim = c(0, limit), ylim = c(0, limit))

# Compute fold-change (biological significance)
fold = bct.mean - vrl.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : ncol(bct.r)) { # For each gene : 
  x = bct.r[,i] 
  y = vrl[,i]
  t = t.test(x, y)
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot(fold, -log10(pvalue), main = "Volcano #1")

fold_cutoff = 1.5
pvalue_cutoff = 0.0001
abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)


# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(e.set.r[filter_by_fold, ])

dim(e.set.r[filter_by_fold,])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(e.set.r[filter_by_pvalue, ])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = e.set.r[filter_combined,]
dim(filtered)

filtered[1:5,1:5]



# Let's generate the volcano plot again,
# highlighting the significantly differential expressed genes
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (fold[filter_combined], -log10(pvalue[filter_combined]),
        pch = 16, col = "red")

# Highlighting up-regulated in red and down-regulated in blue
plot(fold, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (fold[filter_combined & fold < 0],
        -log10(pvalue[filter_combined & fold < 0]),
        pch = 16, col = "red")
points (fold[filter_combined & fold > 0],
        -log10(pvalue[filter_combined & fold > 0]),
        pch = 16, col = "blue")



# saveRDS(, file = 'gset_GSE72809')
#asdf