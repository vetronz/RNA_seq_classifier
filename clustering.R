library(cluster)
library(factoextra)
library(gridExtra)
library(tidyverse)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')

dim(e.set) # megaExp data
dim(status) # megaExp labels

rm(e.set.i, status.iris, targets, targets.iris)
# dim(e.set.i) # iris data
# dim(status.iris) # iris labels

# transpose
e.set.t <- t(e.set)
# e.set.t[1:4,1:4]

label <- as.character(status$most_general)

### DGE
bct <- e.set.t[label == 'bacterial',]
vrl <- e.set.t[label == 'viral',]
dim(bct)
dim(vrl)

e.set <- rbind(bct, vrl)
dim(e.set)

bct.mean <- apply(bct, 2, mean)
vrl.mean <- apply(vrl, 2, mean)

head(bct.mean)
head(vrl.mean)

# Just get the maximum of all the means
limit = max(bct.mean, vrl.mean)

# Scatter plot
plot(vrl.mean ~ bct.mean, xlab = "bct", ylab = "vrl",
     main = "Bct vs Vrl - Scatter", xlim = c(3.5, limit), ylim = c(0, limit))
cor(bct.mean, vrl.mean)

# Compute fold-change (biological significance)
fold = bct.mean - vrl.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : ncol(bct)) { # For each gene : 
  x = bct[,i] 
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
pvalue_cutoff = 1e-14

abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)


# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(e.set[,filter_by_fold])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(e.set[,filter_by_pvalue])
# e.set[,filter_by_pvalue]

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue
filtered = e.set[,filter_combined]
dim(filtered)



############################## CLUSTERING ##############################

unique(label)

idx <- (label == 'bacterial' | label =='viral' |
          label == 'greyb' | label =='greyv' | label == 'greyu')

e.set <- data.frame(e.set.t[idx, filter_combined])
label.i <- label[idx]
status.i <- status[idx,]

# e.set <- scale(e.set) # worse performance with scaling
fviz_nbclust(e.set, kmeans, method = "wss")
fviz_nbclust(e.set, kmeans, method = "silhouette")

gap_stat <- clusGap(e.set, FUN = kmeans, nstart = 25,
                    K.max = 20, B = 50)
fviz_gap_stat(gap_stat)

# distance <- get_dist(e.set)
# fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))

## K2
dim(e.set)

k2 <- kmeans(e.set, centers = 2, nstart = 25)
str(k2)
attributes(k2)
dim(status)


k2$cluster <- as.factor(k2$cluster)

sig.1 <- as.numeric(e.set[,colnames(e.set)[1]])
sig.2 <- as.numeric(e.set[,colnames(e.set)[2]])

ggplot(e.set, aes(sig.1, sig.2, color = k2$cluster, shape=label.i)) +
  geom_point(size=2) +
  xlab(paste('transcript',colnames(e.set))[1]) +
  ylab(paste('transcript',colnames(e.set))[2]) +
  ggtitle("Pairwise Scatter Plot of Most Significantly
          Differentially Expressed Transcripts")

e.set.pca <- prcomp(e.set, scale = TRUE)
summary(e.set.pca)
plot(e.set.pca, type = 'l')

pair1 <- as.data.frame(e.set.pca$x[,1:2])
pair2 <- as.data.frame(e.set.pca$x[,3:4])

fviz_cluster(k2, geom = c("point"),  data = e.set, axes = c(1,2)) +
  ggtitle("PCA Cluster Assignment k = 2")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k2$cluster, shape=label.i), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of First Two Principal Components")

# ggplot(pair2, aes(PC3, PC4)) + geom_point(aes(color=k2$cluster, shape=label.i), size=2) +
#   xlab("Third Principal Component") +
#   ylab("Fourth Principal Component") +
#   ggtitle("Cluster Assignment of Third and Fourth Principal Components")

table(k2$cluster, label.i)

attributes(status)$names
clus1 <- status.i[k2$cluster == 1,c('my_category_2', 'most_general', 'more_general',
                                    'Age..months.', 'Sex', 'ethnicity', 'WBC', 'array.contemporary.CRP')]
clus2 <- status.i[k2$cluster == 2,c('my_category_2', 'most_general', 'more_general',
                                    'Age..months.', 'Sex', 'ethnicity', 'WBC', 'array.contemporary.CRP')]
dim(clus1)
dim(clus2)
summary(clus1)
summary(clus2)

# ix the bcts sneaking into vrl cluster
clus1[clus1$most_general == 'bacterial',]
# https://flowingdata.com/2012/05/15/how-to-visualize-and-compare-distributions/


## K3
k3 <- kmeans(e.set, centers = 3, nstart = 25)
str(k3)
k3$cluster <- as.factor(k3$cluster)

ggplot(as.data.frame(e.set), aes(sig.1, sig.2, color = k3$cluster, shape=label.i)) +
  geom_point(size=2) +
  xlab(paste('transcript',colnames(e.set))[1]) +
  ylab(paste('transcript',colnames(e.set))[2]) +
  ggtitle(
    "Pairwise Scatter Plot of Most Significantly
    Differentially Expressed Transcripts")

table(k3$cluster, label.i)

fviz_cluster(k3, geom = c("point"),  data = e.set, axes = c(1,2)) +
  ggtitle("PCA Cluster Assignment k = 3")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k3$cluster, shape=label.i), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of First Two Principal Components")


## K4
k4 <- kmeans(e.set, centers = 4, nstart = 25)
str(k4)
k4$cluster <- as.factor(k4$cluster)

ggplot(as.data.frame(e.set), aes(sig.1, sig.2, color = k4$cluster, shape=label.i), size=2) +
  geom_point(size=2) +
  xlab(paste('transcript',colnames(e.set))[1]) +
  ylab(paste('transcript',colnames(e.set))[2]) +
  ggtitle(
    "Pairwise Scatter Plot of Most Significantly
    Differentially Expressed Transcripts")

table(k4$cluster, label.i)

# fviz_cluster(k4, geom = c("point"),  data = e.set, axes = c(1,2)) +
  # ggtitle("PCA Cluster Assignment K = 4")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k4$cluster, shape=label.i),size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of First Two Principal Components")



## K10
k10 <- kmeans(e.set, centers = 10, nstart = 25)
# str(k10)
k10$cluster <- as.factor(k10$cluster)

ggplot(as.data.frame(e.set), aes(sig.1, sig.2, color = k10$cluster, shape=label.i), size=2) +
  geom_point(size=2) +
  xlab(paste('transcript',colnames(e.set))[1]) +
  ylab(paste('transcript',colnames(e.set))[2]) +
  ggtitle(
    "Pairwise Scatter Plot of Most Significantly
    Differentially Expressed Transcripts")

table(k10$cluster, label.i)

fviz_cluster(k10, geom = c("point"),  data = e.set, axes = c(1,2)) +
  ggtitle("PCA Cluster Assignment K = 10")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k10$cluster, shape=label.i),size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of First Two Principal Components")






############################################################
## to do
# need to add in phenotypic breakdown to each of the clusters to see how
# the clusters are partitioning them

# please include all genes in an analysis.
# will need to use limma to adjust for confounders age, sex 










#############################################################
# hierarchical

dim(e.set)
e.set[1:5,1:5]

e.set.s <- scale(e.set)
apply(e.set.s, 2, mean)
apply(e.set.s, 2, var)

e.set.s[1:5,1:4]

cor(e.set.s)
cor(e.set.s, method = 'pearson')

dist(1-cor(e.set.s, method="pearson"))

hr <- hclust(as.dist(1-cor(e.set.s, method="pearson")), method = 'complete')
hc <- hclust(as.dist(1-cor(e.set.s, method="spearman")), method = 'complete')


library(gplots)
heatmap.2(as.matrix(e.set),
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none")



TreeC = as.dendrogram(hc, method="average")
plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")


TreeR = as.dendrogram(hr, method="average")
plot(TreeR,
     # leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")


hclusth1.5 = cutree(hr, h=1.5) #cut tree at height of 1.5
hclusth1.0 = cutree(hr, h=1.0) #cut tree at height of 1.0
hclusth0.5 = cutree(hr, h=0.5) #cut tree at height of 0.5

library(dendextend)
#plot the tree
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")

#add the three cluster vectors
the_bars <- cbind(hclusth0.5, hclusth1.0, hclusth1.5)
#this makes the bar
colored_bars(the_bars, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("h=0.5","h=1.0","h=1.5"),cex.rowLabels=0.7)
#this will add lines showing the cut heights
abline(h=1.5, lty = 2, col="grey")
abline(h=1.0, lty = 2, col="grey")
abline(h=0.5, lty = 2, col="grey")


hclustk4 = cutree(hr, k=4)
plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(hclustk4, TreeR, sort_by_labels_order = T, y_shift=-0.1, rowLabels = c("k=4"),cex.rowLabels=0.7)




##############################################################
# behind the shortcuts


# # function to compute total within-cluster sum of square 
# wss <- function(k) {
#   kmeans(e.set, k, nstart = 10 )$tot.withinss
# }
# 
# # Compute and plot wss for k = 1 to k = 15
# k.values <- 1:15
# 
# # extract wss for 2-15 clusters
# wss_values <- map_dbl(k.values, wss)
# 
# plot(k.values, wss_values,
#      type="b", pch = 19, frame = FALSE, 
#      xlab="Number of clusters K",
#      ylab="Total within-clusters sum of squares")


# function to compute average silhouette for k clusters
# avg_sil <- function(k) {
#   km.res <- kmeans(e.set, centers = k, nstart = 25)
#   ss <- silhouette(km.res$cluster, dist(e.set))
#   mean(ss[, 3])
# }
# 
# # Compute and plot wss for k = 2 to k = 15
# k.values <- 2:15
# 
# # extract avg silhouette for 2-15 clusters
# avg_sil_values <- map_dbl(k.values, avg_sil)
# 
# plot(k.values, avg_sil_values,
#      type = "b", pch = 19, frame = FALSE, 
#      xlab = "Number of clusters K",
#      ylab = "Average Silhouettes")
# end