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

e.set.t[1:4,1:4]
dim(e.set.t)
length(status$most_general)
label <- as.character(status$most_general)
class(label)
unique(label)

idx <- (label == 'bacterial' | label =='viral')

idx <- (label == 'bacterial' | label =='viral' |
          label == 'greyb' | label =='greyv')

idx <- (label == 'bacterial' | label =='viral' |
          label == 'greyb' | label =='greyv' |
          label == 'OD' | label == 'HC')

e.set <- e.set.t[idx,]
label <- label[idx]
# rm(e.set.t)
dim(e.set)
length(label)

### DGE
bct <- e.set[label == 'bacterial',]
vrl <- e.set[label == 'viral',]
dim(bct)
dim(vrl)

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
pvalue_cutoff = 0.00000000001
abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)


# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(e.set[,filter_by_fold])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(e.set[,filter_by_pvalue])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue
sum(filter_combined)

filtered = e.set[,filter_combined]
dim(filtered)

e.set[,filter_combined][1:10,1:5]
e.set.f <- e.set[,filter_combined]


### CLUSTERING

k2$cluster <- as.factor(k2$cluster)
colnames(e.set.f)[1]
colnames(e.set.f)[2]

gene_2360348<-as.numeric(e.set.f[,colnames(e.set.f)[1]])
gene_7610053<-as.numeric(e.set.f[,colnames(e.set.f)[2]])
ggplot(e.set.f, aes(gene_2360348, gene_7610053, color = k2$cluster, shape=label)) + geom_point()

k1 <- kmeans(e.set.f, centers = 1, nstart = 25)
k2 <- kmeans(e.set.f, centers = 2, nstart = 25)
k3 <- kmeans(e.set.f, centers = 3, nstart = 25)
k4 <- kmeans(e.set.f, centers = 4, nstart = 25)

# plots to compare
p1 <- fviz_cluster(k1, geom = "point",  data = e.set.f) + ggtitle("k = 1")
p2 <- fviz_cluster(k2, geom = "point",  data = e.set.f) + ggtitle("k = 2")
p3 <- fviz_cluster(k3, geom = "point",  data = e.set.f) + ggtitle("k = 3")
p4 <- fviz_cluster(k4, geom = "point",  data = e.set.f) + ggtitle("k = 4")
grid.arrange(p1, p2, p3, p4, nrow = 2)


set.seed(3)
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(e.set.f, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

fviz_nbclust(e.set.f, kmeans, method = "wss")


# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(e.set.f, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(e.set.f))
  mean(ss[, 3])
}

# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15

# extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")

fviz_nbclust(e.set.f, kmeans, method = "silhouette")

gap_stat <- clusGap(e.set.f, FUN = kmeans, nstart = 25,
                    K.max = 10, B = 50)
fviz_gap_stat(gap_stat)


k2 <- kmeans(e.set.f, centers = 2, nstart = 25)
str(k2)
# k2
table(k2$cluster, label)

k3 <- kmeans(e.set.f, centers = 3, nstart = 25)
str(k3)
# k2
table(k3$cluster, label)

# fviz_cluster(k2, data = e.set.f, geom = c("point"))


# end