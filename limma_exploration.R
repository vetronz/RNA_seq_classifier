library(tidyverse)
library(plyr)
library(limma)
library(DESeq2)
library(Biobase)
library(cluster)
library(factoextra)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

dim(e.set)
e.set.t <- t(e.set)
dim(e.set.t)


idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'
sum(idx)


e.set.f <- e.set.t[idx,]
# label.f <- status[idx,c('most_general')]
# label <- as.character(label.f)
X <- as.matrix(t(e.set.f))
dim(X)

e.set.f <- as.data.frame(e.set.f)
e.set.f$label <- as.character(status[idx,c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

### DESIGN MATRIX
# class(label)
# label <- factor(label)

design <- model.matrix(~label + 0, data = e.set.f)
# design <- model.matrix(~label + sex + 0, data = e.set.f)
design[1:5,]
dim(design)

table(e.set.f$label, e.set.f$sex)

colSums(design)
colnames(design)<- c("bct","greyb","greyv", 'vrl')

contrast.matrix<- makeContrasts(c("bct-vrl"), levels=design)
contrast.matrix
# colnames(fit$coefficients)

dim(X)
dim(design)

fit <- lmFit(X, design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH", number=10)
top.hits <- topTable(fit2, adjust="BH", number=nrow(X))

results <- decideTests(fit2, method='global', p.value = 0.05, adjust.method = 'BH')
summary(results)
vennDiagram(results)

myTopHits <- as_tibble(top.hits, rownames = "geneSymbol")
myTopHits[1:5,1:5]

ggplot(myTopHits, aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneSymbol))) +
  geom_point(size=2) +
  ylim(min(-log10(myTopHits$adj.P.Val)), max(-log10(myTopHits$adj.P.Val))) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1)


rownames(top.hits)[1:5]
dim(top.hits)

# limmaTop <- match(a, pval.df[with(pval.df, order(pvalue)),]['transcript'][[1]])
rownames(top.hits)[6]

hit.idx <- match(rownames(top.hits)[1:100], colnames(e.set.f))
hit.idx

e.set.cut <- e.set.f[,hit.idx]

dim(e.set.cut)
e.set.cut[1:10,1:6]

# cluster centers
fviz_nbclust(e.set.cut, kmeans, method = "wss")
fviz_nbclust(e.set.cut, kmeans, method = "silhouette")

gap_stat <- clusGap(e.set.cut, FUN = kmeans, nstart = 25,
                    K.max = 20, B = 50)
fviz_gap_stat(gap_stat)


# K2 Clustering
k2 <- kmeans(e.set.cut, centers = 2, nstart = 25)
str(k2)
k2$cluster <- as.factor(k2$cluster)

table(k2$cluster, e.set.f$label)

chisq.test(table(k2$cluster, e.set.f$sex), correct = TRUE)

status[idx, c('my_category_2', 'most_general')][k2$cluster == 1,]
status[idx, c('my_category_2', 'most_general')]

clus1 <- status[idx, c('my_category_2', 'most_general', 'more_general',
              'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP')][k2$cluster == 1,]

clus2 <- status[idx, c('my_category_2', 'most_general', 'more_general',
              'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP')][k2$cluster == 2,]


clus1.df <- data.frame(cluster='cluster_1', wbc=clus1$WBC, crp=clus1$array.contemporary.CRP, age=clus1$Age..months., sex=clus1$Sex, label=clus1$most_general)
clus2.df <- data.frame(cluster='cluster_2', wbc=clus2$WBC, crp=clus2$array.contemporary.CRP, age=clus2$Age..months., sex=clus2$Sex, label=clus2$most_general)
df <- rbind(clus1.df, clus2.df)

head(df)
dim(df)

ggplot(df, aes(cluster, fill=sex)) + geom_bar(position="dodge")+
  xlab('Cluster') +
  ylab('Counts') +
  ggtitle("Bar Plot of Gender Split Between Clusters")



# subgroup analysis
e.set.f$label[k2$cluster == 1]
e.set.f$sex[k2$cluster == 1]

clus1.tab <- table(e.set.f$label[k2$cluster == 1], e.set.f$sex[k2$cluster == 1])
addmargins(clus1.tab)

chisq.test(clus1.tab, correct = FALSE)

clus2.tab <- table(e.set.f$label[k2$cluster == 2], e.set.f$sex[k2$cluster == 2])
addmargins(clus2.tab)

chisq.test(clus2.tab, correct = FALSE)


ggplot(clus1.df, aes(label, fill=sex)) + geom_bar(position="dodge")+
  xlab('Diagnosis') +
  ylab('Counts') +
  ggtitle("Cluster One Diagnosis by Gender")

ggplot(clus2.df, aes(label, fill=sex)) + geom_bar(position="dodge")+
  xlab('Diagnosis') +
  ylab('Counts') +
  ggtitle("Cluster Two Diagnosis by Gender")



# PCA
e.set.pca <- prcomp(e.set.cut, scale = TRUE)
# summary(e.set.pca)
plot(e.set.pca, type = 'l')

pair1 <- as.data.frame(e.set.pca$x[,1:2])
pair2 <- as.data.frame(e.set.pca$x[,3:4])

fviz_cluster(k2, geom = c("point"),  data = e.set.cut, axes = c(1,2)) +
  ggtitle("PCA Cluster Assignment k = 2")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k2$cluster, shape=e.set.f$label), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of First Two Principal Components")

ggplot(pair2, aes(PC3, PC4)) + geom_point(aes(color=k2$cluster, shape=e.set.f$label), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("Cluster Assignment of Third-Fourth Principal Components")











# expression vs variability
hist(fit2$Amean)
plotSA(fit2)
# rule of thumb is to keep genes with exp > 5-10 in > k arays
exprs.thresh <- 4
keep <- fit2$Amean > exprs.thresh

fit2 <- eBayes(fit2[keep,], trend=TRUE)
plotSA(fit2)

results <- decideTests(fit2, method="global", adjust.method="BH", p.value=0.05)
dim(results)



colnames(fit$coefficients)
tab.results <- topTable(fit, coef = colnames(fit$coefficients)[4], adjust.method = 'BH', number = ncol(e.set))
# tab.results <- topTable(fit2, adjust.method = 'BH', number = 5)
dim(tab.results)
head(tab.results)
class(tab.results)

p.threshold <- 0.05
tab.results$pval.threshold <- as.logical(tab.results$adj.P.Val < p.threshold)
head(tab.results)







# HAND CRAFTED CHI SQR
# addmargins(table(k2$cluster, e.set.f$sex))
# 
# f.p <- 76/(76+115)
# m.p <- 1 - f.p
# exp.f.1 <- rowSums(table(k2$cluster, e.set.f$sex))[[1]] * f.p
# exp.f.2 <- rowSums(table(k2$cluster, e.set.f$sex))[[2]] * f.p
# exp.m.1 <- rowSums(table(k2$cluster, e.set.f$sex))[[1]] * m.p
# exp.m.2 <- rowSums(table(k2$cluster, e.set.f$sex))[[2]] * m.p
# sum((42 - exp.f.1)^2 / exp.f.1, (34 - exp.f.2)^2 / exp.f.2, (49 - exp.m.1)^2 / exp.m.1 ,(66 - exp.m.2)^2 / exp.m.2)

# tab = matrix(c(menYes,womenYes,menNo,womenNo),2,2)

# f.p <- 42/91
# m.p <- 49/91
# 
# b.f <- 45* f.p
# b.m <- 45 * m.p
# gryb.f <- 29 * f.p
# gryb.m <- 29 * m.p
# vrl.f <- 17 * f.p
# vrl.m <- 17 * m.p
# 
# sum(
#   (26-b.f)^2/b.f, 
#   (19-b.m)^2/b.m, 
#   (13-gryb.f)^2/gryb.f,
#   (16-gryb.m)^2/gryb.m, 
#   (3-vrl.f)^2/vrl.f,
#   (14-vrl.m)^2/vrl.m
# )


## multiple groups
set.seed(42)
sd<- 0.3*sqrt(4/rchisq(100,df=4))
sd
y2<- matrix(rnorm(100*9,sd=sd),100,9)
dim(y2)

rownames(y2)<- paste("Gene",1:100)

y2[1:2,4:6]<- y2[1:2,4:6] + 2
y2[1:2,7:9]<- y2[1:2,7:9] + 2

head(y2)
dim(y2)
class(y2)

y2.t <- as.data.frame(t(y2))
dim(y2.t)
y2.t[,98:101]

y2.t$label <- c('A','A','A','B','B','B','C','C','C')

f<- factor(c(rep("A",3),rep("B",3),rep("C",3)),levels=c("A","B","C"))
f


design <- model.matrix(~label, data = y2.t)
design[,1][4:9]=0
# design<- model.matrix(~0+f)

colnames(design)<- c("A","B","C")

# contrast.matrix<- makeContrasts("B-A", "C-A", "C-B", levels=design)
contrast.matrix<- makeContrasts("C-A", levels=design)
contrast.matrix

fit1<- lmFit(y2,design)
# colnames(fit1$coefficients)
fit2<- contrasts.fit(fit1, contrast.matrix)
fit2<- eBayes(fit2)
topTable(fit2, adjust="BH")





# end