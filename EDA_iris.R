library(tidyverse)
# library(plyr)
# library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
# require(reshape) # for melt()
# require(scales) # for percent
library(gridExtra)
library(dplyr)
library(plotly)

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

dim(e.set.i)

status.iris['most_general'] == 'bacterial'

# cols
dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', "#165bfc")
dx.cols.f <- c("#ed0404", "#fc5716",'#16fc31', '#16e1fc', '#165bfc', "#7a16fc", '#fc16f4')
sex.cols <- c('#fc1676', '#16acfc')

positions <- c('bacterial', 'greyb', 'greyv', 'viral')
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother')

# discrepancy in dimension of transcript and label matrix
dim(e.set.i)
dim(status.iris)

# extract the overlap
v <- intersect(colnames(e.set.i), status.iris$My_code)
length(v)

common <- match(v, status.iris$My_code)

status.i.idx <- status.iris[common,]
dim(status.i.idx)
dim(e.set.i)

as.character(colnames(e.set.i)) == as.character(status.i.idx$My_code) # double check correctly filtered

status.iris.c <- status.iris[common,]
idx <- status.iris.c['most_general'] == 'bacterial' |
  status.iris.c['most_general'] == 'viral' |
  status.iris.c['most_general'] == 'greyb' |
  status.iris.c['most_general'] == 'greyv' |
  status.iris.c['most_general'] == 'HC'
sum(idx)
length(idx)

dim(e.set.i)
dim(e.set.i[,idx])
dim(status.iris.c[idx,])

as.character(status.iris.c[idx,]$My_code)
colnames(e.set.i[,idx])

as.character(status.iris.c[idx,]$My_code) == colnames(e.set.i[,idx])



############ LIMMA ############
e.set.f <- as.data.frame(e.set.i[,idx])
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)
length(as.character(status.iris.c[idx, c('most_general')]))

e.set.f$label <- as.character(status.iris.c[idx, c('most_general')])
e.set.f$sex <- status.iris.c[idx,c('Sex')]
e.set.f$age <- status.iris.c[idx, c('Age')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

# mean/variance calculations
x_var <- apply(e.set.i[,idx], 1, var)
x_mean <- apply(e.set.i[,idx], 1, mean)
plot(log2(x_mean), log2(x_var), pch='.')
abline(v=log2(5), col='red')


dim(e.set.i[,idx])
e.set.i[,idx][which(x_mean > 5)]
dim(t(e.set.i[,idx][which(x_mean > 5),]))
X <- e.set.i[,idx][which(x_mean > 5),]



### DESIGN MATRIX
# site 
design <- model.matrix(~label + sex + age + 0, data = e.set.f)
colnames(design)<- c("bct","greyb","greyv", 'HC', 'vrl', 'sexF', 'sexM', 'age')
design <- design[,-6]

design[1:10,]
dim(design)
colSums(design)

contrast.matrix<- makeContrasts("bct-HC", 'vrl-HC', levels=design)
contrast.matrix
# colnames(fit$coefficients)

fit <- lmFit(X, design)

hist(fit$Amean)
plotSA(fit)
abline(v=5)

keep <- fit$Amean > 5
sum(keep)
fit2<- contrasts.fit(fit, contrast.matrix)
dim(fit2)
fit2 <- eBayes(fit2[keep,], trend = TRUE)
dim(fit2)
plotSA(fit2)

lfc <- 1
pval <- 0.05

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'bct-vrl')
dim(results)
head(results)
summary(results)
vennDiagram(results, include = 'both')
# vennCounts(results, include = 'both')


results.hits <- union(rownames(X[keep,])[results[,1] == 1]
                      ,rownames(X[keep,])[results[,1] == -1])
results.hits
length(results.hits)
results.hits.2 <- union(rownames(X[keep,])[results[,2] == 1]
                        ,rownames(X[keep,])[results[,2] == -1])
results.tot <- union(results.hits, results.hits.2)
length(results.tot)

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'bct-HC')
dim(top.hits)
head(top.hits)
all.hits <- topTable(fit2, number=nrow(fit2), coef = 'bct-HC', 'vrl-HC')
dim(top.hits)
dim(all.hits)

intersect(results.tot, rownames(top.hits))

ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  ggtitle("Volcano Plot of Log Fold Change Against -log10 P Value
          Cutoff - Fold Change:1, P Val:0.05")


results.tot
dim(X)
X[results.tot,]
dim(X[results.tot,])



# PCA
full.pca <- prcomp(t(X[results.tot,]), scale=TRUE)

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

status[idx,]$most_general
dim(pair1)

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', '#464647', "#165bfc")
dx.cols.f <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#bd35fc')
sex.cols <- c('#fc1676', '#16acfc')

positions <- c('bacterial', 'greyb', 'greyv', 'viral', 'HC')
positions.f <- c('bacterialpos', 'bacterialneg', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother', 'HC')

# most general
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$most_general), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general[dx.def]), size=2) +  
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  # labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")


# gender
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$Sex), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Gender')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=sex.cols)+
  ggtitle("Male Female Split against PC1 PC2")

# category
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$Category), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Gender')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  ggtitle("Category")





### Clustering
X.t <- t(X[results.tot,])
dim(X[results.tot,])

fviz_nbclust(X.t, kmeans, method = "wss")
fviz_nbclust(X.t, kmeans, method = "silhouette")
fviz_nbclust(X.t, kmeans, method = "gap_stat")





### K4
k4.pal <- c("#f70d09", "#f76409", '#09f70d', '#2909f7')
site.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
cat.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#16fc31', '#464647', "#165bfc")
# k4<-cmeans(X.t, 4, iter.max = 50, verbose = FALSE,
#            dist = "euclidean", method = "cmeans", m = 2,
#            rate.par = NULL, weights = 1, control = list())
set.seed(2)
k4 <- kmeans(X.t, centers = 4, nstart = 100)
k4$cluster <- as.factor(k4$cluster)


table(k4$cluster, droplevels(status.iris.c[idx,]$most_general))
table(k4$cluster, droplevels(status.iris.c[idx,]$general))

p<-ggplot(status.iris.c[idx,], aes(most_general, fill=k4$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  geom_bar()
ggplotly(p)


positions.more <- c('bacterialgpos', 'bacterialgneg', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status.iris.c[idx,], aes(general, fill=k4$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k4$cluster, size = status.iris.c[idx,]$Age,
        colors=k4.pal, text= ~paste0('age: ', status.iris.c[idx,]$Age, '<br>label:', status.iris.c[idx,]$My_code, '<br>Diagnosis: ',status.iris.c[idx,]$Diagnosis)) %>%
  layout(title = 'Diagnostic Groups by PCA 1-2')

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k4$cluster, size = status.iris.c[idx,]$Age,
        colors = c(k4.pal), text= ~paste0('age: ', status.iris.c[idx,]$Age, '<br>label:', status.iris.c[idx,]$Category, '<br>Diagnosis: ', status.iris.c[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.iris.c[idx,]$Category, size = status.iris.c[idx,]$Age,
        colors=c(cat.pal), text= ~paste0('category: ', status.iris.c[idx,]$Category, '<br>age: ', status.iris.c[idx,]$Age, '<br>label:',status.iris.c[idx,]$My_code, '<br>Diagnosis: ',status.iris.c[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Category Groups by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.iris.c[idx,]$Sex, size = ~status.iris.c[idx,]$Age,
        colors = c(sex.cols), text= ~paste0('age: ', status.iris.c[idx,]$Age, '<br>label:',status.iris.c[idx,]$My_code, '<br>Diagnosis: ',status.iris.c[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Gender by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

ggplot(status.iris.c[idx,], aes(most_general, Age, fill=k4$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k4.pal, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")



# K4 analysis
k4.df <- status.iris.c[idx,][, c('Category', 'My_code', 'most_general', 'general',
                          'Age', 'Sex', 'Diagnosis')]
k4.df$cluster <- k4$cluster
dim(k4.df)

table(k4$cluster, droplevels(status.iris.c[idx,]$most_general))

# k4.df[1:5,(ncol(k4.df)-3): ncol(k4.df)]
# k4.df$array.contemporary.CRP <- as.numeric(as.character(k4.df$array.contemporary.CRP))
sum(k4.df$cluster == 1)
sum(k4.df[k4.df$cluster == 1,]$Category == 'E')
sum(k4.df[k4.df$cluster == 1,]$Category == 'F')

sum(k4.df$cluster == 2 & k4.df$most_general == 'bacterial')
sum(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$Category == 'E')/length(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$Category)
sum(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$Category == 'E')/length(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$Category)
sum(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$Category == 'E')/length(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$Category)
sum(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$Category == 'E')/length(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$Category)

sum(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$Category == 'E')
length(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$Category)

sum(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$Category == 'E')
length(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$Category)

sum(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$Category == 'E')
length(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$Category)

sum(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$Category == 'E')
length(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$Category)

prop.test(x = c(7, 13), n = c(9, 33),
          alternative = "two.sided")

a <- prop.test(x = c(14, 13), n = c(19, 33),
               alternative = "two.sided", correct = FALSE)
a$estimate[1]+a$estimate[2]


