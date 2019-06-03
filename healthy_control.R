
# library(tidyverse)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(dplyr)
library(plotly)
library(e1071)
library(neuralnet)
library(ROCR)

# note checkout lumi package

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')


Sys.setenv("plotly_username"="vetronz1992")
Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', '#384ebc', "#bc38ab")
dx.cols.f <- c("#ed0404", "#fc5716", '#16fc31','#38bc9d', '#55f1fc', '#0934f4', '#384ebc', "#bc38ab")
sex.cols <- c('#fc1676', '#16acfc')


positions <- c('bacterial', 'greyb', 'greyv', 'viral', 'HC')
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother', 'HC')


site.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
cat.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#16fc31', '#464647', "#165bfc")


idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'HC'
sum(idx)


which(as.character(status[idx,]$array.contemporary.CRP) == 'na')
length(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'))

which(status[idx,]$WBC == 0)
clean<-union(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'), which(status[idx,]$WBC == 0))
clean

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]

status[idx,]$array.contemporary.CRP[clean.idx]
wbc <- status[idx,]$WBC[clean.idx]
crp <- as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))



############ LIMMA ############
e.set.f <- as.data.frame(e.set[,idx])
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)
class(e.set.f)
length(as.character(status[idx, c('most_general')]))

e.set.f$label <- as.character(status[idx, c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]
e.set.f$site <- status[idx,c('site')]
e.set.f$age <- status[idx,c('Age..months.')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

# mean/variance calculations
x_var <- apply(e.set[,idx], 1, var)
x_mean <- apply(e.set[,idx], 1, mean)
plot(log2(x_mean), log2(x_var), pch='.')
abline(v=log2(5), col='red')

dim(e.set[,idx])
e.set[,idx][which(x_mean > 5)]
dim(t(e.set[,idx][which(x_mean > 5),]))
X <- e.set[,idx][which(x_mean > 5),]


### DESIGN MATRIX
# site 
design <- model.matrix(~label + sex + age + site + 0, data = e.set.f)
colnames(design)<- c("bct","greyb","greyv", 'HC', 'vrl', 'sexM', 'age', 'CHW', 'EUC', 'FED', 'FPIES', 'KEN', 'OXF','SMH','SOT','UCSD')

design[1:10,]
dim(design)
colSums(design)

# contrast.matrix<- makeContrasts("bct-HC", 'vrl-HC', 'greyb-HC', levels=design)
contrast.matrix<- makeContrasts("HC-bct", 'HC-vrl', levels=design)
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

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc)
dim(results)
# head(results)
# summary(results)
vennDiagram(results, include = 'both')
# vennCounts(results, include = 'both')


results.bct <- union(rownames(X[keep,])[results[,1] == 1]
                      ,rownames(X[keep,])[results[,1] == -1])
length(results.bct)
results.vrl <- union(rownames(X[keep,])[results[,2] == 1]
                      ,rownames(X[keep,])[results[,2] == -1])
results.tot <- union(results.bct, results.vrl)
length(results.tot)

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'HC-bct')
# top.hits
all.hits <- topTable(fit2, number=nrow(fit2), coef = 'HC-bct')
dim(top.hits)
dim(all.hits)

intersect(results.tot, rownames(top.hits))

p<-ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  ggtitle("Volcano Plot of Log Fold Change Against -log10 P Value")
p
# ggplotly(p)


# results.tot
dim(X)
# X[results.tot,]
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

# most_gen 2D
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status[idx,]$most_general), size = status[idx,]$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC1: (", round(pve[2],2), '%)')))
p
# api_create(p, filename = "2d_pca_filt")

# most_gen 3D
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~droplevels(status[idx,]$most_general), size = status[idx,]$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
p
# api_create(p, filename = "3d_pca_filt")

# more_gen
dx.cols.f <- c('#bd35fc', "#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$more_general,
        colors = c(dx.cols.f), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

# category
p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$category, size = status[idx,]$Age..months.,
        colors=c(cat.pal), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Category Groups by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
api_create(p, filename = "3d_pca_cat")

# sex
p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$Sex, size = ~status[idx,]$Age..months.,
        colors = c(sex.cols), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Gender by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
api_create(p, filename = "3d_pca_sex")

# site
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$site, size = ~status[idx,]$Age..months.,
        colors = c(site.pal), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Site Recruitment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


### Clustering
X.t <- t(X[results.tot,])

fviz_nbclust(X.t, kmeans, method = "wss")
# fviz_nbclust(X.t, kmeans, method = "silhouette")
# fviz_nbclust(X.t, kmeans, method = "gap_stat", nboot = 10)

# fviz_nbclust(X.t, cluster::fanny, method = "wss")
# fviz_nbclust(X.t, cluster::fanny, method = "silhouette")
# fviz_nbclust(X.t, cluster::fanny, method = "gap_stat")


### K2
set.seed(42)
k2 <- kmeans(X.t, centers = 2, nstart = 40)
k2$cluster <- as.factor(k2$cluster)

table(k2$cluster, droplevels(status[idx,]$most_general))
table(k2$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols.2, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=dx.cols.2, 'Cluster')+
  geom_bar()
ggplotly(p)


plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k2$cluster[clean.idx],
        colors = c(dx.cols.2), size = crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

ggplot(status[idx,], aes(most_general, Age..months., fill=k2$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.2)+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k2$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.2, 'Cluster')+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k2$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.2, 'Cluster')
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)


bct <- status[idx,]$most_general == 'bacterial'
# sex
plot_ly(pair3D[bct,], x = ~PC1, y = ~PC2, z = ~PC3, color = k2$cluster[bct], size = as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)),
        symbol = droplevels(status[idx,][bct,]$Sex), symbols = c('circle', 'cross'),
        colors = c(dx.cols.2), text= ~paste0('Age: ',status[idx,][bct,]$Age..months., '<br>Sex: ',status[idx,][bct,]$Sex, '<br>WBC: ', status[idx,][bct,]$WBC,'<br>CRP: ', as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)), '<br>label:',status[idx,][bct,]$my_category_2, '<br>Diagnosis: ',status[idx,][bct,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

# microbiology
plot_ly(pair3D[bct,], x = ~PC1, y = ~PC2, z = ~PC3, color = k2$cluster[bct], size = as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)),
        symbol = droplevels(status[idx,][bct,]$category), symbols = c('circle', 'x'),
        colors = c(dx.cols.2), text= ~paste0('Age: ',status[idx,][bct,]$Age..months., '<br>Sex: ',status[idx,][bct,]$Sex, '<br>WBC: ', status[idx,][bct,]$WBC,'<br>CRP: ', as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)), '<br>label:',status[idx,][bct,]$my_category_2, '<br>Diagnosis: ',status[idx,][bct,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

# K2 analysis
k2.df <- status[idx,][, c('category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster
dim(k2.df)

p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(cluster, fill=category)) +
  labs(title = "Barplot of Microbiology Results by Cluster", x = "Diagnosis", y = "Counts")+
  scale_fill_manual(values=c('#6c5ddd', '#dd5dbb'), 'Micro')+
  geom_bar()
ggplotly(p)

table(droplevels(k2.df$category[bct]), droplevels(k2.df$cluster[bct]))
chisq.test(table(droplevels(k2.df$category[bct]), droplevels(k2.df$cluster[bct])))
# not significant


### K4
k4.pal <- c('#09f70d', "#f76409", "#f70d09", '#2909f7')
set.seed(44)
k4 <- kmeans(X.t, centers = 4, nstart = 100)
k4$cluster <- as.factor(k4$cluster)

table(k4$cluster, droplevels(status[idx,]$most_general))
table(k4$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k4$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  geom_bar()
p
api_create(p, filename = "barplot_k4_dx")


positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k4$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k4$cluster, size = status[idx,]$Age..months.,
        colors=k4.pal, text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Cluster Assignment on PC 1-2')

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k4$cluster, size = status[idx,]$Age..months.,
        colors = c(k4.pal), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Cluster Assignment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k4$cluster[clean.idx],
        colors = c(k4.pal), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2[clean.idx], '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k4$cluster[clean.idx])) + geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('WBC Count') +
  xlab('') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  ggtitle("Box Plot of Cluster WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k4$cluster[clean.idx]))+
  geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k4.pal, 'Cluster')
p <- gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)


ggplot(status[idx,], aes(most_general, Age..months., fill=k4$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k4.pal, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")



# K4 analysis
k4.df <- status[idx,][, c('category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k4.df$cluster <- k4$cluster
dim(k4.df)


p<-ggplot(k4.df[k4.df$most_general == 'bacterial',], aes(cluster, fill=category)) +
  labs(title = "Barplot of Microbiology Results by Cluster", x = "Diagnosis", y = "Counts")+
  scale_fill_manual(values=c('#6c5ddd', '#dd5dbb'), 'Micro')+
  geom_bar()
p
api_create(p, filename = "barplot_k4_micro")


table(droplevels(k4.df$category[k4.df$most_general=='bacterial']), droplevels(k4.df$cluster[k4.df$most_general=='bacterial']))
chisq.test(droplevels(k4.df$category[k4.df$most_general=='bacterial']), droplevels(k4.df$cluster[k4.df$most_general=='bacterial']))

addmargins(table(droplevels(k4.df$category[k4.df$most_general=='bacterial']), droplevels(k4.df$cluster[k4.df$most_general=='bacterial'])))
prop.test(x = c(2, 2, 20, 1), n = c(5, 9, 33, 5),
          alternative = "two.sided")

# microbiology
plot_ly(pair3D[bct,], x = ~PC1, y = ~PC2, z = ~PC3, color = k4$cluster[bct], size = as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)),
        symbol = droplevels(status[idx,][bct,]$category), symbols = c('circle', 'x'),
        colors = c(k4.pal), text= ~paste0('Age: ',status[idx,][bct,]$Age..months., '<br>Sex: ',status[idx,][bct,]$Sex, '<br>WBC: ', status[idx,][bct,]$WBC,'<br>CRP: ', as.numeric(as.character(status[idx,][bct,]$array.contemporary.CRP)), '<br>label:',status[idx,][bct,]$my_category_2, '<br>Diagnosis: ',status[idx,][bct,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


### K5
k5.pal <- c("#f70d09", "#f76409", '#09f70d', '#2afcfc', '#2909f7')
set.seed(42)
k5 <- kmeans(X.t, centers = 5, nstart = 100)
k5$cluster <- as.factor(k5$cluster)


table(k5$cluster, droplevels(status[idx,]$most_general))
table(k5$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k5$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k5.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k5$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k5.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k5$cluster, size = status[idx,]$Age..months.,
        colors=k5.pal, text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Cluster Assignment on PC 1-2')

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k5$cluster, size = status[idx,]$Age..months.,
        colors = c(k5.pal), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Cluster Assignment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k5$cluster[clean.idx],
        colors = c(k5.pal), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2[clean.idx], '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k5$cluster[clean.idx])) + geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('WBC Count') +
  xlab('') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k5.pal, 'Cluster')+
  ggtitle("Box Plot of Cluster WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k5$cluster[clean.idx]))+
  geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k5.pal, 'Cluster')
p <- gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

ggplot(status[idx,], aes(most_general, Age..months., fill=k5$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k5.pal, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")



# k5 analysis
k5.df <- status[idx,][, c('category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k5.df$cluster <- k5$cluster
dim(k5.df)

p<-ggplot(k5.df[k5.df$most_general == 'bacterial',], aes(cluster, fill=category)) +
  labs(title = "Barplot of Microbiology Results by Cluster", x = "Diagnosis", y = "Counts")+
  scale_fill_manual(values=c('#6c5ddd', '#dd5dbb'), 'Micro')+
  geom_bar()
ggplotly(p)

table(droplevels(k5.df$category[bct]), droplevels(k5.df$cluster[bct]))
chisq.test(droplevels(k5.df$category[bct]), droplevels(k5.df$cluster[bct]))

prop.test(x = c(0, 19, 2, 4), n = c(3, 33, 5, 11),
          alternative = "two.sided")




############# limma differential expression gram +ve gram -ve #############
dim(X.t)

X.g <- as.data.frame(X.t[status[idx,]$most_general == 'bacterial',])
dim(X.g)
class(X.g)
a <- status[idx,]$category == 'E' | status[idx,]$category == 'F'

sum(status[idx,]$category == 'F')

b <- status[idx,]$most_general == 'bacterial'
sum(a == b)

X.g$gram <- droplevels(status[idx,][a,]$category)


### DESIGN MATRIX
# site
design <- model.matrix(~gram + 0, data = X.g)

colnames(design)<- c('gram.pos', 'gram.neg')

design[1:10,]
dim(design)
colSums(design)

contrast.matrix<- makeContrasts("gram.pos-gram.neg", levels=design)
contrast.matrix
# colnames(fit$coefficients)


dim(t(X.g[-ncol(X.g)]))
dim(design)

fit <- lmFit(t(X.g[-ncol(X.g)]), design)

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

lfc <- 0.5
pval <- 0.1

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'gram.pos-gram.neg')
dim(results)
head(results)
summary(results)
vennDiagram(results, include = 'both')
# vennCounts(results, include = 'both')

dim(X.g[-nrow(X.g)])
colnames(X.g[-ncol(X.g)])

length(colnames(X.g[-ncol(X.g)])[keep])
dim(results)
colnames(X.g[-ncol(X.g)])[keep][results == 1]
colnames(X.g[-ncol(X.g)])[keep][results == -1]
results.tot <- union(colnames(X.g[-ncol(X.g)])[keep][results == 1], colnames(X.g[-ncol(X.g)])[keep][results == -1])

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'gram.pos-gram.neg', number = 19)
head(top.hits)
all.hits <- topTable(fit2, number=nrow(fit2), coef = 'gram.pos-gram.neg')
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

dim(X.t)
dim(X.t[,results.tot])
colnames(X.t[,results.tot])
gram.hits<-colnames(X.t[,results.tot])
gram.hits
# View(X.t[,results.tot])
X.t[,results.tot][a,]



####### hierarchecal clustering #######
# Dissimilarity matrix
d <- dist(t(X.t[,results.tot][a,]), method = "euclidean")
d <- dist(X.t[,results.tot][a,], method = "euclidean")

dim(t(X.t[,results.tot][a,]))

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "average" )
# hc1 <- hclust(d, method = "ward.D" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.7, hang = -1)
rect.hclust(hc1, k = 2, border = c(2,4))
rect.hclust(hc1, k = 3, border = c(2,4))

# outlier
dim(X.t[,results.tot][a,][rownames(X.t[,results.tot][a,]) != 'bacterialgpos_19_SMH',])
e.set.g <- X.t[,results.tot][a,][rownames(X.t[,results.tot][a,]) != 'bacterialgpos_19_SMH',]
status.g <- status[idx,][a,][status[idx,][a,]$my_category_2 != 'bacterialgpos_19_SMH',]
dim(e.set.g)
dim(status.g)


fviz_nbclust(e.set.g, FUN = hcut, method = "wss")
fviz_nbclust(e.set.g, FUN = hcut, method = "silhouette")
gap_stat <- clusGap(e.set.g, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# d <- dist(t(e.set.g), method = "euclidean")
d <- dist(e.set.g, method = "euclidean")
hc1 <- hclust(d, method = "average" )
# hc1 <- hclust(d, method = "ward.D" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.7, hang = -1)
rect.hclust(hc1, k = 2, border = c(2,4))
rect.hclust(hc1, k = 3, border = c(2,4))

# Cut tree into 2 groups
sub_grp <- cutree(hc1, k = 3)
sub_grp
table(droplevels(status.g$category), sub_grp)
chisq.test(table(droplevels(status.g$category), sub_grp))

fviz_cluster(list(data = e.set.g, cluster = sub_grp))

# Hierarchical clustering using Complete Linkage
d1 <- as.dist(1-cor(t(e.set.g), method="pearson"))
d2 <- as.dist(1-cor(e.set.g, method="spearman"))


hc2 <- hclust(d1, method = "average" )
hc3 <- hclust(d2, method = "average" )

# Plot the obtained dendrogram
plot(hc2, cex = 0.7, hang = -1)
rect.hclust(hc2, k = 2, border = c(2,4))
rect.hclust(hc2, k = 3, border = c(2,4))

# Plot the obtained dendrogram
plot(hc3, cex = 0.7, hang = -1)
rect.hclust(hc3, k = 2, border = c(2,4))
rect.hclust(hc3, k = 3, border = c(2,4))

clustRows <- hclust(d1, method="average")
clustColumns <- hclust(d2, method="average")

# module.assign <- cutree(clustColumns, k=3)
# module.assign == 2
# module.assign[module.assign == 2]

module.assign <- cutree(clustRows, k=2)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
myheatcolors2 <- colorRampPalette(colors=c("blue", "yellow","red"))(100)
# produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(e.set.g,
          Rowv=as.dendrogram(clustRows), 
          Colv=NA,
          RowSideColors=module.color,
          col=myheatcolors2, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))


heatmap.2(e.set.g,
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none")


library(heatmaply) #for making interactive heatmaps using plotly
p<-heatmaply(e.set.g,
          colors = myheatcolors2,
          Rowv=as.dendrogram(clustRows),
          RowSideColors=module.color,
          # showticklabels=c(FALSE,FALSE),
          scale='row')
p
api_create(p, filename = "heatmap_micro")

p<-heatmaply(mtcars, k_col = 2, k_row = 3) %>% layout(margin = list(l = 130, b = 40))
api_create(p, filename = "heatmap_test")

library("illuminaHumanv4.db")

getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)

head(illumina)
nrow(illumina)

module.assign[module.assign == 2]
# 5720482 2570300 2100196  990768 3360343
trans <- c(3180392, 2570300, 2100196, 5720482, 7650358, 4780075, 7650433)
# trans <- c(7650358, 5720482, 2570300, 3180392, 5090754)
# trans <- c(3180392)
# trans <- c(2100196)

which(illumina$Array_Address_Id %in% trans)

probeID <- illumina$Probe_Id[which(illumina$Array_Address_Id %in% trans)]
probeID

x <- illuminaHumanv4CHR
x <- illuminaHumanv4NUID
x <- illuminaHumanv4ALIAS2PROBE
x <- illuminaHumanv4ENSEMBL
x <- illuminaHumanv4GENENAME
x <- illuminaHumanv4GO
x <- illuminaHumanv4MAP
x <- illuminaHumanv4REFSEQ
data.frame(trans, Gene=unlist(mget(x = probeID, envir = x)))




# fuzzy clustering allocation
k4.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, most_general) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))


which(k4.df$most_general == 'bacterial' & k4.df$cluster == 2)
which(k4.df$most_general == 'bacterial' & k4.df$cluster == 1)

View(k4.df[which(k4.df$most_general == 'bacterial' & k4.df$cluster == 2),])
View(k4.df[which(k4.df$most_general == 'bacterial' & k4.df$cluster == 1),])

k4$membership[which(status[idx,]$most_general == 'greyb'),]

which.max(k4$membership[which(status[idx,]$most_general == 'greyb'),][,2])

which(status[idx,]$most_general == 'greyb')[24]

status[idx,]$my_category_2[131]
pair3D[131,]
k4$cluster

fcm_centroids <- k4$centers
cor(t(fcm_centroids))


### K7
k7.pal <- c('#09f70d', '#2909f7',"#f70d09", '#5dddbb', '#fc37f5', '#ccf246', "#f76409")
set.seed(18)
k7 <- kmeans(X.t, centers = 7, nstart = 100)
k7$cluster <- as.factor(k7$cluster)

table(k7$cluster, droplevels(status[idx,]$most_general))
table(k7$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k7$cluster, size = status[idx,]$Age..months.,
        colors=k7.pal, text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Cluster Assignment on PC 1-2', xaxis=x, yaxis=y)

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster, size = status[idx,]$Age..months.,
        colors = c(k7.pal), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Cluster Assignment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster[clean.idx],
        colors = c(k7.pal), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2[clean.idx], '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k7$cluster[clean.idx])) + geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('WBC Count') +
  xlab('') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  ggtitle("Box Plot of Cluster WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k7$cluster[clean.idx]))+
  geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k7.pal, 'Cluster')
p <- gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

ggplot(status[idx,], aes(most_general, Age..months., fill=k7$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k7.pal, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")


# K7 analysis
k7.df <- status[idx,][, c('category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k7.df$cluster <- k4$cluster
dim(k7.df)
View(k7.df)

dim(k7.df[k7.df$most_general == 'bacterial',])
k7.df[k7.df$most_general == 'bacterial',]

ggplot(k7.df[k7.df$most_general == 'bacterial',], aes(cluster, fill=category)) +
  labs(title = "Barplot of Microbiology Results by Cluster", x = "Diagnosis", y = "Counts")+
  scale_fill_manual(values=c('#6c5ddd', '#dd5dbb'), 'Micro')+
  geom_bar()







###########################################################################################
### prediction
dim(X.t)

scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

X.s<-data.frame(apply(X.t, 2, scale01))
dim(X.s)

X.s$bacterial <- status[idx,]$most_general == 'bacterial'
dim(X.s)

# train test split
train <- sample(seq(1, dim(X.s)[1]), dim(X.s)[1]*0.7, replace = FALSE)
test <- seq(1, dim(X.s)[1])[-train]

train.x <- X.s[train,]
test.x <- X.s[test,][-ncol(X.s)]
y <- ifelse(status[idx,]$most_general[test] == 'bacterial', TRUE, FALSE)
dim(train.x)
dim(test.x)
y

## fit model using `f`
# colnames(X.s)
# model <- neuralnet(bacterial ~ X1740360 + X450189, data = train.x, hidden=c(4,3), threshold=0.01)

# hi hi hi git
roc.a <- NULL
roc.t <- NULL
h.n <- 9
rep <- 10

for(j in 1:rep){
for(i in 1 : h.n) {
  model <- neuralnet(bacterial ~ ., data = train.x, hidden=i, threshold=0.01)
  pred <- predict(model, test.x, type="class")
  roc.a[i] <- prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
    performance(measure = "auc") %>%
    .@y.values
}
roc.t <- append(roc.t, roc.a)
}


roc.t

df <- data.frame(matrix(unlist(roc.t), nrow=length(roc.t), byrow=T))
colnames(df) <- 'a'
df
as.data.frame(split(df, 1:h.n))
df <- as.data.frame(split(df, 1:h.n))
df.t <- as.data.frame(t(df))
df.t$h.n <- seq(1:nrow(df.t))

boxplot(df, use.cols=TRUE)

# # boxplot reshaping

# colnames(df.t) <- c('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'rep8', 'rep9', 'rep10', 'rep11', 'rep12', 'h.n')
# df.t
# 
# dat.m <- melt(df.t,id.vars='h.n', measure.vars=c('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6'))
# dat.m
# 
# ggplot(dat.m) + geom_boxplot(aes(x=h.n, y=value, color=variable))
# 
# # melt(df, )
# dat.m <- melt(dat,id.vars='ID', measure.vars=c('Freq','Freq.1','Freq.2'))
# 
# roc.m <- as.data.frame(apply(df, 2, mean))
# roc.m$sd <- apply(df, 2, sd)
# roc.m$hidden.nodes <- seq.int(nrow(roc.m))
# roc.m
# colnames(roc.m)[1] <- 'mean.ROC'
# 
# roc.m
# p<-ggplot(roc.m, aes(x=hidden.nodes, y=mean.ROC))+geom_jitter()
# p<-ggplot(roc.m, aes(x=hidden.nodes, y=mean.ROC))+geom_point()



model <- neuralnet(bacterial ~ ., data = train.x, hidden=3, threshold=0.01)
pred <- predict(model, test.x, type="class")

prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "auc") %>%
  .@y.values



# end

