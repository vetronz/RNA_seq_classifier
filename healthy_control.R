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

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

dim(e.set)

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
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother', 'HC')

# most general
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general[dx.def]), size=2) +  
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  # labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")

# site recruitment
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$site), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general[dx.def]), size=2) +  
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols.f)+
  # scale_color_manual(values=dx.cols.2)+
  ggtitle("Site Recruitment Breakdown of based on PC1-PC2")

# gender
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Gender')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=sex.cols)+
  ggtitle("Male Female Split against PC1 PC2")

# category
ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$category), size=2) +
  # ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Gender')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  ggtitle("Category")


library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
wbc.col <- scale_colour_gradientn(colours = myPalette(10), limits=c(min(wbc), max(wbc)))
crp.col <- scale_colour_gradientn(colours = myPalette(10), limits=c(min(crp), max(crp)))


f <- list(
  size = 14,
  color = "#3c3b3d"
)
x <- list(
  title = paste0("PC1: (", round(pve[1],2), '%)'),
  titlefont = f
)
y <- list(
  title = paste0("PC2: (", round(pve[2],2), '%)'),
  titlefont = f
)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~status[idx,]$most_general,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Diagnostic Groups by PCA 1-2', xaxis=x, yaxis=y)


# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$most_general, size = 1,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$most_general[clean.idx],
        colors = c(dx.cols), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))



### Clustering
X.t <- t(X[results.tot,])
dim(X[results.tot,])
 
fviz_nbclust(X.t, kmeans, method = "wss")
# fviz_nbclust(X.t, kmeans, method = "silhouette")
# fviz_nbclust(X.t, kmeans, method = "gap_stat")

# fviz_nbclust(X.t, cluster::fanny, method = "wss")
# fviz_nbclust(X.t, cluster::fanny, method = "silhouette")
# fviz_nbclust(X.t, cluster::fanny, method = "gap_stat")
### fuzzy clustering wants 2 cluster uniformly

# gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 10, K.max = 15, B = 10)
# fviz_gap_stat(gap_stat)
# optimal without filtering = 9
# optimal with filter = 8

### K2
k2 <- kmeans(X.t, centers = 2, nstart = 20)
k2$cluster <- as.factor(k2$cluster)
k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster[clean.idx]
dim(k2.df)
k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))

k2.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, most_general) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

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

dx.cols.r <- c("#165bfc", "#ed0404")
# mode = 'markers', symbol = ~Species, symbols = c('circle','x','o')
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k2$cluster,
        colors = c(dx.cols.2), marker = list(size = 3), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

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




### K4
k4.pal <- c("#f70d09", "#f76409", '#09f70d', '#2909f7')
site.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
cat.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#16fc31', '#464647', "#165bfc")
# k4<-cmeans(X.t, 4, iter.max = 50, verbose = FALSE,
#            dist = "euclidean", method = "cmeans", m = 2,
#            rate.par = NULL, weights = 1, control = list())
k4 <- kmeans(X.t, centers = 4, nstart = 100)
k4$cluster <- as.factor(k4$cluster)



table(k4$cluster, droplevels(status[idx,]$most_general))
table(k4$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k4$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

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
  layout(title = 'Diagnostic Groups by PCA 1-2', xaxis=x, yaxis=y)

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k4$cluster, size = status[idx,]$Age..months.,
        colors = c(k4.pal), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$category, size = status[idx,]$Age..months.,
        colors=c(cat.pal), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$Sex, size = ~status[idx,]$Age..months.,
        colors = c(sex.cols), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$site, size = ~status[idx,]$Age..months.,
        colors = c(site.pal), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
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
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k4$cluster[clean.idx]))+
  geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k4.pal, 'Cluster')
p <- gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, color=k4$cluster[clean.idx])) + geom_jitter(position=position_dodge(width=0.8)) +
  ylab('WBC Count') +
  xlab('') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k4.pal, 'Cluster')+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k4$cluster[clean.idx]))+
  geom_jitter(position=position_dodge(width=0.8)) +
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
# k4.df[1:5,(ncol(k4.df)-3): ncol(k4.df)]
# k4.df$array.contemporary.CRP <- as.numeric(as.character(k4.df$array.contemporary.CRP))
sum(k4.df$cluster == 1)
sum(k4.df[k4.df$cluster == 1,]$category == 'E')
sum(k4.df[k4.df$cluster == 1,]$category == 'F')

sum(k4.df$cluster == 1 & k4.df$most_general == 'bacterial')
sum(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$category == 'E')/length(k4.df[k4.df$cluster == 1 & k4.df$most_general == 'bacterial',]$category)
sum(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$category == 'E')/length(k4.df[k4.df$cluster == 2 & k4.df$most_general == 'bacterial',]$category)
sum(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$category == 'E')/length(k4.df[k4.df$cluster == 3 & k4.df$most_general == 'bacterial',]$category)
sum(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$category == 'E')/length(k4.df[k4.df$cluster == 4 & k4.df$most_general == 'bacterial',]$category)



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


### K6
k6 <- kmeans(X.t, centers = 6, nstart = 20)
k6$cluster <- as.factor(k6$cluster)

k6.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k6.df$cluster <- k6$cluster[clean.idx]
dim(k6.df)
k6.df[1:5,(ncol(k6.df)-3): ncol(k6.df)]
k6.df$array.contemporary.CRP <- as.numeric(as.character(k6.df$array.contemporary.CRP))

k6.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

table(k6$cluster, droplevels(status[idx,]$most_general))
table(k6$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k6$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k6$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)


# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k6$cluster,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Diagnostic Groups by PCA 1-2', xaxis=x, yaxis=y)


# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k6$cluster, size = 1,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k6$cluster[clean.idx],
        colors = c(dx.cols), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

ggplot(status[idx,], aes(most_general, Age..months., fill=k6$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k6$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k6$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)



### k7
dim(X.t)
dim(X[results.tot,])

k7<-cmeans(X.t, 7, iter.max = 50, verbose = FALSE,
           dist = "euclidean", method = "cmeans", m = 2,
           rate.par = NULL, weights = 1, control = list())
k7 <- kmeans(X.t, centers = 7, nstart = 20)
# k7 <- kmeans(X[results.tot,], centers = 7, nstart = 20)

fcm_centroids <- k7$centers
fcm_centroids_df <- as.data.frame(fcm_centroids)
fcm_centroids_df$cluster <- row.names(fcm_centroids_df)
fcm_centroids_df[,1:5]
cor(t(fcm_centroids))

k7$cluster <- as.factor(k7$cluster)
k7.df <- status[idx,][, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k7.df$cluster <- k7$cluster
dim(k7.df)
k7.df[1:5,(ncol(k7.df)-3): ncol(k7.df)]
k7.df$array.contemporary.CRP <- as.numeric(as.character(k7.df$array.contemporary.CRP))
View(k7.df)

k7.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

table(k7$cluster, droplevels(status[idx,]$most_general))
table(k7$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)


# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k7$cluster,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label: ',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Diagnostic Groups by PCA 1-2', xaxis=x, yaxis=y)


# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster, size = 1,
        colors = c(dx.cols), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label: ',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster[clean.idx],
        colors = c(dx.cols), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label: ',status[idx,]$my_category_2[clean.idx], '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

ggplot(status[idx,], aes(most_general, Age..months., fill=k7$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k7$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k7$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.f, 'Cluster')
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)





###########################################################################################
### prediction
dim(X.t)
scale()
?apply()


X.t[1:10,1:5]

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
rep <- 20

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



pred <- predict(model, test.x, type="class")
# pred
# hist(pred)

table(y,pred[,1]>0.5)

prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

# model 1 AUC
prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "auc") %>%
  .@y.values


# end