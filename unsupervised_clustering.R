library(tidyverse)
# library(plyr)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(dplyr)
library(plotly)


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


which(status[idx,]$WBC == 0)
clean<-union(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'), which(status[idx,]$WBC == 0))
clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]

wbc <- status[idx,]$WBC[clean.idx]
crp <- as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))


dim(e.set)
dim(e.set[,idx])

X <- e.set[,idx]
dim(X)

e.set.f <- t(e.set)[idx,]
dim(e.set.f)

# mean/variance calculations
x_var <- apply(X, 1, var)
x_mean <- apply(X, 1, mean)
plot(log2(x_mean), log2(x_var), pch='.')
abline(v=log2(5), col='red')

dim(X)
X[which(x_mean > 5)]
dim(t(X[which(x_mean > 5),]))
X.f <- X[which(x_mean > 5),]




# PCA
full.pca <- prcomp(e.set.f, scale=TRUE)
full.pca <- prcomp(e.set.f[,which(x_mean > 5)], scale=TRUE)

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
  labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols)+
  # scale_color_manual(values=dx.cols.2)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")

# most general
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

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$most_general,
        colors = c(dx.cols), marker = list(size = 3), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3',
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

# plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, symbol = ~droplevels(status[idx,]$most_general), symbols = c('circle','x', 'cross', '0', 'square','o'),
#         colors = c(dx.cols), marker = list(size = 3),  text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))



### Unsupervised Clustering
X.t <- t(X)
X.t <- t(X[which(x_mean > 5),])
dim(X.t)
# 
# fviz_nbclust(X.t, kmeans, method = "wss")
# fviz_nbclust(X.t, kmeans, method = "silhouette")
# 
# gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 5, K.max = 10, B = 5)
# fviz_gap_stat(gap_stat)
# optimal without filtering = 9
# optimal with filter = 8

### K2
k2 <- kmeans(X.t, centers = 2, nstart = 10)
k2$cluster <- as.factor(k2$cluster)
k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster[clean.idx]
dim(k2.df)
k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))

k2.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(most_general, cluster) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))
  
table(k2$cluster, droplevels(status[idx,]$most_general))
table(k2$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols.2, 'Cluster', labels = c("One", "Two"))+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=dx.cols.2, 'Cluster', labels = c("One", "Two"))+
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

ggplot(status[idx,], aes(most_general, Age..months., fill=k2$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.2)+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k2$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.2)+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k2$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.2)
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

### K7
k7 <- kmeans(X.t, centers = 7, nstart = 10)
k7$cluster <- as.factor(k7$cluster)
k7.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general', 'site',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k7.df$cluster <- k7$cluster[clean.idx]
dim(k7.df)
k7.df[1:5,(ncol(k2.df)-3): ncol(k7.df)]
k7.df$array.contemporary.CRP <- as.numeric(as.character(k7.df$array.contemporary.CRP))

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


dx.cols.r <- c("#165bfc", "#ed0404")
# mode = 'markers', symbol = ~Species, symbols = c('circle','x','o')
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster,
        colors = c(dx.cols.f), marker = list(size = 3), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster[clean.idx],
        colors = c(dx.cols.f), size = crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


ggplot(status[idx,], aes(most_general, Age..months., fill=k7$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.f)+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k7$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.f)+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k7$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.f)
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)


### K8
k8 <- kmeans(X.t, centers = 8, nstart = 10)
k8$cluster <- as.factor(k8$cluster)
k8.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general', 'site',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k8.df$cluster <- k8$cluster[clean.idx]
dim(k8.df)
k8.df[1:5,(ncol(k2.df)-3): ncol(k8.df)]
k8.df$array.contemporary.CRP <- as.numeric(as.character(k8.df$array.contemporary.CRP))

k8.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

table(k8$cluster, droplevels(status[idx,]$most_general))
table(k8$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k8$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)

p<-ggplot(status[idx,], aes(more_general, fill=k8$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=dx.cols.f, 'Cluster')+
  geom_bar()
ggplotly(p)


dx.cols.r <- c("#165bfc", "#ed0404")
# mode = 'markers', symbol = ~Species, symbols = c('circle','x','o')
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k8$cluster,
        colors = c(dx.cols.f), marker = list(size = 3), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k8$cluster[clean.idx],
        colors = c(dx.cols.f), size = crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


ggplot(status[idx,], aes(most_general, Age..months., fill=k8$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=dx.cols.f)+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k8$cluster[clean.idx])) + geom_boxplot() +
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=dx.cols.f)+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k8$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=dx.cols.f)
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

# end