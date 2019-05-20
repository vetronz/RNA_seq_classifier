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


# PCA
dim(e.set)
dim(e.set[,idx])

X <- e.set[,idx]
dim(X)

dim(t(e.set)[idx,])

e.set.f <- t(e.set)[idx,]
dim(e.set.f)

full.pca <- prcomp(e.set.f, scale=TRUE)
pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

status[idx,]$most_general
dim(pair1)

# DEF DX PC1 PC2
# fviz_pca_ind(full.pca)

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', '#464647', "#165bfc")
dx.cols.f <- c("#ed0404", "#fc5716", '#16fc31', '#464647', "#165bfc", "#ed0404", "#fc5716", '#16fc31', '#464647', "#165bfc")
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



plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$most_general[clean.idx],
        colors = c(dx.cols), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, symbol = ~droplevels(status[idx,]$most_general), symbols = c('circle','x', 'cross', '0', 'square','o'),
        colors = c(dx.cols), marker = list(size = 3),  text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


pair1[pair1$PC1 < x.pos & pair1$PC2 > y.pos,]
status[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,]
View(status[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,])
status[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,]$Diagnosis
# bacterialgpos_18_SMH    -76.54628 112.37286
# bacterialgpos_23_SOT   -113.17935 117.63147
# bacterialgpos_2_EUC101 -135.65958 107.72649

# a<-status[idx,]$category == 'E' | status[idx,]$category == 'F'
# a
# b <- status$my_category_2 == 'bacterialgpos_18_SMH' | status$my_category_2 == 'bacterialgpos_23_SOT' | status$my_category_2 == 'bacterialgpos_2_EUC101'
# status[b,]$Diagnosis
# 
# dim(pair1[a,])
# View(status[idx,][a,])
# gb.strep <- c(27,31,32,49)
# View(status[idx,][a,][gb.strep,])
# library(ggrepel)
# ggplot(pair1[a,][gb.strep,], aes(PC1, PC2)) + geom_text_repel(label=status[idx,][a,][gb.strep,]$Diagnosis) +
#   xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
#   ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
#   geom_hline(yintercept = y.pos, linetype="longdash", colour="red", size=.5) +
#   geom_vline(xintercept = x.pos, linetype="longdash", colour="red", size=.5) +
#   # scale_color_manual(values=dx.cols)+
#   scale_x_continuous(limits=c(-200,150))+
#   scale_y_continuous(limits=c(-100,200))+
#   ggtitle("Bacterial Viral Split on PC 1 and PC2")








### Unsupervised Clustering
X.t <- t(X)
fviz_nbclust(X.t, kmeans, method = "wss")
fviz_nbclust(X.t, kmeans, method = "silhouette")

# gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 10,
# K.max = 20, B = 20)
# fviz_gap_stat(gap_stat)

### K2
dim(e.set.f)

k2 <- kmeans(X.t, centers = 2, nstart = 10)
k2$cluster <- as.factor(k2$cluster)
table(k2$cluster, droplevels(status[idx,]$most_general))
table(k2$cluster, droplevels(status[idx,]$more_general))
k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster[clean.idx]
dim(k2.df)
k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))

dx.cols.r <- c("#165bfc", "#ed0404")
# mode = 'markers', symbol = ~Species, symbols = c('circle','x','o')
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k2$cluster, size = 1,
        colors = c(dx.cols.r), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k2$cluster[clean.idx],
        colors = c(dx.cols.r), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


status[idx,]
dim(status[idx,])
dim(status[clean.idx,])

status[clean.idx,] %>%
  select(WBC, most_general, array.contemporary.CRP, Age..months., Sex)
group_by(k2$cluster[clean.idx]) %>%
  summarise(wbc.m = mean(WBC), age.m = mean(Age..months.))

k2.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, Sex) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))


### K7
k7 <- kmeans(X.t, centers = 7, nstart = 15)

k7$cluster <- as.factor(k7$cluster)

table(k7$cluster, droplevels(status[idx,]$most_general))



plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster[clean.idx],
        colors = c(dx.cols.f), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))





dx.def <- status[idx,]$most_general == 'bacterial' | status[idx,]$most_general == 'viral'
dx.grey <- status[idx,]$most_general == 'greyb' | status[idx,]$most_general == 'greyv'
dx.def

dim(pair1[dx.def,])
length(status[idx,]$most_general[dx.def])