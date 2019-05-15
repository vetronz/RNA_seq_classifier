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
# library(plotly)

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


v <- intersect(colnames(e.set.i), status.iris$My_code)
common <- match(v, status.iris$My_code)
status.iris[common,]
dim(status.iris[common,])
View(status.iris[common,])

as.character(colnames(e.set.i)) == as.character(status.iris[common,]$My_code) # double check correctly filtered

status.iris.c <- status.iris[common,]
idx <- status.iris.c['most_general'] == 'bacterial' |
  status.iris.c['most_general'] == 'viral' |
  status.iris.c['most_general'] == 'greyb' |
  status.iris.c['most_general'] == 'greyv'
# status['most_general'] == 'OD'
sum(idx)
length(idx)


# PCA
dim(e.set.i)
dim(e.set.i[,idx])

dim(t(e.set.i[,idx]))
full.pca <- prcomp(t(e.set.i[,idx]), scale=TRUE)
pair1 <- as.data.frame(full.pca$x[,1:2])
pair3D <- as.data.frame(full.pca$x[,1:3])
pair2 <- as.data.frame(full.pca$x[,3:4])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]


# DEF DX PC1 PC2
# fviz_pca_ind(full.pca)
dim(pair1)
dim(status.iris[idx,])

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$most_general), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$Sex), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=sex.cols)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")


# top left corner looks like an interesting subset of samples
# want to check if they reemerge following clustering
x.pos <- 10
y.pos <- 100

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status.iris.c[idx,]$most_general, shape=status.iris.c[idx,]$Sex), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Diagnosis')+
  geom_hline(yintercept = y.pos, linetype="longdash", colour="red", size=.5) +
  geom_vline(xintercept = x.pos, linetype="longdash", colour="red", size=.5) +
  scale_color_manual(values=dx.cols)+
  ggtitle("Diagnostic Group Breakdown of based on PC1-PC2")

library(plotly)
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.iris.c[idx,]$most_general, colors = c(dx.cols)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))



pair1[pair1$PC1 < x.pos & pair1$PC2 > y.pos,]
status.iris.c[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,]$Diagnosis
View(status.iris.c[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,])



