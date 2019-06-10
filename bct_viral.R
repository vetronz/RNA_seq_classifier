############ BCT - VIRAL ############

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
library(tidyr)


getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')

# Sys.setenv("plotly_username"="vetronz1992")
# Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- cols[1:5]
dx.cols.6 <- cols

dx.cols.f <- c("#ed0404", "#fc5716", '#16fc31','#38bc9d', '#55f1fc', '#0934f4', '#ac04f9', "#bc38ab")
sex.cols <- c('#fc1676', '#16acfc')
k7.pal <- c("#ed0404", "#fc5716", '#16fc31', '#55f1fc', '#0452f9', '#7304f9', '#f904e9')

positions <- c('bacterial', 'greyb', 'greyv', 'viral', 'greyu')
positions.f <- c('bacterial', 'greyb', 'greyv', 'greyu', 'adeno', 'flu', 'RSV', 'viralother')

site.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
cat.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#16fc31', '#464647', "#165bfc")


idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu' |
  status['most_general'] == 'HC'
sum(idx)

# outlier
which(status$my_category_2 == 'bacterialgpos_19_SMH')
status$my_category_2[28]
idx[28] <- FALSE
sum(idx)
length(idx)



############ LIMMA ############
e.set.f <- as.data.frame(e.set[,idx])
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)

e.set.f$label <- as.character(status[idx, c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]
e.set.f$age <- status[idx,c('Age..months.')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

# mean/variance calculations
x_var <- apply(e.set[,idx], 1, var)
x_mean <- apply(e.set[,idx], 1, mean)
df <- data.frame(log2(x_var), log2(x_mean))

ggplot(df, aes(log2.x_mean., log2.x_var.)) +
  geom_vline(xintercept=log2(5))+
  geom_point(size = 0.2, stroke = 0, shape = 16)+
  labs(title="Mean Variance Scatter Plot",
       x ="log2 Mean Expressioin", y = "log2 Variance")

dim(e.set[,idx])
dim(t(e.set[,idx][which(x_mean > 5),]))
X <- e.set[,idx][which(x_mean > 5),]
dim(X)

### DESIGN MATRIX
design <- model.matrix(~label + sex + age + 0, data = e.set.f)
colnames(design)<- c("bct","greyb",'greyu', "greyv", 'HC', 'vrl', 'sexM', 'age')

design[1:5,]
dim(design)
colSums(design)

# contrast.matrix<- makeContrasts("bct-HC", 'vrl-HC', 'greyb-HC', levels=design)
contrast.matrix<- makeContrasts('bct-vrl', levels=design)
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

bootstraps <- list(c(0, 1), # 1
                   c(0.25, 0.1), # 2
                   c(0.375, 0.1), # 3
                   c(0.5, 0.05), # 4
                   c(0.75, 0.05), # 5
                   c(1, 0.05)) # 6

lfc <- bootstraps[[4]][1]
pval <- bootstraps[[4]][2]
lfc
pval

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc)
dim(results)
head(results)
# summary(results)
vennDiagram(results, include = 'both')

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc)

all.hits <- topTable(fit2, number=nrow(fit2))
dim(all.hits)

p<-ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size = 1, stroke = 0, shape = 16) +
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  labs(title="Volcano Plot of Log Fold Change Against -log10 P Value",
       x ="Log Fold Change", y = "log10 P-value")
p


idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu'

# outlier
idx[28] <- FALSE
idx[28]
sum(idx)
length(idx)

X <- e.set[,idx][which(x_mean > 5),]
dim(X)
dim(results)

results.tot <- ifelse(results[,1] == 0, FALSE, TRUE)

dim(X[results.tot,])

# # PCA
full.pca <- prcomp(t(X[results.tot,]), scale=TRUE)

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

wbc<-bv.k2.df$WBC
crp<-as.numeric(as.character(bv.k2.df$array.contemporary.CRP))

# # most_gen 2D Age
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(bv.k2.df$most_general), size = bv.k2.df$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', bv.k2.df$category, '<br>age: ', bv.k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',bv.k2.df$my_category_2, '<br>Diagnosis: ',bv.k2.df$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups, Age Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p

# most_gen 2D CRP
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(bv.k2.df$most_general), size = as.numeric(as.character(bv.k2.df$array.contemporary.CRP)),
             colors=c(dx.cols), text= ~paste0('category: ', bv.k2.df$category, '<br>age: ', bv.k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',bv.k2.df$my_category_2, '<br>Diagnosis: ',bv.k2.df$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups, CRP Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p

# # api_create(p, filename = "2d_pca_filt")



### Clustering
X.t <- t(X[results.tot,])
dim(X.t)
fviz_nbclust(X.t, kmeans, method = "wss")
# fviz_nbclust(X.t, kmeans, method = "silhouette")
# fviz_nbclust(X.t, kmeans, method = "gap_stat", nboot = 12)

### K2
n = 5
cols = gg_color_hue(n)

set.seed(47)
bv.k2 <- kmeans(X.t, centers = 2, nstart = 50)
bv.k2$cluster <- as.factor(bv.k2$cluster)


bv.df <- status[idx,][, c('barcode_megaexp', 'category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
dim(bv.df)

bv.df$clus.2.4 <- bv.k2$cluster
# View(bv.k2.df)

setwd('/home/patrick/Documents/RNA_seq_classifier/Data')
clin <- read.table('Mega_sub1_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = TRUE)

# Construction: get the full demographics into dataframe
# clin2 <- read.table('Mega_sub2_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
# dim(clin2)
# a <- c(1, 19, 27, 28, 33, 34, 84, 85) # barcode, Dx,  wbc, abs_neut, perc_neut, perc_lymp, dx_1, dx_2
# clin2[clin2$V1 == '9423641116_B', a] # search clin2 using barcode
# View(clin)


bv.df$most_general <- as.character(bv.df$most_general)
bv.df$most_general[bv.df$most_general == 'greyb'] <- 'probable bacterial'
bv.df$most_general[bv.df$most_general == 'greyu'] <- 'unknown'
bv.df$most_general[bv.df$most_general == 'greyv'] <- 'probable viral'
bv.df$most_general <- as.factor(bv.df$most_general)

bv.df$system <- clin$system
bv.df$system.spec <- clin$system_spec
bv.df$micro <- clin$micro
bv.df$sepsis <- clin$sepsis

bv.df$my_wbc <- clin$wbc
bv.df$abs_neut <- clin$abs_neut
bv.df$perc_neut <- clin$perc_neut
bv.df$perc_lymph <- clin$perc_lymph
bv.df$Path_1 <- clin$Path_1
bv.df$Path_2 <- clin$Path_2
bv.df$Path_3 <- clin$Path_3

# have to use this ordering to ensure labels corespond to same color code in alphabetical ggplot
table(droplevels(bv.df$most_general), bv.df$clus.2.4)
dx <- c('viral', 'unknown', 'probable viral', 'probable bacterial', 'bacterial', 'viral', 'unknown', 'probable viral', 'probable bacterial', 'bacterial')

clus1 <- table(bv.df$clus.2.4, droplevels(bv.df$most_general))[1,]
clus2 <- table(bv.df$clus.2.4, droplevels(bv.df$most_general))[2,]
df.1 <- data.frame(clus1, clus2)
# df.1

df.2 <- gather(df.1, key='cluster', value = 'count')
# df.2
df.3 <- mutate(df.2, Diagnosis = dx)
# df.3
library(plyr)
df.4 <- ddply(df.3, "cluster", transform, tot=cumsum(count))
# df.4

ggplot(data=df.4, aes(x=cluster, y=count, fill=Diagnosis)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=tot, label=count), vjust=1.2, color='white', position = position_dodge(width = 0.9))


ggplot(bv.df, aes(bv.k2$cluster, Age..months., fill=bv.k2$cluster)) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(bv.df, aes(x = bv.k2$cluster, y = bv.df$WBC, fill = bv.k2$cluster)) +
  scale_fill_manual(values=cols[c(1,4)])+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(bv.df, aes(x = bv.k2$cluster, y = as.numeric(as.character(bv.df$array.contemporary.CRP)), fill = bv.k2$cluster)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

p1<-ggplot(bv.df, aes(x = bv.k2$cluster, y = bv.df$WBC, fill = bv.df$most_general)) +
  scale_fill_manual(values=dx.cols, name = 'Diagnosis')+
  # guides(fill=FALSE)+
  labs(title="Boxplot WBC Distributions within Clusters",
       x ="Cluster", y = "WBC Count") +
  geom_boxplot()
p2<-ggplot(bv.df, aes(x = bv.k2$cluster, y = as.numeric(as.character(bv.df$array.contemporary.CRP)), fill = bv.df$most_general)) +
  scale_fill_manual(values=dx.cols, name = 'Diagnosis')+
  labs(title="Boxplot CRP Distributions within Clusters",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, nrow=2)

ggplot(bv.df[bv.df$most_general == 'bacterial',], aes(clus.2.4, as.numeric(bv.df[bv.df$most_general == 'bacterial',]$abs_neut), fill=clus.2.4)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Dx')+
  labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Absolute Neutrophil Count") +
  geom_boxplot()

p1<-ggplot(bv.df[bv.df$most_general == 'bacterial',], aes(clus.2.4, as.numeric(bv.df[bv.df$most_general == 'bacterial',]$perc_neut), fill=clus.2.4)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  guides(fill=FALSE)+
  labs(title="Boxplot WBC Percent Neutrophil Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(bv.df[bv.df$most_general == 'bacterial',], aes(clus.2.4, as.numeric(bv.df[bv.df$most_general == 'bacterial',]$perc_lymph), fill=clus.2.4)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  labs(title="Boxplot WBC Percent Lymphocyte Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

# p<-plot_ly(bv.df, x = k2$cluster, y = ~wbc, type = "box",
#            color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'WBC Count'))
# p
# 
# p<-plot_ly(bv.df, x = k2$cluster, y = ~crp,
#            type = "box", color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'CRP Count'))
# p


# K2 analysis
# system
ggplot(bv.df[bv.df$most_general == 'bacterial',], aes(clus.2.4, fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

# micro
ggplot(bv.df[bv.df$most_general == 'bacterial',], aes(clus.2.4, fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

# lookup bacterial outliers consistently assigned clus2
bv.df$most_general == 'bacterial' & bv.df$clus.2.4 == 2
# View(bv.df[bv.df$most_general == 'bacterial' & bv.df$clus.2.4 == 2,])

plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~bv.df$clus.2.4,
        colors=cols, text= ~paste0('category: ', bv.df$category, '<br>age: ', bv.df$Age..months., '<br>WBC: ', bv.df$WBC, '<br>CRP: ', as.numeric(as.character(bv.df$array.contemporary.CRP)), '<br>label:',bv.df$my_category_2, '<br>Micro: ', bv.df$Path_1, '<br>Diagnosis: ',bv.df$Diagnosis),
        symbol = ~ifelse(bv.df$most_general == 'bacterial', 'bct', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))





################## BOOT STRAPPING K2 ################## 
tot <- 239
sum(bv.df$clus.2.3 != bv.df$clus.2.6) # check cluster assignment of 2.2 relative to 2.3 (base cluster)
tot - sum(bv.df$clus.2.3 != bv.df$clus.2.6)

boot.1 <- c(195, 44) # no filtering, full set
boot.2 <- c(220, 19) # 7000
boot.3 <- c(238, 1) # 3000
boot.4 <- c(239, 0) # 2200
boot.5 <- c(213, 26) # 800
boot.6 <- c(192, 47) # 300


df <- data.frame(boot.1, boot.2, boot.3, boot.4, boot.5, boot.6)
df.2<-gather(df)
df.3<-mutate(df.2, cluster = factor(c('Same','Different', 'Same','Different', 'Same','Different', 'Same','Different', 'Same','Different', 'Same','Different')))

ggplot(df.3, aes(x = key, y = value, fill = cluster)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=sex.cols)+
  labs(title = "Bootstrap Cluster Consistency K=2", x = "Bootstrap Sample", y = "Proportion")

View(bv.df)



