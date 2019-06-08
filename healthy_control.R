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

# note checkout lumi package

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')


Sys.setenv("plotly_username"="vetronz1992")
Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

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



# 
# ############ EDA ############
# #Dx
# ggplot(status[idx,], aes(status[idx,]$most_general, fill = status[idx,]$most_general)) +
#   scale_color_manual(values=dx.cols.6) +
#   labs(title="Barplot Diagnostics Groups",
#        x ="", y = "Count") +
#   scale_fill_discrete(name = "Dx") +
#   geom_bar()
# 
# addmargins(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))
# 
# #sex
# sex <- c('F', 'M')
# 
# bacterial <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[1,])
# greyb <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[2,])
# greyu <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[3,])
# greyv <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[4,])
# viral <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[5,])
# HC <- c(table(droplevels(status[idx,]$most_general), status[idx,]$Sex)[6,])
# 
# df <- data.frame(bacterial, greyb, greyu, greyv, HC, viral)
# df.2 <- mutate(df, sex = factor(c('F','M')))
# df.3 <- gather(df.2, dx, count, -sex)
# df.3
# ggplot(df.3, aes(x = dx, y = count, fill = sex)) +
#   geom_bar(position = "fill",stat = "identity")+
#   scale_fill_manual(values=sex.cols)+
#   labs(title = "Barplot Gender Proportions Within Diagnostic Groups", x = "Diagnosis", y = "Proportion")
# ggplotly(p)
# 
# # age
# ggplot(status[idx,], aes(x = status[idx,]$most_general, y = status[idx,]$Age..months., fill = status[idx,]$most_general)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot Age (months) Distributions by Diagnostic Groups",
#        x ="", y = "Age") +
#   scale_fill_discrete(name = "Dx") +
#   geom_boxplot()
# 
# ggplot(status[idx,], aes(x = status[idx,]$most_general, y = status[idx,]$Age..months., fill = status[idx,]$Sex)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot Age (months) Distributions by Gender",
#        x ="", y = "Age") +
#   scale_fill_manual(values=sex.cols, name = "Dx")+
#   geom_boxplot()
# 
# # inflam
# p1<-ggplot(status[idx,], aes(x = status[idx,]$most_general, y = status[idx,]$WBC, fill = status[idx,]$most_general)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot WBC Distributions by Diagnostic Group",
#        x ="", y = "WBC Count") +
#   scale_fill_discrete(name = "Dx") +
#   guides(fill=FALSE)+
#   geom_boxplot()
# p2<-ggplot(status[idx,], aes(x = status[idx,]$most_general, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = status[idx,]$most_general)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot CRP Distributions by Diagnostic Group",
#        x ="", y = "CRP Count") +
#   scale_fill_discrete(name = "Dx") +
#   geom_boxplot()
# grid.arrange(p1, p2, ncol=2)

############ LIMMA ############
e.set.f <- as.data.frame(e.set[,idx])
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)
class(e.set.f)

e.set.f$label <- as.character(status[idx, c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]
# e.set.f$site <- status[idx,c('site')]
e.set.f$age <- status[idx,c('Age..months.')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

# mean/variance calculations
x_var <- apply(e.set[,idx], 1, var)
x_mean <- apply(e.set[,idx], 1, mean)
df <- data.frame(log2(x_var), log2(x_mean))

# plot(log2(x_mean), log2(x_var), pch='.')
# abline(v=log2(5), col='red')

ggplot(df, aes(log2.x_mean., log2.x_var.)) +
  geom_vline(xintercept=log2(5))+
  geom_point(size = 0.2, stroke = 0, shape = 16)+
  labs(title="Mean Variance Scatter Plot",
       x ="log2 Mean Expressioin", y = "log2 Variance")

dim(e.set[,idx])
e.set[,idx][which(x_mean > 5)]
dim(t(e.set[,idx][which(x_mean > 5),]))
X <- e.set[,idx][which(x_mean > 5),]
dim(X)

### DESIGN MATRIX
# site
# design <- model.matrix(~label + sex + age + site + 0, data = e.set.f)
# colnames(design)<- c("bct","greyb","greyv", 'HC', 'vrl', 'sexM', 'age', 'CHW', 'EUC', 'FED', 'FPIES', 'KEN', 'OXF','SMH','SOT','UCSD')

design <- model.matrix(~label + sex + age + 0, data = e.set.f)
colnames(design)<- c("bct","greyb",'greyu', "greyv", 'HC', 'vrl', 'sexM', 'age')

design[1:5,]
dim(design)
colSums(design)

# contrast.matrix<- makeContrasts("bct-HC", 'vrl-HC', 'greyb-HC', levels=design)
# contrast.matrix<- makeContrasts("HC-bct", 'HC-vrl', levels=design)
contrast.matrix<- makeContrasts("((bct+vrl+greyb+greyv+greyu)/5)-HC", levels=design)
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

lfc <- 0.5
pval <- 0.05

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc)
dim(results)
head(results)
# summary(results)
vennDiagram(results, include = 'both')
# vennCounts(results, include = 'both')


# results.bct <- union(rownames(X[keep,])[results[,1] == 1]
#                       ,rownames(X[keep,])[results[,1] == -1])
# length(results.bct)
# results.vrl <- union(rownames(X[keep,])[results[,2] == 1]
#                       ,rownames(X[keep,])[results[,2] == -1])
# results.tot <- union(results.bct, results.vrl)
# length(results.tot)

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc)
top.hits

all.hits <- topTable(fit2, number=nrow(fit2))
# dim(top.hits)
dim(all.hits)

# intersect(results.tot, rownames(top.hits))

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
which(status$my_category_2 == 'bacterialgpos_19_SMH')
status$my_category_2[28]
idx[28]
idx[28] <- FALSE
idx[28]
sum(idx)
length(idx)

X <- e.set[,idx][which(x_mean > 5),]
dim(X)

# getwd()
# setwd('/home/patrick/Documents/RNA_seq_classifier/Data')
# write.csv(labels, file = "Mega_sub1_Demographic.csv", row.names=TRUE, na="")
# write.csv(status[idx,]$barcode_megaexp, file = "barcodes.csv", row.names=TRUE, na="")



status[idx,]$barcode_megaexp
results.tot <- ifelse(results[,1] == 0, FALSE, TRUE)
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

wbc<-status[idx,]$WBC
crp<-as.numeric(as.character(status[idx,]$array.contemporary.CRP))
bct <- status[idx,]$most_general == 'bacterial'

# most_gen 2D Age
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status[idx,]$most_general), size = status[idx,]$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups, Age Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p

# most_gen 2D CRP
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status[idx,]$most_general), size = as.numeric(as.character(status[idx,]$array.contemporary.CRP)),
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups, CRP Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p

# # api_create(p, filename = "2d_pca_filt")

# most_gen 3D
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~droplevels(status[idx,]$most_general), size = as.numeric(as.character(status[idx,]$array.contemporary.CRP)),
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
p
# api_create(p, filename = "3d_pca_filt")
# 
# # more_gen
# dx.cols.f <- c('#bd35fc', "#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
# plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$more_general, size = status[idx,]$Age..months.,
#         colors = c(dx.cols.f), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# 
# # category
# p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$category, size = status[idx,]$Age..months.,
#         colors=c(cat.pal), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Category Groups by PCA 1-2-3, Age Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# api_create(p, filename = "3d_pca_cat")
# 

# sex
# p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$Sex, size = ~status[idx,]$Age..months.,
#         colors = c(sex.cols), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Gender by PCA 1-2-3, Age Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# p
# api_create(p, filename = "3d_pca_sex")

# 
# # site
# plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$site, size = ~status[idx,]$Age..months.,
#         colors = c(site.pal), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Site Recruitment by PCA 1-2-3, Age Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))




### Clustering
X.t <- t(X[results.tot,])
dim(X.t)
# fviz_nbclust(X.t, kmeans, method = "wss")
# fviz_nbclust(X.t, kmeans, method = "silhouette")
# fviz_nbclust(X.t, kmeans, method = "gap_stat", nboot = 10)

# fviz_nbclust(X.t, cluster::fanny, method = "wss")
# fviz_nbclust(X.t, cluster::fanny, method = "silhouette")
# fviz_nbclust(X.t, cluster::fanny, method = "gap_stat")


### K2
n = 5
cols = gg_color_hue(n)
k2.pal <- c(cols[c(1,4)])
set.seed(47)
k2 <- kmeans(X.t, centers = 2, nstart = 100)
k2$cluster <- as.factor(k2$cluster)


setwd('/home/patrick/Documents/RNA_seq_classifier/Data')
clin <- read.table('Mega_sub1_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = TRUE)

# Construction: get the full demographics into dataframe
clin2 <- read.table('Mega_sub2_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
dim(clin2)
# a <- c(1, 19, 27, 28, 33, 34, 84, 85) # barcode, Dx,  wbc, abs_neut, perc_neut, perc_lymp, dx_1, dx_2
# clin2[clin2$V1 == '9423641116_B', a] # search clin2 using barcode
# View(clin)

k2.df <- status[idx,][, c('barcode_megaexp', 'category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]

k2.df$system <- clin$system
k2.df$system.spec <- clin$system_spec
k2.df$micro <- clin$micro
k2.df$sepsis <- clin$sepsis

k2.df$my_wbc <- clin$wbc
k2.df$abs_neut <- clin$abs_neut
k2.df$perc_neut <- clin$perc_neut
k2.df$perc_lymph <- clin$perc_lymph
k2.df$Path_1 <- clin$Path_1
k2.df$Path_2 <- clin$Path_2
k2.df$Path_3 <- clin$Path_3

k2.df$clus.2 <- k2$cluster

# barplot_dx_k2 <- table(k2$cluster, droplevels(status[idx,]$most_general))
# class(barplot_dx_k2)
# barplot_dx_k2$cluster <- seq.int(nrow(barplot_dx_k2))
# barplot_dx_k2
# data$ID <- seq.int(nrow(data))

# setwd('/home/patrick/Documents/RNA_seq_classifier/Data')
# write.csv(barplot_dx_k2, file = "barplot_dx_k2.csv",row.names=TRUE)

table(k2$cluster, droplevels(status[idx,]$more_general))

ggplot(status[idx,], aes(k2$cluster, fill=most_general)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  scale_fill_manual(values=dx.cols, 'Dx')+
  geom_bar()

# # positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother')
# p<-ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
#   labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
#   scale_x_discrete(limits = positions.more)+
#   scale_fill_manual(values=dx.cols.2, 'Cluster')+
#   geom_bar()
# ggplotly(p)

ggplot(status[idx,], aes(k2$cluster, Age..months., fill=k2$cluster)) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=cols[c(1,4)])+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(status[idx,], aes(x = k2$cluster, y = status[idx,]$WBC, fill = k2$cluster)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Dx')+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()

p2<-ggplot(status[idx,], aes(x = k2$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = k2$cluster)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Dx')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

ggplot(status[idx,], aes(x = k2$cluster, y = status[idx,]$WBC, fill = status[idx,]$most_general)) +
  scale_fill_manual(values=dx.cols, name = 'Dx')+
  labs(title="Boxplot WBC Distributions within Clusters",
       x ="Cluster", y = "WBC Count") +
  geom_boxplot()

ggplot(status[idx,], aes(x = k2$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = status[idx,]$most_general)) +
  scale_fill_manual(values=dx.cols, name = 'Dx')+
  labs(title="Boxplot CRP Distributions within Clusters",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()

ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.2, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=clus.2)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Dx')+
  labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Absolute Neutrophil Count") +
  geom_boxplot()

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.2, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=clus.2)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  labs(title="Boxplot WBC Percent Neutrophil Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.2, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=clus.2)) +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  labs(title="Boxplot WBC Percent Lymphocyte Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

# p<-plot_ly(status[idx,], x = k2$cluster, y = ~wbc, type = "box",
#            color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'WBC Count'))
# p
# 
# p<-plot_ly(status[idx,], x = k2$cluster, y = ~crp,
#            type = "box", color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'CRP Count'))
# p


# K2 analysis
table(k2$cluster, droplevels(status[idx,]$most_general))
table(k2$cluster, k2.df$system)
table(k2$cluster, k2.df$system.spec)
table(k2$cluster, k2.df$micro)
table(k2$cluster, k2.df$septic)

# system
ggplot(k2.df[bct,], aes(clus.2, fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

cluster <- c(1, 2)
clus1 <- table(k2.df$clus.2[bct], k2.df$system[bct])[1,]
clus2 <- table(k2.df$clus.2[bct], k2.df$system[bct])[2,]
df.1 <- data.frame(clus1, clus2)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Infection System Proportions by Cluster", x = "Cluster", y = "Proportion")

table(k2.df$clus.2[bct], k2.df$system[bct])
chisq.test(table(k2.df$clus.2[bct], k2.df$system[bct]))

# micro
ggplot(k2.df[bct,], aes(clus.2, fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

cluster <- c(1, 2)
clus1 <- table(k2.df$clus.2[bct], k2.df$Path_1[bct])[1,]
clus2 <- table(k2.df$clus.2[bct], k2.df$Path_1[bct])[2,]
df.1 <- data.frame(clus1, clus2)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")

table(k2.df$clus.2[bct], k2.df$Path_1[bct])
chisq.test(table(k2.df$clus.2[bct], k2.df$Path_1[bct]))


# # pca bacterial cluster 1
# plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df$clus.2,
#              colors=cols[c(1,4)], text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis),
#               symbol = ~k2.df$clus.2 == 2 & k2.df$most_general == 'bacterial', symbols = c('circle','x')) %>%
#   add_markers() %>%
#   layout(title = 'PCA of Diagnostic Groups',
#          xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#          yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
# 
# 
# plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df$clus.2,
#         colors=cols[c(1,4)], text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis),
#         symbol = ~k2.df$clus.2 == 2 & k2.df$most_general == 'bacterial', symbols = c('circle','x')) %>%
#   add_markers() %>%
#   layout(title = 'PCA of Diagnostic Groups',
#          xaxis = list(title = paste0("PC1: (", round(pve[3],2), '%)')),
#          yaxis = list(title = paste0("PC2: (", round(pve[4],2), '%)')))
# 







### K4
# k4.pal <- c('#09f70d', "#f76409", "#f70d09", '#2909f7')
set.seed(47)
k4 <- kmeans(X.t, centers = 4, nstart = 100)
k4$cluster <- as.factor(k4$cluster)
k2.df$clus.4 <- k4$cluster


table(k2.df$clus.4, droplevels(status[idx,]$most_general))

ggplot(status[idx,], aes(k4$cluster, fill=most_general)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  scale_fill_manual(values=dx.cols, 'Dx')+
  geom_bar()

ggplot(status[idx,], aes(k4$cluster, Age..months., fill=k4$cluster)) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=cols)+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(status[idx,], aes(x = k4$cluster, y = status[idx,]$WBC, fill = k4$cluster)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()

p2<-ggplot(status[idx,], aes(x = k4$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = k4$cluster)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)


ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.4, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=clus.4)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Absolute Neutrophil Count") +
  geom_boxplot()

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.4, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=clus.4)) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot WBC Percent Neutrophil Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.4, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=clus.4)) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot WBC Percent Lymphocyte Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

# system
ggplot(k2.df[bct,], aes(clus.4, fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

cluster <- c(1, 2, 3, 4)
clus1 <- table(k2.df$clus.4[bct], k2.df$system[bct])[1,]
clus2 <- table(k2.df$clus.4[bct], k2.df$system[bct])[2,]
clus3 <- table(k2.df$clus.4[bct], k2.df$system[bct])[3,]
clus4 <- table(k2.df$clus.4[bct], k2.df$system[bct])[4,]
df.1 <- data.frame(clus1, clus2, clus3, clus4)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Infection System Proportions by Cluster", x = "Cluster", y = "Proportion")


# micro
ggplot(k2.df[bct,], aes(clus.4, fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

clus1 <- table(k2.df$clus.4[bct], k2.df$Path_1[bct])[1,]
clus2 <- table(k2.df$clus.4[bct], k2.df$Path_1[bct])[2,]
clus3 <- table(k2.df$clus.4[bct], k2.df$Path_1[bct])[3,]
clus4 <- table(k2.df$clus.4[bct], k2.df$Path_1[bct])[4,]
df.1 <- data.frame(clus1, clus2, clus3, clus4)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3



ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")

plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df$clus.4,
             colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
              symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))


plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df$clus.2,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))














### K9
set.seed(47)
n = 9
cols = gg_color_hue(n)

k9 <- kmeans(X.t, centers = 9, nstart = 100)
k9$cluster <- as.factor(k9$cluster)
k2.df$clus.9 <- k9$cluster

table(k2.df$clus.9, droplevels(status[idx,]$most_general))

ggplot(status[idx,], aes(k9$cluster, fill=most_general)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=9", x = "", y = "Counts")+
  scale_fill_manual(values=dx.cols, 'Dx')+
  geom_bar()

ggplot(status[idx,], aes(k9$cluster, Age..months., fill=k9$cluster)) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  # scale_fill_manual(values=cols)+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(status[idx,], aes(x = k9$cluster, y = status[idx,]$WBC, fill = k9$cluster)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(status[idx,], aes(x = k9$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = k9$cluster)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.9, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=clus.9)) +
  scale_fill_manual(values=cols, name = 'Dx')+
  labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Absolute Neutrophil Count") +
  geom_boxplot()

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.9, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=clus.9)) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot WBC Percent Neutrophil Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.9, as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=clus.9)) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot WBC Percent Lymphocyte Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)

# system
ggplot(k2.df[bct,], aes(clus.9, fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

cluster <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
clus1 <- table(k2.df$clus.9[bct], k2.df$system[bct])[1,]
clus2 <- table(k2.df$clus.9[bct], k2.df$system[bct])[2,]
clus3 <- table(k2.df$clus.9[bct], k2.df$system[bct])[3,]
clus4 <- table(k2.df$clus.9[bct], k2.df$system[bct])[4,]
clus5 <- table(k2.df$clus.9[bct], k2.df$system[bct])[5,]
clus6 <- table(k2.df$clus.9[bct], k2.df$system[bct])[6,]
clus7 <- table(k2.df$clus.9[bct], k2.df$system[bct])[7,]
clus8 <- table(k2.df$clus.9[bct], k2.df$system[bct])[8,]
clus9 <- table(k2.df$clus.9[bct], k2.df$system[bct])[9,]

df.1 <- data.frame(clus1, clus2, clus3, clus4, clus5, clus6, clus7, clus8, clus9)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Infection System Proportions by Cluster", x = "Cluster", y = "Proportion")


# micro
p<-ggplot(k2.df[bct,], aes(clus.9, fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()
ggplotly(p)

clus1 <- table(k2.df$clus.9[bct], k2.df$Path_1[bct])[1,]
clus2 <- table(k2.df$clus.9[bct], k2.df$Path_1[bct])[2,]
clus3 <- table(k2.df$clus.9[bct], k2.df$Path_1[bct])[3,]
clus4 <- table(k2.df$clus.9[bct], k2.df$Path_1[bct])[4,]
df.1 <- data.frame(clus1, clus2, clus3, clus4)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")



k2.df$clus.9
plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df$clus.9,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))


plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df$clus.2,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))






# 
# table(k9$cluster, droplevels(status[idx,]$most_general))
# # table(k9$cluster, droplevels(status[idx,]$more_general))
# barplot_dx_k9 <- table(k9$cluster, droplevels(status[idx,]$most_general))
# 
# write.csv(barplot_dx_k4, file = "barplot_dx_k4.csv", row.names=TRUE)
# 
# ggplot(status[idx,], aes(k4$cluster, fill=most_general)) +
#   labs(title = "Barplot Diagnostic Groups by Cluster K=4", x = "Cluster", y = "Counts")+
#   scale_fill_manual(values=dx.cols, 'Dx')+
#   geom_bar()
# 
# # ggplot(status[idx,], aes(k4$cluster, Age..months., fill=most_general)) + geom_boxplot()+
# #   # scale_x_discrete(limits = positions)+
# #   xlab('Diagnostic Group') +
# #   ylab('Age') +
# #   scale_fill_manual(values=dx.cols, name='Dx')+
# #   ggtitle("Age (months) by Diagnostic Group, Split by Cluster")
# 
# #inflam
# p1<-ggplot(status[idx,], aes(x = k4$cluster, y = status[idx,]$WBC, fill = k4$cluster)) +
#   scale_fill_manual(values=dx.cols, name = 'Cluster')+
#   labs(title="Boxplot WBC Distributions by Cluster",
#        x ="Cluster", y = "WBC Count") +
#   guides(fill=FALSE)+
#   geom_boxplot()
# p2<-ggplot(status[idx,], aes(x = k4$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = k4$cluster)) +
#   scale_fill_manual(values=dx.cols, name = 'Cluster')+
#   labs(title="Boxplot CRP Distributions by Cluster",
#        x ="Cluster", y = "CRP Count") +
#   geom_boxplot()
# grid.arrange(p1, p2, ncol=2)
# 
# # k4 analysis
# 
# 
# table(k4$cluster, droplevels(status[idx,]$most_general))
# table(k4$cluster, k2.df$system)
# table(k4$cluster, k2.df$system.spec)
# table(k4$cluster, k2.df$micro)
# table(k4$cluster, k2.df$septic)
# 
# 
# # system
# cluster <- c(1, 2, 3, 4)
# clus1 <- table(k4$cluster, k2.df$system)[1,]
# clus2 <- table(k4$cluster, k2.df$system)[2,]
# clus3 <- table(k4$cluster, k2.df$system)[3,]
# clus4 <- table(k4$cluster, k2.df$system)[4,]
# 
# df.1 <- data.frame(clus1, clus2, clus3, clus4)
# df.1
# df.2 <- mutate(df.1, system=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -system)
# df.3
# 
# ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot Infection System by Cluster", x = "Cluster", y = "Proportion")
# 
# table(k4$cluster, k2.df$system)
# chisq.test(table(k4$cluster, k2.df$system))
# 
# # micro
# table(k4$cluster[micro], k2.df$micro[micro])
# chisq.test(k4$cluster[micro], k2.df$micro[micro])
# 
# cluster <- c(1, 2, 3, 4)
# clus1 <- table(k4$cluster[micro], k2.df$micro[micro])[1,]
# clus2 <- table(k4$cluster[micro], k2.df$micro[micro])[2,]
# clus3 <- table(k4$cluster[micro], k2.df$micro[micro])[3,]
# clus4 <- table(k4$cluster[micro], k2.df$micro[micro])[4,]
# 
# df.1 <- data.frame(clus1, clus2, clus3, clus4)
# df.1
# df.2 <- mutate(df.1, micro=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -micro)
# df.3
# 
# ggplot(df.3, aes(x = cluster, y = count, fill = micro)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")
# 
# ggplot(k2.df[micro,], aes(k2.df$clus.9[micro], fill=k2.df$micro[micro]))+
#   labs(title = "Barplot Microbiology by Cluster", x = "Diagnosis", y = "Counts")+
#   geom_bar()
# 
# ggplot(k2.df[bct,], aes(k2.df$clus.9[bct], fill=droplevels(k2.df$category[bct])))+
#   labs(title = "Barplot Microbiology by Cluster", x = "Diagnosis", y = "Counts")+
#   scale_fill_manual(values=cols[c(2,4)], name = 'Gram')+
#   geom_bar()



# 
# ### K7
# # k4.pal <- c('#09f70d', "#f76409", "#f70d09", '#2909f7')
# set.seed(47)
# k7 <- kmeans(X.t, centers = 7, nstart = 100)
# k7$cluster <- as.factor(k7$cluster)
# 
# table(k7$cluster, droplevels(status[idx,]$most_general))
# # table(k7$cluster, droplevels(status[idx,]$more_general))
# 
# ggplot(status[idx,], aes(k7$cluster, fill=most_general)) +
#   labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
#   scale_fill_manual(values=dx.cols, 'Dx')+
#   geom_bar()
# 
# ggplot(status[idx,], aes(k7$cluster, Age..months., fill=most_general)) + geom_boxplot()+
#   # scale_x_discrete(limits = positions)+
#   xlab('Diagnostic Group') +
#   ylab('Age') +
#   scale_fill_manual(values=dx.cols, name='Dx')+
#   ggtitle("Age (months) by Diagnostic Group, Split by Cluster")
# 
# #inflam
# ggplot(status[idx,], aes(x = k7$cluster, y = status[idx,]$WBC, fill = status[idx,]$most_general)) +
#   scale_fill_manual(values=dx.cols, name = 'Dx')+
#   labs(title="Boxplot Diagnosis Groups WBC Distributions",
#        x ="Cluster", y = "WBC Count") +
#   geom_boxplot()
# 
# ggplot(status[idx,], aes(x = k7$cluster, y = as.numeric(as.character(status[idx,]$array.contemporary.CRP)), fill = status[idx,]$most_general)) +
#   scale_fill_manual(values=dx.cols, name = 'Dx')+
#   labs(title="Boxplot CRP Distributions within Clusters",
#        x ="Cluster", y = "CRP Count") +
#   geom_boxplot()
# 
# 
# # k7 analysis
# k2.df$clus.7 <- k7$cluster
# 
# table(k7$cluster, droplevels(status[idx,]$most_general))
# table(k7$cluster, k2.df$system)
# table(k7$cluster, k2.df$system.spec)
# table(k7$cluster, k2.df$micro)
# table(k7$cluster, k2.df$septic)
# 
# 
# # system
# cluster <- c(1, 2, 3, 4, 5, 6, 7)
# clus1 <- table(k7$cluster, k2.df$system)[1,]
# clus2 <- table(k7$cluster, k2.df$system)[2,]
# clus3 <- table(k7$cluster, k2.df$system)[3,]
# clus4 <- table(k7$cluster, k2.df$system)[4,]
# clus5 <- table(k7$cluster, k2.df$system)[5,]
# clus6 <- table(k7$cluster, k2.df$system)[6,]
# clus7 <- table(k7$cluster, k2.df$system)[7,]
# 
# df.1 <- data.frame(clus1, clus2, clus3, clus4, clus5, clus6, clus7)
# df.1
# df.2 <- mutate(df.1, system=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -system)
# df.3
# 
# ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot of Gender Proportions Within Diagnostic Groups", x = "Diagnosis", y = "Proportion")
# 
# table(k7$cluster, k2.df$system)
# chisq.test(table(k7$cluster, k2.df$system))
# 
# # micro
# table(k7$cluster[micro], k2.df$micro[micro])
# chisq.test(k7$cluster[micro], k2.df$micro[micro])
# 
# clus1 <- table(k7$cluster[micro], k2.df$micro[micro])[1,]
# clus2 <- table(k7$cluster[micro], k2.df$micro[micro])[2,]
# clus3 <- table(k7$cluster[micro], k2.df$micro[micro])[3,]
# clus4 <- table(k7$cluster[micro], k2.df$micro[micro])[4,]
# clus5 <- table(k7$cluster[micro], k2.df$micro[micro])[5,]
# clus6 <- table(k7$cluster[micro], k2.df$micro[micro])[6,]
# clus7 <- table(k7$cluster[micro], k2.df$micro[micro])[7,]
# 
# df.1 <- data.frame(clus1, clus2, clus3, clus4, clus5, clus6, clus7)
# df.2 <- mutate(df.1, micro=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -micro)
# df.3
# 
# ggplot(df.3, aes(x = cluster, y = count, fill = micro)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")
# 
# ggplot(k2.df[micro,], aes(k2.df$clus.7[micro], fill=k2.df$micro[micro]))+
#   labs(title = "Barplot of Microbiology by Cluster", x = "Diagnosis", y = "Counts")+
#   geom_bar()















# end

