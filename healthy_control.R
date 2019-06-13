# library(tidyverse)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(dplyr)
library(plyr)
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


Sys.setenv("plotly_username"="vetronz1992")
Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

n = 10
cols.10 = gg_color_hue(n)

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- cols[1:5]
dx.cols.6 <- cols

dx.cols.f <- c("#ed0404", "#fc5716", '#16fc31','#38bc9d', '#55f1fc', '#0934f4', '#ac04f9', "#bc38ab")
sex.cols <- c('#fc1676', '#16acfc')
k7.pal <- c("#ed0404", "#fc5716", '#16fc31', '#55f1fc', '#0452f9', '#7304f9', '#f904e9')

positions <- c('bacterial', 'greyb', 'greyv', 'viral', 'greyu')
positions.f <- c('bacterial', 'greyb', 'greyv', 'greyu', 'adeno', 'flu', 'RSV', 'viralother')




# ### supervised
idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu' |
  status['most_general'] == 'HC'
sum(idx)
dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral', 'healthy_control') # supervised

### unsup
# idx <- status['most_general'] == 'bacterial' |
#   status['most_general'] == 'viral' |
#   status['most_general'] == 'greyb' |
#   status['most_general'] == 'greyv'|
#   status['most_general'] == 'greyu'
# sum(idx)
# dx <- c('bacterial', 'probable bacterial', 'unknown', 'probable viral', 'viral')


### outlier
# which(status[idx,]$my_category_2 == 'bacterialgpos_19_SMH')
# status[idx,]$my_category_2[which(status[idx,]$my_category_2 == 'bacterialgpos_19_SMH')]
# idx[which(status[idx,]$my_category_2 == 'bacterialgpos_19_SMH')]
# idx[which(status.idx$my_category_2 == 'bacterialgpos_19_SMH')] <- FALSE
# idx[which(status[idx,]$my_category_2 == 'bacterialgpos_19_SMH')]
# sum(idx)
# length(idx)

status.idx <- status[idx,]
dim(status.idx)

# rename most_general
status.idx$most_general <- as.character(status.idx$most_general)
status.idx$most_general[status.idx$most_general == 'greyb'] <- 'probable_bacterial'
status.idx$most_general[status.idx$most_general == 'greyu'] <- 'unknown'
status.idx$most_general[status.idx$most_general == 'greyv'] <- 'probable_viral'
status.idx$most_general[status.idx$most_general == 'HC'] <- 'healthy_control' # toggle for unsupervised
status.idx$most_general <- as.factor(status.idx$most_general)

levels(status.idx$most_general)
status.idx$most_general <- factor(status.idx$most_general, levels = dx)
levels(status.idx$most_general)

status.idx$array.contemporary.CRP <- as.numeric(as.character(status.idx$array.contemporary.CRP))

status.idx$most_general
# # ############ EDA ############
p<-ggplot(status.idx, aes(most_general, fill = most_general)) +
  scale_color_manual(values=dx.cols.6) +
  labs(title="Barplot Diagnostics Groups",
       x ="", y = "Count") +
  scale_fill_discrete(name = "Diagnostic Group") +
  geom_bar()
p
p+theme(axis.text=element_text(size=12))
# api_create(p, filename = "barplot_dx_groups")
addmargins(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))

#sex
sex <- c('F', 'M')
bacterial <- c(table(status.idx$most_general, status.idx$Sex)[1,])
probable_bacterial <- c(table(status.idx$most_general, status.idx$Sex)[2,])
unknown <- c(table(status.idx$most_general, status.idx$Sex)[3,])
probable_viral <- c(table(status.idx$most_general, status.idx$Sex)[4,])
viral <- c(table(status.idx$most_general, status.idx$Sex)[5,])
healthy_control <- c(table(status.idx$most_general, status.idx$Sex)[6,])

df <- data.frame(bacterial, probable_bacterial, unknown, probable_viral, viral, healthy_control)
df
df.2 <- mutate(df, sex = factor(c('F','M')))
df.2
df.3 <- gather(df.2, dx, count, -sex)
df.3$dx <- factor(df.3$dx, levels=dx)
levels(df.3$dx)

status.idx$most_general <- factor(status.idx$most_general, levels = dx)
p<-ggplot(df.3, aes(x = dx, y = count, fill = sex)) +
  geom_bar(position = "fill",stat = "identity")+
  scale_fill_manual(values=sex.cols)+
  labs(title = "Barplot Gender Proportions Within Diagnostic Groups", x = "Diagnosis", y = "Proportion")
p+theme(axis.text=element_text(size=12))

# # age
p<-ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$Age..months., fill = status.idx$most_general)) +
  # scale_fill_manual(values=dx.cols)+
  labs(title="Boxplot of Age (months) by Diagnostic Groups",
       x ="", y = "Age") +
  scale_fill_discrete(name = "Diagnostic Group") +
  geom_boxplot()
p+theme(axis.text=element_text(size=12))

ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$Age..months., fill = status.idx$Sex)) +
  # scale_fill_manual(values=dx.cols)+
  labs(title="Boxplot Age (months) Distributions by Gender",
       x ="", y = "Age") +
  scale_fill_manual(values=sex.cols, name = "Dx")+
  geom_boxplot()
#
# # inflam
p1<-ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$WBC, fill = status.idx$most_general)) +
  # scale_fill_manual(values=dx.cols)+
  labs(title="Boxplot WBC Distributions by Diagnostic Group",
       x ="", y = "WBC Count") +
  scale_fill_discrete(name = "Diagnostic Groups") +
  guides(fill=FALSE)+
  geom_boxplot()
p1 <- p1+theme(axis.text=element_text(size=11))
p2<-ggplot(status.idx, aes(x = status.idx$most_general, y = as.numeric(as.character(status.idx$array.contemporary.CRP)), fill = status.idx$most_general)) +
  # scale_fill_manual(values=dx.cols)+
  labs(title="Boxplot CRP Distributions by Diagnostic Group",
       x ="", y = "CRP Count") +
  scale_fill_discrete(name = "Diagnostic Groups") +
  geom_boxplot()
p2<-p2+theme(axis.text=element_text(size=8))
grid.arrange(p1, p2, ncol=2)

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
# e.set[,idx][which(x_mean > 5)]
dim(t(e.set[,idx][which(x_mean > 5),]))
X <- e.set[,idx][which(x_mean > 5),]
X.t <- t(X)


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

bootstraps <- list(c(0, 1), # 1 full
                   c(0.25, 0.1), # 2 7141
                   c(0.375, 0.1), # 3
                   c(0.5, 0.1), # 4 gap stat 6 (9 with p val 0.05)
                   c(0.75, 0.05), # 5 gap stat of three
                   c(1, 0.05)) # 6


boot <- 1
lfc <- bootstraps[[boot]][1]
pval <- bootstraps[[boot]][2]
lfc
pval

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc)
dim(results)
# head(results)
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

dim(status.idx)

# idx <- status['most_general'] == 'bacterial' |
#   status['most_general'] == 'probable bacterial' |
#   status['most_general'] == 'unknown' |
#   status['most_general'] == 'probable viral'|
#   status['most_general'] == 'viral'
# sum(idx)

# which(status$my_category_2 == 'bacterialgpos_19_SMH')
# # outlier
# idx[28] <- FALSE
# idx[28]
# sum(idx)
# length(idx)

dim(results)

results.tot <- ifelse(results[,1] == 0, FALSE, TRUE)

dim(X[results.tot,])

# # PCA
full.pca <- prcomp(X.t, scale=TRUE) # unsupervised
# full.pca <- prcomp(t(X[results.tot,]), scale=TRUE) # supervised

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]
#
# wbc<-status.idx$WBC
# crp<-as.numeric(as.character(status.idx$array.contemporary.CRP))

# # most_gen 2D Age
# p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status.idx$most_general), size = status.idx$Age..months.,
#              colors=c(dx.cols), text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'PCA of Diagnostic Groups, Age Size Mapping',
#          xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#          yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
# p

# # most_gen 2D
p <- plot_ly(pair1, x = ~PC1, y = ~PC2, color = status.idx$most_general, size = status.idx$array.contemporary.CRP,
             colors=cols, text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PC 1-2 of Diagnostic Groups, CRP Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p
p <- plot_ly(pair2, x = ~PC3, y = ~PC4, color = status.idx$most_general, size = status.idx$Age..months.,
             colors=cols, text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PC 3-4 of Diagnostic Groups, Age Size Mapping',
         xaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC4: (", round(pve[4],2), '%)')))
p

# api_create(p, filename = "2d_pca_filt")

## most_gen 3D
View(status.idx)
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$most_general, size = ~as.numeric(as.character(status.idx$array.contemporary.CRP)),
             colors=c(dx.cols), text= ~paste0('<br>age: ', status.idx$Age..months., '<br>Sex:', status.idx$Sex, '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
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
# plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$more_general, size = status.idx$Age..months.,
#         colors = c(dx.cols.f), text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
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
# View(status.idx)
# p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$Sex, size = ~status.idx$Age..months.,
#              colors=c(sex.cols), text= ~paste0('<br>age: ', status.idx$Age..months., '<br>Sex:', status.idx$Sex, '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
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
p<-fviz_nbclust(X.t, kmeans, method = "wss")
p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_tss_boot1")

p<-fviz_nbclust(X.t, kmeans, method = "silhouette")
p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_silhouette_boot1")

# p<-fviz_nbclust(X.t, kmeans, method = "gap_stat", nboot = 10)
# p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_gap_boot1")

# fviz_nbclust(X.t, cluster::fanny, method = "wss")
# fviz_nbclust(X.t, cluster::fanny, method = "silhouette")
# fviz_nbclust(X.t, cluster::fanny, method = "gap_stat")

# k2.df.old <- status[idx,][, c('barcode_megaexp', 'category', 'my_category_2', 'most_general', 'more_general', 'site',
                          # 'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]

k2.df <- status.idx[c('barcode_megaexp', 'category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]

# write.csv(k2.df, file = "k2.df.csv", row.names=TRUE)

setwd('/home/patrick/Documents/RNA_seq_classifier/Data')
clin <- read.table('Mega_sub1_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = TRUE)


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

# View(k2.df)

### K2
# n = 5
# cols = gg_color_hue(n)
k2.pal <- c(cols[c(1,4)])
set.seed(47)
k2 <- kmeans(X.t, centers = 2, nstart = 50)
k2$cluster <- as.factor(k2$cluster)

# Construction: get the full demographics into dataframe
# clin2 <- read.table('Mega_sub2_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
# dim(clin2)
# a <- c(1, 19, 27, 28, 33, 34, 84, 85) # barcode, Dx,  wbc, abs_neut, perc_neut, perc_lymp, dx_1, dx_2
# clin2[clin2$V1 == '9423641116_B', a] # search clin2 using barcode
# View(clin)


clus.boot <-paste0('clus.', boot, '.2')
clus.boot

k2.df$clus <- k2$cluster # have to assign using clus then rename it
colnames(k2.df)[ncol(k2.df)] <- clus.boot
k2.df[clus.boot]
colnames(k2.df)

table(k2$cluster, k2.df$most_general) # sanity check
table(k2.df[[clus.boot]], k2.df$most_general) # sanity check

cluster <- c(1, 2)
clus1<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[1,]
clus2<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[2,]

df.1 <- data.frame(clus1, clus2)
df.1

df.2 <- mutate(df.1, Diagnosis=factor(levels(k2.df$most_general)))
df.2
df.3 <- gather(df.2, cluster, count, -Diagnosis)
df.3

p<-ggplot(df.3, aes(x = cluster, y = count, fill = Diagnosis)) +
  geom_bar(position = "fill",stat = "identity")+
  labs(title = "Barplot of Diagnostic Group Proportions by Cluster", x = "Cluster", y = "Proportion")
ggplotly(p)
# api_create(p, filename = "barplot_dx_clus.1.2")

ggplot(k2.df, aes(k2$cluster, Age..months., fill=k2$cluster)) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=cols[c(1,4)], name = 'Cluster')+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = k2.df$WBC, fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols.10[c(2,7)], name = 'Cluster')+
  labs(title="Boxplot of WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  geom_boxplot()
p2<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = as.numeric(as.character(k2.df$array.contemporary.CRP)), fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols.10[c(2,7)], name = 'Cluster')+
  labs(title="Boxplot of CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
# grid.arrange(p1, p2, ncol=2)
ggplotly(p1)
ggplotly(p2)
# api_create(p1, filename = "boxplot_wbc_clus.1.2")
# api_create(p2, filename = "barplot_crp_clus.1.2")

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], y = k2.df[k2.df$most_general == 'bacterial',]$WBC, fill = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols.10[c(1,4)], name = 'Diagnostic Group')+
  labs(title="Boxplot of WBC Distributions for Definite Bacterial Cases by Cluster",
       x ="Cluster", y = "WBC Count") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], y = k2.df[k2.df$most_general == 'bacterial',]$array.contemporary.CRP, fill = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols.10[c(1,4)], name = 'Diagnostic Group')+
  labs(title="Boxplot of CRP Distributions for Definite Bacterial Cases by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
# grid.arrange(p1, p2, ncol=2)
ggplotly(p1)
ggplotly(p2)
# api_create(p1, filename = "boxplot_wbc_bct_clus.1.2")
# api_create(p2, filename = "barplot_crp_bct_clus.1.2")

# p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
#   scale_fill_manual(values=cols[c(2,7)], name = 'Cluster')+
#   labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
#        x ="Cluster", y = "Absolute Neutrophil Count") +
#   geom_boxplot()
# ggplotly(p)

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols.10[c(6,9)], name = 'Cluster')+
  guides(fill=FALSE)+
  labs(title="Boxplot of Percent Neutrophil Count for Definite Bacterials Cases",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols.10[c(6,9)], name = 'Cluster')+
  labs(title="Boxplot of Percent Lymphocyte Count for Definite Bacterials Cases",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
# grid.arrange(p1, p2, ncol=2)
ggplotly(p1)
ggplotly(p2)
# api_create(p1, filename = "boxplot_neut_bct_clus.1.2")
# api_create(p2, filename = "boxplot_lymp_bct_clus.1.2")


p<-plot_ly(status.idx, x = k2$cluster, y = ~status.idx$WBC, type = "box",
           color = ~k2$cluster, colors = c(dx.cols.2)) %>%
  layout(boxmode = "group",
         title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender',
         xaxis = list(title = 'Diagnosis'),
         yaxis = list(title = 'WBC Count'))
p

# 
# p<-plot_ly(status[idx,], x = k2$cluster, y = ~crp,
#            type = "box", color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'CRP Count'))
# p


# K2 analysis
table(k2.df[[clus.boot]], droplevels(k2.df$most_general))
table(k2.df[[clus.boot]], k2.df$system)
table(k2.df[[clus.boot]], k2.df$system.spec)
table(k2.df[[clus.boot]], k2.df$micro)
table(k2.df[[clus.boot]], k2.df$sepsis)

# system
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()
ggplotly(p)

cluster <- c(1, 2)
clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[2,]
df.1 <- data.frame(clus1, clus2)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

p<-ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Infection System Proportions by Cluster", x = "Cluster", y = "Proportion")
ggplotly(p)
# api_create(p, filename = "barplot_system_prop_clus.1.2")


# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial']))

# micro
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()
ggplotly(p)

cluster <- c(1, 2)
clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[2,]
df.1 <- data.frame(clus1, clus2)
df.1
df.2 <- mutate(df.1, system=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -system)
df.3

p<-ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")
ggplotly(p)
# api_create(p, filename = "barplot_micro_prop_clus.1.2")

# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial']))


# lookup bacterial outliers consistently assigned clus2
k2.df$most_general == 'bacterial' & k2.df[[clus.boot]] == 2
View(k2.df[k2.df$most_general == 'bacterial' & k2.df[[clus.boot]] == 2,])

plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df[[clus.boot]],
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$most_general == 'bacterial', 'bct', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))

plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df[[clus.boot]],
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$most_general == 'bacterial', 'bct', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[4],2), '%)')))


### TABLES FOR PRES ###
clus.1.2 <- k2.df[k2.df$most_general == 'bacterial' & k2.df$clus.1.2 == 2,][c('Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis', 'system', 'Path_1', 'Path_2' , 'clus.1.2')]
clus.1.2 <- clus.1.2[order(clus.1.2$system),]

rownames(clus.1.2) <- seq(1, nrow(clus.1.2))
plotly.table <- clus.1.2
# View(plotly.table)
colnames(plotly.table)
colnames(plotly.table) <- c('Age', 'Sex', 'WBC', 'CRP', 'Presentation', 'System', 'Micro_1', 'Micro_2', 'Cluster')

p <- plot_ly(
  type = 'table',
  columnorder = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9),
  columnwidth = c(25, 20, 20, 20, 20, 80, 30, 35, 30, 20),
  header = list(
    values = c("<b>Patients</b>", names(plotly.table)),
    align = c('left', rep('center', ncol(plotly.table))),
    line = list(width = 1, color = 'black'),
    fill = list(color = '#444444'),
    font = list(family = "Arial", size = 14, color = "white")
  ),
  cells = list(
    values = rbind(
      rownames(plotly.table), 
      t(as.matrix(unname(plotly.table)))
    ),
    align = c('left', rep('center', ncol(plotly.table))),
    line = list(color = "black", width = 1),
    fill = list(color = c('#9a9e9d')),
    font = list(family = "Arial", size = 12, color = c("black"))
  ))
p
api_create(p, filename = "table_clus.1.2")
# write.csv(clus.1.2, file = "clus.1.2.csv", row.names=TRUE)




############ K4 ############
set.seed(47)
k4 <- kmeans(X.t, centers = 4, nstart = 50)
k4$cluster <- as.factor(k4$cluster)

clus.boot <-paste0('clus.', boot, '.4')
clus.boot

k2.df$clus <- k4$cluster # have to assign using clus then rename it
colnames(k2.df)[ncol(k2.df)] <- clus.boot
k2.df[clus.boot]
colnames(k2.df)

table(k4$cluster, droplevels(k2.df$most_general)) # sanity check
table(k2.df[[clus.boot]], droplevels(k2.df$most_general)) # sanity check

View(k2.df)
# getwd()
# setwd('/home/patrick/Documents/RNA_seq_classifier/Data')

# system
p<-ggplot(k2.df, aes(k2.df[[clus.boot]], fill=most_general)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()
ggplotly(p)

cluster <- c(1, 2, 3, 4)
clus1<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[1,]
clus2<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[2,]
clus3<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[3,]
clus4<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[4,]

df.1 <- data.frame(clus1, clus2, clus3, clus4)
df.1
df.2 <- mutate(df.1, Diagnosis=factor(levels(k2.df$most_general)))
df.2
df.3 <- gather(df.2, cluster, count, -Diagnosis)
df.3

p<-ggplot(df.3, aes(x = cluster, y = count, fill = Diagnosis)) +
  geom_bar(position = "fill",stat = "identity")+
  labs(title = "Barplot of Diagnostic Group Proportions by Cluster", x = "Cluster", y = "Proportion")
ggplotly(p)
# api_create(p, filename = "barplot_dx_clus.1.4")

ggplot(k2.df, aes(k2.df[[clus.boot]], k2.df$Age..months., fill=k2.df[[clus.boot]])) + geom_boxplot()+
  # scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=cols, name = 'Cluster')+
  ggtitle("Age (months) by Cluster")

#inflam
p1<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = k2.df$WBC, fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols)+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = as.numeric(as.character(k2.df$array.contemporary.CRP)), fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()

# grid.arrange(p1, p2, ncol=2)
# api_create(p1, filename = "boxplot_wbc_clus.1.4")
# api_create(p2, filename = "barplot_crp_clus.1.4")

# p1<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = k2.df$WBC, fill = k2.df$most_general)) +
#   scale_fill_manual(values=cols, name = 'Cluster')+
#   labs(title="Boxplot WBC Distributions by Cluster",
#        x ="Cluster", y = "WBC Count") +
#   # guides(fill=FALSE)+
#   geom_boxplot()
# p2<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = as.numeric(as.character(k2.df$array.contemporary.CRP)), fill = k2.df$most_general)) +
#   scale_fill_manual(values=cols, name = 'Cluster')+
#   labs(title="Boxplot CRP Distributions by Cluster",
#        x ="Cluster", y = "CRP Count") +
#   geom_boxplot()
# grid.arrange(p1, p2, nrow=2)

ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols, name = 'Diagnosis')+
  labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Absolute Neutrophil Count") +
  geom_boxplot()

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  # guides(fill=FALSE)+
  labs(title="Boxplot WBC Percent Neutrophil Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Boxplot WBC Percent Lymphocyte Distributions of Definite Bacterials within Clusters",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
p<-grid.arrange(p1, p2, ncol=2)
# api_create(p1, filename = "barplot_neut_bct_clus.1.4")
# api_create(p2, filename = "barplot_lymph_bct_clus.1.4")
ggplotly(p)


p1<-plot_ly(k2.df[k2.df$most_general  == 'bacterial',], x = k2.df[[clus.boot]][k2.df$most_general  == 'bacterial'], y = as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), type = "box",
            color = k2.df[[clus.boot]][k2.df$most_general  == 'bacterial'], colors = cols) %>%
  layout(title = '',
         xaxis = list(title = ''),
         yaxis = list(title = 'Perc Lymp'))
p1
p2<-plot_ly(k2.df[k2.df$most_general  == 'bacterial',], x = k2.df[[clus.boot]][k2.df$most_general  == 'bacterial'], y = as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), type = "box",
            color = k2.df[[clus.boot]][k2.df$most_general  == 'bacterial'], colors = cols) %>%
  layout(title = '',
         xaxis = list(title = ''),
         yaxis = list(title = 'Perc Lymp'))
p2
p<- subplot(p1, p2)
p
# api_create(p, filename = "boxplot_differential_clus.1.4") # style > axis > lines > show

# 
# p<-plot_ly(status[idx,], x = k2$cluster, y = ~crp,
#            type = "box", color = ~k2$cluster, colors = c(dx.cols.2)) %>%
#   layout(boxmode = "group",
#          title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
#          xaxis = list(title = 'Diagnosis'),
#          yaxis = list(title = 'CRP Count'))
# p

# K4 analysis
table(k2.df[[clus.boot]], droplevels(k2.df$most_general))
table(k2.df[[clus.boot]], k2.df$system)
table(k2.df[[clus.boot]], k2.df$system.spec)
table(k2.df[[clus.boot]], k2.df$micro)
table(k2.df[[clus.boot]], k2.df$sepsis)

# system
ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=ifelse(k2.df[k2.df$most_general == 'bacterial',]$sepsis, 'septic', 'non-septic'))) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  scale_fill_manual(values=cols[c(4,1)], name = 'Sepsis')+
  geom_bar()

clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[2,]
clus3 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[3,]
clus4 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[4,]
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

# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial']))

# micro
ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=ifelse(k2.df[k2.df$most_general == 'bacterial',]$category == 'E', 'Gram +ve', 'Gram -ve'))) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  scale_fill_manual(values=cols[c(5,4)], name = 'Micro')+
  geom_bar()

ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=micro)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[2,]
clus3 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[3,]
clus4 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[4,]
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

# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial']))


# MENINGOCOCCAL ANALYSIS
plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df[[clus.boot]],
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$micro == 'meningococcal', 'meningococcal', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))

plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df[[clus.boot]],
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$micro == 'meningococcal', 'meningococcal', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[4],2), '%)')))


### TABLES FOR PRES ###
View(k2.df[k2.df$micro == 'meningococcal',])
clus.1.4 <- k2.df[k2.df$micro == 'meningococcal',][c('my_category_2', 'Diagnosis', 'system', 'Path_1', 'clus.1.2', 'clus.1.4')]
clus.1.4 <- clus.1.4[order(clus.1.4$clus.1.2),]
View(clus.1.4)
# write.csv(clus.1.4, file = "clus.1.4.csv", row.names=TRUE)




############ bootstraping ############
# write.csv(k2.df, file = "k2.df_bootstrapping.csv", row.names=TRUE) # original bootstrap sample file
# constructed conversion between clusters manually
k2.df$clus.3.4_con <- ifelse(k2.df$clus.3.4 == 1, 4, ifelse(k2.df$clus.3.4 == 2, 3, ifelse(k2.df$clus.3.4 == 3, 2, ifelse(k2.df$clus.3.4 == 4, 1, 0))))
k2.df$clus.4.4_con <- ifelse(k2.df$clus.4.4 == 1, 3, ifelse(k2.df$clus.4.4 == 2, 4, ifelse(k2.df$clus.4.4 == 3, 1, ifelse(k2.df$clus.4.4 == 4, 2, 0))))
k2.df$clus.5.4_con <- ifelse(k2.df$clus.5.4 == 1, 3, ifelse(k2.df$clus.5.4 == 2, 2, ifelse(k2.df$clus.5.4 == 3, 4, ifelse(k2.df$clus.5.4 == 4, 1, 0))))
k2.df$clus.6.4_con <- ifelse(k2.df$clus.6.4 == 1, 2, ifelse(k2.df$clus.6.4 == 2, 3, ifelse(k2.df$clus.6.4 == 3, 4, ifelse(k2.df$clus.6.4 == 4, 1, 0))))

tot <- 239
sum(k2.df$clus.1.4 == k2.df$clus.2.4)
tot - sum(k2.df$clus.1.4 == k2.df$clus.2.4)

sum(k2.df$clus.1.4 == k2.df$clus.3.4_con)
tot - sum(k2.df$clus.1.4 == k2.df$clus.3.4_con)

sum(k2.df$clus.1.4 == k2.df$clus.4.4_con)
tot - sum(k2.df$clus.1.4 == k2.df$clus.4.4_con)

sum(k2.df$clus.1.4 == k2.df$clus.5.4_con)
tot - sum(k2.df$clus.1.4 == k2.df$clus.5.4_con)

sum(k2.df$clus.1.4 == k2.df$clus.6.4_con)
tot - sum(k2.df$clus.1.4 == k2.df$clus.6.4_con)

b2 <- c(198, 41)
b3 <- c(176, 63)
b4 <- c(178, 61)
b5 <- c(161, 78)
b6 <- c(155, 84)
df.1 <- data.frame(b2, b3, b4, b5, b6)
df.1
df.2 <- gather(df.1, boot, count)
df.2

df.3 <- mutate(df.2, cluster.assigned=rep(c('same', 'different'), times = 5))

ggplot(df.3, aes(x = boot, y = count, fill = cluster.assigned)) +
  geom_bar(position = "fill",stat = "identity")+
  scale_fill_manual(values=cols[c(1,4)])+
  labs(title = "Barplot Proportions of Bootstrap Samples with Same Cluster Assignment", x = "Bootstrap Sample", y = "Proportion")









############ K9 ############
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
ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.9, fill=system)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()

cluster <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
clus1 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[2,]
clus3 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[3,]
clus4 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[4,]
clus5 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[5,]
clus6 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[6,]
clus7 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[7,]
clus8 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[8,]
clus9 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[9,]

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
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(clus.9, fill=Path_1)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
  geom_bar()
ggplotly(p)

clus1 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[2,]
clus3 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[3,]
clus4 <- table(k2.df$clus.9[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[4,]
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



# pca vis of clustering
plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df$clus.9,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$most_general == 'bacterial', 'bct', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))

plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df$clus.9,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~ifelse(k2.df$most_general == 'bacterial', 'bct', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[4],2), '%)')))


plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df$clus.9,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))

plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df$clus.9,
        colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
        symbol = ~k2.df$micro == 'meningococcal', symbols = c('circle','x')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC1: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[4],2), '%)')))


# lookup bacterial outliers
sum(k2.df$most_general == 'bacterial' & k2.df$clus.4 == 1) # bacterials assigned clus.4 == 1
# View(k2.df[k2.df$most_general == 'bacterial' & k2.df$clus.4 == 1,])
View(k2.df[k2.df$most_general == 'bacterial',])









### BARPLOT WITH COUNT VALUES GRAVEYARD ###
# have to use this ordering to ensure labels corespond to same color code in alphabetical ggplot
dx <- c('viral', 'unknown', 'probable viral', 'probable bacterial', 'bacterial', 'viral', 'unknown', 'probable viral', 'probable bacterial', 'bacterial')
dx <- rep(levels(k2.df$most_general), times=2)

clus1<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[1,]
clus2<-table(k2.df[[clus.boot]], droplevels(k2.df$most_general))[2,]

df.1 <- data.frame(clus1, clus2)
df.1
df.2 <- gather(df.1, key='cluster', value = 'count')
df.2
df.3 <- mutate(df.2, Diagnosis = dx)
df.3

df.4 <- ddply(df.3, "cluster",
              transform, tot=cumsum(count))
df.4
ggplot(data=df.4, aes(x=cluster, y=count, fill=Diagnosis)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=tot, label=count))








# end

