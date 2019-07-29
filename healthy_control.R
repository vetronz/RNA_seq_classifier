library(tidyverse)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(plyr)
# library(plotly)
library(e1071)
library(neuralnet)
library(ROCR)
library(tidyr)
library("illuminaHumanv4.db")
library(Mfuzz)
library(caret)
library(class)
library(randomForest)
library(dplyr)

getwd()
setwd('/home/patrick/Code/R')
rm(list=setdiff(ls(), 'all'))
load('esets.RData')

setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')
cyber.s <- read.table('CIBERSORT.Output_Job14.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data')
clin <- read.table('Mega_sub1_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = TRUE)


# Sys.setenv("plotly_username"="vetronz1992")
# Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

# Sys.setenv("plotly_username"="vetronz")
# Sys.setenv("plotly_api_key"="OhacJkwCAaZOcC0wHPhp")



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

n = 8
cols.8 = gg_color_hue(n)

n = 10
cols.10 = gg_color_hue(n)

n = 14
cols.14 = gg_color_hue(n)

dx.cols <- c('#D81D22', '#FF3A3A', '#AD4187' , '#776BB9' , '#4A8AC3', '#32A46A')
sex.cols <- c('#fc1676', '#16acfc')

clus.cols <- c('#FFDD38' , '#56DD5F', '#6763CF', '#FF5338')


###### Cyber Sort ######
dim(status)
dim(cyber.s)
status[1:5,1:8]
cyber.s[1:5,1:5]
colnames(cyber.s)[1] <- 'my_category_2'

# join the status and cibersort dataframes together on the common column
# join rather than merge preserves the row order
status <- join(status, cyber.s, by='my_category_2')
dim(status)

# select data from clin cases where we have neutrophil perc and compare to ciber sort pred
df.2 <- clin[!is.na(clin$perc_neut),c('category', 'perc_neut')]
colnames(df.2)[1] <- 'my_category_2'
match(df.2$my_category_2, status$my_category_2)
df.3 <- status[match(df.2$my_category_2, status$my_category_2),c('my_category_2', 'Neutrophils')]
df.4 <- merge(df.3, df.2)

ggplot(df.4, aes(x=Neutrophils, y=perc_neut)) + 
  geom_point()+
  geom_smooth(method=lm, level=0.95)

cor(df.4$Neutrophils, df.4$perc_neut)

# ### unsup
# idx <- status['most_general'] == 'bacterial' |
#   status['most_general'] == 'viral' |
#   status['most_general'] == 'greyb' |
#   status['most_general'] == 'greyv'|
#   status['most_general'] == 'greyu'
# sum(idx)
# dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral')

### supervised
idx <- status$most_general == 'bacterial' |
  status$most_general == 'viral' |
  status$most_general == 'greyb' |
  status$most_general == 'greyv'|
  status$most_general == 'greyu' |
  status$most_general == 'HC'
sum(idx)
class(idx)
dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral', 'healthy_control') # supervised

### outlier
which(status$my_category_2 == 'bacterialgpos_19_SMH')
idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')]
idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')] <- FALSE
idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')]

status.idx <- status[idx,]
status.idx$most_general[1:10]

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
# status.idx$most_general

status.idx$array.contemporary.CRP <- as.numeric(as.character(status.idx$array.contemporary.CRP))

dim(status.idx)
# 
# # # ############ EDA ############
# p<-ggplot(status.idx, aes(most_general, fill = most_general)) +
#   scale_fill_manual(values=dx.cols, 'Diagnostic Groups')+
#   labs(title="Barplot Diagnostics Groups",
#   x ="", y = "Count") +
#   geom_bar()
# p
# p+theme(axis.text=element_text(size=12))
# p<-ggplotly(p)
# # api_create(p, filename = "barplot_dx_breakdown")
# addmargins(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))
# 
# #sex
# sex <- c('F', 'M')
# bacterial <- c(table(status.idx$most_general, status.idx$Sex)[1,])
# probable_bacterial <- c(table(status.idx$most_general, status.idx$Sex)[2,])
# unknown <- c(table(status.idx$most_general, status.idx$Sex)[3,])
# probable_viral <- c(table(status.idx$most_general, status.idx$Sex)[4,])
# viral <- c(table(status.idx$most_general, status.idx$Sex)[5,])
# # healthy_control <- c(table(status.idx$most_general, status.idx$Sex)[6,])
# 
# df <- data.frame(bacterial, probable_bacterial, unknown, probable_viral, viral)
# # df <- data.frame(bacterial, probable_bacterial, unknown, probable_viral, viral, healthy_control)
# df
# df.2 <- mutate(df, sex = factor(c('F','M')))
# df.2
# df.3 <- gather(df.2, dx, count, -sex)
# df.3$dx <- factor(df.3$dx, levels=dx)
# levels(df.3$dx)
# 
# status.idx$most_general <- factor(status.idx$most_general, levels = dx)
# p<-ggplot(df.3, aes(x = dx, y = count, fill = sex)) +
#   geom_bar(position = "fill",stat = "identity")+
#   scale_fill_manual(values=sex.cols)+
#   labs(title = "Barplot Gender Proportions Within Diagnostic Groups", x = "Diagnosis", y = "Proportion")
# p+theme(axis.text=element_text(size=12))
# p <- ggplotly(p)
# p
# # api_create(p, filename = "barplot_dx_sex_dist")
# 
# # # age
# p<-ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$Age..months., fill = status.idx$most_general)) +
#   scale_fill_manual(values=dx.cols, name = "Diagnostic Group")+
#   labs(title="Boxplot of Age (months) by Diagnostic Groups",
#        x ="", y = "Age") +
#   geom_boxplot()
# p<-p+theme(axis.text=element_text(size=12))
# p
# api_create(p, filename = "box_whisker_age")
# 
# ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$Age..months., fill = status.idx$Sex)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot Age (months) Distributions by Gender",
#        x ="", y = "Age") +
#   scale_fill_manual(values=sex.cols, name = "Dx")+
#   geom_boxplot()
# #
# # # inflam
# p1<-ggplot(status.idx, aes(x = status.idx$most_general, y = status.idx$WBC, fill = status.idx$most_general)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot WBC Distributions by Diagnostic Group",
#        x ="", y = "WBC Count") +
#   scale_fill_manual(values=dx.cols, name = "Diagnostic Groups") +
#   geom_boxplot()
# p1 <- p1+theme(axis.text=element_text(size=11))
# p2<-ggplot(status.idx, aes(x = status.idx$most_general, y = as.numeric(as.character(status.idx$array.contemporary.CRP)), fill = status.idx$most_general)) +
#   # scale_fill_manual(values=dx.cols)+
#   labs(title="Boxplot CRP Distributions by Diagnostic Group",
#        x ="", y = "CRP Count") +
#   scale_fill_manual(values=dx.cols, name = "Diagnostic Groups") +
#   geom_boxplot()
# p2<-p2+theme(axis.text=element_text(size=8))
# grid.arrange(p1, p2, ncol=2)
# # api_create(p1, filename = "box_whisker_wbc")
# # api_create(p2, filename = "box_whisker_crp")





############ LIMMA ############
dim(e.set)
# View(head(e.set))
e.set.f <- e.set[,idx]
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)

# check that the patient order has not been changed by the join step
# e.set.f defined using the raw idx whereas we are appending columns from status.idx which have been processed
status[idx,]$my_category_2[1:10]
status.idx$my_category_2[1:10]

e.set.f$label <- status.idx$most_general
e.set.f$sex <- status.idx$Sex
e.set.f$age <- status.idx$Age..months.
e.set.f$Neutrophils <- status.idx$Neutrophils
e.set.f$Monocytes <- status.idx$Monocytes
e.set.f$NK.cells.resting <- status.idx$NK.cells.resting
e.set.f$NK.cells.activated <- status.idx$NK.cells.activated
e.set.f$T.cells.CD8 <- status.idx$T.cells.CD8
# e.set.f$site <- status[idx,c('site')]

# check that we have added 8 columns
dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-9):ncol(e.set.f)]

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
# X <- e.set[,idx] # full unfiltered 47323 genes
X <- e.set[,idx][which(x_mean > 5),]
X.t <- t(X)
dim(X.t)

### DESIGN MATRIX
# colnames(design)<- c("bct","greyb","greyv", 'HC', 'vrl', 'sexM', 'age', 'CHW', 'EUC', 'FED', 'FPIES', 'KEN', 'OXF','SMH','SOT','UCSD')

# design <- model.matrix(~label + sex + age + 0, data = e.set.f)
# colnames(design)<- c("bct","greyb",'greyu', "greyv", 'vrl', 'HC', 'sexM', 'age')

design <- model.matrix(~label + sex + age + round(Neutrophils, 3) +
                         round(Monocytes, 3) + round(NK.cells.resting, 3) +
                         round(NK.cells.activated, 3) + round(T.cells.CD8, 3) + 0,
                       data = e.set.f)
colnames(design)<- c("bct","greyb",'greyu', "greyv", 'vrl', 'HC', 'sexM', 'age', 'neut', 'mono', 'nk.rest', 'nk.act', 'CD8')

# double check design with status to ensure correct cibersort values in design
design[1:5,]
status.idx[1:5,c('Neutrophils', 'Monocytes')]

# check sums of design
dim(design)
colSums(design)

contrast.matrix <- makeContrasts("bct-vrl", levels=design)
# contrast.matrix<- makeContrasts("HC-bct", 'HC-vrl', levels=design)
# contrast.matrix<- makeContrasts("((bct+vrl+greyb+greyv+greyu)/5)-HC", levels=design)
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
fit2 <- eBayes(fit2[keep,], trend = TRUE) # same result with or without the keep filter as we have pre-processed
dim(fit2)
plotSA(fit2)

bootstraps <- list(c(0, 1), # 1 full
                   c(0.15, 0.2), # 2 6628
                   c(0.25, 0.1), # 3 3054
                   c(0.375, 0.1), # 4
                   c(0.5, 0.1), # 5 gap stat 6 (9 with p val 0.05)
                   c(0.75, 0.05), # 6 gap stat of three
                   c(1, 0.05), # 7 gap stat of three
                   c(1.25, 0.001), # 8 gap stat of three
                   c(2, 0.0001)) # 9


boot <- 7
lfc <- bootstraps[[boot]][1]
pval <- bootstraps[[boot]][2]
# pval <- 5.e-5
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
top.hits[1:5,]

all.hits <- topTable(fit2, number=nrow(fit2))
# dim(top.hits)
dim(all.hits)
all.hits[1:5,]

all.filt <- all.hits[abs(all.hits$logFC) > lfc & all.hits$adj.P.Val < pval,]
# all.filt[1:5,]
# all.filt <- all.hits
dim(all.filt)

p<-ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size = 1, stroke = 0, shape = 16) +
  # text = ~paste("<br>Ensembl: ", ensemb, '<br>Gene: ', gene)
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  labs(title="Volcano Plot of Log Fold Change Against -log10 P Value Boot=5",
       x ="Log Fold Change", y = "log10 P-value")
p

# p<-plot_ly(all.filt, x=~logFC, y=~-log10(adj.P.Val),
#            text = ~paste('<br>Gene: ', gene, '<br>Ensembl: ', ensemb),
#            type='scatter', mode = "markers", color = ~logThresh)
# add_segments(x = lfc, xend = lfc, y = 0, yend = 25, name = paste0('lfc >:', lfc)) %>%
# add_segments(x = -lfc, xend = -lfc, y = 0, yend = 25, name = paste0('lfc <:', -1*lfc))
# p
# api_create(p, filename = "volcano_b5")


####### SELECTION OF TRANSCRIPTS #######
# load('esets.RData')
# saveRDS(X.t, "X.t.rds")
# rm(X.t)
# X.t <- readRDS("X.t.rds")

dim(results)
results.tot <- ifelse(results[,1] == 0, FALSE, TRUE)
dim(X[results.tot,])
X.diff <- X[results.tot,]
# filter out the healthy controls
idx.d <- status.idx$most_general != 'healthy_control'
X.r <- t(X.diff[,idx.d])
dim(X.r)


###### PCA ######
full.pca <- prcomp(X.r, scale=TRUE) # unsupervised
# full.pca <- prcomp(t(X[results.tot,]), scale=TRUE) # supervised

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]


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
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$most_general, size = status.idx$Age..months.,
             colors=c(dx.cols), text= ~paste0('<br>age: ', status.idx$Age..months., '<br>Sex:', status.idx$Sex, '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
p
# api_create(p, filename = "3d_pca_dx_unfilt")

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




############ Clustering ############
# p<-fviz_nbclust(X.r, kmeans, method = "wss", k.max = 15)
# p
# p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_tss_boot1")

# p<-fviz_nbclust(X.r, kmeans, method = "silhouette", k.max = 15)
# p
# p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_silhouette_boot1")

# p<-fviz_nbclust(X.r, kmeans, method = "gap_stat", nboot = 10)
# p
# p<-ggplotly(p)
# api_create(p, filename = "opt_cluster_gap_boot1")


# 
# ###### FUZZY CLUSTERING ######
# dim(X.r)
# df.1 <- X.r
# 
# ### MFUZZ
# #save it to a temp file so ti doesnt clutter up my blog directory
# tmp <- tempfile()
# write.table(df.1, file=tmp, sep='\t', quote = F, col.names=NA)
# 
# #read it back in as an expression set
# data <- table2eset(file=tmp)
# 
# data.s <- standardise(data)
# 
# m1 <- mestimate(data.s)
# m1
# 
# # built in estimation of opt clust number
# # Dmin(data.s, m=m1, crange=seq(2,22,1), repeats=5, visu=TRUE) # reccomend ~ 10 clusters
# 
# k.max <- 15
# x <- as.matrix(df.1)
# diss <- stats::dist(x)
# v <- rep(0, k.max)
# v1 <- rep(0, k.max)
# 
# for(i in 2:k.max){
#   print(paste0('iter: ', i))
#   c<-mfuzz(data.s, c = i, m=m1)
#   v[i] <- c$withinerror
#   ss <- silhouette(c$cluster, diss)
#   v1[i] <- summary(ss)[4][[1]]
# }
# plot(v[-1])
# v1
# plot(v1[-1])
# 
# 
# gap_stat <- clusGap(X.r, FUN = fanny, nstart = 10,
#                     K.max = 10, B = 10)
# 
# # mfuzz uses cmeans
# # validation using cmeans directly
# 
# # Get total within sum of square
# .get_withinSS <- function(d, cluster){
#   d <- stats::as.dist(d)
#   cn <- max(cluster)
#   clusterf <- as.factor(cluster)
#   clusterl <- levels(clusterf)
#   cnn <- length(clusterl)
#   
#   if (cn != cnn) {
#     warning("cluster renumbered because maximum != number of clusters")
#     for (i in 1:cnn) cluster[clusterf == clusterl[i]] <- i
#     cn <- cnn
#   }
#   cwn <- cn
#   # Compute total within sum of square
#   dmat <- as.matrix(d)
#   within.cluster.ss <- 0
#   for (i in 1:cn) {
#     cluster.size <- sum(cluster == i)
#     di <- as.dist(dmat[cluster == i, cluster == i])
#     within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size
#   }
#   within.cluster.ss
# }
# 
# k.max <- 15
# v <- rep(0, k.max)
# for (i in 2:k.max) {
#   print(paste0('iter: ', i))
#   clust <- cmeans(df.1, centers = i, iter.max=30, verbose=FALSE, dist="euclidean",
#                   method="cmeans", m=m1, rate.par = NULL)
#   v[i] <- .get_withinSS(diss, clust$cluster)
# }
# plot(v[-1])
# 
# 
# 
# # distance metrix investigations
# x <- mtcars["Honda Civic",] 
# y <- mtcars["Camaro Z28",] 
# dist(rbind(x, y))
# 
# # custom functions to evaluate the cluster assignments
# # norm vec to calc euclidian dist and then 2BitBios version which i think is actually for sum sq error
# norm_vec <- function(x) sqrt(sum(x^2))
# norm_vec(as.numeric(y))-norm_vec(as.numeric(x)) # think the dist function performs euclid dist on matrix
# 
# 
# 
# 
# 
# 
# 
# 
# 
# #### clinical df construction
# k2.df <- status.idx[status.idx$most_general != 'healthy_control',c('barcode_megaexp', 'category', 'my_category_2', 'most_general', 'more_general', 'site',
#                                                                    'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
# dim(k2.df)
# dim(clin)
# # write.csv(k2.df, file = "k2.df.csv", row.names=TRUE)
# 
# k2.df$system <- clin$system
# k2.df$system.spec <- clin$system_spec
# k2.df$micro <- clin$micro
# k2.df$sepsis <- clin$sepsis
# 
# k2.df$my_wbc <- clin$wbc
# k2.df$abs_neut <- clin$abs_neut
# k2.df$perc_neut <- clin$perc_neut
# k2.df$perc_lymph <- clin$perc_lymph
# k2.df$Path_1 <- clin$Path_1
# k2.df$Path_2 <- clin$Path_2
# k2.df$Path_3 <- clin$Path_3
# 
# k2.df$Path_1[is.na(k2.df$Path_1)] <- 'unknown'
# k2.df$Path_1[k2.df$Path_1 == ''] <- 'unknown'
# # View(k2.df)
# 
# 
# 
# ### K2
# # n = 5
# # cols = gg_color_hue(n)
# k2.pal <- c(cols[c(1,4)])
# set.seed(47)
# # k2 <- kmeans(X.r, centers = 2, nstart = 25)
# # k2$tot.withinss
# 
# k2 <- cmeans(X.r, centers = 2, iter.max=30, verbose=FALSE, dist="euclidean",
#              method="cmeans", m=1.3, rate.par = NULL)
# 
# k2$membership
# k2$cluster <- as.factor(k2$cluster)
# 
# clus.boot <-paste0('clus.', boot, '.2')
# clus.boot
# 
# k2.df$clus <- k2$cluster # have to assign using clus then rename it
# colnames(k2.df)[ncol(k2.df)] <- clus.boot
# k2.df[clus.boot]
# colnames(k2.df)
# 
# table(k2$cluster, k2.df$most_general) # sanity check
# dx_clus.1.2 <- addmargins(table(k2.df[[clus.boot]], k2.df$most_general))
# dx_clus.1.2
# # write.csv(dx_clus.1.2, file = "dx_clus.1.2.csv", row.names=TRUE)
# 
# 
# p<-ggplot(k2.df, aes(k2.df[[clus.boot]], fill=most_general)) +
#   labs(title = "Barplot of Diagnostic Groups by Cluster K=2", x = "", y = "Counts")+
#   scale_fill_manual(values=dx.cols, 'Diagnostic Groups')+
#   geom_bar()
# p
# # p + guides(fill=guide_legend(title="Diagnostic Groups"))






###### PREDICTION ######
scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

# setup the disease only demographic matrix
idx.d <- status.idx$most_general != 'healthy_control'
status.idx.d <- status.idx[idx.d,]
table(status.idx.d$most_general)


X.r <- t(X.diff[,idx.d])
X.r <- as.data.frame(X.r)
dim(X.r)

# scale the transcript data
# X.s <- data.frame(X.r)
X.s<-data.frame(apply(X.r, 2, scale01))
dim(X.s)


# add the class labels to the transcript data
# X.s$bacterial <- status.idx.d$most_general == 'bacterial'
# X.s$p.b <- status.idx.d$most_general == 'probable_bacterial'
X.s$bct <- NULL
table(status.idx.d$most_general)
X.s$bct <- status.idx.d$most_general == 'bacterial'
sum(X.s$bct)

prop1 <- 0.8

logistic.m <- NULL
knn.m <- NULL
randForrest.m <- NULL
nn.m <- NULL
svm.m <- NULL
for (i in 1:16){
  print(paste0('iter: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  
  # factor labels required for some ml models 
  train.f <- train
  test.f <- test
  train.f$bct <- factor(train.f$bct)
  test.f$bct <- factor(test.f$bct)
  
  # dim(train)
  # dim(test)
  
  ### LOGISTIC REGRESSION
  model <- glm(bct~ ., data=train, family=binomial(link='logit'), maxit = 500)
  # summary(model)
  # anova(model, test="Chisq")
  
  p <- predict(model, test[-ncol(test)])
  pr <- prediction(p, test$bct)
  
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  
  logistic.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  # pr <- ifelse(p > 0.5,1,0)
  # misClasificError <- mean(pr != test$bct)
  # print(paste('Accuracy',round(1-misClasificError, 3)))

  #
  ### KNN
  # opt neighbours
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 10)
  model <- train(bct ~., data = train.f, method = "knn",
                 trControl=trctrl,
                 tuneLength = 10)

  # ggplot(model$results, aes(k, Accuracy)) + geom_point()

  knn.opt <- model$results$k[which.max(model$results$Accuracy)]
  # train
  p <- knn(train[-ncol(train)], test[-ncol(test)], train$bct,  k=knn.opt, prob=TRUE)

  # display the confusion matrix
  table(p, test$bct)
  # attributes(pred)
  p<-attr(p, "prob")
  p<-1-p

  # plot(model, print.thres = 0.5, type="S")
  # confusionMatrix(p, test$bct)
  pr <- prediction(p, test$bct)

  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)

  knn.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  
  ### RANDOM FORREST
  model <- randomForest(bct ~ . , data = train.f)
  # plot(model)
  # attributes(model)
  # model$mtry # number of variables considered by each tree
  pred<-predict(model , test.f[-ncol(test)])
  # attributes(pred)
  # model=randomForest(x,y,xtest=x,ytest=y,keep.forest=TRUE)
  model.prob <- predict(model, test.f, type="prob")
  p <- 1-model.prob[,1]
  pr <- prediction(p, test$bct)
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  randForrest.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  ### Neural Net
  nn <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic")
  p <- predict(nn, test[-ncol(test)])
  p
  pr <- prediction(p, test$bct)
  
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  
  nn.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  
  ### SVM
  model <- svm(bct ~ . , train.f, probability = TRUE)
  pred <- predict(model, test.f, probability = TRUE)
  p<-attr(pred, "prob")
  pr <- prediction(1-p[,1], test$bct)
  # plot(prf)
  svm.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

# df.1 <- as.data.frame(cbind(logistic.m, knn.m, randForrest.m, nn.m, svm.m))
df.1 <- as.data.frame(cbind(logistic.m, knn.m, randForrest.m, nn.m, svm.m))

df.2 <- gather(df.1, learner, roc)
df.2$roc <- unlist(df.2$roc)
df.2$learner <- factor(df.2$learner)
ggplot(df.2, aes(learner, roc, color = learner)) + geom_boxplot()
ggplot(df.2, aes(roc, fill= learner, color = learner)) + geom_density( alpha=0.1)+
  labs(title=paste0('Roc Area Density with pval: ', pval, ' and lfc: ', lfc),
     x ="density", y = "roc area")


detach("package:plyr", unload=TRUE)
df.2 %>%
  group_by(learner) %>%
  summarise(roc.m = mean(roc), roc.med = median(roc), roc.sd = sd(roc))

# ggplot(df.2, aes(learner, roc, color = learner)) + geom_jitter()


###### HYPER PARAM OPTIMIZATION FOR PROMISING MODELS ######
# remove pseudo labels
status.idx.d <- status.idx[idx.d,]
X.s$bct <- status.idx.d$most_general == 'bacterial'
table(status.idx.d$most_general)
sum(X.s$bct)

prop1 <- 4/5
prop2 <- 3/4
# subdivide into train and test set

set.seed(44)
index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
train <- X.s[index, ]
test <- X.s[-index, ]

dim(train)
dim(test)


### train, val, test split function 
# split01 <- function(df, train.prop, val.prop, test.prop){
#   spec = c(train = train.prop, validate = val.prop, test = test.prop)
#   
#   split <- cut(
#     seq(nrow(df)), 
#     nrow(df)*cumsum(c(0,spec)),
#     labels = names(spec)
#   )
#   
#   g <- sample(split, replace = FALSE)
#   split.res <- split(df, g) # split returns a list of dataframes
# }
# res <- split01(X.s, 0.6, 0.2, 0.2)
# sapply(res, nrow)/nrow(X.s) # check props


### NEURAL NET
# optimize hidden nodes, activation function (logistic over tanh), cost function (sse over ce)
h.n <- 5
boot <- 32
roc.a <- NULL
roc.t <- NULL
for(i in 1:h.n){
  print(paste0('hidden Nodes: ', i))
  for(j in 1:boot){
    index <- sample(nrow(train), round(prop2*nrow(train)))
    train.cv <- train[index,]
    test.cv <- train[-index,]
    
    nn1 <- neuralnet(bct~ ., train.cv, linear.output = FALSE, act.fct = "logistic",
                     hidden = c(1), rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
    pred <- predict(nn1, test.cv[-ncol(test.cv)])
    
    roc.a[j] <- prediction(pred, test.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
  }
  roc.t <- append(roc.t, roc.a)
}

df <- data.frame(matrix(unlist(roc.t), nrow=length(roc.t), byrow=T))
colnames(df) <- 'roc.A'
df$h.n <- as.factor(sort(rep(seq(1:h.n), boot)))
# df

df.1 <- df %>%
  group_by(h.n) %>%
  summarise(roc.m = mean(roc.A), roc.med = median(roc.A), roc.sd = sd(roc.A))

ggplot(df, aes(roc.A, color=h.n, fill = h.n)) + geom_density(alpha=0.1)
df.1

mixed.max <- which.max(df.1$roc.med + df.1$roc.m/2)
med.max <- which.max(df.1$roc.med)

if(mixed.max == med.max){
  opt.h.n <- which.max(df.1$roc.med)
  print(paste0('optimal hidden nodes set to: ', opt.h.n))
} else {
  print(paste0('mixed max: ', mixed.max, ' median max: ', med.max))
  opt.h.n <- mixed.max
}


# weight initialization. Cant get neuralnet package to accept so just using random init
# xavier.w <- rnorm(dim(train)[1], mean = 0, sd = sqrt(1/dim(train)[1]))
# hist(xavier.w)

### TEST SET NEURAL
# opt.h.n <- 2
boot <- 128
nn.test <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  nn1 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn1, test[-ncol(test)])
  
  nn.test[i] <- prediction(pred, test$bct) %>%
    performance(measure = "auc") %>%
    .@y.values
}

pr <- prediction(pred, test$bct)
prf <- performance(pr, measure = "tpr", x.measure = "fpr")
plot(prf)
# plot(nn1)

mean(unlist(nn.test))
median(unlist(nn.test))
sd(unlist(nn.test))



### RANDOM FORREST
# change the train.f and test sets to factors
train.f$bct <- factor(train.f$bct)
test$bct <- factor(test$bct)

# index.val <- sample(nrow(train.f), round(prop2*nrow(train.f)))
# train <- train.f[index.val,]
# val <- train.f[-index.val,]

model <- randomForest(bct ~ . , data = train.f, nodesize = 1, maxnodes = NULL)
plot(model)
attributes(model)

model$mtry
# [1] 12

opt.tree <- which.min(model$err.rate[,1])
opt.tree
# hyper.tree <- opt.tree * 2
hyper.tree <- round(opt.tree * 1.5)

model$err.rate[which.min(model$err.rate[,1])]

hyper_grid <- NULL

hyper_grid <- expand.grid(
  mtry       = seq(7, 14, by = 1),
  node_size  = seq(1, 9, by = 2),
  sampe_size = c(.54, .57, .6, .632, .66, .69, .72)
)

head(hyper_grid)
dim(hyper_grid)

for(i in 1:nrow(hyper_grid)) {
  # train model
  model <- randomForest(
    formula = bct ~ ., 
    data = train.f, 
    ntree = hyper.tree,
    mtry = hyper_grid$mtry[i],
    nodesize = hyper_grid$node_size[i],
    samplesize = hyper_grid$sampe_size[i]
    # seed = 123
  )
  # extract error
  hyper_grid$OOB_MSE[i] <- model$err.rate[which.min(model$err.rate[,1])]
}

hyper_top <- hyper_grid %>% 
  dplyr::arrange(OOB_MSE) %>%
  head(10)
hyper_top
mtry.opt <- hyper_top[1,][[1]]
node.opt <- hyper_top[1,][[2]]
sample.opt <- hyper_top[1,][[3]]

### TEST SET EVAL FORREST
boot <- 32
rf.test <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  model <- randomForest(formula = bct ~ .,
                        data = train.f,
                        ntree = hyper.tree,
                        mtry = mtry.opt,
                        nodesize = node.opt,
                        samplesize = sample.opt,
                        maxnodes = NULL)
  p <- predict(model, test, type="prob")
  pr <- prediction(p[,2], test$bct)
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  
  rf.test[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

mean(unlist(rf.test))
median(unlist(rf.test))

mean(unlist(nn.test))
median(unlist(nn.test))

df.1 <- as.data.frame(matrix(unlist(c(nn.test, rf.test))), nrow=2, byrow=T)
colnames(df.1) <- 'roc.a'
df.1$algo <- as.factor(sort(rep(seq(1:2), boot)))
ggplot(df.1, aes(roc.a, color = algo, fill=algo)) + geom_density(alpha=.1)





####### PSEUDO LABELING #######
# remove pseudo labels
status.idx.d <- status.idx[idx.d,]
X.s$bct <- status.idx.d$most_general == 'bacterial'
table(status.idx.d$most_general)
print(paste0('bacterial cases: ', sum(X.s$bct)))


# record nn.test with full dataset rather than just train set
nn.b.test <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  
  nn1 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn1, test[-ncol(test)])
  pr <- prediction(pred, test$bct)
  nn.b.test[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

roc.h <- NULL
iters <- 60
for (i in 1:iters){
  print(paste0('iter: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  nn1 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn1, test[-ncol(test)])
  pr <- prediction(pred, test$bct)
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")

  roc.h[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  ### checks passing test index to demographic data selects same pts as pred
  # as.character(status.idx.d[-index,]$my_category_2) == rownames(as.data.frame(pred))
  
  # construct pb filter, select pb from status
  pb.filter <- status.idx.d[-index,]$most_general == 'probable_bacterial'
  pred[pb.filter,]
  
  if(max(pred[pb.filter,]) > 0.99){
    # select most probable bacterial case and store as ppv
    which.max(pred[pb.filter,])
    ppb <- names(which.max(pred[pb.filter,]))
    
    # select ppb form demographic subset using ppb
    # print(paste0('PPB case: ', status.idx.d[-index,]$my_category_2[status.idx.d[-index,]$my_category_2 == ppb], ', ', status.idx.d[-index,]$most_general[status.idx.d[-index,]$my_category_2 == ppb]))
    status.idx.d[-index,]$most_general[status.idx.d[-index,]$my_category_2 == ppb] <- 'bacterial'
    # print(paste0('PPB case: ', status.idx.d[-index,]$my_category_2[status.idx.d[-index,]$my_category_2 == ppb], ', ', status.idx.d[-index,]$most_general[status.idx.d[-index,]$my_category_2 == ppb]))
    Sys.sleep(1)
  }
}

table(status.idx.d$most_general)
# add these pseudo labeled ppb cases to the X.s
X.s$bct <- status.idx.d$most_general == 'bacterial'
print(paste0('bacterial cases: ', sum(X.s$bct)))

### PSEUDO NETWORK TEST
nn.psd.test <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  
  nn1 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn1, test[-ncol(test)])
  
  pr <- prediction(pred, test$bct)
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  nn.psd.test[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

df.1 <- as.data.frame(matrix(unlist(c(nn.test, nn.b.test, nn.psd.test))), nrow=3, byrow=T)
colnames(df.1) <- 'roc.a'
df.1$algo <- as.factor(sort(rep(seq(1:3), boot)))
df.1$algo <- ifelse(df.1$algo == 1, 'normal', ifelse(df.1$algo == 2, 'normal_boot', 'pseudo_label'))
ggplot(df.1, aes(roc.a, color = algo, fill=algo)) + geom_density(alpha=.1)


mean(unlist(nn.b.test))
median(unlist(nn.b.test))
sd(unlist(nn.b.test))

mean(unlist(nn.psd.test))
median(unlist(nn.psd.test))
sd(unlist(nn.psd.test))


###### IRIS VALIDATION ######

# discrepancy in dimension of transcript and label matrix
dim(e.set.i)
dim(status.iris)

X.i <- t(e.set.i)
dim(X.i)

# extract the overlap
common.my_cat <- intersect(rownames(X.i), status.iris$My_code)
length(common.my_cat)

# find position of common.my_cat in status.iris
com.idx <- match(common.my_cat, status.iris$My_code)

# pass com.idx to filter the status matrix
dim(status.iris[com.idx,])

status.i.idx <- status.iris[com.idx,]

dim(status.i.idx)
dim(X.i)

# remove non relevent patients
status.i.idx$most_general
i.idx<- status.i.idx$most_general == 'bacterial' |
  status.i.idx$most_general == 'viral' |
  status.i.idx$most_general == 'greyb' |
  status.i.idx$most_general == 'greyv'|
  status.i.idx$most_general == 'greyu'
sum(i.idx)

status.i.idx <- status.i.idx[i.idx,]
X.i <- X.i[i.idx,]

dim(status.i.idx)
dim(X.i)

# rename most_general
status.i.idx$most_general <- as.character(status.i.idx$most_general)
status.i.idx$most_general[status.i.idx$most_general == 'greyb'] <- 'probable_bacterial'
status.i.idx$most_general[status.i.idx$most_general == 'greyu'] <- 'unknown'
status.i.idx$most_general[status.i.idx$most_general == 'greyv'] <- 'probable_viral'
status.i.idx$most_general[status.i.idx$most_general == 'HC'] <- 'healthy_control'

status.i.idx$most_general <- as.factor(status.i.idx$most_general)
dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral')
levels(status.i.idx$most_general)
status.i.idx$most_general <- factor(status.i.idx$most_general, levels = dx)
levels(status.i.idx$most_general)

# transcript processing
dim(X.i)
x.i.var <- apply(X.i, 2, var)
x.i.mean <- apply(X.i, 2, mean)
df <- data.frame(log2(x.i.var), log2(x.i.mean))
colnames(df) <- c('V1', 'V2')

ggplot(df, aes(V1, V2)) +
  geom_vline(xintercept=log2(5))+
  geom_point(size = 0.2, stroke = 0, shape = 16)+
  labs(title="Mean Variance Scatter Plot",
       x ="log2 Mean Expressioin", y = "log2 Variance")









###### LEARNING CURVE ######
j_train <- NULL
j_test <- NULL
roc.train <- NULL
roc.test <- NULL
roc.train.me <- NULL
roc.test.me <- NULL
learning_curve.df <- NULL
p.h <- NULL
k <- 30

for(p in 15:85){
  prop <- p/100
  print(prop)
  for(i in 1:k){
    index <- sample(nrow(X.s), round(prop*nrow(X.s)))
    train_cv <- X.s[index, ]
    test_cv <- X.s[-index, ]
    
    nn_cv <- neuralnet(bct~ ., train_cv, linear.output = FALSE, act.fct = "logistic", hidden = opt.h.n)
    pred_train <- predict(nn_cv, train_cv[-ncol(train_cv)])
    pred_test <- predict(nn_cv, test_cv[-ncol(test_cv)])
    
    j_train[i] <- prediction(pred_train[,1], train_cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
    j_test[i] <- prediction(pred_test[,1], test_cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values  
  }
  
  # class.calls <- c(j_test)
  class.calls <- c(j_train, j_test)
  full.list <- class.calls
  full.df <- data.frame(matrix(unlist(full.list), nrow=length(full.list), byrow=T))
  colnames(full.df) <- 'roc.A'
  
  full.df$class <- as.factor(sort(rep(seq(1:(length(class.calls)/k)), k)))
  
  roc.stats <- full.df %>%
    group_by(class) %>%
    summarise(roc.m = mean(roc.A), roc.v = var(roc.A), roc.sd = sd(roc.A))
  
  roc.stats <- roc.stats %>% mutate(
    roc.se = roc.sd/sqrt(k),
    z.stat = qnorm(0.975),
    roc.me = z.stat * roc.se
  )
  p.h[p] <- p
  roc.train[p] <- roc.stats$roc.m[1]
  roc.train.me[p] <- roc.stats$roc.me[1]
  roc.test[p] <- roc.stats$roc.m[2]
  roc.test.me[p] <- roc.stats$roc.me[2]
  
}

learning_curve.df <- as.data.frame(cbind(roc.train, roc.test, roc.train.me, roc.test.me, p.h))
colnames(learning_curve.df) <- c('train', 'test', 'train.me', 'test.me', 'perc')
learning_curve.df

ggplot(learning_curve.df, aes(x=p.h, y=train)) +
  geom_line(aes(y=train, color='train'))+
  geom_errorbar(aes(ymin=train-train.me, ymax=train+train.me), width=0.1)+
  geom_line(aes(y=test, color='test'))+
  geom_errorbar(aes(ymin=test-test.me, ymax=test+test.me), width=0.1)+
  labs(title=paste0('Learning Curve with ', opt.h.n, ' hidden nodes'), x ="training Data Percentage", y = "ROCA")




########################END#####################################

#inflam
p1<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = k2.df$WBC, fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols.10[c(6,9)], name = 'Cluster')+
  labs(title="Boxplot of WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = as.numeric(as.character(k2.df$array.contemporary.CRP)), fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols.10[c(6,9)], name = 'Cluster')+
  labs(title="Boxplot of CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)
# ggplotly(p1)
# ggplotly(p2)
# api_create(p1, filename = "boxplot_wbc_clus.1.2")
# api_create(p2, filename = "barplot_crp_clus.1.2")

p1<-ggplot(k2.df[k2.df$more_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$more_general == 'bacterial'], y = k2.df$WBC[k2.df$more_general == 'bacterial'], fill=k2.df$more_general[k2.df$more_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Diagnostic Group')+
  labs(title="Boxplot of WBC Counts in Bacterial Cases",
       x ="Cluster", y = "WBC Count") +
  guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$more_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$more_general == 'bacterial'], y = k2.df$array.contemporary.CRP[k2.df$more_general == 'bacterial'], fill=k2.df$more_general[k2.df$more_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Diagnostic Group')+
  labs(title="Boxplot of CRP Counts in Bacterial Cases",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)
# api_create(p1, filename = "boxplot_wbc_bct_clus.1.2")
# api_create(p2, filename = "barplot_crp_bct_clus.1.2")

# p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'],
#                                                           y = k2.df[k2.df$most_general == 'bacterial',]$WBC, fill = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
#   scale_fill_manual(values=cols.10[c(1,4)], name = 'Diagnostic Group')+
#   labs(title="Boxplot of WBC Distributions for Definite Bacterial Cases by Cluster",
#        x ="Cluster", y = "WBC Count") +
#   geom_boxplot()
# p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], y = k2.df[k2.df$most_general == 'bacterial',]$array.contemporary.CRP, fill = k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
#   scale_fill_manual(values=cols.10[c(1,4)], name = 'Diagnostic Group')+
#   labs(title="Boxplot of CRP Distributions for Definite Bacterial Cases by Cluster",
#        x ="Cluster", y = "CRP Count") +
#   geom_boxplot()
# grid.arrange(p1, p2, ncol=2)
# ggplotly(p1)
# ggplotly(p2)
# api_create(p1, filename = "boxplot_wbc_bct_clus.1.2")
# api_create(p2, filename = "barplot_crp_bct_clus.1.2")

# p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$abs_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
#   scale_fill_manual(values=cols[c(2,7)], name = 'Cluster')+
#   labs(title="Boxplot Absolute Neutrophil Count Distributions of Definite Bacterials within Clusters",
#        x ="Cluster", y = "Absolute Neutrophil Count") +
#   geom_boxplot()
# ggplotly(p)

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=k2.df$most_general[k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Cluster')+
  guides(fill=FALSE)+
  labs(title="Percent Neutrophil Count for Definite Bacterials Cases",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=k2.df$most_general[k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Cluster')+
  labs(title="Percent Lymphocyte Count for Definite Bacterials Cases",
       x ="Cluster", y = "Lymphocyte Percent") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)
# ggplotly(p1)
# ggplotly(p2)
# api_create(p1, filename = "boxplot_neut_bct_clus.1.2")
# api_create(p2, filename = "boxplot_lymp_bct_clus.1.2")



### system
table(k2.df[k2.df$most_general == 'bacterial',][[clus.boot]], k2.df[k2.df$most_general == 'bacterial',]$system)
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=system)) +
  labs(title = "Barplot of Microbioloty by Cluster in Bacterial Cases K=2", x = "", y = "Counts")+
  scale_fill_manual(values = cols.14[c(2:14)], name = "System")+
  geom_bar()
p<-ggplotly(p)
p

# cluster <- c(1, 2)
# clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[1,]
# clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])[2,]
# df.1 <- data.frame(clus1, clus2)
# df.1
# df.2 <- mutate(df.1, system=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -system)
# df.3
# 
# p<-ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot of Infection System in Bacterial Cases K=2", x = "Cluster", y = "Proportion")

# ggplotly(p)
# api_create(p, filename = "barplot_system_prop_clus.1.2")


# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial']))

# micro
table(k2.df[k2.df$most_general == 'bacterial',][[clus.boot]], k2.df[k2.df$most_general == 'bacterial',]$Path_1)
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=Path_1)) +
  labs(title = "Barplot of Microbioloty by Cluster in Bacterial Cases K=2", x = "", y = "Counts")+
  scale_fill_manual(values = cols.14[c(2:14)], name = "Microbe")+
  geom_bar()
p
# api_create(p, filename = "barplot_micro_clus.1.2")


# 
# cluster <- c(1, 2)
# clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[1,]
# clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])[2,]
# df.1 <- data.frame(clus1, clus2)
# df.1
# df.2 <- mutate(df.1, system=factor(rownames(df.1)))
# df.2
# df.3 <- gather(df.2, cluster, count, -system)
# df.3

# p<-ggplot(df.3, aes(x = cluster, y = count, fill = system)) +
#   geom_bar(position = "fill",stat = "identity")+
#   # scale_fill_manual(values=dx.cols.f)+
#   labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")
# p
# ggplotly(p)
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
# clus.1.2 <- k2.df[k2.df$most_general == 'bacterial' & k2.df$clus.1.2 == 2,][c('Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis', 'system', 'Path_1', 'Path_2' , 'clus.1.2')]
clus.1.2 <- k2.df[k2.df$most_general == 'bacterial',][c('Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis', 'system', 'Path_1', 'Path_2' , 'clus.1.2')]
clus.1.2 <- clus.1.2[order(clus.1.2$Cluster),]


rownames(clus.1.2) <- seq(1, nrow(clus.1.2))
colnames(clus.1.2) <- c('Age', 'Sex', 'WBC', 'CRP', 'Presentation', 'System', 'Path_1', 'Path_2', 'Cluster')
# write.csv(clus.1.2, file = "clus.1.2.csv", row.names=TRUE)
plotly.table <- clus.1.2
View(plotly.table)
colnames(plotly.table)

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
# api_create(p, filename = "table_clus.1.2")









############ K4 ############
set.seed(47)
k4 <- kmeans(X.r, centers = 4, nstart = 25)
k4$cluster <- as.factor(k4$cluster)

clus.boot <-paste0('clus.', boot, '.4')
clus.boot

k2.df$clus <- k4$cluster # have to assign using clus then rename it
colnames(k2.df)[ncol(k2.df)] <- clus.boot
k2.df[clus.boot]
colnames(k2.df)

table(k4$cluster, droplevels(k2.df$most_general)) # sanity check
table(k2.df[[clus.boot]], droplevels(k2.df$most_general)) # sanity check
addmargins(table(k2.df[[clus.boot]], droplevels(k2.df$most_general)))
dx_clus.1.4 <- addmargins(table(k2.df[[clus.boot]], droplevels(k2.df$most_general)))
# write.csv(dx_clus.1.4, file = "dx_clus.1.4.csv", row.names=TRUE)

# View(k2.df)
# getwd()
# setwd('/home/patrick/Documents/RNA_seq_classifier/Data')

p<-ggplot(k2.df, aes(k2.df[[clus.boot]], fill=most_general)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster K=4", x = "", y = "Counts")+
  scale_fill_manual(values = dx.cols, name = 'Diagnostic_Groups')+
  geom_bar()
p<-ggplotly(p)
p
# api_create(p, filename = "barplot_dx_clus.1.4")

# ggplotly(p)
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
p
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
  scale_fill_manual(values=cols.10[c(6,9,3,7)], name='Cluster')+
  labs(title="Boxplot WBC Distributions by Cluster",
       x ="Cluster", y = "WBC Count") +
  # guides(fill=FALSE)+
  geom_boxplot()
p2<-ggplot(k2.df, aes(x = k2.df[[clus.boot]], y = as.numeric(as.character(k2.df$array.contemporary.CRP)), fill = k2.df[[clus.boot]])) +
  scale_fill_manual(values=cols.10[c(6,9,3,7)], name='Cluster')+
  labs(title="Boxplot CRP Distributions by Cluster",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)
# api_create(p1, filename = "boxplot_wbc_clus.1.4")
# api_create(p2, filename = "barplot_crp_clus.1.4")


p1<-ggplot(k2.df[k2.df$more_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$more_general == 'bacterial'], y = k2.df$WBC[k2.df$more_general == 'bacterial'], fill=k2.df$more_general[k2.df$more_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Diagnostic Group')+
  labs(title="Boxplot of WBC Counts in Bacterial Cases",
       x ="Cluster", y = "WBC Count") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$more_general == 'bacterial',], aes(x = k2.df[[clus.boot]][k2.df$more_general == 'bacterial'], y = k2.df$array.contemporary.CRP[k2.df$more_general == 'bacterial'], fill=k2.df$more_general[k2.df$more_general == 'bacterial'])) +
  scale_fill_manual(values=dx.cols, name = 'Diagnostic Group')+
  labs(title="Boxplot of CRP Counts in Bacterial Cases",
       x ="Cluster", y = "CRP Count") +
  geom_boxplot()
grid.arrange(p1, p2, ncol=2)
# api_create(p1, filename = "boxplot_wbc_bct_clus.1.4")
# api_create(p2, filename = "barplot_crp_bct_clus.1.4")

p1<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_neut), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  # guides(fill=FALSE)+
  labs(title="Percent Neutrophil Count for Definite Bacterials Cases",
       x ="Cluster", y = "Neutrophil Percentage") +
  geom_boxplot()
p2<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], as.numeric(k2.df[k2.df$most_general == 'bacterial',]$perc_lymph), fill=k2.df[[clus.boot]][k2.df$most_general == 'bacterial'])) +
  scale_fill_manual(values=cols, name = 'Cluster')+
  labs(title="Percent Lymphocyte Count for Definite Bacterials Cases",
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

# system
addmargins(table(k2.df[k2.df$most_general == 'bacterial',][[clus.boot]], k2.df[k2.df$most_general == 'bacterial',]$system))
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=system)) +
  labs(title = "Barplot of Infection System by Cluster in Bacterial Cases K=4", x = "", y = "Counts")+
  scale_fill_manual(values = cols.10[c(4:10)])+
  geom_bar()
p<-ggplotly(p)
p
# api_create(p, filename = "barplot_system_clus.1.4")

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
  labs(title = "Barplot of Organ System Infection Proportions by Cluster", x = "Cluster", y = "Proportion")

# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$system[k2.df$most_general == 'bacterial']))

# micro
ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=ifelse(k2.df[k2.df$most_general == 'bacterial',]$category == 'E', 'Gram +ve', 'Gram -ve'))) +
  labs(title = "Barplot of Microbioloty by Cluster in Bacterial Cases K=4", x = "", y = "Counts")+
  scale_fill_manual(values = cols.14[c(2,10)], name = "Microbe")+
  geom_bar()

addmargins(table(k2.df[k2.df$most_general == 'bacterial',][[clus.boot]], k2.df[k2.df$most_general == 'bacterial',]$Path_1))
p<-ggplot(k2.df[k2.df$most_general == 'bacterial',], aes(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], fill=Path_1)) +
  labs(title = "Barplot of Microbioloty by Cluster in Bacterial Cases K=4", x = "", y = "Counts")+
  scale_fill_manual(values = cols.14[c(2:14)], name = "Microbe")+
  geom_bar()
p<-ggplotly(p)
p
# api_create(p, filename = "barplot_micro_clus.1.4")

clus1 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[1,]
clus2 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[2,]
clus3 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[3,]
clus4 <- table(k2.df[[clus.boot]][k2.df$most_general == 'bacterial'], k2.df$micro[k2.df$most_general == 'bacterial'])[4,]
df.1 <- data.frame(clus1, clus2, clus3, clus4)
df.1
df.2 <- mutate(df.1, micro=factor(rownames(df.1)))
df.2
df.3 <- gather(df.2, cluster, count, -micro)
df.3

ggplot(df.3, aes(x = cluster, y = count, fill = micro)) +
  geom_bar(position = "fill",stat = "identity")+
  # scale_fill_manual(values=dx.cols.f)+
  labs(title = "Barplot of Microbiology Proportions by Cluster", x = "Diagnosis", y = "Proportion")

# table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial'])
# chisq.test(table(k2.df$clus.2[k2.df$most_general == 'bacterial'], k2.df$Path_1[k2.df$most_general == 'bacterial']))

sum(k2.df$Path_1 == 'meningococcus')

# MENINGOCOCCAL ANALYSIS
p<-plot_ly(pair1, x = ~PC1, y = ~PC2, color = ~k2.df[[clus.boot]],
           colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
           symbol = ~ifelse(k2.df$micro == 'meningococcal', 'meningococcal', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PC 1-2 Meningococcal Cluster Assignment',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
p
# api_create(p, filename = "2d_pca_mening_clus.1.4")

p<-plot_ly(pair2, x = ~PC3, y = ~PC4, color = ~k2.df[[clus.boot]],
           colors=cols, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$Age..months., '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis),
           symbol = ~ifelse(k2.df$micro == 'meningococcal', 'meningococcal', 'other'), symbols = c('x','circle')) %>%
  add_markers() %>%
  layout(title = 'PCA of Diagnostic Groups',
         xaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)')),
         yaxis = list(title = paste0("PC4: (", round(pve[4],2), '%)')))
p
# api_create(p, filename = "2d_pca_mening_2_clus.1.4")

p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z=~PC3, color = ~k2.df[[clus.boot]],
           colors=cols, size = k2.df$array.contemporary.CRP, text= ~paste0('category: ', k2.df$category, '<br>age: ', k2.df$array.contemporary.CRP, '<br>WBC: ', k2.df$WBC, '<br>CRP: ', as.numeric(as.character(k2.df$array.contemporary.CRP)), '<br>label:',k2.df$my_category_2, '<br>Micro: ', k2.df$Path_1, '<br>Diagnosis: ',k2.df$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'PC Plot of Cluster Assignment, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
p
# api_create(p, filename = "3d_pca_mening_clus.1.4")




### TABLES FOR PRES ###
# View(k2.df[k2.df$micro == 'meningococcal',])
sum(k2.df$Path_1 == 'meningococcus')
clus.1.4 <- k2.df[k2.df$Path_1 == 'meningococcus',]
colnames(clus.1.4)
clus.1.4[c('Age..months.', 'WBC', 'array.contemporary.CRP', 'Path_1', 'Path_2', 'clus.1.2')]
clus.1.4 <- clus.1.4[order(clus.1.4$clus.1.2),]
# View(clus.1.4)
write.csv(clus.1.4, file = "mening_clus.1.4.csv", row.names=TRUE)



mening <- k2.df[k2.df$Path_1 == 'meningococcus',][c('Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis', 'system', 'Path_1', 'Path_2' , 'clus.3.2', 'clus.3.4')]
dim(mening)
mening <- mening[order(mening$clus.3.4),]
rownames(mening) <- seq(1, nrow(mening))
colnames(mening) <- c('Age', 'Sex', 'WBC', 'CRP', 'Presentation', 'System', 'Path_1', 'Path_2', 'K2_Cluster', 'K4_Cluster')

plotly.table <- mening
View(plotly.table)
# colnames(plotly.table)

p <- plot_ly(
  type = 'table',
  columnorder = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  columnwidth = c(25, 20, 20, 20, 20, 80, 30, 35, 30, 28, 28),
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
# api_create(p, filename = "table_mening")


















############ bootstraping ############
# write.csv(k2.df, file = "k2.df_bootstrapping.csv", row.names=TRUE) # original bootstrap sample file

bs <- read.table('k2.df_bootstrapping.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = TRUE)
# View(bs)
# constructed conversion between clusters manually

bs$clus.3.4_con <- ifelse(bs$clus.3.4 == 1, 4, ifelse(bs$clus.3.4 == 2, 3, ifelse(bs$clus.3.4 == 3, 2, ifelse(bs$clus.3.4 == 4, 1, 0))))
bs$clus.4.4_con <- ifelse(bs$clus.4.4 == 1, 3, ifelse(bs$clus.4.4 == 2, 4, ifelse(bs$clus.4.4 == 3, 1, ifelse(bs$clus.4.4 == 4, 2, 0))))
bs$clus.5.4_con <- ifelse(bs$clus.5.4 == 1, 3, ifelse(bs$clus.5.4 == 2, 2, ifelse(bs$clus.5.4 == 3, 4, ifelse(bs$clus.5.4 == 4, 1, 0))))
bs$clus.6.4_con <- ifelse(bs$clus.6.4 == 1, 2, ifelse(bs$clus.6.4 == 2, 3, ifelse(bs$clus.6.4 == 3, 4, ifelse(bs$clus.6.4 == 4, 1, 0))))

tot <- 239
sum(bs$clus.1.4 == bs$clus.2.4)
tot - sum(bs$clus.1.4 == bs$clus.2.4)

sum(bs$clus.1.4 == bs$clus.3.4_con)
tot - sum(bs$clus.1.4 == bs$clus.3.4_con)

sum(bs$clus.1.4 == bs$clus.4.4_con)
tot - sum(bs$clus.1.4 == bs$clus.4.4_con)

sum(bs$clus.1.4 == bs$clus.5.4_con)
tot - sum(bs$clus.1.4 == bs$clus.5.4_con)

sum(bs$clus.1.4 == bs$clus.6.4_con)
tot - sum(bs$clus.1.4 == bs$clus.6.4_con)

b1 <- c(198, 41)
b2 <- c(176, 63)
b3 <- c(178, 61)
b4 <- c(161, 78)
b5 <- c(155, 84)
df.1 <- data.frame(b1, b2, b3, b4, b5)
df.1
df.2 <- gather(df.1, boot, count)
df.2

df.3 <- mutate(df.2, cluster.assigned=rep(c('same', 'different'), times = 5))

p<-ggplot(df.3, aes(x = boot, y = count, fill = cluster.assigned)) +
  geom_bar(position = "fill",stat = "identity")+
  scale_fill_manual(values=c('#0CE60C','#E10BA0'), name='Assignment')+
  labs(title = "Barplot of Bootstrap Samples Assigned to Same Cluster", x = "Bootstrap Sample", y = "Proportion")
p<-ggplotly(p)
p
api_create(p, filename = "barplot_bootstrap")
















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





# Construction: get the full demographics into dataframe
# clin2 <- read.table('Mega_sub2_Demographic.csv', sep = ',', stringsAsFactors = FALSE, fill = TRUE, header = FALSE)
# dim(clin2)
# a <- c(1, 19, 27, 28, 33, 34, 84, 85) # barcode, Dx,  wbc, abs_neut, perc_neut, perc_lymp, dx_1, dx_2
# clin2[clin2$V1 == '9423641116_B', a] # search clin2 using barcode
# View(clin)


# end

