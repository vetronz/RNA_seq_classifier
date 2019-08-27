###### PACKAGES ######
library(limma)
library(cluster)
# library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
# library(gridExtra)
library(e1071)
library(tidyr)
# library("illuminaHumanv4.db")
# library(Mfuzz)
library(caret)
library(class)
library(plyr)
library(dplyr)

library(tidyverse)
library(neuralnet)
library(ROCR)
library(randomForest)
library(sva)
library(splitstackshape) # stratified sampling
library(plotly)
library(glmnet)
library(ggrepel)

options(scipen=999)

# install.packages("XYZ")

# library('lumi')
# library("illuminaHumanv4.db")
# library(lumiHumanIDMapping)

ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
ip[which(ip$Package == 'sva'),]
ip[which(ip$Package == 'limma'),]

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

table(status$most_general)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)

scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

# confusiton matrix tips
# when passing vectors to table method, pass the ground truth, then predicted
# rows are ground truth values
# top left true neg, top right false pos
# bottom left false neg, bottom right true pos

f1.score <- function(x){
  tpr <- x[4]/(x[4]+x[2]) # tpr == sensitivity == recall
  ppv <- x[4]/(x[4]+x[3]) # ppv == precision
  2*((tpr*ppv) / (tpr+ppv))
}

tpr <- function(x){
  tpr <- x[4]/(x[2]+x[4])
  return(tpr)
}

tnr <- function(x){
  tnr <- x[1]/(x[1]+x[3])
  return(tnr)
}

euclidean_distance <- function(p,q){
  sqrt(sum((p - q)^2))
}

logistic.func <- function(x){
  1/(1+(exp(-(x))))
}

norm_vec <- function(x) sqrt(sum(x^2))



###### COLORS ######

cohort.cols <- c('#ADB814', '#9A0794')
bv.cols <- c('#B43C22', '#283EA9') # red, blue

bv.cols.drawio <- c('#F8A9A4', '#7E9ABF')

# sex.cols <- c('#fc1676', '#16acfc')
# clus.cols <- c('#FFDD38' , '#56DD5F', '#6763CF', '#FF5338')


###### Merge ######
# discovery prep
idx <- status$most_general == 'bacterial' |
  status$most_general == 'viral' |
  status$most_general == 'greyb' |
  status$most_general == 'greyv'|
  status$most_general == 'greyu' |
  status$most_general == 'HC'
dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral', 'healthy_control') 

### remove outlier
idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')] <- FALSE
sum(idx) # 301

X.d <- e.set[,idx]
status.idx <- status[idx,]

# rename most_general
status.idx$most_general <- as.character(status.idx$most_general)
status.idx$most_general[status.idx$most_general == 'greyb'] <- 'probable_bacterial'
status.idx$most_general[status.idx$most_general == 'greyu'] <- 'unknown'
status.idx$most_general[status.idx$most_general == 'greyv'] <- 'probable_viral'
status.idx$most_general[status.idx$most_general == 'HC'] <- 'healthy_control'

status.idx$most_general <- as.factor(status.idx$most_general)
levels(status.idx$most_general)
status.idx$most_general <- factor(status.idx$most_general, levels = dx)
status.idx$most_general
# status.idx$array.contemporary.CRP <- as.numeric(as.character(status.idx$array.contemporary.CRP))

# discovery data
dim(status.idx)
dim(X.d)

## IRIS VALIDATION
# discrepancy in dimension of transcript and label matrix
dim(e.set.i)
dim(status.iris)

X.i.t <- t(e.set.i)
# dim(X.i)

# extract the overlap
common.my_cat <- intersect(rownames(X.i.t), status.iris$My_code)
length(common.my_cat)

# find position of common.my_cat in status.iris
com.idx <- match(common.my_cat, status.iris$My_code)

# pass com.idx to filter the status matrix
status.i <- status.iris[com.idx,]

dim(status.i)
dim(X.i.t)

# create disease index
status.i$most_general
i.idx<- status.i$most_general == 'bacterial' |
  status.i$most_general == 'viral' |
  status.i$most_general == 'greyb' |
  status.i$most_general == 'greyv'|
  status.i$most_general == 'greyu'|
  status.i$most_general == 'HC'
sum(i.idx)

# single kawasaki case removed
which(status.i$most_general == 'KD')


# filter status.i and eset.i (stored as X.i.t) by the i.idx to remove KD
status.i.idx <- status.i[i.idx,]
X.i <- t(X.i.t[i.idx,])
# X.i <- t(X.i) # revert the transpose

# rename most_general
status.i.idx$most_general <- as.character(status.i.idx$most_general)
status.i.idx$most_general[status.i.idx$most_general == 'greyb'] <- 'probable_bacterial'
status.i.idx$most_general[status.i.idx$most_general == 'greyu'] <- 'unknown'
status.i.idx$most_general[status.i.idx$most_general == 'greyv'] <- 'probable_viral'
status.i.idx$most_general[status.i.idx$most_general == 'HC'] <- 'healthy_control'

status.i.idx$most_general <- as.factor(status.i.idx$most_general)
levels(status.i.idx$most_general)
status.i.idx$most_general <- factor(status.i.idx$most_general, levels = dx)
levels(status.i.idx$most_general)


## PREPED DATA
# discovery
dim(status.idx)
dim(X.d)

# validation
dim(status.i.idx)
dim(X.i)

X.d[1:5,1:5]
X.i[1:5,1:5]

# select the intersection, transcripts that are contained within both datasets
int <- intersect(rownames(X.d), rownames(X.i))
length(int)

# pass the int vector of common genes to match to pull rows
dim(X.d[match(int, rownames(X.d)),])
dim(X.i[match(int, rownames(X.i)),])

# check that all rows (transcripts) match between the 2 matrices
sum(rownames(X.i[match(int, rownames(X.i)),]) != rownames(X.d[match(int, rownames(X.d)),]))

X.c <- cbind(X.d[match(int, rownames(X.d)),], X.i[match(int, rownames(X.i)),])
X.c.t <- (t(X.c))

dim(X.c.t)
class(X.c.t)

dim(status.idx)[1]
dim(status.i.idx)[1]

# construct batch and response vectors
cohort <- as.factor(ifelse(c(rep(1, dim(status.idx)[1]),
                            rep(2, dim(status.i.idx)[1])) == 1, 'dis', 'val'))
bacterial <- as.factor(c(ifelse(status.idx$most_general == 'bacterial', 'pos', 'neg'), ifelse(status.i.idx$most_general == 'bacterial', 'pos', 'neg')))

### PCA
# full.pca <- prcomp(X.c.t, scale=TRUE)

pair1 <- as.data.frame(full.pca$x[,1:2])
pcs <- as.data.frame(full.pca$x[,1:4])
pcs$cohort <- cohort
pcs$bacterial <- bct.vec
pcs[1:5,]

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

# PCA Non Combat
ggplot(data = pcs, aes(PC1, PC2, color=cohort))+geom_point(alpha=0.8)+
scale_color_manual(values=cohort.cols)+
  labs(x =paste0('PC1: ', round(pve[1]), ' % Variance '),
       y = paste0('PC2: ', round(pve[2]), ' % Variance '))+
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))+
  guides(color = guide_legend(override.aes = list(size=4)))


                           
ggplot(data = pcs, aes(PC1, PC2, color=bacterial))+geom_point(alpha=0.8)+
  scale_color_manual(values=bv.cols)+
  labs(x =paste0('PC1: ', round(pve[1]), ' % Variance '),
       y = paste0('PC2: ', round(pve[2]), ' % Variance '))+
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))+
  guides(color = guide_legend(override.aes = list(size=4)))


###### ComBat ######
# mod <- model.matrix(~label, data = X.c.t)
# mod0 <- model.matrix(~1, data=X.c.t)
modcombat <- model.matrix(~1, data=as.data.frame(X.c.t))
# modcombat
class(modcombat)

dim(X.c)
length(cohort)
# needed to transpose the matrix to work

X.comb <- ComBat(X.c, batch=cohort, mod=NULL)
# ComBat(X.c.t, batch=cohort, mod=mod0, par.prior=TRUE, prior.plots=FALSE)

dim(X.comb)

# subtle adjustment between the original and the combat matrix
X.comb[1:5,1:5]
X.c[1:5,1:5]

# transpose for PCA
X.comb <- as.data.frame(X.comb)
X.comb.t <- t(X.comb)
# pca.comb <- prcomp(X.comb.t, scale=TRUE)

pcs.comb <- as.data.frame(pca.comb$x[,1:4])

# add partition vectors to the PCs
pcs.comb$cohort <- cohort
pcs.comb$bacterial <- bct.vec

ve <- pca.comb$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

# holy hell its worked
ggplot(data = pcs.comb, aes(PC1, PC2, color=cohort))+geom_point(alpha=0.8) + 
  scale_color_manual(values=cohort.cols)+
  labs(x =paste0('PC1: ', round(pve[1]), ' % Variance '),
       y = paste0('PC2: ', round(pve[2]), ' % Variance '))+
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))+
  guides(color = guide_legend(override.aes = list(size=4)))

ggplot(data = pcs.comb, aes(PC1, PC2, color=bacterial))+geom_point(alpha=0.8)+
  scale_color_manual(values=bv.cols)+
  labs(x =paste0('PC1: ', round(pve[1]), ' % Variance '),
       y = paste0('PC2: ', round(pve[2]), ' % Variance '))+
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))+
  guides(color = guide_legend(override.aes = list(size=4)))



# split the combat normalized matrix back into discovery and val datasets
dim(status.idx)[1]
dim(status.i.idx)[1]
dim(X.comb.t)
X.dis <- X.c.t[cohort == 'dis',]
X.val <- X.c.t[cohort == 'val',]
dim(X.dis)
dim(X.val)


###### Cyber Sort ######
dim(status)
dim(cyber.s)
status[1:5,1:8]
cyber.s[1:5,1:5]
colnames(cyber.s)[1] <- 'my_category_2'

# join the status and cibersort dataframes together on the common column
# join rather than merge preserves the row order
status.cyber <- join(status, cyber.s, by='my_category_2')
dim(status.cyber) # (23+26) = 49 - 1 for the common merge column = 48 columns

# use idx to subset the cybersort dataframe to select out classes of interest
length(idx)
status.cyber.idx <- status.cyber[idx,]
dim(status.cyber.idx)
dim(status.idx)

# subset the cybersort matrix selecting b.v cases
dim(status.cyber.idx[status.idx$most_general == 'bacterial' | status.idx$most_general == 'viral',])

# select all cybersort cell lines for b.v cases
names(status.cyber)[24:45]
df.1 <- status.cyber.idx[status.idx$most_general == 'bacterial' |
                           status.idx$most_general == 'viral',
                         names(status.cyber)[24:45]]

# add class labels
dim(df.1)
df.1$most_general <- status.cyber.idx[status.idx$most_general == 'bacterial' | status.idx$most_general == 'viral',]$most_general
df.1$most_general <- droplevels(df.1$most_general)

# split the b.v cases into 2 separate dataframes to calc sig diff cell lines
split.df <- split(df.1, df.1$most_general)
class(split.df)
dim(split.df[[1]])
dim(split.df[[2]])

test.stat <- NULL
p.val <- NULL
cell <- NULL
for(i in 1:(ncol(split.df[[1]])-1)){
  print(i)
  cell[i] <- names(split.df[[1]][i])
  x <- split.df[[1]][[i]]
  y <- split.df[[2]][[i]]
  
  a<-t.test(x, y,
            alternative = c("two.sided"), paired = FALSE, var.equal = FALSE,
            conf.level = 0.95)
  test.stat[i] <- a$statistic
  p.val[i] <- a$p.value
}

cell.lines<- as.data.frame(cbind(cell, test.stat, p.val))

# remove undetecteble cell lines
cell.lines <- cell.lines[-c(which(as.character(cell.lines$p.val)=='NaN')),]
cell.lines$p.val <- as.numeric(as.character(cell.lines$p.val))

# select significantly differentially expressed cell lines between b.v
cell.thresh <- 0.05
cell.lines[cell.lines$p.val < cell.thresh ,]
sig.cells <- as.character(cell.lines[cell.lines$p.val < cell.thresh ,]$cell)
sig.cells

# filter df.1 by sig genes
df.1.sig <- df.1[,match(sig.cells, colnames(df.1))]

# reshape df for plotting
a<-gather(df.1.sig)
dim(a)

# create a label colup to split boxplot on
label.col <- rep(c(rep('bacterial', dim(split.df[[1]])[1]), rep('viral', dim(split.df[[2]])[1])), length(sig.cells))
length(label.col)
a$label <- label.col
a$key[a$key == 'T.cells.CD4.memory.activated'] <- 'T.cells.CD4.activated' # rename for space
# cibersort boxplot
ggplot(a, aes(key, value, fill=label.col))+geom_boxplot()+
  scale_fill_manual(values=bv.cols)+
  labs(x='Cell Subset', y='Proportion')+
  guides(fill=guide_legend(title="Diagnosis"),
         color = guide_legend(override.aes = list(size=4))) + # legend title
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 18))+


# remove the renamed CD4 cells
a<-gather(df.1.sig)
a$label <- label.col


### check correlation
# select data from clin cases where we have neutrophil perc and compare to ciber sort pred
df.2 <- clin[!is.na(clin$perc_neut),c('category', 'perc_neut')]
colnames(df.2)[1] <- 'my_category_2'

# find index positions of these cases in the status matrix or status.cyber matrix, they are same
match(df.2$my_category_2, status.cyber$my_category_2)

# pass index to status.cyber to filter out the cases with associated blood count predictions
df.3 <- status.cyber[match(df.2$my_category_2, status$my_category_2), c('my_category_2', 'Neutrophils')]

# merge dataframes so we have actual and predicted counts
df.4 <- merge(df.3, df.2)
dim(df.4)

ggplot(df.4, aes(x=Neutrophils*100, y=perc_neut)) + 
  scale_y_continuous(limits = c(1,100))+
  geom_point()+
  labs(x='Cybersort Predicted Neutrophil Proportion', y='Clinical Data Neutrophil Proportion')+
  geom_smooth(method=lm, level=0.95)+
  theme(axis.title=element_text(size=21),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))


lin.mod <- lm(Neutrophils ~ perc_neut, data = df.4)
attributes(lin.mod)

cor(df.4$Neutrophils, df.4$perc_neut)
cor(df.4$Neutrophils, df.4$perc_neut)^2

# as.character(status.cyber.idx$my_category_2) == rownames(X.dis)


############ LIMMA ############
# mean/variance calculations
x_var <- apply(X.dis, 2, var)
x_mean <- apply(X.dis, 2, mean)
df <- data.frame(log2(x_var), log2(x_mean))
colnames(df) <- c('V1', 'V2')

ggplot(df, aes(V2, V1)) +
  geom_vline(xintercept=log2(5))+
  geom_point(size = 0.4, stroke = 0, shape = 16)+
  labs(x ="log2 Mean Expressioin", y = "log2 Variance")+
  theme(axis.title=element_text(size=21),
      # legend.title=element_text(size=21),
      # legend.text=element_text(size=20),
      axis.text.x = element_text(size = 18),
      axis.text.y = element_text(size = 18))

X.dis.fit <- as.data.frame(X.dis[,x_mean > 5])
dim(X.dis.fit)


# # adding CRP to limma
# status.cyber.idx$array.contemporary.CRP <- as.numeric(as.character(status.cyber.idx$array.contemporary.CRP))
# a <- as.numeric(as.character(status.cyber.idx$array.contemporary.CRP))
# b<-is.na(a)
# b
# 
# # empty crp cases, HC and Viral
# status.cyber.idx$most_general[b]
# 
# # empty viral crp index, reverse it to select vrl with crp values
# is.na(a[status.cyber.idx$most_general=='viral'])==FALSE
# 
# # define median viral crp
# vrl.crp.med <- median(a[status.cyber.idx$most_general=='viral'][is.na(a[status.cyber.idx$most_general=='viral'])==FALSE])
# vrl.crp.med
# 
# # set empty viral crp values to vrl.crp.med
# is.na(a[status.cyber.idx$most_general=='viral'])
# status.cyber.idx$array.contemporary.CRP[status.cyber.idx$most_general=='viral'][is.na(a[status.cyber.idx$most_general=='viral'])] <- vrl.crp.med
# 
# # define HC NA index, pass status.cyber.idx and set na values to 1
# is.na(status.cyber.idx$array.contemporary.CRP[status.cyber.idx$most_general == 'HC'])
# status.cyber.idx$array.contemporary.CRP[is.na(status.cyber.idx$array.contemporary.CRP)] <- 1
# 
# # check no NA values
# sum(is.na(status.cyber.idx$array.contemporary.CRP))


# initialize d.dis.lim to add the cybersort covariates
X.dis.lim <- X.dis.fit
sig.cells

X.dis.lim$label <- status.idx$most_general
X.dis.lim$sex <- status.idx$Sex
X.dis.lim$age <- status.idx$Age..months.
X.dis.lim$B.cells.naive <- status.cyber.idx$B.cells.naive
X.dis.lim$B.cells.memory <- status.cyber.idx$B.cells.memory
X.dis.lim$Plasma.cells <- status.cyber.idx$Plasma.cells
X.dis.lim$T.cells.CD8 <- status.cyber.idx$T.cells.CD8
X.dis.lim$T.cells.CD4.naive <- status.cyber.idx$T.cells.CD4.naive
X.dis.lim$T.cells.CD4.memory.activated <- status.cyber.idx$T.cells.CD4.memory.activated
X.dis.lim$NK.cells.resting <- status.cyber.idx$NK.cells.resting
X.dis.lim$Monocytes <- status.cyber.idx$Monocytes
X.dis.lim$Macrophages.M0 <- status.cyber.idx$Macrophages.M0
# X.dis.lim$Neutrophils <- status.cyber.idx$Neutrophils

# X.dis.lim$crp <- status.cyber.idx$array.contemporary.CRP


# check 8 columns added
dim(X.dis.lim)

### DESIGN MATRIX
design <- model.matrix(~label + sex + age + B.cells.naive + B.cells.memory +
                         Plasma.cells + T.cells.CD8 + T.cells.CD4.naive +
                         T.cells.CD4.memory.activated + NK.cells.resting +
                         Monocytes + Macrophages.M0 + 0,
                         # Monocytes + Macrophages.M0 + Neutrophils + 0,
                         # Monocytes + Macrophages.M0 + Neutrophils + crp + 0, # crp line
                       data = X.dis.lim)
colnames(design)<- c("bct","greyb",'greyu', "greyv", 'vrl', 'HC', 'sexM', 'age',
                     'b.naive', 'b.mem', 'plasma', 'CD8', 'CD4.naive', 'CD4.mem',
                     'nk.rest', 'mono', 'macrophage')
                     # 'nk.rest', 'mono', 'macrophage', 'neut')
                     # 'nk.rest', 'mono', 'macrophage', 'neut', 'crp') # crp line

# check preserved order between cyber and X.dis.lim
design[1:5,]
round(status.cyber.idx[1:5,c('Neutrophils', 'Monocytes')],3)

# check sums of design
dim(design)
colSums(design)

# contrast.matrix <- makeContrasts("bct-vrl", levels=design)
contrast.matrix<- makeContrasts("vrl-bct", 'vrl-greyb', levels=design)
# contrast.matrix<- makeContrasts("((bct+vrl+greyb+greyv+greyu)/5)-HC", levels=design)
contrast.matrix
# colnames(fit$coefficients)


# fit <- lmFit(X.dis.lim, design)
fit <- lmFit(t(X.dis.fit), design)

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
                   c(0.5, 0.05), # 6 gap stat of three
                   c(1, 0.05), # 7 gap stat of three
                   c(1.25, 0.001), # 8 gap stat of three
                   c(2, 0.0001)) # 9


boot <- 6
lfc <- bootstraps[[boot]][1]
pval <- bootstraps[[boot]][2]
# pval <- 5.e-5
lfc
pval




all.hits <- topTable(fit2, adjust.method = 'BH', number=nrow(fit2))
all.hits[1:5,]
colnames(all.hits)[1]


results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc)
dim(results)
head(results)
summary(results)
colnames(results) <- c('DV-DB', 'DV-PB')
vennDiagram(results, include = 'both')

# peek at the top 10 in each group
# vrl.bct, brl.greyb appear to be to lfc for each gene with respect to contrast
top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc)
top.hits[1:5,]

# pull all genes passing number=nrow(fit2)
all.hits <- topTable(fit2, adjust.method = 'BH', number=nrow(fit2))
dim(all.hits)
all.hits[1:5,]

# we set max.lfc to the vrl.bct value unless abs value vrl.greyb is bigger
all.hits$max.lfc <- all.hits$vrl.bct
for(i in 1:nrow(all.hits)){
  b.v <- all.hits$vrl.bct[i]
  gb.v <- all.hits$vrl.greyb[i]
  if(abs(gb.v) > abs(b.v)){
    all.hits$max.lfc[i] <- gb.v
  }
}

# filter out the sig.transcripts
all.filt <- all.hits[abs(all.hits$max.lfc) > lfc & all.hits$adj.P.Val < pval,]
dim(all.filt)

# volcano
p<-ggplot(all.hits, aes(y=-log10(adj.P.Val), x=max.lfc)) +
  geom_point(size = 2, stroke = 0, shape = 16) +
  scale_x_continuous(limits = c(-1.5,1.5))+
  scale_y_continuous(limits = c(0, 10))+
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  labs(x ="Log Fold Change", y = "log10 P-value")+
  theme(axis.title=element_text(size=30),
        axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 30))
p









# subset the disc and validation matrices by the sig genes
dim(results)

# get sig.transcripts
attributes(results)
trans <- dimnames(results)[[1]]
length(trans)

b.v.res <- ifelse(results[,1] == 0, FALSE, TRUE)
pb.v.res <- ifelse(results[,2] == 0, FALSE, TRUE)

# creates vector of TRUE FALSE, TRUE if either b.v or pb.v is significant
length(ifelse(b.v.res == TRUE  | pb.v.res == TRUE, 1, 0) == 1)
sum(ifelse(b.v.res == TRUE  | pb.v.res == TRUE, 1, 0) == 1)

# pass this vector to names to extract the sig genes
sig.trans <- trans[ifelse(b.v.res == TRUE  | pb.v.res == TRUE, 1, 0) == 1]

# select sig transcripts from the dis and val dataset
X.diff <- X.dis.fit[,match(sig.trans, colnames(X.dis.fit))]
X.diff.val <- as.data.frame(X.val[,match(sig.trans, colnames(X.val))])

# filter out the healthy controls
X.diff <- X.diff[status.idx$most_general != 'healthy_control',]
X.diff.val <- X.diff.val[status.i.idx$most_general != 'healthy_control',]

status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
status.i.idx.d <- status.i.idx[status.i.idx$most_general != 'healthy_control',]

dim(X.diff)
dim(status.idx.d)

dim(X.diff.val)
dim(status.i.idx.d)


# scale data
X.s<-data.frame(apply(X.diff, 2, scale01))
X.s.val <-data.frame(apply(X.diff.val, 2, scale01))

X.s$bct <- status.idx.d$most_general=='bacterial'
X.s.val$bct <- status.i.idx.d$most_general=='bacterial'

dim(X.s)
dim(X.s.val)

# filter X.s.val to select only the definite cases for testing
bv.filt <- status.i.idx.d$most_general=='viral' | status.i.idx.d$most_general=='bacterial'
X.s.val.bv <- X.s.val[bv.filt,]
dim(X.s.val.bv)


## STRATIFIED SAMPLING
# set.seed(46)
# test.set.df <- stratified(X.s, c('bct'), (1-prop1), select = NULL, replace = FALSE,
#                           keep.rownames = TRUE, bothSets = FALSE)
# 
# print(paste0('test set proportion: ', round(dim(test.set.df)[1]/239, 2)))
# print(paste0('test set bacterial proportion: ', round(sum(test.set.df$bct)/dim(test.set.df)[1], 2))) # preserves overal bct prop in train
# sum(test.set.df$bct)
# 
# index <- match(setdiff(rownames(X.s), test.set.df$rn), rownames(X.s))
# train <- X.s[index, ]
# test <- X.s[-index, ]
# 
# train.fac <- train
# test.fac <- test
# train.fac$bct <- as.factor(train.fac$bct)
# test.fac$bct <- as.factor(test.fac$bct)


###### GLM FEATURE SELECTION ######
dim(X.s)
X <- as.matrix(X.s[-ncol(X.s)])

y <- X.s$bct # feature selection for b.v split
# y <- status.idx.d$most_general == 'bacterial' | status.idx.d$most_general == 'probable_bacterial' # feature selection for b.v pb.v split

# Elastic Net Selection of Alpha
fold_id <- sample(1:10, size = dim(X)[1], replace=TRUE)

# search across a range of alphas
tuning_grid <- tibble::tibble(
  alpha      = seq(0, 1, by = 0.05),
  mse_min    = NA,
  mse_1se    = NA,
  lambda_min = NA,
  lambda_1se = NA
)

for(i in seq_along(tuning_grid$alpha)) {
  print(paste0('alpha parameter search iter: ', i))
  # fit CV model for each alpha value
  fit <- cv.glmnet(X, y, alpha = tuning_grid$alpha[i], foldid = fold_id)
  
  # extract MSE and lambda values
  tuning_grid$mse_min[i]    <- fit$cvm[fit$lambda == fit$lambda.min]
  tuning_grid$mse_1se[i]    <- fit$cvm[fit$lambda == fit$lambda.1se]
  tuning_grid$lambda_min[i] <- fit$lambda.min
  tuning_grid$lambda_1se[i] <- fit$lambda.1se
  Sys.sleep(0.25)
}

# select the alpha with the lowest associated mse
opt.alpha <- tuning_grid$alpha[which.min(tuning_grid$mse_min)]
# opt.alpha <- 0.5
print(paste0('optimal alpha for elastic-net: ', opt.alpha))

tuning_grid$mse_1se # gives the upper bound (mse estimate + 1 SE)
tuning_grid$mse_min # gives the mse estimate

tuning_grid %>%
  mutate(se = mse_1se - mse_min) %>%
  ggplot(aes(alpha, mse_min)) +
  geom_line(size = 1) +
  scale_y_continuous(limits = c(0.07, 0.11))+
  # geom_vline(xintercept = opt.alpha)+
  geom_ribbon(aes(ymax = mse_min + se, ymin = mse_min - se), alpha = .25)+
  theme(axis.title=element_text(size=30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30))


# fit a cv model with opt alpha to assess lambda
fit.cv <- cv.glmnet(X, y, alpha = opt.alpha, foldid = fold_id)
plot(fit.cv)

alpha.df <- as.data.frame(cbind(fit.cv$lambda, fit.cv$cvm, fit.cv$cvsd, fit.cv$cvup, fit.cv$cvlo))
colnames(alpha.df) <- c('lambda', 'mse', 'mse.sd', 'cvup', 'cvlo')
alpha.df$log.lambda <- log(alpha.df$lambda)
alpha.df[1:5,]

ggplot(alpha.df, aes(log.lambda, mse))+
  geom_line(size = 1)+
  geom_ribbon(aes(ymax = mse + mse.sd, ymin = mse - mse.sd), alpha = .25)+
theme(axis.title=element_text(size=30),
      axis.text.x = element_text(size = 30),
      axis.text.y = element_text(size = 30))



# opt lambda set using the lambda.min which is calc for us
opt.lambda <- fit.cv$lambda.min
opt.lambda

# log and add to plot to make sure it looks reasonable
log(opt.lambda)
abline(v=log(opt.lambda))

## opt elastic
elastic <- glmnet(X, y, alpha = opt.alpha)

plot(elastic, xvar = "lambda")
abline(v=log(opt.lambda), col="black", lwd = 1.5)

elastic

coefs <- coef(elastic, s = opt.lambda)
length(coefs)
# coefs1 <- coef(fit.cv, s = "lambda.min")
# sum(coefs == coefs1) # 532

# remove the intercept term
coefs <- coefs[,1][-1]
coefs <- coefs[coefs!=0]

# filter the training and validation data by the selected genes
names(coefs)
X.s.e <- X.s[,match(names(coefs), colnames(X.s))]
X.s.e$bct <- X.s$bct

X.s.e.val <- X.s.val[, match(names(coefs), colnames(X.s.val))]
X.s.e.val$bct <- X.s.val$bct

dim(X.s.e)
dim(X.s.e.val)

# training setup either full or elasticnet features
# X.s<-data.frame(apply(X.dis, 2, scale01))
X.s <- X.s.e
dim(X.s)


###### CROSS VALIDATION MODEL EVALUATION ######
print(paste0('bacterial cases: ', sum(X.s$bct)))

prop1 <- 0.75
prop2 <- 0.7
boot <- 1

n_folds <- 5
n.train <- round(nrow(X.s))

set.seed(2)
folds.i <- sample(rep(1:n_folds, length.out = n.train))

logistic.m <- NULL
knn.m <- NULL
randForrest.m <- NULL
nn.m <- NULL
svm.m <- NULL
df.2.list <- NULL
df.3 <- NULL
# 
# for(j in 1:boot){
#   for (k in 1:n_folds) {
#     print(paste0('boot: ', j, ', fold: ', k))
#     
#     test.i <- which(folds.i == k)
#     
#     train.cv <- X.s[-test.i, ]
#     test.cv <- X.s[test.i, ]
#     
#     # factoring cv data
#     train.cv.fac <- train.cv
#     train.cv.fac$bct <- as.factor(train.cv$bct)
#     
#     test.cv.fac <- test.cv
#     test.cv.fac$bct <- as.factor(test.cv$bct)
#     
#     ### LOGISTIC REGRESSION
#     model <- glm(bct~ ., data=train.cv, family=binomial(link='logit'), maxit = 128)
#     # summary(model)
#     # anova(model, test="Chisq")
#     pred.test <- predict(model, test.cv[-ncol(test.cv)])
#     pr.test <- prediction(pred.test, test.cv$bct)
# 
#     logistic.m[k] <- pr.test %>%
#       performance(measure = "auc") %>%
#       .@y.values
# 
# 
#     ### KNN
#     trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
#     model <- train(bct ~., data = train.cv.fac, method = "knn",
#                    trControl=trctrl,
#                    tuneLength = 10)
#     knn.opt <- model$results$k[which.max(model$results$Accuracy)]
#     p <- knn(train.cv[-ncol(train.cv)], test.cv[-ncol(test.cv)], train.cv$bct,  k=knn.opt, prob=TRUE)
#     p<-attr(p, "prob")
#     p<-1-p
#     pr <- prediction(p, test.cv$bct)
#     knn.m[k] <- pr %>%
#       performance(measure = "auc") %>%
#       .@y.values
# 
#     ### RANDOM FORREST
#     model <- randomForest(bct ~ . , data = train.cv.fac,
#                           nodesize = 1)
#     pred<-predict(model , test.cv.fac[-ncol(test.cv.fac)])
#     model.prob <- predict(model, test.cv.fac, type="prob")
#     p <- model.prob[,2]
#     pr <- prediction(p, test.cv$bct)
#     randForrest.m[k] <- pr %>%
#       performance(measure = "auc") %>%
#       .@y.values
#     
#     ### Neural Net
#     nn <- neuralnet(bct ~ . , train.cv, linear.output = FALSE, act.fct = "logistic",
#                        hidden = 1, rep = 5, stepmax = 1e+07, startweights = NULL, err.fct = "sse")
#     
#     p <- predict(nn, test.cv[-ncol(test.cv)])
#     pr <- prediction(p, test.cv$bct)
#     nn.m[k] <- pr %>%
#       performance(measure = "auc") %>%
#       .@y.values
#     
#     ### SVM
#     model <- svm(bct ~ . , train.cv.fac, kernel = "linear", probability = TRUE)
#     pred <- predict(model, test.cv.fac, probability = TRUE)
#     p <- attr(pred, "prob")[,2]
#     pr <- prediction(p, test.cv$bct)
#     svm.m[k] <- pr %>%
#       performance(measure = "auc") %>%
#       .@y.values
#     Sys.sleep(1)
#   }
# 
# df.1 <- as.data.frame(cbind(logistic.m, knn.m, randForrest.m, nn.m, svm.m))
# # df.1 <- as.data.frame(cbind(randForrest.m, nn.m, svm.m))
# df.2 <- gather(df.1, learner, roc)
# df.2$roc <- unlist(df.2$roc)
# df.2$learner <- factor(df.2$learner)
# df.2.list[[j]] <- df.2
# }
# 
# print(paste0('boots: ', boot))
# 
# df.3 <- rbind(df.2.list[[1]])
# 
# dim(df.3)
# 5*n_folds*boot
# 
# ggplot(df.3, aes(roc, fill= learner, color = learner)) + geom_density( alpha=0.1)+
#   labs(title=paste0('Roc Area Density with pval: ', pval, ' and lfc: ', lfc),
#      x ="density", y = "roc area")
# ggplot(df.3, aes(learner, roc, color = learner)) + geom_boxplot()
# 
# # detach("package:plyr", unload=TRUE)
# df.3 %>%
#   group_by(learner) %>%
#   summarise(roc.m = mean(roc), roc.med = median(roc), roc.sd = sd(roc))



# df.3 <- rbind(df.2.list[[1]], df.2.list[[2]], df.2.list[[3]],
# df.2.list[[4]])

# df.3 <- rbind(df.2.list[[1]], df.2.list[[2]], df.2.list[[3]],
# df.2.list[[4]], df.2.list[[5]], df.2.list[[6]],
# df.2.list[[7]], df.2.list[[8]])

# df.3 <- rbind(df.2.list[[1]], df.2.list[[2]], df.2.list[[3]],
#               df.2.list[[4]], df.2.list[[5]], df.2.list[[6]],
#               df.2.list[[7]], df.2.list[[8]], df.2.list[[9]],
#               df.2.list[[10]], df.2.list[[11]], df.2.list[[12]],
#               df.2.list[[13]], df.2.list[[14]], df.2.list[[15]],
#               df.2.list[[16]])


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


# coef(test_lasso, s = "lambda.1se") %>%
#   tidy() %>%
#   filter(row != "(Intercept)") %>%
#   ggplot(aes(value, reorder(row, value), color = value > 0)) +
#   geom_point(show.legend = FALSE) +
#   ggtitle("Influential variables") +
#   xlab("Coefficient") +
#   ylab(NULL)


###### NEURAL ######
### NEURAL OPTIMIZATION
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)
X.s$bct <- status.idx.d$most_general == 'bacterial'
print(paste0('bacterial cases: ', sum(X.s$bct==TRUE)))

### neural net grid search
boot <- 16
h.n <- 30
hyper_grid <- NULL
hyper_grid <- expand.grid(
  h.n = seq(1:h.n),
  activation = c('logistic', 'tanh'),
  error = c('sse', 'ce')
)

# had to manually select each chunk of the grid to search separately
# as a full grid search iterating over entire grid in one loop took too long

hyper_grid
# hyper_grid <- hyper_grid[1:(h.n * 3),]
# hyper_grid <- hyper_grid[(h.n+1):(h.n*2),]
hyper_grid <- hyper_grid[((h.n*2)+1):(h.n*3),]

dim(hyper_grid)

prop1 <- 0.75
roc.a <- NULL
roc.t <- NULL
j.train <- NULL
j.test <- NULL
h.n.hx <- NULL
roc.train <- NULL
roc.train.me <- NULL
roc.test <- NULL
roc.test.me <- NULL

for(i in 1:nrow(hyper_grid)) {
  # for (k in 1:n_folds) {
  for (k in 1:boot) {
    print(paste0('hyper: ', i, ', fold/boot: ', k))
    
    # boot
    index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
    train.cv <- X.s[index,]
    test.cv <- X.s[-index,]
    
    # CV
    # test.i <- which(folds.i == k)
    # train.cv <- X.s[-test.i, ]
    # test.cv <- X.s[test.i, ]
    
    # dim(train.cv)
    
    # train model
    nn1 <- neuralnet(bct~ ., train.cv, linear.output = FALSE,
                     act.fct = hyper_grid$activation[i],
                     err.fct = hyper_grid$error[i],
                     hidden = hyper_grid$h.n[i],
                     rep = 3, stepmax = 1e+06, startweights = NULL)
    
    pred.train <- predict(nn1, train.cv[-ncol(train.cv)])
    pred.test <- predict(nn1, test.cv[-ncol(test.cv)])
  
    # extract error
    j.train[k] <- prediction(pred.train[,1], train.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
    j.test[k] <- prediction(pred.test[,1], test.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    Sys.sleep(0.5)
  }
  
  full.list <- c(j.train, j.test)
  full.df <- data.frame(matrix(unlist(full.list), nrow=length(full.list), byrow=T))
  colnames(full.df) <- 'roc.A'
  
  full.df$class <- as.factor(sort(rep(seq(1:(length(full.list)/k)), k)))
  full.df$class <- ifelse(full.df$class == 1, 'j.train', 'j.test')
  
  roc.stats <- full.df %>%
    group_by(class) %>%
    summarise(roc.m = mean(roc.A), roc.med = median(roc.A), roc.sd = sd(roc.A))
  
  roc.stats <- roc.stats %>% mutate(
    roc.se = roc.sd/sqrt(k),
    z.stat = qnorm(0.975),
    roc.me = z.stat * roc.se
  )
  h.n.hx[i] <- i
  roc.test[i] <- roc.stats$roc.med[1]
  roc.test.me[i] <- roc.stats$roc.me[1]
  roc.train[i] <- roc.stats$roc.med[2]
  roc.train.me[i] <- roc.stats$roc.me[2]
  Sys.sleep(2)
}

hyper_grid$roc.train <- roc.train
hyper_grid$roc.train.me <- roc.train.me

hyper_grid$roc.test <- roc.test
hyper_grid$roc.test.me <- roc.test.me


hyper_grid$funcs <- ifelse(hyper_grid$activation == 'tanh', 'tanh_sse',
                           ifelse(hyper_grid$error == 'ce',
                                  'logistic_ce', 'logistic_sse'))

hyper_grid
### saving the line search outputs to rds objects
# logistic_sse_lineSearch
# tanh_sse_lineSearch
# logistic_ce_lineSearch

# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/')
# saveRDS(logistic_ce_lineSearch, "logistic_ce_lineSearch.rds")

### load previously saved values
# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/')
logistic_sse_lineSearch <- readRDS("logistic_sse_lineSearch.rds")
tanh_sse_lineSearch <- readRDS("tanh_sse_lineSearch.rds")
logistic_ce_lineSearch <- readRDS("logistic_ce_lineSearch.rds")

# reconstruct the grid from the line searches
grid_search<-as.data.frame(rbind(logistic_sse_lineSearch, tanh_sse_lineSearch, logistic_ce_lineSearch))
rownames(grid_search) <- seq(1:dim(grid_search)[1])

grid_search
grid_search[order(grid_search$AUC.test, decreasing = TRUE),][1:10,]

opt.h.n <- grid_search[order(grid_search$roc.test, decreasing = TRUE),][1,][[1]]
opt.error <- as.character(grid_search[order(grid_search$roc.test, decreasing = TRUE),][1,][[3]])

# rename the columns for plotting
colnames(logistic_sse_lineSearch)[4:7] <- c('AUC.train', 'AUC.train.me', 'AUC.test', 'AUC.test.me')
colnames(tanh_sse_lineSearch)[4:7] <- c('AUC.train', 'AUC.train.me', 'AUC.test', 'AUC.test.me')
colnames(logistic_ce_lineSearch)[4:7] <- c('AUC.train', 'AUC.train.me', 'AUC.test', 'AUC.test.me')

# complexity plots
ggplot(logistic_sse_lineSearch, aes(x=h.n, y=AUC.train, color='AUC.train')) +
  geom_line(aes(y=AUC.train))+
  scale_y_continuous(limits = c(.65,1))+
  geom_errorbar(aes(ymin=AUC.train-AUC.train.me, ymax=AUC.train+AUC.train.me), width=.4, position=pd)+
  geom_line(aes(y=AUC.test, color='AUC.test'))+
  geom_errorbar(aes(ymin=AUC.test-AUC.test.me, ymax=AUC.test+AUC.test.me, color='AUC.test'), width=0.4)+
  scale_colour_manual(values = train.test.cols)+
  labs(x =paste0('1 - ', h.n, ' hidden nodes'), y = "AUC", color='Cohort')+
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))
 

ggplot(tanh_sse_lineSearch, aes(x=h.n, y=AUC.train, color='AUC.train')) +
  geom_line(aes(y=AUC.train))+
  scale_y_continuous(limits = c(.65,1))+
  geom_errorbar(aes(ymin=AUC.train-AUC.train.me, ymax=AUC.train+AUC.train.me), width=.4, position=pd)+
  geom_line(aes(y=AUC.test, color='AUC.test'))+
  geom_errorbar(aes(ymin=AUC.test-AUC.test.me, ymax=AUC.test+AUC.test.me, color='AUC.test'), width=0.4)+
  scale_colour_manual(values = train.test.cols)+
  labs(x =paste0('1 - ', h.n, ' hidden nodes'), y = "AUC", color='Cohort')+
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))


ggplot(logistic_ce_lineSearch, aes(x=h.n, y=AUC.train, color='AUC.train')) +
  geom_line(aes(y=AUC.train))+
  scale_y_continuous(limits = c(.65,1))+
  geom_errorbar(aes(ymin=AUC.train-AUC.train.me, ymax=AUC.train+AUC.train.me), width=.4, position=pd)+
  geom_line(aes(y=AUC.test, color='AUC.test'))+
  geom_errorbar(aes(ymin=AUC.test-AUC.test.me, ymax=AUC.test+AUC.test.me, color='AUC.test'), width=0.4)+
  scale_colour_manual(values = train.test.cols)+
  labs(x =paste0('1 - ', h.n, ' hidden nodes'), y = "AUC", color='Cohort')+
  guides(color = guide_legend(override.aes = list(size=4))) +
  theme(axis.title=element_text(size=21),
        legend.title=element_text(size=21),
        legend.text=element_text(size=20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20))




# neural opt, val
f1.opt <- NULL
roc.opt <- NULL
for (i in 1:boot){
  print(i)
  
  nn.opt <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                      hidden = opt.h.n, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = opt.error)
  pred.opt.val <- predict(nn.opt, X.s.e.val[bv.filt,][-ncol(X.s.e.val)])
  pred.opt.val <- ifelse(pred.opt.val > 0.5, TRUE, FALSE)
  f1.opt[i] <- f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
  
  prob.opt.val <- predict(nn.opt, X.s.e.val[bv.filt,][-ncol(X.s.e.val)])
  pr <- prediction(prob.opt.val, status.i.idx.d[bv.filt,]$most_general=='bacterial')
  roc.opt[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  Sys.sleep(0.3)
}
table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)
f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
mean(unlist(roc.opt))
median(unlist(roc.opt))



# tpr tnr threshold
nn.opt <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                    hidden = opt.h.n, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
prob.opt.val <- predict(nn.opt, X.s.val.bv)

p.scale <- 1:99
tpr.h <- NULL
tnr.h <- NULL
for(i in p.scale){
  p.thresh <- i/100
  print(p.thresh)
  pred.opt.val <- ifelse(prob.opt.val > p.thresh, TRUE, FALSE)
  table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)
  tpr.h[i] <- tpr(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
  tnr.h[i] <- tnr(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
}
df.1 <- as.data.frame(cbind(tpr.h, tnr.h))
df.2 <- gather(df.1, 'metric', 'result')
df.2$P_threshold <- rep(seq(1:length(p.scale))/100,2)

ggplot(df.2, aes(x=P_threshold, y=result, group=metric, color=metric))+geom_line()

opt.thresh <- max(which(df.1$tpr.h > df.1$tnr.h))/100

pred.opt.val <- ifelse(prob.opt.val > opt.thresh, TRUE, FALSE)
table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)
f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))


###### PSEUDO LABELING ######
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)
X.s$bct <- status.idx.d$most_general == 'bacterial'
# X.s$bct <- as.factor(X.s$bct)
print(paste0('bacterial cases: ', sum(X.s$bct==TRUE)))

# Pseudo Labels
X.psd <- X.s
X.psd$most_general <- status.idx.d$most_general
dim(X.psd)


b.thresh <- 0.999
nn.psd.df <- NULL
ppb.h <- NULL
ppb.prob.h <- NULL
roc.h.psd <- NULL
b.thresh.h <- NULL
for(i in 1:250){
  if (sum(X.psd$most_general == 'probable_bacterial' & X.psd$bct==FALSE) == 0) {
    print(paste0('Breaking at iteration: ' , i, ', PB cases: ', sum(X.psd$most_general == 'probable_bacterial' & X.psd$bct==FALSE)))
    break
  }
  index <- sample(nrow(X.psd), round(prop1*nrow(X.psd)))
  train.cv <- X.psd[index, ]
  test.cv <- X.psd[-index, ]
  
  # dim(train.cv)
  # dim(test.cv)
  
  pb.psd <- which(test.cv$most_general == 'probable_bacterial' & test.cv$bct == TRUE)
  # dim(test.cv[pb.psd,])
  
  if(identical(pb.psd, integer(0))){
    print('pass')
  } else{
    print('some pseudo labels')
    
    # add the pseudo-labaled PB cases to train set
    train.cv <- rbind(train.cv, test.cv[pb.psd,])
    
    # select all rows except the matched pb.psd case and remove from test.cv
    test.cv <- test.cv[-c(match(rownames(test.cv[pb.psd,]), rownames(test.cv))),]
  }
  
  model <- neuralnet(bct ~ . , train.cv[-ncol(train.cv)], linear.output = FALSE, act.fct = "logistic",
                     hidden = opt.h.n, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = opt.error)
  
  pred <- predict(model, test.cv[-((ncol(test.cv)-1):ncol(test.cv))])
  
  # create a filter to extract the pb cases from both
  pb.filt <- test.cv$most_general == 'probable_bacterial'
  
  # which.max(pred[pb.filt,])
  ppb <- names(which.max(pred[pb.filt,]))
  ppb.prob <- max(pred[pb.filt])
  
  if(ppb.prob > b.thresh){
    # check ROC performance BEFORE changing PB label
    # this is despite the fact we trained the network outside of if statement
    # this prevents the PB pseudolabels in to train set arrtificially elevating ROC
    # by improving detection of PB case.
    pr <- prediction(pred, test.cv$bct)
    roc.h.psd[i] <- pr %>%
      performance(measure = "auc") %>%
      .@y.values
    
    # pseudo label the ppb case
    X.psd$most_general[match(ppb, rownames(X.psd))]
    X.psd$bct[match(ppb, rownames(X.psd))]
    X.psd$bct[match(ppb, rownames(X.psd))] <- TRUE
    
    # pseudo label the ppb case
    # status.idx.d$most_general[which(status.idx.d$my_category_2 == ppb)] <- 'bacterial'
    # X.s$bct <- status.idx.d$most_general == 'bacterial'
    # X.s$bct <- as.factor(X.s$bct)
    print(paste0('PB Case: ', ppb, ' added with P: ', round(ppb.prob, 5)))
    print(paste0('iteration: ', i, ', bacterial cases: ', sum(X.psd$bct==TRUE), ', bac.threshold: ', b.thresh))
    ppb.h[i] <- ppb
    ppb.prob.h[i] <- ppb.prob
    b.thresh.h[i] <- b.thresh
  }  else {
    if (b.thresh <= 0) {print(paste0('passing - B.Thresh: ', b.thresh)) }
    else {
      if(i < 90){
        b.thresh = b.thresh-0.001
      } else if(i < 150){
        b.thresh = b.thresh-0.001
      } else {
        b.thresh = b.thresh-0.02
      } 
    }
  }
  Sys.sleep(0.25)
}


print(paste0('bacterial cases: ', sum(X.psd$bct==TRUE),
             ', PB Cases: ', sum(X.psd$most_general == 'probable_bacterial' & X.psd$bct == FALSE),
             ', b.threshold: ', b.thresh))


# # create ppb.h.df
# ppb.h
# b.thresh.h
# roc.h.psd
# ppb.prob.h
# 
# # remove empty values
# ppb.h.f <- ppb.h[!is.na(ppb.h)]
# b.thresh.h.f <- b.thresh.h[!is.na(b.thresh.h)]
# roc.h.psd.f <- unlist(roc.h.psd)
# ppb.prob.h.f <- ppb.prob.h[!is.na(ppb.prob.h)]
# 
# a <- as.data.frame(cbind(ppb.h.f, b.thresh.h.f, roc.h.psd.f, ppb.prob.h.f))
# colnames(a) <- c('pb.case', 'threshold', 'roc.a', 'prob')
# 
# # make numeric and round
# a$threshold <- as.numeric(as.character(a$threshold))
# a$roc.a <- round(as.numeric(as.character(a$roc.a)), 5)
# a$prob <- round(as.numeric(as.character(a$prob)), 5)
# 
# # set the below 0 thresholds to 0
# a$threshold[a$threshold < 0] = 0
# 
# # set to ppb dataframe
# ppb.h.df <- a
# 
# ppb.h.df$index <- 1:nrow(ppb.h.df)


### load previously saved values
# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/')
# saveRDS(ppb.h.df, "ppb.h.df.rds")
ppb.h.df <- readRDS("ppb.h.df.rds")
ppb.h.df

# select ROC derivative as query pseudo labeling cut off
a <- lm(formula = roc.a ~ splines::bs(index, 3), data=ppb.h.df)

ppb.opt <- which.max(a$fitted.values)[[1]]
ppb.opt
ggplot(ppb.h.df, aes(index, roc.a))+
  geom_text(aes(label=pb.case), check_overlap = TRUE, size=2.7)+
  labs(title='Iterative Pseudo-labeling Error', x='PB', y='ROCA')+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE)+
  geom_vline(xintercept = ppb.opt)

# select the abs max jump to set as a query pseudo labeling cut off
psd.opt <- which.max(abs(diff(ppb.h.df$prob)))
ggplot(ppb.h.df, aes(index, prob))+
  geom_vline(xintercept = ppb.opt + 0.5)+
  # geom_text(aes(label=ppb.h), check_overlap = TRUE, size=2)
  geom_text_repel(aes(label=pb.case), size=2, segment.colour=NA)


ppb.opt
psd.opt


###  ADD PSEUDO LABELS FOR TRAINING
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)
X.s$bct <- status.idx.d$most_general == 'bacterial'
print(paste0('bacterial cases: ', sum(X.s$bct==TRUE)))

# ppb.opt <- psd.opt # if we use the P drop cutoff
ppb.h.df[1:ppb.opt,]$pb.case

# find optimal pb labels in status.idx.d
ppb.pos <- match(ppb.h.df[1:ppb.opt,]$pb.case, status.idx.d$my_category_2)

# check they are the same as ppb.h and that they are pb
status.idx.d$my_category_2[ppb.pos]
status.idx.d$most_general[ppb.pos]

# change the labels and add to X.s
status.idx.d$most_general[ppb.pos] <- 'bacterial' # add all pseudo labeled bct cases
# status.idx.d$most_general[ppb.pos[1:22]] <- 'bacterial' # add only top 20 pbs to bct labels
X.s$bct <- status.idx.d$most_general == 'bacterial'
table(status.idx.d$most_general)
print(paste0('bacterial cases: ', sum(X.s$bct==TRUE)))





# neural opt, val, pseud bct-vrl
f1.opt.psd <- NULL
roc.opt.psd <- NULL
for (i in 1:boot){
  print(i)
  nn.opt.psd <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                          hidden = opt.h.n, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = opt.error)
  
  pred.opt.val.psd <- predict(nn.opt.psd, X.s.e.val[bv.filt,][-ncol(X.s.e.val)])
  pred.opt.val <- ifelse(pred.opt.val.psd > 0.5, TRUE, FALSE)
  f1.opt.psd[i] <- f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
  
  prob.opt.val.psd <- predict(nn.opt.psd, X.s.e.val[bv.filt,][-ncol(X.s.e.val)])
  pr <- prediction(prob.opt.val.psd, status.i.idx.d[bv.filt,]$most_general=='bacterial')
  roc.opt.psd[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  Sys.sleep(0.2)
}

table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)
f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))

nn.opt.psd <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                        hidden = opt.h.n.psd, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")

prob.opt.val.psd <- predict(nn.opt.psd, X.s.e.val[-ncol(X.s.e.val)])
prob.opt.val.b.v.psd <- predict(nn.opt.psd, X.s.e.val[bv.filt,][-ncol(X.s.e.val)])

pr.all <- prediction(prob.opt.val.psd, status.i.idx.d$most_general=='bacterial')
pr.b.v <- prediction(prob.opt.val.b.v.psd, status.i.idx.d[bv.filt,]$most_general=='bacterial')

roc.perf.all <- performance(pr.all, measure = "tpr", x.measure = "fpr")
pr.all %>%
  performance(measure = "auc") %>%
  .@y.values

roc.perf.b.v <- performance(pr.b.v, measure = "tpr", x.measure = "fpr")
pr.b.v %>%
  performance(measure = "auc") %>%
  .@y.values

roc.perf.all
x <- attr(roc.perf.all, 'x.values')
x <- x[[1]]
y <- attr(roc.perf.all, 'y.values')
y <- y[[1]]

df.1 <- as.data.frame(cbind(x,y))

x <- attr(roc.perf.b.v, 'x.values')
x <- x[[1]]
y <- attr(roc.perf.b.v, 'y.values')
y <- y[[1]]

df.2 <- as.data.frame(cbind(x,y))

ggplot(df.1, aes(x,y))+
  geom_line()+
  geom_line(data = df.2, aes(x,y), color='red')+
  labs(title='ROC Curve for Pseudo-labeled Network', x='1 - TNR', y='TPR')


# compare F1 scores between normal and psed bootstraps
a <- as.data.frame(cbind(f1.opt, f1.opt.psd))
colnames(a) <- c('f1_normal', 'f1_pseudo_labeled')
b<-gather(a, 'model', 'result')
ggplot(b, aes(result, colour=model, fill=model))+geom_density(alpha=0.2)+
  labs(title = 'F1 Score of Normal and Pseudo-Labeled Models in Iris Validation Cohort', x='F1 Score', y='Density')
ggplot(b, aes(model, result, colour=model, fill=model))+geom_boxplot(alpha=0.2)+
  labs(title = 'F1 Score of Normal and Pseudo-Labeled Models in Iris Validation Cohort', x='Model', y='F1 Score')

# compare ROC.A scores between normal and psed bootstraps
a <- as.data.frame(cbind(unlist(roc.opt), unlist(roc.opt.psd)))
colnames(a) <- c('ROC_normal', 'ROC_pseudo_labeled')
b<-gather(a, 'model', 'result')
b %>% 
  group_by(model)%>%
  summarise(x.mean = mean(result), x.med = median(result))

ggplot(b, aes(result, colour=model, fill=model))+geom_density(alpha=0.2)
ggplot(b, aes(model, result, colour=model, fill=model))+geom_boxplot(alpha=0.2)+
  scale_y_continuous(limits = c(0.8,1))+
  labs(title = 'ROC Area Score of Normal and Pseudo-Labeled Models in Iris Validation Cohort', x='Model', y='ROC Area')


# PB and Unknown distributions
pb.dist.df <- as.data.frame(cbind(predict(nn.opt, X.s.val[status.i.idx.d$most_general=='probable_bacterial',]),
                                  predict(nn.opt.psd, X.s.val[status.i.idx.d$most_general=='probable_bacterial',])))
colnames(pb.dist.df) <- c('normal', 'pseudo_labeled')
pb.dist.df <- melt(pb.dist.df)
colnames(pb.dist.df)[1] <- 'model'
hist1 <- ggplot(pb.dist.df, aes(value))+geom_histogram(bins = 10)
# hist1 + facet_grid(model ~ .)
hist1 + facet_wrap(model ~ .)

unknown.dist.df <- as.data.frame(cbind(predict(nn.opt, X.s.val[status.i.idx.d$most_general=='unknown',]),
                                       predict(nn.opt.psd, X.s.val[status.i.idx.d$most_general=='unknown',])))
colnames(unknown.dist.df) <- c('normal', 'pseudo_labeled')
unknown.dist.df <- melt(unknown.dist.df)
colnames(unknown.dist.df)[1] <- 'model'
hist2 <- ggplot(unknown.dist.df, aes(value))+geom_histogram(bins = 10)
# hist1 + facet_grid(model ~ .)
hist2 + facet_wrap(model ~ .)


# tpr tnr threshold
nn.opt.psd <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                        hidden = opt.h.n.psd, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
prob.opt.val <- predict(nn.opt.psd, X.s.val.bv)

p.scale <- 1:99
tpr.h <- NULL
tnr.h <- NULL
f1.h <- NULL
for(i in p.scale){
  p.thresh <- i/100
  print(p.thresh)
  pred.opt.val <- ifelse(prob.opt.val > p.thresh, TRUE, FALSE)
  table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)
  tpr.h[i] <- tpr(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
  tnr.h[i] <- tnr(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
  f1.h[i] <-  f1.score(table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val))
}
df.0 <- as.data.frame(cbind(tpr.h, tnr.h, (p.scale/100), f1.h))
opt.thresh <- which.max(df.0$f1.h)
df.0[opt.thresh,]

df.1 <- as.data.frame(cbind(tpr.h, tnr.h))
df.2 <- gather(df.1, 'metric', 'result')
df.2$P_threshold <- rep(seq(1:length(p.scale))/100,2)
ggplot(df.2, aes(x=P_threshold, y=result, group=metric, color=metric))+geom_line()+
  labs(title = 'TRP (Sensitivity - Blue) and TNR (Specificity - RED) With Varying Classifier Cutoffs',
       x='Probability Threshold', y='TPR & TNR')+
  geom_vline(xintercept = opt.thresh/100, linetype="dashed", 
             color = "black", size=0.2)

# comparison of confusion matrix with 0.5 and optimal cutoff
pred.opt.val <- ifelse(prob.opt.val > 0.5, TRUE, FALSE)
table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)

pred.opt.val <- ifelse(prob.opt.val > (opt.thresh/100), TRUE, FALSE)
table(status.i.idx.d[bv.filt,]$most_general == 'bacterial', pred.opt.val)


### learning curves
j_train <- NULL
j_test <- NULL
roc.train <- NULL
roc.test <- NULL
roc.train.me <- NULL
roc.test.me <- NULL
learning_curve.df <- NULL
p.h <- NULL
boot<- 16
props <- seq(from=10, to=90, by=5)/100

for (j in 1:length(props)){
  for(i in 1:boot){
    prop <- props[j]
    print(paste0('proportion: ', prop, ', bootstrap: ', i))
    index <- sample(nrow(X.s), round(prop*nrow(X.s)))
    train.cv <- X.s[index, ]
    test.cv <- X.s[-index, ]
    
    model <- neuralnet(bct ~ . , train.cv, linear.output = FALSE, act.fct = "logistic",
                       hidden = opt.h.n, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
    
    pred_train <- predict(model, train.cv[-ncol(train.cv)])
    pred_test <- predict(model, test.cv[-ncol(test.cv)])
    
    j_train[i] <- prediction(pred_train[,1], train.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
    j_test[i] <- prediction(pred_test[,1], test.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values  
  }
  # class.calls <- c(j_test)
  class.calls <- c(j_train, j_test)
  full.list <- class.calls
  full.df <- data.frame(matrix(unlist(full.list), nrow=length(full.list), byrow=T))
  colnames(full.df) <- 'roc.A'
  
  full.df$class <- as.factor(sort(rep(seq(1:(length(class.calls)/boot)), boot)))
  
  roc.stats <- full.df %>%
    group_by(class) %>%
    summarise(roc.m = mean(roc.A), roc.v = var(roc.A), roc.sd = sd(roc.A))
  
  roc.stats <- roc.stats %>% mutate(
    roc.se = roc.sd/sqrt(k),
    z.stat = qnorm(0.975),
    roc.me = z.stat * roc.se
  )
  roc.train[j] <- roc.stats$roc.m[1]
  roc.train.me[j] <- roc.stats$roc.me[1]
  roc.test[j] <- roc.stats$roc.m[2]
  roc.test.me[j] <- roc.stats$roc.me[2]
}

learning_curve.df <- as.data.frame(cbind(roc.train, roc.test, roc.train.me, roc.test.me))
learning_curve.df$prop <- props
colnames(learning_curve.df) <- c('train', 'test', 'train.me', 'test.me', 'prop')

ggplot(learning_curve.df, aes(x=prop, y=train)) +
  # scale_y_continuous(limits = c(0.5,1))+
  geom_point(aes(y=train, color='train'))+
  geom_errorbar(aes(ymin=train-train.me, ymax=train+train.me, color='train'), width=.02, alpha=0.75)+
  geom_point(aes(y=test, color='test'))+
  geom_errorbar(aes(ymin=test-test.me, ymax=test+test.me, color='test'), width=0.02, alpha=0.75)+
  labs(title=paste0('Learning Curve with ', opt.h.n, ' hidden nodes'), x ="training Data Percentage", y = "ROCA")





###### 2 TRANSCRIPT COMPARISON ######
# install.packages("pROC")
library(pROC)
library(gmodels)
setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/DRS')

drs.val = read.table("DRS_results_on_validation-2.6266572_probes.txt", 
               sep="\t",
               fill=FALSE, 
               strip.white=TRUE)

class(drs.val)
colnames(drs.val) <- c('DRS', 'most_general')
drs.val$most_general <- ifelse(drs.val$most_general==1, 'viral',
                               ifelse(drs.val$most_general==2, 'probable_viral',
                                      ifelse(drs.val$most_general==3, 'unknown',
                                             ifelse(drs.val$most_general==4, 'probable_bacterial', 'bacterial'))))

rownames(drs.val)

dim(drs.val)
dim(X.s.e.val)

sum(rownames(drs.val) == rownames(X.s.e.val))

# add bct vector
drs.val$bct <- drs.val$most_general == 'bacterial'


nn.opt.psd <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                        hidden = opt.h.n.psd, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")

prob.opt.val.psd <- predict(nn.opt.psd, X.s.e.val[-ncol(X.s.e.val)])
pr <- prediction(prob.opt.val.psd, status.i.idx.d$most_general=='bacterial')
pr %>%
  performance(measure = "auc") %>%
  .@y.values

nn.psd.val <- as.data.frame(prob.opt.val.psd)
colnames(nn.psd.val) <- 'prob'
nn.psd.val$most_general <- status.i.idx.d$most_general

# p <- ggplot(nn.psd.val, aes(most_general, prob, fill=most_general)) +geom_boxplot()
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, alpha=0.5)

# p <- ggplot(drs.val, aes(most_general, DRS, fill=most_general))+geom_boxplot()
# p + geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6, alpha=0.5)



ggplot(nn.psd.val, aes(most_general, prob, color=most_general)) +
  geom_violin(alpha=0.2)+
  geom_jitter(width = 0.25, height = 0.002)+
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "black" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "black" )

ggplot(drs.val, aes(most_general, DRS, color=most_general)) +
  geom_violin(alpha=0.2)+
  geom_jitter(width = 0.2, height = 0.002)+
  stat_summary(fun.y = "mean", geom = "point", 
               shape = 8, size = 3, color = "black" ) +
  stat_summary(fun.y = "median", geom = "point", 
               shape = 2, size = 3, color = "black" )

a <- cbind(nn.psd.val, drs.val)

dim(a)

a <- cbind(nn.psd.val, drs.val)

a[,c(2,4,5)] <- NULL

range01(a$DRS)
logistic.func(a$DRS)
a$scale.drs <- scale01(a$DRS)

a
b <- a[status.i.idx.d$most_general=='probable_viral',]

b$DRS <- NULL
c<-gather(b, 'model', 'result')
ggplot(c, aes(result, fill=model))+geom_density(alpha=0.2)





# create pb and unknown filter
pb.filt <- status.i.idx.d$most_general=='probable_bacterial'
u.filt <- status.i.idx.d$most_general=='unknown'

# extract pos and neg nn predictions from these
pb.pos <- rownames(nn.psd.val[pb.filt,][nn.psd.val[pb.filt,]$prob > 0.5,])
pb.neg <- rownames(nn.psd.val[pb.filt,][nn.psd.val[pb.filt,]$prob < 0.5,])

u.pos <- rownames(nn.psd.val[u.filt,][nn.psd.val[u.filt,]$prob > 0.5,])
u.neg <- rownames(nn.psd.val[u.filt,][nn.psd.val[u.filt,]$prob < 0.5,])

# index to pass to View
match(u.pos, status.i.idx.d$My_code)

# PB partitions
View(status.i.idx.d[match(pb.pos, status.i.idx.d$My_code),])
View(status.i.idx.d[match(pb.neg, status.i.idx.d$My_code),])

# Unknown Partitions
View(status.i.idx.d[match(u.pos, status.i.idx.d$My_code),])
View(status.i.idx.d[match(u.neg, status.i.idx.d$My_code),])

# lots of unknows in neg group could well have had abx therapy
# perhaps not required?



# cannot see any easy partition of these cases
# all PBs look sick from the limited clinical information in the iris dataset


# roc 2 transcript signature
a <- roc(drs.val, bct, DRS, ci=TRUE, conf.level=0.95, boot.n=100)
attributes(a)
a
plot(a)



# coefs are the same as colnames in validation iris set
names(coefs) == colnames(X.s.e.val[-ncol(X.s.e.val)])

# remove the leading X
coefs.array <- substring(names(coefs),2)


setwd('~/Documents/Masters/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
dim(illumina)


coefs.probe <- illumina$Probe_Id[match(colnames(X.diff), illumina$Array_Address_Id)]
# coefs.probe <- illumina$Probe_Id[match(coefs.array, illumina$Array_Address_Id)]
coefs.probe


library('biomaRt')

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'illumina_humanht_12_v4', 
              values = coefs.probe, 
              mart = ensembl)

df.1

myrsini.sig <- c('ENSG00000182118', 'ENSG00000137959') # ifi and fama looked up on gene cards
match(myrsini.sig, df.1$ensembl_gene_id)

df.1[match(myrsini.sig, df.1$ensembl_gene_id),]



###### CLUSTERING ######
dim(e.set.i)
dim(X.s.e.val)

e.set.i.d <- t(e.set.i)[status.i$most_general != 'HC',]
dim(e.set.i.d)

k.test <- cmeans(e.set.i.d, centers = 2, iter.max=10, verbose=FALSE, dist="euclidean",
       method="cmeans", rate.par = NULL)

wss <- function(d, cluster){
  # d <- dist(d, method='euclidean', diag=TRUE)
  dmat <- as.matrix(d)
  
  cn <- max(cluster)
  clusterf <- as.factor(cluster)
  clusterl <- levels(clusterf)
  cnn <- length(clusterl)
  
  within.cluster.ss <- 0
  for (i in 1:cn) {
    cluster.size <- sum(cluster == i)
    di <- as.dist(dmat[cluster == i, cluster == i])
    within.cluster.ss <- within.cluster.ss + sum(di^2)/cluster.size
  }
  within.cluster.ss
}


k.max <- 8
v <- rep(0, k.max)
# x <- iris[-ncol(iris)]
# x <- e.set.i.d
x <- X.s.e.val
dim(x)
x <- dist(x, method='euclidean', diag=TRUE)


for (i in 2:k.max) {
  print(paste0('iter: ', i))
  clust <- cmeans(x, centers = i, iter.max=300, verbose=FALSE, dist="euclidean",
                  method="cmeans", rate.par = NULL)
  v[i] <- wss(x, clust$cluster)
}
plot(v[-1])
v

clust <- cmeans(x, centers = 2, iter.max=300, verbose=FALSE, dist="euclidean",
                method="cmeans", m = 1.2, rate.par = NULL)

clust$membership
df.1 <- clust$centers

df.1 <- as.data.frame(clust$membership)
df.1$cluster <- clust$cluster

ggplot(df.1, aes(cluster, fill=status.i.idx.d$most_general)) + 
  geom_bar()

# compare the fuzzy cluster results with the nn results
# create PB filter
pb.filt <- status.i.idx.d$most_general=='probable_bacterial'


nn.opt.psd <- neuralnet(bct ~ . , X.s, linear.output = FALSE, act.fct = "logistic",
                        hidden = opt.h.n.psd, rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
prob.opt.val.psd <- predict(nn.opt.psd, X.s.e.val[-ncol(X.s.e.val)])
pr <- prediction(prob.opt.val.psd, status.i.idx.d$most_general=='bacterial')
pr %>%
  performance(measure = "auc") %>%
  .@y.values


dim(prob.opt.val.psd)

# nn pb probabilities
prob.opt.val.psd[pb.filt,]

# fuzzy pb probabilities
clust$membership[pb.filt,][,2]


full.pca <- prcomp(x[pb.filt,])

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

# fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]


dim(pair1)
# separation based on batch effect
ggplot(data = pair1, aes(PC1, PC2, color = clust$membership[pb.filt,][,1]>0.5))+geom_point() + 
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color=cohort))+geom_point()+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[3]), ' %'), y = paste0('variance: ', round(pve[4]), ' %'))
ggplot(data = pair1, aes(PC1, PC2, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))




























### end

