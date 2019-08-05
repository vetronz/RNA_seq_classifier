###### PACKAGES ######
library(limma)
library(cluster)
# library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
# library(gridExtra)
# library(plotly)
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
library(splitstackshape)

# install.packages("XYZ")

ip <- as.data.frame(installed.packages()[,c(1,3:4)])
rownames(ip) <- NULL
ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
ip[which(ip$Package == 'sva'),]
ip[which(ip$Package == 'splitstackshape'),]

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


# 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 5
cols = gg_color_hue(n)
# 
# n = 8
# cols.8 = gg_color_hue(n)
# 
# n = 10
# cols.10 = gg_color_hue(n)
# 
# n = 14
# cols.14 = gg_color_hue(n)
# 
# dx.cols <- c('#D81D22', '#FF3A3A', '#AD4187' , '#776BB9' , '#4A8AC3', '#32A46A')
# sex.cols <- c('#fc1676', '#16acfc')
# clus.cols <- c('#FFDD38' , '#56DD5F', '#6763CF', '#FF5338')


###### COMBI ######
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
batch <- as.factor(ifelse(c(rep(1, dim(status.idx)[1]),
                            rep(2, dim(status.i.idx)[1])) == 1, 'discovery', 'iris'))
bct.vec <- as.factor(c(ifelse(status.idx$most_general == 'bacterial', 'pos', 'neg'), ifelse(status.i.idx$most_general == 'bacterial', 'pos', 'neg')))

### PCA
# full.pca <- prcomp(X.c.t, scale=TRUE)

pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])
pair3D <- as.data.frame(full.pca$x[,1:3])

# fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

# separation based on batch effect
ggplot(data = pair1, aes(PC1, PC2, color=batch))+geom_point() + 
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color=batch))+geom_point()+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[3]), ' %'), y = paste0('variance: ', round(pve[4]), ' %'))
ggplot(data = pair1, aes(PC1, PC2, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))


# COMBAT
# mod <- model.matrix(~label, data = X.c.t)
# mod0 <- model.matrix(~1, data=X.c.t)
modcombat <- model.matrix(~1, data=as.data.frame(X.c.t))
# modcombat
class(modcombat)

dim(X.c)
length(batch)
# needed to transpose the matrix to work

X.comb <- ComBat(X.c, batch=batch, mod=NULL)
# ComBat(X.c.t, batch=batch, mod=mod0, par.prior=TRUE, prior.plots=FALSE)

dim(X.comb)

# subtle adjustment between the original and the combat matrix
X.comb[1:5,1:5]
X.c[1:5,1:5]

# transpose for PCA
X.comb <- as.data.frame(X.comb)
X.comb.t <- t(X.comb)
# pca.comb <- prcomp(X.comb.t, scale=TRUE)

pair1 <- as.data.frame(pca.comb$x[,1:2])
pair2 <- as.data.frame(pca.comb$x[,3:4])

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

# holy hell its worked
ggplot(data = pair1, aes(PC1, PC2, color=batch))+geom_point() + 
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color=batch))+geom_point()+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[3]), ' %'), y = paste0('variance: ', round(pve[4]), ' %'))
ggplot(data = pair1, aes(PC1, PC2, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 1 - 2 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))
ggplot(data = pair2, aes(PC3, PC4, color = bct.vec))+geom_point()+
  scale_color_manual(values=cols[c(3,5)])+
  labs(title="PCA 3 - 4 Discovery & Iris Dataset", x =paste0('variance: ', round(pve[1]), ' %'), y = paste0('variance: ', round(pve[2]), ' %'))


# split the combat normalized matrix back into discovery and val datasets
dim(status.idx)[1]
dim(status.i.idx)[1]
dim(X.comb.t)
X.dis <- X.c.t[batch == 'discovery',]
X.val <- X.c.t[batch == 'iris',]
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
  labs('Scatter Plot of Cybersort Neutrophil % Prediction against Clinical Data', x='Cybersort Prediction', y='Clinical Neut Count')+
  geom_smooth(method=lm, level=0.95)

cor(df.4$Neutrophils, df.4$perc_neut)

# as.character(status.cyber.idx$my_category_2) == rownames(X.dis)


############ LIMMA ############
# mean/variance calculations
x_var <- apply(X.dis, 2, var)
x_mean <- apply(X.dis, 2, mean)
df <- data.frame(log2(x_var), log2(x_mean))
colnames(df) <- c('V1', 'V2')

ggplot(df, aes(V2, V1)) +
  geom_vline(xintercept=log2(5))+
  geom_point(size = 0.2, stroke = 0, shape = 16)+
  labs(title="Mean Variance Scatter Plot",
       x ="log2 Mean Expressioin", y = "log2 Variance")

X.dis.fit <- as.data.frame(X.dis[,x_mean > 5])
dim(X.dis.fit)

# initialize d.dis.lim to add the cybersort covariates
X.dis.lim <- X.dis.fit

X.dis.lim$label <- status.idx$most_general
X.dis.lim$sex <- status.idx$Sex
X.dis.lim$age <- status.idx$Age..months.
X.dis.lim$Neutrophils <- status.cyber.idx$Neutrophils
X.dis.lim$Monocytes <- status.cyber.idx$Monocytes
X.dis.lim$NK.cells.resting <- status.cyber.idx$NK.cells.resting
X.dis.lim$NK.cells.activated <- status.cyber.idx$NK.cells.activated
X.dis.lim$T.cells.CD8 <- status.cyber.idx$T.cells.CD8

# check 8 columns added
dim(X.dis.lim)

### DESIGN MATRIX
design <- model.matrix(~label + sex + age + round(Neutrophils, 3) +
                         round(Monocytes, 3) + round(NK.cells.resting, 3) +
                         round(NK.cells.activated, 3) + round(T.cells.CD8, 3) + 0,
                       data = X.dis.lim)
colnames(design)<- c("bct","greyb",'greyu', "greyv", 'vrl', 'HC', 'sexM', 'age', 'neut', 'mono', 'nk.rest', 'nk.act', 'CD8')

design[1:5,]
round(status.cyber.idx[1:5,c('Neutrophils', 'Monocytes')],3)

# check sums of design
dim(design)
colSums(design)

contrast.matrix <- makeContrasts("bct-vrl", levels=design)
# contrast.matrix<- makeContrasts("HC-bct", 'HC-vrl', levels=design)
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


####### SELECTION OF TRANSCRIPTS #######
# subset the disc and validation matrices by the sig genes
dim(results)
results.tot <- ifelse(results[,1] == 0, FALSE, TRUE)
dim(X.dis.fit[,results.tot])
X.diff <- X.dis.fit[,results.tot]
X.diff.val <- X.val[,match(names(results[results.tot,]), colnames(X.val))]

dim(X.diff)
dim(X.diff.val)

# filter out the healthy controls
X.dis <- X.diff[status.idx$most_general != 'healthy_control',]
X.val <- X.diff.val[status.i.idx$most_general != 'healthy_control',]
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
status.i.idx.d <- status.i.idx[status.i.idx$most_general != 'healthy_control',]

dim(X.dis)
dim(status.idx.d)

dim(X.val)
dim(status.i.idx.d)



###### PREDICTION ######
# reset pseudo labels


scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

# scale the transcript data
# X.s <- data.frame(X.r)
X.s<-data.frame(apply(X.dis, 2, scale01))
dim(X.s)


# add the class labels to the transcript data
# X.s$bacterial <- status.idx.d$most_general == 'bacterial'
# X.s$p.b <- status.idx.d$most_general == 'probable_bacterial'
X.s$bct <- NULL
table(status.idx$most_general)
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
X.s$bct <- status.idx.d$most_general == 'bacterial'
sum(X.s$bct)


## STRATIFIED SAMPLING

prop1 <- 0.75
prop2 <- 0.7
prop3 <- 0.9

boot <- 32
 
# # STRAT
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


## CROSS VALIDATION
logistic.m <- NULL
knn.m <- NULL
randForrest.m <- NULL
nn.m <- NULL
svm.m <- NULL

n_folds <- 10
n.train <- round(nrow(X.s))
folds.i <- sample(rep(1:n_folds, length.out = n.train))
cv_tmp <- matrix(NA, nrow = n_folds, ncol = length(df))
for (k in 1:n_folds) {
  test.i <- which(folds.i == k)
  train.cv <- X.s[-test.i, ]
  test.cv <- X.s[test.i, ]
  
  # factoring cv data
  train.cv.fac <- train.cv
  train.cv.fac$bct <- as.factor(train.cv$bct)
  test.cv.fac <- test.cv
  test.cv.fac$bct <- as.factor(test.cv$bct)
  
  # dim(train.cv)[1] +dim(test.cv)[1]
  
  print(paste0('fold: ', k))
  
  ### LOGISTIC REGRESSION
  model <- glm(bct~ ., data=train.cv, family=binomial(link='logit'), maxit = 128)
  # summary(model)
  # anova(model, test="Chisq")
  pred.test <- predict(model, test.cv[-ncol(test.cv)])
  pr.test <- prediction(pred.test, test.cv$bct)
  
  logistic.m[k] <- pr.test %>%
    performance(measure = "auc") %>%
    .@y.values
  
  
  ### KNN
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 5)
  model <- train(bct ~., data = train.cv.fac, method = "knn",
                 trControl=trctrl,
                 tuneLength = 10)
  
  knn.opt <- model$results$k[which.max(model$results$Accuracy)]
  # train
  p <- knn(train.cv[-ncol(train.cv)], test.cv[-ncol(test.cv)], train.cv$bct,  k=knn.opt, prob=TRUE)
  
  # display the confusion matrix
  # table(p, test$bct)
  # attributes(pred)
  p<-attr(p, "prob")
  p<-1-p
  
  # plot(model, print.thres = 0.5, type="S")
  # confusionMatrix(p, test$bct)
  pr <- prediction(p, test.cv$bct)
  
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  
  knn.m[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  

  ### RANDOM FORREST
  model <- randomForest(bct ~ . , data = train.cv.fac)
  # plot(model)
  # attributes(model)
  # model$mtry # number of variables considered by each tree
  pred<-predict(model , test.cv.fac[-ncol(test.cv.fac)])
  # attributes(pred)
  # model=randomForest(x,y,xtest=x,ytest=y,keep.forest=TRUE)
  model.prob <- predict(model, test.cv.fac, type="prob")
  p <- 1-model.prob[,1]
  pr <- prediction(p, test.cv$bct)
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)
  randForrest.m[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  ### Neural Net
  nn <- neuralnet(bct~ ., train.cv, linear.output = FALSE, act.fct = "logistic")
  p <- predict(nn, test.cv[-ncol(test.cv)])
  pr <- prediction(p, test.cv$bct)
  nn.m[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  
  ### SVM
  model <- svm(bct ~ . , train.cv.fac, probability = TRUE)
  pred <- predict(model, test.cv.fac, probability = TRUE)
  p<-attr(pred, "prob")
  pr <- prediction(1-p[,1], test.cv$bct)
  # plot(prf)
  svm.m[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}


logistic.m <- NULL
knn.m <- NULL
randForrest.m <- NULL
nn.m <- NULL
svm.m <- NULL
for (i in 1:8){
  print(paste0('iter: ', i))
  
  ### LOGISTIC REGRESSION
  model <- glm(bct~ ., data=train, family=binomial(link='logit'), maxit = 128)
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
  model <- train(bct ~., data = train.fac, method = "knn",
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
  model <- randomForest(bct ~ . , data = train.fac)
  # plot(model)
  # attributes(model)
  # model$mtry # number of variables considered by each tree
  pred<-predict(model , test.fac[-ncol(test)])
  # attributes(pred)
  # model=randomForest(x,y,xtest=x,ytest=y,keep.forest=TRUE)
  model.prob <- predict(model, test.fac, type="prob")
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
  model <- svm(bct ~ . , train.fac, probability = TRUE)
  pred <- predict(model, test.fac, probability = TRUE)
  p<-attr(pred, "prob")
  pr <- prediction(1-p[,1], test$bct)
  # plot(prf)
  svm.m[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

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


###### HYPER PARAM OPTIMIZATION FOR PROMISING MODELS ######

# remove pseudo labels

# X.dis <- X.diff[status.idx$most_general != 'healthy_control',]
# X.val <- X.diff.val[status.i.idx$most_general != 'healthy_control',]
# status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
# status.i.idx.d <- status.i.idx[status.i.idx$most_general != 'healthy_control',]
# 
# dim(status.idx)
# 
# status.idx.d <- status.idx[status.idx.d$most_general == 'bacterial',]
# X.s$bct <- status.idx.d$most_general == 'bacterial'
# table(status.idx.d$most_general)
# sum(X.s$bct)


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
h.n <- 40
roc.a <- NULL
roc.t <- NULL
j.train <- NULL
j.test <- NULL
h.n.hx <- NULL
roc.train <- NULL
roc.train.me <- NULL
roc.test <- NULL
roc.test.me <- NULL
boot <- 16
for(i in 1:h.n){
  for(j in 1:boot){
    print(paste0('hidden Nodes: ', i, ', bootstrap: ', j))
    index <- sample(nrow(train), round(prop2*nrow(train)))
    train.cv <- train[index,]
    test.cv <- train[-index,]
    # dim(train.cv)
    # dim(test.cv)
    nn1 <- neuralnet(bct~ ., train.cv, linear.output = FALSE, act.fct = "logistic",
                     hidden = c(i), rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
    
    pred.train <- predict(nn1, train.cv[-ncol(train.cv)])
    pred.test <- predict(nn1, test.cv[-ncol(test.cv)])
    
    j.train[j] <- prediction(pred.train[,1], train.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
    j.test[j] <- prediction(pred.test[,1], test.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
  }
    full.list <- c(j.train, j.test)
    full.df <- data.frame(matrix(unlist(full.list), nrow=length(full.list), byrow=T))
    colnames(full.df) <- 'roc.A'
    
    full.df$class <- as.factor(sort(rep(seq(1:(length(full.list)/boot)), boot)))
    
    roc.stats <- full.df %>%
      group_by(class) %>%
      summarise(roc.m = mean(roc.A), roc.med = median(roc.A), roc.sd = sd(roc.A))
    
    roc.stats <- roc.stats %>% mutate(
      roc.se = roc.sd/sqrt(boot),
      z.stat = qnorm(0.975),
      roc.me = z.stat * roc.se
    )
    h.n.hx[i] <- i
    roc.train[i] <- roc.stats$roc.m[1]
    roc.train.me[i] <- roc.stats$roc.me[1]
    roc.test[i] <- roc.stats$roc.m[2]
    roc.test.me[i] <- roc.stats$roc.me[2]
}


mod.complexity.df <- as.data.frame(cbind(roc.train, roc.test, roc.train.me, roc.test.me, h.n.hx))
colnames(mod.complexity.df) <- c('train', 'test', 'train.me', 'test.me', 'h.n')

pd <- position_dodge(0.2)
ggplot(mod.complexity.df, aes(x=h.n, y=train, color='train')) +
  # scale_y_continuous(limits = c(,1))+
  geom_line(aes(y=train))+
  geom_errorbar(aes(ymin=train-train.me, ymax=train+train.me), width=.4, position=pd)+
  geom_line(aes(y=test, color='test'))+
  geom_errorbar(aes(ymin=test-test.me, ymax=test+test.me), width=0.4, color='red')+
  labs(title=paste0('Bias-Variance Trade Off'), x =paste0('1 - ', h.n, ' hidden nodes'), y = "ROCA")

which.max(mod.complexity.df$test)

### NEURAL HYPER GRID OPTIMIZATION
# hyper_grid <- NULL
# hyper_grid <- expand.grid(
#   # h.n       = seq(1, 100, by = 3), # beyond 40 overfitting
#   h.n       = seq(5, 40, by = .2), # maxed between 10 - 20 roc.a ~ 90
#   train.cv.prop = seq(from=60, to=70, by=5)/100
# )
# 
# head(hyper_grid)
# dim(hyper_grid)
# 
# roc.a <- NULL
# for(i in 1:nrow(hyper_grid)) {
#   print(paste0('iter: ', i))
#   # train model
#   index <- sample(nrow(train), round(hyper_grid$train.cv.prop[i] * nrow(train)))
#   train.cv <- train[index,]
#   test.cv <- train[-index,]
#   nn1 <- neuralnet(bct~ ., train.cv, linear.output = FALSE, act.fct = "logistic",
#                    hidden = c(hyper_grid$h.n[i]), rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
#   pred <- predict(nn1, test.cv[-ncol(test.cv)])
#   # extract error
#   roc.a[i] <- prediction(pred, test.cv$bct) %>%
#     performance(measure = "auc") %>%
#     .@y.values
# }
# 
# hyper_grid$roc.a <- unlist(roc.a)
# top.n <- nrow(hyper_grid)
# 
# # arranges the highest top.n rows of hypergrid roc.a values in descending order
# hyper_top <- hyper_grid %>% 
#   dplyr::arrange(desc(unlist(roc.a)))%>%
#   head(top.n)
# 
# dim(hyper_top)
# 
# # h.n to factor for plotting
# hyper_top$h.n <- as.factor(hyper_top$h.n)
# 
# # compute mean, median, sdev values for top.n rows in hyper_grid
# hyp.df <- hyper_top %>%
#   group_by(h.n) %>%
#   summarise(roc.m = mean(roc.a), roc.med = median(roc.a), roc.sd = sd(roc.a))
# 
# # ggplot(hyper_top, aes(roc.a, fill=h.n, color=h.n)) + geom_density(alpha=0.2)
# ggplot(hyper_top, aes(h.n, roc.a, fill=h.n, color=h.n)) + geom_boxplot(alpha=0.7)
# 
# print(paste0('combined mean and median maxed at: ', hyp.df$h.n[which.max(hyp.df$roc.m+hyp.df$roc.med/2)], ' hidden node(s)'))
# 
# hyp.df
# hyp.df[which.max(hyp.df$roc.m+hyp.df$roc.med/2),]
# 
# # check the number of times h.n occured in hyper_top
# # ensures that we do not set h.n.opt to a freak one off occurence of high roc
# sort(table(hyper_top$h.n), decreasing = TRUE)
# 
# 
# # top 10 rows
# hyper_top[1:20,]




# weight initialization. Cant get neuralnet package to accept so just using random init
# xavier.w <- rnorm(dim(train)[1], mean = 0, sd = sqrt(1/dim(train)[1]))
# hist(xavier.w)

### TEST SET NEURAL

opt.h.n <- which.max(mod.complexity.df$test)
# opt.h.n <- 1
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
model <- randomForest(bct ~ . , data = train.fac, nodesize = 1, maxnodes = NULL)
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
    data = train.fac, 
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
rf.test <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  model <- randomForest(formula = bct ~ .,
                        data = train.fac,
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
df.1$algo <- as.factor(ifelse(df.1$algo == 1, 'neural_net', 'randomForest'))
df.1$algo
ggplot(df.1, aes(roc.a, color = algo, fill=algo)) + geom_density(alpha=.1)





####### PSEUDO LABELING #######
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
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

# train iter neural nets on randomly drawn train sets
# evaluate on test set. Any PB with P bct > bct.thresh is labeled as bct but the dx label is not changed until outside of the loop

roc.h <- NULL
bct.thresh <- 0.99
iters <- 50
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
  
  if(max(pred[pb.filter,]) > bct.thresh){
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
  
  nn2 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn2, test[-ncol(test)])
  
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


#################################################
nn.psd.val <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  
  
  nn2 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  
  pred.2 <- predict(nn2, X.val)
  pr <- prediction(pred.2, status.i.idx$most_general == 'bacterial')
  
  nn.psd.val[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}


# remove pseudo labels
status.idx.d <- status.idx[idx.d,]
X.s$bct <- status.idx.d$most_general == 'bacterial'
table(status.idx.d$most_general)
print(paste0('bacterial cases: ', sum(X.s$bct)))

nn.norm.val <- NULL
for (i in 1:boot){
  print(paste0('boot: ', i))
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train <- X.s[index, ]
  test <- X.s[-index, ]
  
  nn1 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
                   hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
  pred <- predict(nn1, X.val)
  pr <- prediction(pred, status.i.idx$most_general == 'bacterial')
  
  nn.norm.val[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

mean(unlist(nn.psd.val))
median(unlist(nn.psd.val))
sd(unlist(nn.psd.val))

mean(unlist(nn.norm.val))
median(unlist(nn.norm.val))
sd(unlist(nn.norm.val))

# index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
# train <- X.s[index, ]
# test <- X.s[-index, ]
# 
# nn2 <- neuralnet(bct~ ., train, linear.output = FALSE, act.fct = "logistic",
#                  hidden = c(opt.h.n), rep = 3, stepmax = 1e+06, startweights = NULL)
# 
# pred.2 <- predict(nn2, X.val)
# pr <- prediction(pred.2, status.i.idx$most_general == 'bacterial')
# pr %>%
#   performance(measure = "auc") %>%
#   .@y.values



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















###### exploratory data analysis ######
### unsup
# idx <- status['most_general'] == 'bacterial' |
#   status['most_general'] == 'viral' |
#   status['most_general'] == 'greyb' |
#   status['most_general'] == 'greyv'|
#   status['most_general'] == 'greyu'
# sum(idx)
# dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral')

### supervised
# idx <- status$most_general == 'bacterial' |
#   status$most_general == 'viral' |
#   status$most_general == 'greyb' |
#   status$most_general == 'greyv'|
#   status$most_general == 'greyu' |
#   status$most_general == 'HC'
# sum(idx)
# class(idx)
# dx <- c('bacterial', 'probable_bacterial', 'unknown', 'probable_viral', 'viral', 'healthy_control') # supervised
# 
# ### outlier
# which(status$my_category_2 == 'bacterialgpos_19_SMH')
# idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')]
# idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')] <- FALSE
# idx[which(status$my_category_2 == 'bacterialgpos_19_SMH')]
# 
# status.idx <- status[idx,]
# status.idx$most_general[1:10]
# 
# # rename most_general
# status.idx$most_general <- as.character(status.idx$most_general)
# status.idx$most_general[status.idx$most_general == 'greyb'] <- 'probable_bacterial'
# status.idx$most_general[status.idx$most_general == 'greyu'] <- 'unknown'
# status.idx$most_general[status.idx$most_general == 'greyv'] <- 'probable_viral'
# status.idx$most_general[status.idx$most_general == 'HC'] <- 'healthy_control' # toggle for unsupervised
# 
# status.idx$most_general <- as.factor(status.idx$most_general)
# 
# levels(status.idx$most_general)
# status.idx$most_general <- factor(status.idx$most_general, levels = dx)
# levels(status.idx$most_general)
# # status.idx$most_general
# 
# status.idx$array.contemporary.CRP <- as.numeric(as.character(status.idx$array.contemporary.CRP))
# 
# dim(status.idx)
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



# ###### PCA ######
# full.pca <- prcomp(X.r, scale=TRUE) # unsupervised
# # full.pca <- prcomp(t(X[results.tot,]), scale=TRUE) # supervised
# 
# pair1 <- as.data.frame(full.pca$x[,1:2])
# pair2 <- as.data.frame(full.pca$x[,3:4])
# pair3D <- as.data.frame(full.pca$x[,1:3])
# 
# fviz_eig(full.pca)
# 
# ve <- full.pca$sdev^2
# pve <- ve/sum(ve)*100
# pve[1:5]
# 
# 
# # # most_gen 2D Age
# # p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status.idx$most_general), size = status.idx$Age..months.,
# #              colors=c(dx.cols), text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', wbc, '<br>CRP: ', crp, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
# #   add_markers() %>%
# #   layout(title = 'PCA of Diagnostic Groups, Age Size Mapping',
# #          xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
# #          yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
# # p
# 
# # # most_gen 2D
# p <- plot_ly(pair1, x = ~PC1, y = ~PC2, color = status.idx$most_general, size = status.idx$array.contemporary.CRP,
#              colors=cols, text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'PC 1-2 of Diagnostic Groups, CRP Size Mapping',
#          xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#          yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')))
# p
# p <- plot_ly(pair2, x = ~PC3, y = ~PC4, color = status.idx$most_general, size = status.idx$Age..months.,
#              colors=cols, text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'PC 3-4 of Diagnostic Groups, Age Size Mapping',
#          xaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)')),
#          yaxis = list(title = paste0("PC4: (", round(pve[4],2), '%)')))
# p
# 
# # api_create(p, filename = "2d_pca_filt")
# 
# ## most_gen 3D
# p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$most_general, size = status.idx$Age..months.,
#              colors=c(dx.cols), text= ~paste0('<br>age: ', status.idx$Age..months., '<br>Sex:', status.idx$Sex, '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# p
# # api_create(p, filename = "3d_pca_dx_unfilt")
# 
# # # more_gen
# # dx.cols.f <- c('#bd35fc', "#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
# # plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$more_general, size = status.idx$Age..months.,
# #         colors = c(dx.cols.f), text= ~paste0('category: ', status.idx$category, '<br>age: ', status.idx$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status.idx$my_category_2, '<br>Diagnosis: ',status.idx$Diagnosis)) %>%
# #   add_markers() %>%
# #   layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
# #          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
# #                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
# #                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# # 
# # # category
# # p<-plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$category, size = status[idx,]$Age..months.,
# #         colors=c(cat.pal), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
# #   add_markers() %>%
# #   layout(title = 'Category Groups by PCA 1-2-3, Age Size Mapping',
# #          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
# #                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
# #                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# # api_create(p, filename = "3d_pca_cat")
# # 
# 
# # sex
# # View(status.idx)
# # p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status.idx$Sex, size = ~status.idx$Age..months.,
# #              colors=c(sex.cols), text= ~paste0('<br>age: ', status.idx$Age..months., '<br>Sex:', status.idx$Sex, '<br>WBC: ', status.idx$WBC, '<br>CRP: ', status.idx$array.contemporary.CRP, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
# #   add_markers() %>%
# #   layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
# #          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
# #                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
# #                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# # p
# # api_create(p, filename = "3d_pca_sex")
# 
# # 
# # # site
# # plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$site, size = ~status[idx,]$Age..months.,
# #         colors = c(site.pal), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
# #   add_markers() %>%
# #   layout(title = 'Site Recruitment by PCA 1-2-3, Age Size Mapping',
# #          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
# #                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
# #                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
# 



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


