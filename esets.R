library(dplyr)
library(ggfortify)
library(glmnet)
library(ROCR)
library(dplyr)
library(ggplot2)
library(neuralnet)

# library(GGally)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()

rm(list=setdiff(ls(), "x"))
load('esets.RData')


dim(e.set) # megaExp data
dim(status) # megaExp labels

dim(e.set.i) # iris data
dim(status.iris) # iris labels

# attributes(status)
# status$Diagnosis
# status$category

# transpose
e.set.t <- t(e.set)
e.set.i.t <- t(e.set.i)

dim(e.set.t)
dim(e.set.i.t)

# rm(targets, targets.iris, e.set, e.set.i)
# e.set.t[1:5,(ncol(e.set.t)-2):ncol(e.set.t)]

# scale
e.set.s <- scale(e.set.t)
e.set.i.s <- scale(e.set.i.t)

# e.set.i.s <- as.data.frame(scale(e.set.i.t))
# e.set.s <- as.data.frame(scale(e.set.t))
# e.set.i.s <- as.data.frame(scale(e.set.i.t))

# apply(e.set.s, 2, var)
# apply(e.set.s, 2, mean)


# labels
label.i <- rownames(e.set.i.t)
label.s.i <- status.iris$My_code
common <- match(label.i, label.s.i)

e.set.l <- as.character(status$most_general)
e.set.i.l <- as.character(status.iris$most_general[common]) # pass common index to extract common iris labels

e.set.df <- as.data.frame(e.set.s)
e.set.i.df <- as.data.frame(e.set.i.s)

e.set.df$label <- e.set.l
e.set.df$exp <- 0

e.set.i.df$label <- e.set.i.l
e.set.i.df$exp <- 1

e.set.df[1:5, (ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:5, (ncol(e.set.i.df)-4):ncol(e.set.i.df)]


# combine datasets
common.genes <-intersect(colnames(e.set.df), colnames(e.set.i.df))
length(common.genes) # 42004

dim(e.set.df[,common.genes])
dim(e.set.i.df[,common.genes])

e.set.df[1:10,(ncol(e.set.df)-4):ncol(e.set.df)]
e.set.i.df[1:10,(ncol(e.set.i.df)-4):ncol(e.set.i.df)]

e.set.c <- rbind(e.set.df[,common.genes], e.set.i.df[,common.genes])

dim(e.set.c)
e.set.c[1:10,(ncol(e.set.c)-5):ncol(e.set.c)]
e.set.c[(nrow(e.set.c)-10):nrow(e.set.c),(ncol(e.set.c)-5):ncol(e.set.c)]

e.set.r <- e.set.c[,-((ncol(e.set.c)-1):ncol(e.set.c))]
dim(e.set.r)


e.set.c$label == 'bacterial'
e.set.c$label == 'viral'
bac.bin <- ifelse(e.set.c$label == 'bacterial', 1, 0)

mod.label <- ifelse(e.set.c$label == 'bacterial', 'bacterial', ifelse(e.set.c$label == 'viral', 'viral', 'other'))

remove(e.set, e.set.df, e.set.s, e.set.i.s, e.set.i, e.set.t)
# ls()


## PCA
e.set.pca <- prcomp(e.set.r, scale = FALSE)
# e.set.pca <- gset.pca
summary(gset.pca)
plot(e.set.pca, type = 'l')

# autoplot(gset.pca, data = e.set.c, colour = 'exp')

pair1 <- e.set.pca$x[,1:2]
pair2 <- e.set.pca$x[,3:4]
pair3 <- e.set.pca$x[,5:6]
ggplot(pair1, aes(PC1, PC2, color=e.set.c$exp)) + geom_point() +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=e.set.c$exp)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=e.set.c$exp)) + geom_point()

ggplot(pair1, aes(PC1, PC2, color=mod.label)) + geom_point()
ggplot(pair2, aes(PC3, PC4, color=mod.label)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=mod.label)) + geom_point()


## OUTLIERS
a<-pair2[,2]>200
which(a)
# OD_38_KEN KDno_1_UCSD 
# 41         434

b<-pair3[,2]>180
which(b)
# OD_38_KEN HC_53_SMH 
# 41       489

outlier <- c(41, 434, 489)
e.set.clean <- e.set.r[-outlier,]
mod.label.clean <- mod.label[-outlier]
exp.label.clean <- e.set.c$exp[-outlier]
dim(e.set.clean)

e.set.c['OD_38_KEN',1:10]
e.set.clean['OD_38_KEN',1:10]


## PCA 2
e.set.pca <- prcomp(e.set.clean, scale = FALSE)
# e.set.pca <- gset.pca
summary(e.set.pca)
plot(e.set.pca, type = 'l')
# autoplot(gset.pca, data = e.set.c, colour = 'label')
# autoplot(gset.pca, data = e.set.c, colour = 'exp')

pair1 <- e.set.pca$x[,1:2]
pair2 <- e.set.pca$x[,3:4]
pair3 <- e.set.pca$x[,5:6]
ggplot(pair1, aes(PC1, PC2, color=exp.label.clean)) + geom_point() +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=exp.label.clean)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=exp.label.clean)) + geom_point()

ggplot(pair1, aes(PC1, PC2, color=mod.label.clean)) + geom_point()
ggplot(pair2, aes(PC3, PC4, color=mod.label.clean)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=mod.label.clean)) + geom_point()


ve <- e.set.pca$sdev^2
pve <- ve / sum(ve)
round(pve, 5)

remove(e.set.c, e.set.r, e.set.i.df)

## DGE
bct <- e.set.clean[mod.label.clean == 'bacterial',]
vrl <- e.set.clean[mod.label.clean == 'viral',]
dim(bct)
dim(vrl)

bct.mean <- apply(bct, 2, mean)
vrl.mean <- apply(vrl, 2, mean)

head(bct.mean)
head(vrl.mean)

# Just get the maximum of all the means
limit = max(bct.mean, vrl.mean)

# Scatter plot
plot(vrl.mean ~ bct.mean, xlab = "bct", ylab = "vrl",
     main = "Bct vs Vrl - Scatter", xlim = c(0, limit), ylim = c(0, limit))
cor(bct.mean, vrl.mean)

# Compute fold-change (biological significance)
fold = bct.mean - vrl.mean

# Histogram of the fold differences
hist(fold, col = "gray")

# Compute statistical significance (using t-test)
pvalue = NULL # Empty list for the p-values
tstat = NULL # Empty list of the t test statistics

for(i in 1 : ncol(bct)) { # For each gene : 
  x = bct[,i] 
  y = vrl[,i]
  t = t.test(x, y)
  pvalue[i] = t$p.value
  tstat[i] = t$statistic
}

# Histogram of p-values (-log10)
hist(-log10(pvalue), col = "gray")

# Volcano: put the biological significance (fold-change)
# and statistical significance (p-value) in one plot
plot(fold, -log10(pvalue), main = "Volcano #1")

fold_cutoff = 1.5
pvalue_cutoff = 0.000000000000001
abline(v = fold_cutoff, col = "blue", lwd = 2)
abline(v = -fold_cutoff, col = "red", lwd = 2)
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 2)


# Fold-change filter for "biological" significance
filter_by_fold = abs(fold) >= fold_cutoff
dim(e.set.clean[filter_by_fold,])
dim(e.set.clean[filter_by_fold,])

# P-value filter for "statistical" significance
filter_by_pvalue = pvalue <= pvalue_cutoff
dim(e.set.clean[filter_by_pvalue,])

# Combined filter (both biological and statistical)
filter_combined = filter_by_fold & filter_by_pvalue

filtered = e.set.clean[filter_combined,]
dim(filtered)

e.set.clean[,filter_combined][1:5,1:5]


## TEST TRAIN DIFFERENTIAL EXPRESSED TRANSCRIPTS VIRAL BACTERIAL
bac.vrl.index <- mod.label.clean == 'bacterial' | mod.label.clean == 'viral'
x <- e.set.clean[bac.vrl.index,filter_combined]

y <- mod.label.clean[bac.vrl.index]
dim(x)
length(y)

x[1:5,1:5]
y[1:5]

set.seed(3)
n <- nrow(x)
index <- seq(1:n)
train = sample(1:n, round(n*0.8))
test = index[-train]
intersect(train, test)

x_train <- x[train,]
x_train$label <- ifelse(y[train] == 'viral', TRUE, FALSE)
x_test <- x[test,]
dim(x_train)
dim(x_test)
ytest = y[test]

## NEURAL NET
set.seed(123)

nn.1 <- neuralnet(f, data = x_train, hidden=c(10), err.fct = 'ce', likelihood = TRUE)

net.names <- names(x_train)[-(length(x_train))]
f <- as.formula(paste("label ~", paste(net.names[!net.names %in% "label"], collapse = " + ")))
f
nn <- neuralnet(label ~ 3130600, hidden=c(10), data = x_train)

nn <- neuralnet(f,data=x_train, hidden=c(10),linear.output=T)

nn <- neuralnet(
  label ~ .,
  data=x_train, hidden=2, err.fct='ce',
  linear.output=FALSE)



x <- cbind(seq(-10, 40, 1), seq(51, 100, 1))
y <- x[,1]*x[,2]
y <- ifelse(y > 3000, TRUE, FALSE)
colnames(x) <- c('x1', 'x2')
names(y) <- 'y'
dt <- data.frame(x, y)
dt
model <- neuralnet(label~., x_train, hidden=10, threshold=0.01, likelihood = TRUE)

model <- neuralnet(label ~ "3130600", x_train, hidden=10, threshold=0.01, likelihood = TRUE)
colnames(x_train)




model <- neuralnet(label ~ "3130600", x_train, hidden=10, threshold=0.01, likelihood = TRUE)
feats <- names(x_train)[-(ncol(x_train))]
f <- paste(feats,collapse=' + ')
f <- paste('label ~',f)
f <- as.formula(f)
f
library(neuralnet)
nn <- neuralnet(f,x_train,hidden=c(10,10,10),linear.output=FALSE)












# logistic
logistic.mod1 <- glm(label ~., family = "binomial", data = x_train)
summary(logistic.mod1)

log.pred <- predict(logistic.mod1, x_test, type = 'response')

prediction(log.pred, ytest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

# model 1 AUC
prediction(log.pred, ytest) %>%
  performance(measure = "auc") %>%
  .@y.values

# extracts gene expression data from confident model predictions
x[names(which.min(log.pred)),] # looks bacterial
log.pred[which.min(log.pred)] # 5% prob of being viral

x[names(which.max(log.pred)),] # looks viral
log.pred[which.max(log.pred)] # 98% prob of being viral

# manual probability prediction double check
regression <- sum(logistic.mod1$coefficients[-1] * x[names(which.min(log.pred)),]) + logistic.mod1$coefficients[1]
regression <- sum(logistic.mod1$coefficients[-1] * x[names(which.max(log.pred)),][-ncol(x)]) + logistic.mod1$coefficients[1]
1/(1+exp(1)^-(regression))


# find optimal prediction cutoff to maximise f1 score
cutoff <- seq(0.01, 0.99, 0.01)
f1.list <- c()
for (i in 1:length(cutoff)){
  # print(cutoff[i])
  cat.pred <- ifelse(log.pred < cutoff[i], 'bacterial', 'viral')
  table(ytest, cat.pred)
  tpr <- table(ytest, cat.pred)[1] / sum(table(ytest, cat.pred)[1] + table(ytest, cat.pred)[3]) # sensitivity / recall
  tnr <- table(ytest, cat.pred)[4] / sum(table(ytest, cat.pred)[4] + table(ytest, cat.pred)[2]) # specificity
  ppp <- table(ytest, cat.pred)[1] / sum(table(ytest, cat.pred)[1] + table(ytest, cat.pred)[2]) # precision
  f1 <- 2 * (ppp * tpr) / (ppp + tpr)
  f1.list[i] = f1
}
plot(cutoff, f1.list)
opt.cutoff <- cutoff[which.max(f1.list)]
f1.list[which.max(f1.list)]
# cutoff <- 0.6
cat.pred <- ifelse(log.pred < opt.cutoff, 'bacterial', 'viral')

table(ytest, cat.pred)

# saveRDS(, file = 'gset_GSE72809')
#asdf