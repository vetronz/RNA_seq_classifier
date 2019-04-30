library(dplyr)
library(ggfortify)
library(glmnet)
library(ROCR)
library(dplyr)
library(ggplot2)
library(boot)
library(neuralnet)
# library(knitr)
# library(GGally)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

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
# dim(e.set.t)
# dim(e.set.i.t)

# scale
e.set.s <- scale(e.set.t)
e.set.i.s <- scale(e.set.i.t)
# e.set.s <- e.set.t
# e.set.i.s <- e.set.i.t

# labels
label.i <- rownames(e.set.i.t)
label.s.i <- status.iris$My_code
common <- match(label.i, label.s.i)

e.set.l <- as.character(status$most_general)
e.set.i.l <- as.character(status.iris$most_general[common]) # pass common index

e.set.df <- as.data.frame(e.set.s)
e.set.i.df <- as.data.frame(e.set.i.s)

e.set.df$label <- e.set.l
e.set.df$exp <- 'MegaExperiment'

e.set.i.df$label <- e.set.i.l
e.set.i.df$exp <- 'IRIS'

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
# e.set.c[1:10,(ncol(e.set.c)-5):ncol(e.set.c)]
# e.set.c[(nrow(e.set.c)-10):nrow(e.set.c),(ncol(e.set.c)-5):ncol(e.set.c)]

e.set.r <- e.set.c[,-((ncol(e.set.c)-1):ncol(e.set.c))]
dim(e.set.r)


mod.label <- ifelse(e.set.c$label == 'bacterial', 'bacterial',
                    ifelse(e.set.c$label == 'viral', 'viral', 'other'))

remove(e.set, e.set.df, e.set.s, e.set.i.s, e.set.i, e.set.t, common.genes, common)
ls()

## PCA
e.set.pca <- prcomp(e.set.r, scale = FALSE)
# e.set.pca <- gset.pca
summary(gset.pca)
plot(e.set.pca, type = 'l')
# autoplot(gset.pca, data = e.set.c, colour = 'exp')
experiment <- e.set.c$exp
organism <- mod.label

pair1 <- e.set.pca$x[,1:2]
pair2 <- e.set.pca$x[,3:4]
pair3 <- e.set.pca$x[,5:6]
ggplot(pair1, aes(PC1, PC2, color=experiment)) + geom_point() +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=experiment)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=experiment)) + geom_point()

ggplot(pair1, aes(PC1, PC2, color=organism)) + geom_point() +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=organism)) + geom_point()
ggplot(pair3, aes(PC5, PC6, color=organism)) + geom_point()


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
# summary(e.set.pca)
plot(e.set.pca, type = 'l')
experiment <- exp.label.clean
# organism <- mod.label

pair1 <- e.set.pca$x[,1:2]
pair2 <- e.set.pca$x[,3:4]
pair3 <- e.set.pca$x[,5:6]
ggplot(pair1, aes(PC1, PC2, color=experiment)) + geom_point() +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=experiment)) + geom_point() +
  xlab("Third Principal Component") +
  ylab("Fourth Principal Component") +
  ggtitle("Third and Fourth Principal Components of Combined Mega Iris Data")
ggplot(pair3, aes(PC5, PC6, color=experiment)) + geom_point() +
  xlab("Fifth Principal Component") +
  ylab("Sixth Principal Component") +
  ggtitle("Fifth and Sixth Principal Components of Combined Mega Iris Data")


ggplot(pair1, aes(PC1, PC2, color=mod.label.clean)) + geom_point() +
  xlab("First Principal Component") + 
  ylab("Second Principal Component") + 
  ggtitle("First Two Principal Components of Combined Mega Iris Data")
ggplot(pair2, aes(PC3, PC4, color=mod.label.clean)) + geom_point() +
  xlab("Third Principal Component") +
  ylab("Fourth Principal Component") +
  ggtitle("Third and Fourth Principal Components of Combined Mega Iris Data")
ggplot(pair3, aes(PC5, PC6, color=mod.label.clean)) + geom_point() +
  xlab("Fifth Principal Component") +
  ylab("Sixth Principal Component") +
  ggtitle("Fifth and Sixth Principal Components of Combined Mega Iris Data")


ve <- e.set.pca$sdev^2
pve <- ve / sum(ve)
round(pve, 5)[1:10]

remove(e.set.c, e.set.r, e.set.i.df, e.set.i.t, e.set.pca, pair1, pair2, pair3)

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
x <- e.set.clean[bac.vrl.index, filter_combined]
y <- mod.label.clean[bac.vrl.index]

set.seed(3)
n <- nrow(x)
index <- seq(1:n)
train = sample(1:n, round(n*0.7))
test = index[-train]
intersect(train, test)

x_train <- x[train,]
x_test <- x[test,]

y_train <- ifelse(y[train] == 'bacterial', TRUE, FALSE)
y_test = ifelse(y[test] == 'bacterial', TRUE, FALSE)


## NEURAL NET
library(neuralnet)
x.l <- list()
for(i in 1 : ncol(x_train)) { # For each gene : 
  x.l[[i]] <- x_train[,i]
}
length(x.l)

nn_train <-do.call(cbind, x.l) # apply cbind to each item in x.l list
nn_train <- nn_train[,1:5]

nn1 <- neuralnet(y_train~., nn_train,
                 hidden=5,
                 threshold=0.01,
                 err.fct = 'ce',
                 linear.output = FALSE,
                 likelihood = TRUE)

# plot(nn1)
# dev.off()
nn1_train_error <- nn1$result.matrix[1,1]
paste("CE Error: ", round(nn1_train_error, 3)) 
nn1_aic <- nn1$result.matrix[4,1]
paste("AIC: ", round(nn1_aic,3))
nn1_bic <- nn1$result.matrix[5,1]
paste("BIC: ", round(nn1_bic, 3))


nn1.pred <- compute(nn1, x_test[,1:5])
nn1.r <-nn1.pred$net.result

nn1.pred.c <- ifelse(nn1.r > 0.5, 'bacterial', 'viral' )
nn1.pred.c
y[test]

table(y[test], nn1.pred.c)
res <-ifelse(y[test] == 'bacterial',1,0)

# model 1 AUC
detach(package:neuralnet,unload = T)
prediction(nn1.pred$net.result, res) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

prediction(nn1.pred$net.result, res) %>%
  performance(measure = "auc") %>%
  .@y.values


# 2-Hidden Layers, Layer-1 2-neurons, Layer-2, 1-neuron
set.seed(123)
nn1 <- neuralnet(y_train ~., b, 
                 linear.output = FALSE, 
                 err.fct = 'ce', 
                 likelihood = TRUE,
                 hidden = c(2,1))

# 2-Hidden Layers, Layer-1 2-neurons, Layer-2, 2-neurons
set.seed(123)
nn2 <- nn1 <- neuralnet(y_train ~., b, 
                        linear.output = FALSE, 
                        err.fct = 'ce', 
                        likelihood = TRUE, 
                        hidden = c(2,2))

# 2-Hidden Layers, Layer-1 1-neuron, Layer-2, 2-neuron
set.seed(123)
nn3 <- nn1 <- neuralnet(y_train ~., b, 
                        linear.output = FALSE, 
                        err.fct = 'ce', 
                        likelihood = TRUE, 
                        hidden = c(1,2))


# Bar plot of results
nn_comp <- tibble('Network' = rep(c("NN1", "NN2", "NN3"), each = 3), 
                       'Metric' = rep(c('AIC', 'BIC', 'ce Error'), length.out = 9),
                       'Value' = c(nn1$result.matrix[4,1], nn1$result.matrix[5,1], 
                                   nn1$result.matrix[1,1], 
                                   nn2$result.matrix[4,1], nn2$result.matrix[5,1], 
                                   nn2$result.matrix[1,1], nn3$result.matrix[4,1], 
                                   nn3$result.matrix[5,1], nn3$result.matrix[1,1]))


nn_comp %>%
  ggplot(aes(Network, Value, fill = Metric)) +
  geom_col(position = 'dodge')  +
  ggtitle("AIC, BIC, and Cross-Entropy Error of the Classification ANNs")
# AIC measures information lost by model therefore lower better

length(b)
length(x_test[,1])


nn1$model.list

pr.nn <- compute(nn,test_[,1:13])

## LOGISTIC
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
regression <- sum(
  logistic.mod1$coefficients[-1] * x[names(which.min(log.pred)),]) +
  logistic.mod1$coefficients[1]
regression <- sum(
  logistic.mod1$coefficients[-1] * x[names(which.max(log.pred)),][-ncol(x)]) +
  logistic.mod1$coefficients[1]
1/(1+exp(1)^-(regression))


# find optimal prediction cutoff to maximise f1 score
cutoff <- seq(0.01, 0.99, 0.01)
f1.list <- c()
for (i in 1:length(cutoff)){
  # print(cutoff[i])
  cat.pred <- ifelse(log.pred < cutoff[i], 'bacterial', 'viral')
  table(ytest, cat.pred)
  # sensitivity / recall
  tpr <- table(ytest, cat.pred)[1] / sum(table(ytest, cat.pred)[1] + table(ytest, cat.pred)[3])
  # specificity
  tnr <- table(ytest, cat.pred)[4] / sum(table(ytest, cat.pred)[4] + table(ytest, cat.pred)[2])
  # precision
  ppp <- table(ytest, cat.pred)[1] / sum(table(ytest, cat.pred)[1] + table(ytest, cat.pred)[2])
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