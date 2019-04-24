library(Biobase)
library(GEOquery)
library(glmnet)
library(dplyr)
library(ROCR)
library(dplyr)
library(ggplot2)
# library(tidyr)
# library(rsample)
# library(caret)
# library(h2o)
# library(tidyverse)  # data manipulation and visualization
# library(modelr)     # provides easy pipeline modeling functions
# library(broom)

## GET DATA FROM GEO
# gset <- getGEO('GSE72809')
# head(gset)
# saveRDS(gset, file = 'gset_GSE72809')

## LOAD DATA FROM rds object
setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
getwd()
gset <-readRDS(file = "gset_GSE72809")
head(gset)
# exprs(gset_t[[1]])
n <- ncol(gset[[1]])


gset.df <- exprs(gset[[1]])
gset.df[1:5,1:5]

colnames(phenoData(gset[[1]]))
phenoData(gset[[1]])$'category:ch1'
class(phenoData(gset[[1]])$'category:ch1')



## ADD DEF BACT AND DEF VIRAL TO LISTS
bct.index <- c()
vrl.index <- c()
for (i in 1:n){
  if (phenoData(gset[[1]])$'category:ch1'[i] == 'Definite Bacterial'){
    bct.index[[i]] <- i
  }
  if (phenoData(gset[[1]])$'category:ch1'[i] == 'Definite Viral'){
    vrl.index[[i]] <- i
  }
}

vrl.index <- vrl.index[!is.na(vrl.index)]



## SUBSET THE GSET TO GET BACT VIRAL DF
index.binary <- c(bct.index, vrl.index)
gset.binary <- gset_df[,index.binary]
gset.binary.t <- as.data.frame(t(gset.binary))

# dim(gset_binary)
dim(gset.binary.t)

## CONSTRUCTS GROUND TRUTH LIST
truth <- character(0)
index.binary <- list(bct.index, vrl.index)
index.binary
for (i in 1:length(index.binary)){
  for (j in 1:length(index.binary[[i]])){
    print(paste(i, j))
    truth = append(truth,i)
  }
}

label <- factor(truth, levels = c(1,2), labels = c('bacterial', 'viral'))
label

attributes(label)

gset.binary.t$label <- label
gset.binary.t[1:5, 47320:47324]
gset.binary.t[1:5, c(1,2,47323, 47324)]
dim(gset.binary.t)


## DEFINE TRAINING AND TEST SETS
x <- gset.binary.t[, c(1,2,47323, 47324)]
dim(x)
# x$label <- y
y <- gset.binary.t[,ncol(gset.binary.t)]

# n <- nrow(x)
# index <- seq(1:n)
# train = sample(1:n, round(n*0.8))
# test = index[-train]
# 
# # intersect(train, test)
# # integer(0)
# 
x_train <- x[train,]
dim(x_train)
x_test <- x[test,][-ncol(x)]
dim(x_test)
# ytest = y[test]


# Basic scatter plot
# x_test <- x[test,]
p <- ggplot(x_train, aes(x=ILMN_1343291, y=ILMN_3311190))
# p <- ggplot(x_train, aes(x=ILMN_1343291, y=ILMN_1343295))
p + geom_point(aes(colour = factor(label)), size = 2)


# logistic
logistic.mod1 <- glm(label ~., family = "binomial", data = x_train)
summary(logistic.mod1)
# logistic.mod1

log.pred <- predict(logistic.mod1, x_test, type = 'response')
cat.pred <- ifelse(log.pred < 0.5, 'bact', 'viral')
cat.pred

# x_test['GSM1872559',]
# 
# x_test['GSM1872497',]
# 
# x_test['GSM1872495',]




x_train[1,]
logistic.mod1$coefficients
log.pred[1]

(0.5562*4.9273 + 5.898*-1.12124)+4.5931
a <- 0.556119*4.927301
b <- -1.121235 * 5.897797
c <- 4.593186

# classification
list(
  table(ytest, cat.pred) # %>% prop.table() %>% round(3)
  # table(y.cat.test, r.class) %>% prop.table() %>% round(3)
)


numeric.truth <- ifelse(truth == '1', 0, 1)
dev.off(dev.list()["RStudioGD"])

prediction(log.pred, ytest) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()



## HAND MADE PREDICTIONS
logistic.mod1$coefficients
x_test[1,]

logistic.mod1$coefficients[-1] * x_test[1,]
sum(logistic.mod1$coefficients[-1] * x_test[1,])


length(x_test[1,][-length(x_test)])
length(logistic.mod1$coefficients[-1])

x_test[1,][-length(x_test)] * logistic.mod1$coefficients[-1]
sum(x_test[1,][-length(x_test)] * logistic.mod1$coefficients[-1]) + logistic.mod1$coefficients[1]


# model 1 AUC
prediction(cat.pred, ytest) %>%
  performance(measure = "auc") %>%
  .@y.values








## TRAIN TEST SPLIT
set.seed(123)
index <- sample(1:nrow(gset_binary_t), round(nrow(gset_binary_t) * 0.7))
index
train_1 <- gset_binary_t[index, ]
test_1  <- gset_binary_t[-index, ]

dim(train_1)
dim(test_1)


bac_prop <- summary(gset_binary_t$label)[[1]]/nrow(gset_binary_t)
train_bac_prop <- summary(train_1$label)[[1]]/dim(train_1)[1]
test_bac_prop <- summary(test_1$label)[[1]]/dim(test_1)[1]

print(paste('bacterial proportion:', bac_prop))
print(paste('bacterial train proportion:', train_bac_prop))
print(paste('bacterial test proportion:', test_bac_prop))






## NORMALIZATION
# identify only the predictor variables
features <- setdiff(names(train_1), 'truth')

# pre-process estimation based on training features
pre_process <- preProcess(
  x      = train_1[, features],
  method = c("center", "scale")
)

# apply to both training & test
train_x <- predict(pre_process, train_1[, features])
test_x  <- predict(pre_process, test_1[, features])

train_x[1:5, 47320:47323]
test_x[1:5, 1:5]
# remove(train_x, test_x)