library(Biobase)
library(GEOquery)
library(ggplot2)
library(tidyr)
# library(rsample)
library(caret)
# library(h2o)
library(dplyr)
library(tidyverse)  # data manipulation and visualization
library(modelr)     # provides easy pipeline modeling functions
# library(broom)
library(glmnet)   # implementing regularized regression approaches

## GET DATA FROM GEO
gset <- getGEO('GSE72809')
head(gset)
saveRDS(gset, file = 'gset_GSE72809')

## LOAD DATA FROM rds object
# setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# getwd()
# gset <-readRDS(file = "gset_GSE72809")
# head(gset)
# exprs(gset_t[[1]])
n <- ncol(gset[[1]])


gset_df <- exprs(gset[[1]])
gset_df[1:5,1:5]

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


# bct_df <- gset_df[, bct]
# vrl_df <- gset_df[, vrl]
# bct_df <- gset.sig[, bct]
# vrl_df <- gset.sig[, vrl]
# dim(bct_df)
# dim(vrl_df)

## SUBSET THE GSET TO GET BACT VIRAL DF
index_binary <- c(bct.index, vrl.index)
gset_binary <- gset_df[,index_binary]
gset_binary_t <- as.data.frame(t(gset_binary))

# gset_binary[1:5,1:5]
# gset_binary_t[1:5,1:5]

dim(gset_binary)
dim(gset_binary_t)

## CONSTRUCTS GROUND TRUTH LIST 
truth <- character(0)
index_binary <- list(bct.index, vrl.index)
for (i in 1:length(index_binary)){
  for (j in 1:length(index_binary[[i]])){
    print(paste(i, j))
    truth = append(truth,i)
  }
}
truth

# yyz$b <- as.numeric(as.character(yyz$b))

truth_num <- as.numeric(truth)-1
label <- as.factor(truth_num)

gset_binary_t$label <- label
gset_binary_t[1:5, 47320:47324]

dim(gset_binary_t)


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



# Create training and testing feature model matrices and response vectors.
# we use model.matrix(...)[, -1] to discard the intercept
gset_train_x <- model.matrix(label ~ ., train_1)[,-1]
gset_train_y <- train_1$label

gset_test_x <- model.matrix(label ~ ., test_1)[,-1]
gset_test_y <- test_1$label

# What is the dimension of of your feature matrix?
dim(gset_train_x)


gset_ridge <- glmnet(
  x = gset_train_x,
  y = gset_train_y,
  alpha = 0
)

plot(gset_ridge, xvar = "lambda")












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

## LOAD DATA FROM rds object
# saveRDS(train_x, file = 'train_x_GSE72809')
# saveRDS(test_x, file = 'test_x_GSE72809')
train_x <- readRDS(file = 'train_x_GSE72809')
test_x <- readRDS(file = 'test_x_GSE72809')

train_x$truth <- truth

model1 <- glm(truth ~ ILMN_1657871, family = "binomial", data = train_1)
summary(model1)

pred <- predict(model1, test_1, type='response')

pred
pred.class <- ifelse(pred > 0.5, 1, 0)
pred.class

test_1


