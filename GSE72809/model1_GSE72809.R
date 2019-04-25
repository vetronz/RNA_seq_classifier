library(Biobase)
library(GEOquery)
library(glmnet)
library(dplyr)
library(ROCR)
library(dplyr)
library(ggplot2)


## GET DATA FROM GEO
# gset <- getGEO('GSE72809')
# head(gset)
# saveRDS(gset, file = 'gset_GSE72809')


## LOAD DATA FROM rds object
setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
getwd()
gset <-readRDS(file = "gset_GSE72809")
head(gset)
n <- ncol(gset[[1]])

gset.df <- exprs(gset[[1]])
gset.df[1:5,1:5]

colnames(phenoData(gset[[1]]))
phenoData(gset[[1]])$'category:ch1'


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
gset.binary <- gset.df[,index.binary]
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

# creates a two level factor from the ground truth vector
label <- factor(truth, levels = c(1,2), labels = c('bacterial', 'viral'))
attributes(label)


# scale
gset.s <- as.data.frame(scale(gset.binary.t))
# apply(gset.s, 2, mean) # checks the means which should be around 0
# apply(gset.s, 2, sd) # checks the sd which should be around 1

dim(gset.s)
gset.s$label <- label
dim(gset.s)

gset.s[1:5, 47320:47324] # head
gset.s[139:144, 47320:47324] # tail


## DEFINE TRAINING AND TEST SETS
custom.rows <- c(1,2,47323, 47324)
# custom.rows <- c(47300:47324)
x <- gset.s[, custom.rows]
dim(x)
y <- gset.s[,ncol(gset.s)]

# training test split
set.seed(3)
n <- nrow(x)
index <- seq(1:n)
train = sample(1:n, round(n*0.8))
test = index[-train]
intersect(train, test)

x_train <- x[train,]
dim(x_train)
x_test <- x[test,][-ncol(x)] # strip the labels
dim(x_test)
ytest = y[test]

# Basic scatter plot
# x_test <- x[test,]
# p <- ggplot(x_train, aes(x=ILMN_1343291, y=ILMN_3311190))
# p <- ggplot(x_train, aes(x=ILMN_3311190, y=ILMN_1343295))
# p + geom_point(aes(colour = factor(label)), size = 2)

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
regression <- sum(logistic.mod1$coefficients[-1] * x[names(which.min(log.pred)),][-ncol(x)]) + logistic.mod1$coefficients[1]
regression <- sum(logistic.mod1$coefficients[-1] * x[names(which.max(log.pred)),][-ncol(x)]) + logistic.mod1$coefficients[1]
1/(1+exp(1)^-(regression))


# find optimal prediction cutoff to maximise f1 score
cutoff <- seq(0.2, 0.9, 0.01)
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



###  REGULARIZED APPROACHES

# need data matrix and y ground truth vector
x <- data.matrix(gset.s[,-ncol(gset.s)])
y


gset.s[1:5,(ncol(gset.s)-4):ncol(gset.s)]

gset.mat <- data.matrix(gset.s[-ncol(gset.s)])
dim(gset.mat)

lambda <- 10^seq(10, -2, length = 100)
ridge.mod <- glmnet(x[train,], y[train], family = 'binomial', standardize = FALSE)


plot(ridge.mod, xvar = "lambda")


ridge.mod$lambda %>% head()

# coefficients for the largest and smallest lambda parameters
coef(ridge.mod)[c('ILMN_3311175','ILMN_3311180'),1:2]

dim(coef(ridge.mod))

coef(ridge.mod)[c(1,2,3),1:5]
coef(ridge.mod)[c('ILMN_1343291','ILMN_1343295'),1:5]

# highest and lowest coefficients for lambda values of 1 and 100
coef(ridge.mod)[1,1]
coef(ridge.mod)[1,100]

# [c("Gr_Liv_Area", "TotRms_AbvGrd"), 100]
##   Gr_Liv_Area TotRms_AbvGrd 
##  0.0001004011  0.0096383231



# dev.off(dev.list()["RStudioGD"])




