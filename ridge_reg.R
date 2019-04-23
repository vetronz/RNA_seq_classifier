library(glmnet)   # implementing regularized regression approaches
library(dplyr)
library(ROCR)

swiss <- datasets::swiss
swiss
dim(swiss)

x <- model.matrix(Fertility~., swiss)[,-1]
y <- swiss$Fertility
y.cat <- y > 80

lambda <- 10^seq(10, -2, length = 100)

#create test and training sets
set.seed(489)
train = sample(1:nrow(x), round(nrow(x)*0.7))
test = (-train)
ytest = y[test]
y.cat.test = y.cat[test]

#OLS
swisslm <- lm(Fertility~., data = swiss)
coef(swisslm)

#ridge
ridge.mod <- glmnet(x, y, alpha = 0, lambda = lambda)
predict(ridge.mod, s = 0, type = 'coefficients')

swisslm <- lm(Fertility~., data = swiss, subset = train)
coef(swisslm)

ridge.mod <- glmnet(x[train,], y[train], alpha = 0, lambda = lambda)
#find the best lambda from our list via cross-validation
cv.out <- cv.glmnet(x[train,], y[train], alpha = 0)

bestlam <- cv.out$lambda.min

#make predictions
s.pred <- predict(swisslm, newdata = as.data.frame(x)[test,])
ridge.pred <- predict(ridge.mod, s = bestlam, newx = x[test,])



## MODEL EVALUATION
#check MSE
mean((s.pred-ytest)^2)
mean((ridge.pred-ytest)^2)

# classification
s.class <- s.pred > 80
r.class <- ridge.pred > 80

s.class == r.class

list(
  table(y.cat.test, s.class) %>% prop.table() %>% round(3),
  table(y.cat.test, r.class) %>% prop.table() %>% round(3)
)



par(mfrow=c(1, 2))

prediction(s.pred, y.cat.test) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

prediction(ridge.pred, y.cat.test) %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()
caret::varImp(swisslm)



# model 1 AUC
prediction(s.pred, y.cat.test) %>%
  performance(measure = "auc") %>%
  .@y.values

# model 1 AUC
prediction(ridge.pred, y.cat.test) %>%
  performance(measure = "auc") %>%
  .@y.values


# caret::varImp(ridge.mod)

# anova(swisslm, ridge.mod, test = "Chisq")


# lasso.mod <- glmnet(x[train,], y[train], alpha = 1, lambda = lambda)
# lasso.pred <- predict(lasso.mod, s = bestlam, newx = x[test,])
# mean((lasso.pred-ytest)^2)
# 
# lasso.coef  <- predict(lasso.mod, type = 'coefficients', s = bestlam)[1:6,]
# lasso.coef
