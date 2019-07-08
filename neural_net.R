library(neuralnet)

scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}


X <- scale01(iris[-ncol(iris)])
y <- iris$Species

X.train <- cbind(X,y)


set.seed(500)
# 10 fold cross validation
k <- 10
# Results from cv
setosa <- NULL
versicolor <- NULL
virginica <- NULL

# Train test split proportions
proportion <- 0.70
# Set to 0.995 for LOOCV
# not able to do this with iris dataset because of 2 class requirement.
# reducing the proportion in the test set means we get some iters where there are no seratosa
# this breaks the roc calculation

# Crossvalidate
for(i in 1:k){
  index <- sample(nrow(X.train), round(proportion*nrow(X.train)))
  train_cv <- X.train[index, ]
  test_cv <- X.train[-index, ]

  nn_cv <- neuralnet(y~ ., train_cv, linear.output = FALSE, act.fct = "logistic", hidden = c(10, 3))
  pred <- predict(nn_cv, test_cv[-ncol(test_cv)])
  pred
  setosa[i] <- prediction(pred[,1], test_cv$y == 'setosa') %>%
    performance(measure = "auc") %>%
    .@y.values
  versicolor[i] <- prediction(pred[,2], test_cv$y == 'versicolor') %>%
    performance(measure = "auc") %>%
    .@y.values
  virginica[i] <- prediction(pred[,3], test_cv$y == 'virginica') %>%
    performance(measure = "auc") %>%
    .@y.values
}

mean(unlist(setosa))
mean(unlist(versicolor))
mean(unlist(virginica))
