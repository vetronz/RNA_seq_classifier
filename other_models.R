# Colours

# DB
#F8A9A4

# PB
#F5C3B8

# U
#E1D5E7

# PV
#DAE8FC

# V
#7E9ABF

# HC
#D5E8D4


# discovery
# validation


# train
#9A0794

# test
#E4F132


###
# messing around with the DEG

results <- decideTests(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, method = 'global')
dim(results)
head(results)
summary(results)
# colnames(results) <- c('DV-DB', 'DV-PB')
vennDiagram(results, include = 'both') # total sig 462

# check the ven diagram maths
sum(ifelse(results[,1]==0,0,1))
# 248+126 = 374
sum(ifelse(results[,2]==0,0,1))
# 126+88 = 214

sum(ifelse(results[,1]==0,0,1)) + sum(ifelse(results[,2]==0,0,1)) - 126

ifelse(results[,1]==0,0,1)
rownames(results[ifelse(results[,2]==0,0,1),])


results.mat <- as.matrix(results)
results.df <- as.data.frame(results.mat)
colnames(results.df) <- c('v.b', 'v.pb')

v.b.sig <- rownames(results.df[results.df$v.b == 1 | results.df$v.b == -1,])
v.pb.sig <- rownames(results.df[results.df$v.pb == 1 | results.df$v.pb == -1,])

length(union(v.b.sig, v.pb.sig))
length(v.b.sig) + length(v.pb.sig) - length(intersect(v.b.sig, v.pb.sig))


decide.hits <- union(v.b.sig, v.pb.sig)


# pull all genes passing number=nrow(fit2)
toptable.df <- topTable(fit2, adjust.method = 'BH', number=nrow(fit2),
                        p.value = pval, lfc=lfc)

toptable.df <- topTable(fit2, adjust.method = 'BH', number=nrow(fit2),
                        coef = c(1,2))

# not equal to the 462 genes id above
dim(toptable.df)


toptable.hits <- rownames(toptable.df)

# check the intersection to ensure not a completely different set
length(intersect(toptable.hits, decide.hits)) # all decide hits in top table hits
length(intersect(toptable.hits, decide.hits)) / length(decide.hits)

# it appears 96% of the genes ID using decide test are also picked up using filtered top table

# locate the 445 COMMON hits in toptable.hits
match(intersect(decide.hits, toptable.hits), toptable.hits)

# select them from toptable.df
filtered.toptable.df <- toptable.df[match(intersect(decide.hits, toptable.hits), toptable.hits),]
dim(filtered.toptable.df)

filtered.toptable.df[1:5,]

# we set max.lfc to the vrl.bct value 
# we then overwrite this if abs gb.v > abs b.v
filtered.toptable.df$max.lfc <- filtered.toptable.df$vrl.bct
for(i in 1:nrow(filtered.toptable.df)){
  v.b <- filtered.toptable.df$vrl.bct[i]
  v.pb <- filtered.toptable.df$vrl.greyb[i]
  if(abs(v.pb) > abs(v.b)){
    filtered.toptable.df$max.lfc[i] <- v.pb
  }
}






# OPTIMIZATION
boot <- 32
h.n <- 15
roc.a <- NULL
roc.t <- NULL
j.train <- NULL
j.test <- NULL
h.n.hx <- NULL
roc.train <- NULL
roc.train.me <- NULL
roc.test <- NULL
roc.test.me <- NULL

for(i in 1:h.n){
  for (k in 1:boot) {
    # for (k in 1:n_folds) {
    # print(paste0('hidden_nodes: ', i, ', fold: ', k))
    
    # test.i <- which(folds.i == k)
    # train.cv <- X.s[-test.i, ]
    # test.cv <- X.s[test.i, ]
    
    print(paste0('hidden Nodes: ', i, ', bootstrap: ', k))
    index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
    train.cv <- X.s[index,]
    test.cv <- X.s[-index,]
    # dim(train.cv)
    # dim(test.cv)
    nn1 <- neuralnet(bct~ ., train.cv, linear.output = FALSE, act.fct = "logistic",
                     hidden = c(i), rep = 3, stepmax = 1e+06, startweights = NULL, err.fct = "sse")
    
    pred.train <- predict(nn1, train.cv[-ncol(train.cv)])
    pred.test <- predict(nn1, test.cv[-ncol(test.cv)])
    
    j.train[k] <- prediction(pred.train[,1], train.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
    j.test[k] <- prediction(pred.test[,1], test.cv$bct) %>%
      performance(measure = "auc") %>%
      .@y.values
    
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
  Sys.sleep(0.1)
}

mod.complexity.df <- as.data.frame(cbind(roc.train, roc.test, roc.train.me, roc.test.me, h.n.hx))
colnames(mod.complexity.df) <- c('train', 'test', 'train.me', 'test.me', 'h.n')

pd <- position_dodge(0.2)
ggplot(mod.complexity.df, aes(x=h.n, y=train, color='train')) +
  scale_y_continuous(limits = c(0.80,1))+
  geom_line(aes(y=train))+
  geom_errorbar(aes(ymin=train-train.me, ymax=train+train.me), width=.4, position=pd)+
  geom_line(aes(y=test, color='test'))+
  geom_errorbar(aes(ymin=test-test.me, ymax=test+test.me), width=0.4, color='red')+
  labs(title=paste0('Bias-Variance Trade Off - Non Pseudo Labeled'), x =paste0('1 - ', h.n, ' hidden nodes'), y = "ROCA")

which.max(mod.complexity.df$test)
mod.complexity.df[which.max(mod.complexity.df$test),]

mod.complexity.df %>%
  dplyr::arrange(desc(test))%>%
  head(10)

opt.h.n <- which.max(mod.complexity.df$test)

###### RANDOM FORREST OPTIMIZATION ########
model <- randomForest(bct ~ . , data = train.fac, nodesize = 1, maxnodes = NULL)
plot(model)
attributes(model)

model$mtry
# [1] 12

model$err.rate[which.min(model$err.rate[,1])]

hyper_grid <- NULL
hyper_grid <- expand.grid(
  mtry       = seq(7, 14, by = 1),
  node_size  = seq(1, 9, by = 1),
  sampe_size = c(.54, .57, .6, .632, .66, .69, .72)
)

head(hyper_grid)
dim(hyper_grid)

for(i in 1:nrow(hyper_grid)) {
  print(i)
  # train model
  model <- randomForest(
    formula = bct ~ ., 
    data = train.fac,
    mtry = hyper_grid$mtry[i],
    nodesize = hyper_grid$node_size[i],
    samplesize = hyper_grid$sampe_size[i]
    # seed = 123
  )
  # extract error
  hyper_grid$OOB_MSE[i] <- model$err.rate[which.min(model$err.rate[,1])]
}

# plot_ly(hyper_grid, x = ~mtry, y = ~node_size, z = ~OOB_MSE, color = ~OOB_MSE)

hyper_top <- hyper_grid %>% 
  dplyr::arrange(OOB_MSE) %>%
  head(10)
hyper_top
mtry.opt <- hyper_top[1,][[1]]
node.opt <- hyper_top[1,][[2]]
sample.opt <- hyper_top[1,][[3]]


### RF BENCHMARKING
# basic rf, training
rf.mod <- randomForest(formula = bct ~ ., data = train.fac)

table(status.idx.d$most_general=='bacterial', rf.mod$predicted)
f1.score(table(status.idx.d$most_general=='bacterial', rf.mod$predicted))

prob <- rf.mod$votes
pr <- prediction(prob[,2], train.fac$bct)
pr %>%
  performance(measure = "auc") %>%
  .@y.values

# opt rf, training
rf.opt <- randomForest(
  formula = bct ~ .,
  data = train.fac,
  mtry = mtry.opt,
  nodesize = node.opt,
  samplesize = sample.opt
)
table(status.idx.d$most_general=='bacterial', rf.opt$predicted)
f1.score(table(status.idx.d$most_general=='bacterial', rf.opt$predicted))

prob <- rf.opt$votes
pr <- prediction(prob[,2], train.fac$bct)
pr %>%
  performance(measure = "auc") %>%
  .@y.values


# basic rf, val
pred <- predict(rf.mod, X.s.val)
# table(status.i.idx.d$most_general == 'bacterial', pred)
print(paste0('basic rf val F1 Score: ', round(f1.score(table(status.i.idx.d$most_general == 'bacterial', pred)),3)))

prob <- predict(rf.mod, X.s.val, type = 'prob')
pr <- prediction(prob[,2], status.i.idx.d$most_general=='bacterial')

a <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('basic rf, val: ', round(unlist(a),3)))

# opt rf, val
pred.opt.val <- predict(rf.opt, X.s.val)
table(status.i.idx.d$most_general == 'bacterial', pred.opt.val)
print(paste0('basic rf val F1 Score: ', round(f1.score(table(status.i.idx.d$most_general == 'bacterial', pred.opt.val)),3)))

prob.opt.val <- predict(rf.opt, X.s.val, type = 'prob')
pr <- prediction(prob.opt.val[,2], status.i.idx.d$most_general=='bacterial')
a <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('basic rf, val: ', round(unlist(a),3)))


### CROSS VAL OPTIMIZED RANDOM FORREST
j.test <- NULL
for (k in 1:n_folds) {
  print(paste0('fold: ', k))
  
  test.i <- which(folds.i == k)
  train.cv <- X.s[-test.i, ]
  test.cv <- X.s[test.i, ]
  
  train.cv$bct <- as.factor(train.cv$bct)
  test.cv$bct <- as.factor(test.cv$bct)
  
  model <- randomForest(
    formula = bct ~ ., 
    data = train.cv, 
    mtry = mtry.opt,
    nodesize = node.opt,
    samplesize = sample.opt
  )
  # extract error
  pred.train <- predict(model, train.cv[-ncol(train.cv)], type='prob')
  pr <- prediction(train.cv[,2], train.cv$bct)
  j.test[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
  
  pred.test <- predict(model, test.cv[-ncol(test.cv)], type='prob')
  pr <- prediction(pred.test[,2], test.cv$bct)
  j.test[k] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}

mean(unlist(j.test))
median(unlist(j.test))



####### rf.pseud labeling #######
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)

X.s$bct <- status.idx.d$most_general == 'bacterial'
train.fac <- X.s
train.fac$bct <- as.factor(train.fac$bct)
print(paste0('bacterial cases: ', sum(train.fac$bct==TRUE)))

mse <- NULL
for(i in 1:42){
  train.fac <- X.s
  train.fac$bct <- as.factor(train.fac$bct)
  
  sum(status.idx.d$most_general == 'probable_bacterial')
  
  rf.mod <- randomForest(
    formula = bct ~ ., 
    data = train.fac,
    mtry = mtry.opt,
    nodesize = node.opt,
    samplesize = sample.opt
  )
  
  pred <- predict(rf.mod, train.fac, type = 'prob')
  pr <- prediction(pred[,2], status.idx.d$most_general == 'bacterial')
  # pr %>%
  #   performance(measure = "auc") %>%
  #   .@y.values
  mse[i] <- rf.mod$err.rate[which.min(rf.mod$err.rate[,1])]
  
  sum(as.character(status.idx.d$my_category_2) == rownames(pred))
  
  # status.idx.d$most_general == 'probable_bacterial'
  # pred[status.idx.d$most_general == 'probable_bacterial', 2]
  
  # names(which.max(pred[status.idx.d$most_general == 'probable_bacterial', 2]))
  ppb <- which(status.idx.d$my_category_2 == names(which.max(pred[status.idx.d$most_general == 'probable_bacterial', 2])))
  status.idx.d$most_general[ppb] <- 'bacterial'
  X.s$bct <- status.idx.d$most_general == 'bacterial'
  print(paste0('iteration: ', i, ', bacterial cases: ', sum(X.s$bct)))
  Sys.sleep(1)
}

mse.df <- as.data.frame(mse)
mse.df$pb <- seq(1:nrow(mse.df))

ggplot(mse.df, aes(pb, mse))+geom_point()+
  labs(title='Iterative Pseudo-labeling Error', x='Mean Squared Error', y='MSE')+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE)

a <- lm(formula = mse ~ splines::bs(pb, 3), data=mse.df)
ppb.opt <- which.min(a$fitted.values)

x <- seq(1:42)
y <- a$fitted.values
ycs <- cumsum(y)
ycs.prime <- diff(ycs)/diff(x)
plot(ycs.prime)
which.min(ycs.prime)

ppb.opt # confirmed which was the lowest fitted value


# # numerical derivative estimation
# x = seq(2, 5, 0.5)
# y = x^(2)-2
# plot(x, y)
# ycs = cumsum(y)
# ycs.prime <- diff(ycs)/diff(x)
# which.min(ycs.prime)
# plot(ycs.prime)
# ycs.prime[which.min(ycs.prime)]


# ADD THE OPTIMAL NUMBER OF CASES
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)

X.s$bct <- status.idx.d$most_general == 'bacterial'
train.fac <- X.s
train.fac$bct <- as.factor(train.fac$bct)
print(paste0('bacterial cases: ', sum(train.fac$bct==TRUE)))

mse <- NULL
for(i in 1:ppb.opt[[1]]){
  rf.mod <- randomForest(
    formula = bct ~ ., 
    data = train.fac, 
    mtry = mtry.opt,
    nodesize = node.opt,
    samplesize = sample.opt
  )
  
  pred <- predict(rf.mod, train.fac, type = 'prob')
  pr <- prediction(pred[,2], status.idx.d$most_general == 'bacterial')
  mse[i] <- rf.mod$err.rate[which.min(rf.mod$err.rate[,1])]
  
  ppb <- which(status.idx.d$my_category_2 == names(which.max(pred[status.idx.d$most_general == 'probable_bacterial', 2])))
  status.idx.d$most_general[ppb] <- 'bacterial'
  X.s$bct <- status.idx.d$most_general == 'bacterial'
  train.fac <- X.s
  train.fac$bct <- as.factor(train.fac$bct)
  print(paste0('iteration: ', i, ', bacterial cases: ', sum(train.fac$bct==TRUE)))
  Sys.sleep(1)
}

print(paste0('bacterial cases: ', sum(train.fac$bct==TRUE)))


# RF OPTIMIZATION WITH PSEUDOLABALED DATA
model <- randomForest(bct ~ . , data = train.fac, nodesize = 1, maxnodes = NULL)
plot(model)

model$mtry
model$err.rate[which.min(model$err.rate[,1])]

hyper_grid <- NULL
hyper_grid <- expand.grid(
  mtry       = seq(7, 14, by = 1),
  node_size  = seq(1, 9, by = 1),
  sampe_size = c(.54, .57, .6, .632, .66, .69, .72)
)

head(hyper_grid)
dim(hyper_grid)

for(i in 1:nrow(hyper_grid)) {
  print(i)
  # train model
  model <- randomForest(
    formula = bct ~ ., 
    data = train.fac,
    mtry = hyper_grid$mtry[i],
    nodesize = hyper_grid$node_size[i],
    samplesize = hyper_grid$sampe_size[i]
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

rf.opt.psd <- randomForest(
  formula = bct ~ ., 
  data = train.fac, 
  mtry = mtry.opt,
  nodesize = node.opt,
  samplesize = sample.opt
)


### RF validation performance
# rf opt, val
pred.opt.val <- predict(rf.opt, X.s.val)
table(status.i.idx.d$most_general == 'bacterial', pred.opt.val)
print(paste0('Opt RF Validation F1 Score: ', round(f1.score(table(status.i.idx.d$most_general == 'bacterial', pred.opt.val)),3)))

prob.opt.val <- predict(rf.opt, X.s.val, type = 'prob')
pr <- prediction(prob.opt.val[,2], status.i.idx.d$most_general=='bacterial')
a <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('Opt RF Validation ROCA: ', round(unlist(a),3)))

# rf opt val psd
pred.opt.val.psd <- predict(rf.opt.psd, X.s.val)
# table(status.i.idx.d$most_general == 'bacterial', pred.opt.val.psd)
print(paste0('Pseudo Opt RF Validation F1 Score: ',
             round(f1.score(table(status.i.idx.d$most_general == 'bacterial', pred.opt.val.psd)),3)))

prob.opt.val.psd <- predict(rf.opt.psd, X.s.val, type = 'prob')
pr <- prediction(prob.opt.val.psd[,2], status.i.idx.d$most_general=='bacterial')
b <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('Pseudo Opt RF Validation ROCA: ', round(unlist(b),3)))

pb.dist.df <- as.data.frame(cbind(prob.opt.val[status.i.idx.d$most_general=='probable_bacterial'],
                                  prob.opt.val.psd[status.i.idx.d$most_general=='probable_bacterial']))
colnames(pb.dist.df) <- c('normal', 'pseudo-labeled')
hist(pb.dist.df$normal)
hist(pb.dist.df$`pseudo-labeled`)
pb.dist.df <- melt(pb.dist.df)
colnames(pb.dist.df)[1] <- 'model'
ggplot(pb.dist.df, aes(value, color=model, fill=model))+geom_density(alpha=0.2)+
  labs(title = 'Probability Density that Probable Bacterial Cases Represent Genuine Bacterial Infection', 
       x='Probability of Bacterial', y='Density')


# filtering out the non b and v cases
dim(X.s.val)
dim(status.i.idx.d)
bv.filt <- status.i.idx.d$most_general == 'bacterial' | status.i.idx.d$most_general == 'viral'
# status.i.idx.d$most_general == 'probable_viral' | status.i.idx.d$most_general == 'unknown' |
# status.i.idx.d$most_general == 'probable_bacterial'

sum(bv.filt)
X.s.val.bv <- X.s.val[bv.filt,]
pred.opt.val <- predict(rf.opt, X.s.val.bv)
table(status.i.idx.d$most_general[bv.filt] == 'bacterial', pred.opt.val)
print(paste0('basic rf val F1 Score: ', round(f1.score(table(status.i.idx.d$most_general[bv.filt], pred.opt.val)),3)))

prob.opt.val <- predict(rf.opt, X.s.val.bv, type = 'prob')
pr <- prediction(prob.opt.val[,2], status.i.idx.d$most_general[bv.filt]=='bacterial')
a <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('basic rf, val: ', round(unlist(a),3)))

# rf opt val psd
pred.opt.val.psd <- predict(rf.opt.psd, X.s.val.bv)
table(status.i.idx.d$most_general[bv.filt] == 'bacterial', pred.opt.val.psd)
print(paste0('basic rf val F1 Score: ', round(f1.score(table(status.i.idx.d$most_general[bv.filt] == 'bacterial', pred.opt.val.psd)),3)))

prob.opt.val.psd <- predict(rf.opt.psd, X.s.val.bv, type = 'prob')
pr <- prediction(prob.opt.val.psd[,2], status.i.idx.d$most_general[bv.filt]=='bacterial')
b <- pr %>%
  performance(measure = "auc") %>%
  .@y.values
print(paste0('basic rf, val: ', round(unlist(b),3)))


# pseudolabeling process, similar ROCA for non psd vs psd (87-87%)
# however the pattern of this ROC is different with more true positives, fewer false neg
# but compensated for by more false positives

# pred.opt.val
# FALSE TRUE
# FALSE    94   13
# TRUE     10   13

# pred.opt.val.psd
# FALSE TRUE
# FALSE    82   25
# TRUE      5   18

# on removing the pb, pv, other cases leaving b and v in the test iris set we see ROCA 97% in both models



###### SVM ######
# remove pseudo labels
status.idx.d <- status.idx[status.idx$most_general != 'healthy_control', ]
table(status.idx.d$most_general)
X.s$bct <- status.idx.d$most_general == 'bacterial'
X.s$bct <- as.factor(X.s$bct)
print(paste0('bacterial cases: ', sum(X.s$bct==TRUE)))


# cross validation
# set.seed(42)
folds.i <- sample(rep(1:n_folds, length.out = n.train))

# svm.m <- NULL
# for (k in 1:n_folds) {
#   print(paste0('fold: ', k))
# 
#   test.i <- which(folds.i == k)
#   train.cv <- X.s[-folds.i]
#   test.cv <- X.s[folds.i]
# 
#   ### SVM
#   model <- svm(bct ~ . , train.cv, probability = TRUE)
#   pred <- predict(model, train.cv, probability = TRUE)
#   # dim(attr(pred, "probabilities"))
#   pr <- prediction(attr(pred, "probabilities")[,2], test.cv$bct)
#   # plot(prf)
#   svm.m[k] <- pr %>%
#     performance(measure = "auc") %>%
#     .@y.values
# }
# mean(unlist(svm.m))
# median(unlist(svm.m))
# sd(unlist(svm.m))


###### svm.psd ######
attributes(model)

roc.h.psd <- NULL
for (i in 1:42){
  model <- svm(bct ~ . , X.s, probability = TRUE)
  pred <- predict(model, X.s, probability = TRUE)
  pred <- attr(pred, "probabilities")
  dim(pred)
  
  ppb <- names(which.max(pred[status.idx.d$most_general=='probable_bacterial',2]))
  ppb
  status.idx.d$most_general[which(status.idx.d$my_category_2 == ppb)] <- 'bacterial'
  X.s$bct <- status.idx.d$most_general == 'bacterial'
  X.s$bct <- as.factor(X.s$bct)
  print(paste0('iteration: ', i, ', bacterial cases: ', sum(X.s$bct==TRUE)))
  
  index <- sample(nrow(X.s), round(prop1*nrow(X.s)))
  train.cv <- X.s[index, ]
  test.cv <- X.s[-index, ]
  
  model <- svm(bct ~ . , train.cv, probability = TRUE)
  
  pred <- predict(model, test.cv[-ncol(test.cv)], probability = TRUE)
  # dim(attr(pred, "probabilities"))
  pr <- prediction(attr(pred, "probabilities")[,2], test.cv$bct)
  # plot(prf)
  roc.h.psd[i] <- pr %>%
    performance(measure = "auc") %>%
    .@y.values
}


svm.psd.df <- as.data.frame(unlist(roc.h.psd))
colnames(svm.psd.df) <- 'roc.a'
svm.psd.df$pb <- seq(1:nrow(svm.psd.df))
svm.psd.df

ggplot(svm.psd.df, aes(pb, roc.a))+geom_point()+
  labs(title='Iterative Pseudo-labeling Error', x='PB', y='ROCA')+
  geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = TRUE)

svm.psd.df$roc.a
a <- lm(formula = roc.a ~ splines::bs(pb, 3), data=svm.psd.df)

x <- seq(1:42)
y <- a$fitted.values
ycs <- cumsum(y)
ycs.prime <- diff(ycs)/diff(x)
plot(ycs.prime)
which.max(ycs.prime)
ppb.opt <- which.max(a$fitted.values)[[1]]









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
