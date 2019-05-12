library(tidyverse)
library(plyr)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(dplyr)


getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

dim(e.set)
e.set.t <- t(e.set)
dim(e.set.t)


idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'
  # status['most_general'] == 'OD'
sum(idx)

# e.set.f <- e.set[,idx]
# status.f <- status[idx,]
# dim(e.set.f)
# dim(status.f)

attributes(status)$names

which(as.character(status[idx,]$array.contemporary.CRP) == 'na')
length(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'))

which(status[idx,]$WBC == 0)
clean<-union(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'), which(status[idx,]$WBC == 0))
clean
length(clean)

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]
clean.idx

status[idx,]$array.contemporary.CRP[clean.idx]
as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))

cor(status[idx,]$WBC[clean.idx], as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))


ggplot(status[idx,][clean.idx,], aes(x=WBC, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))) +
  geom_point(shape=1)+
  geom_smooth(method=lm)


# Diagnosis Breakdown
table(droplevels(status[idx,]$most_general))
table(droplevels(status[idx,]$more_general))

dx.cols <- c("#ed0404", "#fc5716", '#16fc31', "#165bfc")
dx.cols.f <- c("#ed0404", "#fc5716",'#16fc31', '#16e1fc', '#165bfc', "#7a16fc", '#fc16f4')

positions <- c('bacterial', 'greyb', 'greyv', 'viral')
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother')

ggplot(status[idx,], aes(most_general, fill=most_general, alpha=Sex)) +
  scale_alpha_manual(values=c(0.6, 1)) +
  labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols)+
  geom_bar()

ggplot(status[idx,], aes(most_general, group=Sex)) + 
  geom_bar(aes(y = ..prop.., fill=factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  scale_fill_manual(values=dx.cols)+
  facet_grid(~Sex) +
  coord_flip()+
  scale_y_continuous(labels = scales::percent)+
  labs(title = "Barplot of Percentage Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Percentage")

ggplot(status[idx,], aes(more_general, group=Sex)) + 
  geom_bar(aes(y = ..prop.., fill=factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count") +
  scale_fill_manual(values=dx.cols.f)+
  scale_x_discrete(limits=positions.f)+
  facet_grid(~Sex) +
  coord_flip()+
  scale_y_continuous(labels = scales::percent)+
  labs(title = "Barplot of Percentage Full Diagnostic Groups Breakdown by Gender", x = "Diagnosis", y = "Percentage")

table(status[idx,]$Sex)
chisq.test(table(status[idx,]$Sex))

table(droplevels(status[idx,]$most_general), status[idx,]$Sex)
table(droplevels(status[idx,]$more_general), status[idx,]$Sex)

chisq.test(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))
# chisq.test(table(droplevels(status[idx,]$more_general), status[idx,]$Sex))


# Age Breakdown
ggplot(status[idx,], aes(status[idx,]$most_general, Age..months., fill=Sex)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c('#16acfc', '#fc1676'))+
  ggtitle("Box and Whisker Plot of Age by Diagnostic Group, Split by Gender")

ggplot(status[idx,], aes(status[idx,]$most_general, Age..months., color=Sex)) + geom_jitter(alpha = .9, width = 0.2)+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_x_discrete(limits = positions)+
  scale_color_manual(values=c('#16acfc', '#fc1676'))+
  ggtitle("Jitter Plot of Age by Diagnostic Group, Split by Gender")


ddply(status[idx,], ~most_general, summarise, mean=mean(Age..months.))

a <- as.data.frame(droplevels(status[idx,]$most_general))
colnames(a) <- 'dx'
a$age <-status[idx,]$Age..months.
head(a)
summary(a)

b <- ddply(a,~dx)

anova.res <- aov(age ~ dx, data = b)
summary(anova.res)


# Inflamatory Marker Breakdown
p1 <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=Sex)) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c('#16acfc', '#fc1676'))+
  ggtitle("Box and Whisker Plot of WBC Count by Diagnosis, Split by Gender")

p2 <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,][clean.idx,]$array.contemporary.CRP)), fill=Sex))+
  geom_boxplot()+
  # scale_y_continuous(limits = c(0, 230))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c('#16acfc', '#fc1676'))+
  ggtitle("Box and Whisker Plot of CRP Count by Diagnosis, Split by Gender")
gridExtra::grid.arrange(p1, p2, nrow = 2)

crp <- as.numeric(as.character(status[idx,][clean.idx,]$array.contemporary.CRP))
crp[status[idx,][clean.idx,]$Sex=='M' & status[idx,][clean.idx,]$most_general == 'bacterial']
crp[status[idx,][clean.idx,]$Sex=='F' & status[idx,][clean.idx,]$most_general == 'bacterial']

t.test(crp[status[idx,][clean.idx,]$Sex=='M' & status[idx,][clean.idx,]$most_general == 'bacterial'],
       crp[status[idx,][clean.idx,]$Sex=='F' & status[idx,][clean.idx,]$most_general == 'bacterial'])
# p-value = 0.6789
# no stat sig difference in crp between males and females in bct group



# 
# # WBC
# ggplot(status[idx,], aes(Sex, WBC)) + geom_boxplot() +
#   xlab('Cluster Assignment') +
#   ylab('WBC Count') +
#   ggtitle("Box Plot of WBC Count by Cluster Assignment")
# 
# ggplot(status[idx,], aes(WBC, fill=Sex)) + geom_histogram(bins = 20)
#   xlab('WBC') +
#   ylab('Density') +
#   ggtitle("Histogram of WBC Count by Cluster Assignment")
# 
# ggplot(status[idx,], aes(WBC, fill=Sex)) + geom_density(alpha=.5)+
#   xlab('WBC') +
#   ylab('Density') +
#   ggtitle("Density Plot of WBC Count by Cluster Assignment")
# 
# ps <- (seq(0,499) + 0.5)/500
# qs <- quantile(status[idx,]$WBC, ps)
# normalqs <- qnorm(ps, mean(status[idx,]$WBC), sd(status[idx,]$WBC))
# plot(normalqs,qs,xlab="Normal Percentiles",ylab="WBC Percentiles")
# abline(0,1)
# 
# # bin
# q <- quantile(status[idx,]$WBC, seq(0,1,0.2))
# wbc_bin <- cut(status[idx,]$WBC, breaks=q, include.lowest = T)
# 
# tab <- with(status[idx,], table(wbc_bin, Sex))
# print(prop.table(tab, 2))
# chisq.test(tab)
# 
# df1 <- melt(ddply(status[idx,],.(Sex),function(x){prop.table(table(wbc_bin))}),id.vars = 1)
# 
# ggplot(df1, aes(x = variable,y = value)) +
#   facet_wrap(~Sex, nrow=2, ncol=1) +
#   scale_y_continuous(labels=percent) +
#   geom_bar(stat = "identity") +
#   xlab('Binned WBC Levels') +
#   ylab('Percentage Bin Assignment') +
#   ggtitle("% Assignment to Binned WBC Levels by Cluster")
# 
# 
# # CRP
# ggplot(status[idx,], aes(Sex, as.numeric(array.contemporary.CRP))) + geom_boxplot() +
#   xlab('Gender') +
#   ylab('CRP') +
#   ggtitle("Box Plot of CRP Count by Gender")
# 
# ggplot(status[idx,], aes(as.numeric(array.contemporary.CRP), fill=Sex)) + geom_histogram(bins = 10)
# xlab('CRP') +
#   ylab('Density') +
#   ggtitle("Histogram of CRP by Gender")
# 
# ggplot(status[idx,], aes(as.numeric(array.contemporary.CRP), fill=Sex)) + geom_density(alpha=.5)+
#   xlab('CRP') +
#   ylab('Density') +
#   ggtitle("Density Plot of CRP by Gender")
# 
# 
# ps <- (seq(0,499) + 0.5)/500
# qs <- quantile(as.numeric(status[idx,]$array.contemporary.CRP), ps)
# normalqs <- qnorm(ps, mean(as.numeric(status[idx,]$array.contemporary.CRP)), sd(as.numeric(status[idx,]$array.contemporary.CRP)))
# plot(normalqs,qs,xlab="Normal Percentiles",ylab="CRP Percentiles")
# abline(0,1)
# 
# # bin
# q <- quantile(as.numeric(status[idx,]$array.contemporary.CRP), seq(0,1,0.25))
# crp_bin <- cut(as.numeric(status[idx,]$array.contemporary.CRP), breaks=q, include.lowest = T)
# 
# tab <- with(status[idx,], table(crp_bin, Sex))
# print(prop.table(tab, 2))
# chisq.test(tab)
# 
# df1 <- melt(ddply(status[idx,],.(Sex),function(x){prop.table(table(crp_bin))}),id.vars = 1)
# df2 <- melt(ddply(status[idx,],.(WBC),function(x){prop.table(table(wbc_bin))}),id.vars = 1)
# 
# ggplot(df1, aes(x = variable,y = value, fill=Sex)) +
#   facet_wrap(~Sex, nrow=2, ncol=1) +
#   scale_y_continuous(labels=percent) +
#   geom_bar(stat = "identity") +
#   xlab('Binned CRP Levels') +
#   ylab('Percentage Bin Assignment') +
#   ggtitle("% Assignment to Binned CRP Levels by Cluster")



###### differential gene expression with limma ######
dim(e.set)
e.set.t <- t(e.set)
dim(e.set.t)


e.set.f <- e.set.t[idx,]
dim(e.set.f)
# label.f <- status[idx,c('most_general')]
# label <- as.character(label.f)
X <- as.matrix(t(e.set.f))
dim(X)

e.set.f <- as.data.frame(e.set.f)
e.set.f$label <- as.character(status[idx,c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]
e.set.f$age <- status[idx,c('Age..months.')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

rm(e.set, e.set.i, e.set.t)

### DESIGN MATRIX
design <- model.matrix(~label + sex + age + 0, data = e.set.f)
colnames(design)<- c("bct","greyb","greyv", 'vrl', 'sexM', 'age')

design[1:10,]
dim(design)
colSums(design)

contrast.matrix<- makeContrasts("bct-vrl", levels=design)
contrast.matrix
# colnames(fit$coefficients)

dim(X)
dim(design)

fit <- lmFit(X, design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

dim(fit2$coefficients)
fit2$coefficients[1:10,]

top.hits <- topTable(fit2, number=100, coef = 'bct-vrl', lfc = 1.5, p.value = 0.05,  adjust="BH")
dim(top.hits)
X.t[1:5, rownames(top.hits)]

results <- decideTests(fit2, method='global', p.value = 0.05, adjust.method = 'BH', lfc=1.5)
head(results)
summary(results)

vennDiagram(results, include = 'both')
# heatDiagram(results, include = 'both', coef = 'bct-vrl')

dim(results)
intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))
rownames(results)[intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))]

rownames(top.hits)

intersect(rownames(top.hits), rownames(results)[intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))])

# 
# # proof that in above we are not taking the intersect of genes that are over expressed on one
# # and under expressed in the other
# a<-which(results[,1] == 1)
# b<-which(results[,2] == 1)
# c<-which(results[,1] == -1)
# d<-which(results[,2] == -1)
# 
# rownames(results)[union(intersect(a,b), intersect(c,d))]
# 
# # still get 66 taking the intersect of the two above methods
# intersect(union(intersect(a,b), intersect(c,d)),
#           intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1)))
# 
# 
# rownames(results)
############################## CLUSTERING ##############################
rm(X, fit, fit2)

### PCA
X.t <- t(X)
dim(X.t[,rownames(top.hits)])
X.pca <- prcomp(X.t[,rownames(top.hits)], scale = TRUE)
# summary(e.set.pca)
plot(X.pca, type = 'l')

pair1 <- as.data.frame(X.pca$x[,1:2])
pair2 <- as.data.frame(X.pca$x[,3:4])


fviz_nbclust(X.t[,rownames(top.hits)], kmeans, method = "wss")
fviz_nbclust(X.t[,rownames(top.hits)], kmeans, method = "silhouette")

gap_stat <- clusGap(X.t[,rownames(top.hits)], FUN = kmeans, nstart = 25,
                    K.max = 20, B = 25)
fviz_gap_stat(gap_stat)

### K2
k2 <- kmeans(X.t[,rownames(top.hits)], centers = 2, nstart = 25)
str(k2)
k2$cluster <- as.factor(k2$cluster)


k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster[clean.idx]
dim(k2.df)
k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))

addmargins(table(k2$cluster, droplevels(status[idx,]$most_general)))
addmargins(table(k2$cluster, droplevels(status[idx,]$more_general)))

k2.clus.col <- c("#ed0404", "#165bfc")
positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother')
ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k2.clus.col)+
  geom_bar()
# table(k2$cluster, status[idx, c('Sex')])

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k2$cluster, shape=status[idx,]$most_general), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  scale_color_manual(values=clus.col)+
  ggtitle("Cluster Assignment of First Two Principal Components")

ggplot(pair2, aes(PC3, PC4)) + geom_point(aes(color=k2$cluster, shape=status[idx,]$most_general), size=2) +
  xlab("Third Principal Component") +
  ylab("Fourth Principal Component") +
  scale_color_manual(values=clus.col)+
  ggtitle("Cluster Assignment of Third-Fourth Principal Components")

ggplot(status[idx,], aes(most_general, Age..months., fill=k2$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=clus.col)+
  ggtitle("Age (months) by Diagnostic Group, Split by Gender")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k2$cluster[clean.idx])) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_fill_manual(values=clus.col)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k2$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  scale_fill_manual(values=clus.col)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")

gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)


filter.bct <- status[idx,]$most_general == 'bacterial'
# filter.bct <- status[idx,]$most_general == 'greyb'
filter.clus <- k2$cluster == 2

filter.comb <- filter.bct & filter.clus


sum(filter.comb)
# View(status[idx,])
status[idx,]$Diagnosis[filter.comb]
View(status[idx,][filter.comb, c('my_category_2', 'most_general',
                                 'more_general', 'Age..months.', 'Sex', 'WBC',
                                 'array.contemporary.CRP', 'Diagnosis')])

detach('package:plyr', unload = TRUE, character.only = TRUE)
k2.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, most_general) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

k2.df$WBC[k2.df$cluster == 1 & k2.df$most_general == 'bacterial']
k2.df$WBC[k2.df$cluster == 2 & k2.df$most_general == 'bacterial']
t.test(k2.df$WBC[k2.df$cluster == 1 & k2.df$most_general == 'bacterial'], k2.df$WBC[k2.df$cluster == 2 & k2.df$most_general == 'bacterial'])
# not enough power to detect sig diff within bct group given few bct assigned to class 2
# no difference with crp

k2.df$WBC[k2.df$cluster == 1]
k2.df$WBC[k2.df$cluster == 2]
t.test(k2.df$WBC[k2.df$cluster == 1], k2.df$WBC[k2.df$cluster == 2])
# t = 4.628, df = 138.86, p-value = 8.385e-06

k2.df$array.contemporary.CRP[k2.df$cluster == 1]
k2.df$array.contemporary.CRP[k2.df$cluster == 2]
t.test(k2.df$array.contemporary.CRP[k2.df$cluster == 1]
       ,k2.df$array.contemporary.CRP[k2.df$cluster == 2])
# t = 8.1468, df = 127.89, p-value = 2.92e-13




### K3
k3 <- kmeans(X.t[,rownames(top.hits)], centers = 3, nstart = 25)
# str(k3)
k3$cluster <- as.factor(k3$cluster)

addmargins(table(k3$cluster, droplevels(status[idx,]$most_general)))
addmargins(table(k3$cluster, droplevels(status[idx,]$more_general)))

k3.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k3.df$cluster <- k3$cluster[clean.idx]
dim(k3.df)
k3.df[1:8,(ncol(k2.df)-3): ncol(k2.df)]
k3.df$array.contemporary.CRP <- as.numeric(as.character(k3.df$array.contemporary.CRP))

# had to un attach plyr to get working
k3.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, Sex) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))

k3.clus.col <- c("#fcac16","#ed0404", "#165bfc")
ggplot(status[idx,], aes(more_general, fill=k3$cluster)) +
  labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k3.clus.col)+
  geom_bar()

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k3$cluster, shape=status[idx,]$most_general), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  scale_color_manual(values=k3.clus.col)+
  ggtitle("Cluster Assignment of First Two Principal Components")

ggplot(status[idx,], aes(k3$cluster, Age..months., fill=Sex)) + geom_boxplot()+
  xlab('Diagnostic Group') +
  ylab('Age') +
  ggtitle("Age (months) by Diagnostic Group, Split by Gender")

k3.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k3$cluster[clean.idx])) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k3.clus.col)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")

k3.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k3$cluster[clean.idx]))+
  geom_boxplot()+
  # scale_y_continuous(limits = c(0, 230))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_fill_manual(values=k3.clus.col)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")
gridExtra::grid.arrange(k3.wbc, k3.crp, nrow = 2)

k3.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, color=k3$cluster[clean.idx])) + geom_jitter(alpha = .9, width = 0.2) +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  scale_color_manual(values=k3.clus.col)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")

k3.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
                                               as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k3$cluster[clean.idx]))+
  geom_jitter(alpha = .9, width = 0.2)+
  # scale_y_continuous(limits = c(0, 230))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_color_manual(values=k3.clus.col)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")
gridExtra::grid.arrange(k3.wbc, k3.crp, nrow = 2)


filter.bct <- status[idx,]$most_general == 'bacterial'
# filter.bct <- status[idx,]$most_general == 'greyb'
filter.clus <- k3$cluster == 2
filter.clus <- k3$cluster == 1
filter.clus <- k3$cluster == 3
filter.clus <- k3$cluster == 1 | k3$cluster == 3

filter.comb <- filter.bct & filter.clus

filter.comb
# View(status[idx,])
status[idx,]$Diagnosis[filter.comb]
View(status[idx,][filter.comb, c('my_category_2', 'most_general',
                                 'more_general', 'Age..months.', 'Sex', 'WBC',
                                 'array.contemporary.CRP', 'Diagnosis')])




### K5
k5 <- kmeans(X.t[,rownames(top.hits)], centers = 5, nstart = 25)
k5$cluster <- as.factor(k5$cluster)

addmargins(table(k5$cluster, droplevels(status[idx,]$most_general)))
addmargins(table(k5$cluster, droplevels(status[idx,]$more_general)))

k5.clus.col <- c("#fcac16","#165bfc", '#16fcd2', '#fc16af', "#ed0404")
ggplot(status[idx,], aes(more_general, fill=k5$cluster)) +
  labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k5.clus.col)+
  geom_bar()

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k5$cluster, shape=status[idx,]$most_general), size=2) +
  xlab("First Principal Component") +
  ylab("Second Principal Component") +
  scale_color_manual(values=k5.clus.col)+
  ggtitle("Cluster Assignment of First Two Principal Components")

k5.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k5$cluster[clean.idx])) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k5.clus.col)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")
k5.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k5$cluster[clean.idx]))+
  geom_boxplot()+
  # scale_y_continuous(limits = c(0, 250))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_fill_manual(values=k5.clus.col)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")
gridExtra::grid.arrange(k5.wbc, k5.crp, nrow = 2)

k5.wbc <- ggplot(status[idx,][clean.idx,], aes(more_general, WBC, color=k5$cluster[clean.idx])) +
  geom_jitter(alpha = .9, width = 0.2) +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  # scale_x_discrete(limits = positions)+
  scale_color_manual(values=k5.clus.col)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")
k5.crp <- ggplot(status[idx,][clean.idx,], aes(more_general,
  as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k5$cluster[clean.idx]))+
  geom_jitter(alpha = .9, width = 0.2)+
  # scale_y_continuous(limits = c(0, 250))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_color_manual(values=k5.clus.col)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")
gridExtra::grid.arrange(k5.wbc, k5.crp, nrow = 2)


# 
# ### K7
# k7 <- kmeans(X.t[,rownames(top.hits)], centers = 7, nstart = 27)
# k7$cluster <- as.factor(k7$cluster)
# 
# addmargins(table(k7$cluster, droplevels(status[idx,]$most_general)))
# addmargins(table(k7$cluster, droplevels(status[idx,]$more_general)))
# 
# k7.clus.col <- c("#fcac16","#165bfc", '#16fcd2', '#fc16af', "#ed0404", '#f8fc16', '#16fc1a')
# ggplot(status[idx,], aes(more_general, fill=k7$cluster)) +
#   labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
#   scale_x_discrete(limits = positions.more)+
#   scale_fill_manual(values=k7.clus.col)+
#   geom_bar()
# 
# ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k7$cluster, shape=status[idx,]$most_general), size=2) +
#   xlab("First Principal Component") +
#   ylab("Second Principal Component") +
#   scale_color_manual(values=k7.clus.col)+
#   ggtitle("Cluster Assignment of First Two Principal Components")
# 
# k7.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k7$cluster[clean.idx])) + geom_boxplot() +
#   # scale_y_continuous(limits = c(0, 70))+
#   ylab('WBC Count') +
#   scale_x_discrete(limits = positions)+
#   scale_fill_manual(values=k7.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# k7.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
#                                                as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k7$cluster[clean.idx]))+
#   geom_boxplot()+
#   # scale_y_continuous(limits = c(0, 270))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_fill_manual(values=k7.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k7.wbc, k7.crp, nrow = 2)
# 
# k7.wbc <- ggplot(status[idx,][clean.idx,], aes(more_general, WBC, color=k7$cluster[clean.idx])) +
#   geom_jitter(alpha = .9, width = 0.2) +
#   # scale_y_continuous(limits = c(0, 70))+
#   ylab('WBC Count') +
#   # scale_x_discrete(limits = positions)+
#   scale_color_manual(values=k7.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# k7.crp <- ggplot(status[idx,][clean.idx,], aes(more_general,
#                                                as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k7$cluster[clean.idx]))+
#   geom_jitter(alpha = .9, width = 0.2)+
#   # scale_y_continuous(limits = c(0, 250))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_color_manual(values=k7.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k7.wbc, k7.crp, nrow = 2)








# end