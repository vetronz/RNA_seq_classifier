library(tidyverse)
# library(plyr)
library(limma)
library(cluster)
library(factoextra)
library(ggplot2)
require(reshape) # for melt()
require(scales) # for percent
library(gridExtra)
library(dplyr)
library(plotly)

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

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]

dx <- status[idx,]$most_general
dx.clean <- status[idx,]$most_general[clean.idx]
status[idx,]$array.contemporary.CRP[clean.idx]
wbc <- status[idx,]$WBC[clean.idx]
crp <- as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))

cor(wbc, crp)

ggplot(status[idx,][clean.idx,], aes(x=WBC, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))) +
  geom_point(shape=1)+
  geom_smooth(method=lm)


# Diagnosis Breakdown
table(droplevels(status[idx,]$most_general))
table(droplevels(status[idx,]$more_general))

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', "#165bfc")
dx.cols.f <- c("#ed0404", "#fc5716",'#16fc31', '#16e1fc', '#165bfc', "#7a16fc", '#fc16f4')
sex.cols <- c('#fc1676', '#16acfc')

positions <- c('bacterial', 'greyb', 'greyv', 'viral')
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother')

ggplot(status[idx,], aes(most_general, fill=most_general, alpha=Sex)) +
  scale_alpha_manual(values=c(0.6, 1)) +
  labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=dx.cols)+
  geom_bar()


c$most_general == 'bacterial' & c$Sex == 'M'
sum(c$most_general == 'bacterial' & c$Sex == 'M')
sum(c$most_general == 'bacterial' & c$Sex == 'F')
sum(c$most_general == 'greyb' & c$Sex == 'M')
sum(c$most_general == 'greyb' & c$Sex == 'F')
sum(c$most_general == 'greyv' & c$Sex == 'M')
sum(c$most_general == 'greyv' & c$Sex == 'F')
sum(c$most_general == 'viral' & c$Sex == 'M')
sum(c$most_general == 'viral' & c$Sex == 'F')


sex <- c('M', 'F')
bacterial <- c(22, 30)
greyb <- c(24,18)
greyv <- c(4,1)
viral <- c(65,27)
df <- data.frame(bacterial, greyb, greyv, viral)
df.2 <- mutate(df, sex = factor(c('M','F')))
df.3 <- gather(df.2, dx, count, -sex)
df.3
ggplot(df.3, aes(x = dx, y = count, fill = sex)) + 
  geom_bar(position = "fill",stat = "identity")+
  scale_fill_manual(values=sex.cols)
  # geom_text(df.3, sex, labels=round(sex))


# ggplot(status[idx,], aes(most_general, group=Sex)) +
#   geom_bar(aes(y = ..prop.., fill=factor(..x..)), stat="count") +
#   geom_text(aes( label = scales::percent(..prop..),
#                  y= ..prop.. ), stat= "count", vjust = -.5) +
#   scale_fill_manual(values=dx.cols)+
#   facet_grid(~Sex) +
#   coord_flip()+
#   scale_y_continuous(labels = scales::percent)+
#   labs(title = "Barplot of Percentage Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Percentage")
# 
# ggplot(status[idx,], aes(more_general, group=Sex)) + 
#   geom_bar(aes(y = ..prop.., fill=factor(..x..)), stat="count") +
#   geom_text(aes( label = scales::percent(..prop..),
#                  y= ..prop.. ), stat= "count") +
#   scale_fill_manual(values=dx.cols.f)+
#   scale_x_discrete(limits=positions.f)+
#   facet_grid(~Sex) +
#   coord_flip()+
#   scale_y_continuous(labels = scales::percent)+
#   labs(title = "Barplot of Percentage Full Diagnostic Groups Breakdown by Gender", x = "Diagnosis", y = "Percentage")

table(status[idx,]$Sex)
chisq.test(table(status[idx,]$Sex))

table(droplevels(status[idx,]$most_general), status[idx,]$Sex)
# table(droplevels(status[idx,]$more_general), status[idx,]$Sex)

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

# ggplot(status[idx,], aes(status[idx,]$most_general, Age..months., color=Sex)) + geom_jitter(alpha = .9, width = 0.2)+
#   scale_x_discrete(limits = positions)+
#   xlab('Diagnostic Group') +
#   ylab('Age') +
#   scale_x_discrete(limits = positions)+
#   scale_color_manual(values=c('#16acfc', '#fc1676'))+
#   ggtitle("Jitter Plot of Age by Diagnostic Group, Split by Gender")


# ddply(status[idx,], ~most_general, summarise, mean=mean(Age..months.))
# 
# a <- as.data.frame(droplevels(status[idx,]$most_general))
# colnames(a) <- 'dx'
# a$age <-status[idx,]$Age..months.
# head(a)
# summary(a)
# 
# b <- ddply(a,~dx)
# 
# anova.res <- aov(age ~ dx, data = b)
# summary(anova.res)


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


# PCA
dim(e.set)
e.set[,idx]
dim(e.set[,idx])

X <- e.set[,idx]
dim(X)

e.set.t <- t(e.set)
dim(e.set.t)

e.set.f <- e.set.t[idx,]
dim(e.set.f)


full.pca <- prcomp(e.set.f, scale=TRUE)
pair1 <- as.data.frame(full.pca$x[,1:2])
pair2 <- as.data.frame(full.pca$x[,3:4])

fviz_eig(full.pca)

ve <- full.pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

status[idx,]$most_general
dim(pair1)

dim()
# DEF DX PC1 PC2
fviz_pca_ind(full.pca)

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general), size=2) +
# ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general[dx.def]), size=2) +  
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Diagnosis')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=dx.cols)+
  # scale_color_manual(values=dx.cols.2)+
  ggtitle("Bacterial Viral Split on PC 1 and PC2")

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex), size=2) +
# ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  labs(col='Gender')+
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=sex.cols)+
  ggtitle("Male Female Split against PC1 PC2")


library("RColorBrewer")
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
wbc.col <- scale_colour_gradientn(colours = myPalette(10), limits=c(min(wbc), max(wbc)))
crp.col <- scale_colour_gradientn(colours = myPalette(10), limits=c(min(crp), max(crp)))


ggplot(pair1[clean.idx,], aes(PC1, PC2, color = wbc, shape=dx.clean)) + geom_point(size=2.5)+
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  wbc.col +
  labs(col='WBC', shape='Diagnosis')+
  theme(panel.background = element_rect(fill = '#303030', colour = 'red'))+
  ggtitle("Bacterial Viral Split on PC 1 and PC2")

ggplot(pair1[clean.idx,], aes(PC1, PC2, color = crp, shape=dx.clean)) + geom_point(size=2.5)+
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  crp.col +
  labs(col='CRP', shape='Diagnosis')+
  theme(panel.background = element_rect(fill = '#303030', colour = 'red'))+
  ggtitle("Bacterial Viral Split on PC 1 and PC2")


# top left corner looks like an interesting subset of samples
# want to check if they reemerge following clustering
x.pos <- -50
y.pos <- 80

ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$most_general, shape=status[idx,]$Sex), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  # geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  # geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_hline(yintercept = y.pos, linetype="longdash", colour="red", size=.5) +
  geom_vline(xintercept = x.pos, linetype="longdash", colour="red", size=.5) +
  scale_color_manual(values=dx.cols)+
  ggtitle("Bacterial Viral Split on PC 1 and PC2")

pair1[pair1$PC1 < x.pos & pair1$PC2 > y.pos,]
status[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,]$Diagnosis
View(status[idx,][pair1$PC1 < x.pos & pair1$PC2 > y.pos,])

# plotly play
pal <- c("red", "blue")
pal <- setNames(pal, c("M", "F"))
plot_ly(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=status[idx,]$Sex[dx.def]), size=2)
plot_ly(data = pair1[dx.def,], x = ~PC1, y = ~PC2, color = status[idx,]$Sex[dx.def], colors = pal)

dim(e.set.f)
dim(X.t)

### Unsupervised Clustering
fviz_nbclust(X.t, kmeans, method = "wss")
fviz_nbclust(X.t, kmeans, method = "silhouette")

gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 10,
K.max = 20, B = 20)
fviz_gap_stat(gap_stat)

### K2
dim(e.set.f)

k2 <- kmeans(X, centers = 2, nstart = 10)

k2$cluster <- as.factor(k2$cluster)

table(k2$cluster, droplevels(status[idx,]$most_general))
table(k2$cluster, droplevels(status[idx,]$more_general))

k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
                                   'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k2.df$cluster <- k2$cluster[clean.idx]
dim(k2.df)
k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))


k2.clus.col <- c("#ed0404", "#165bfc")
positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother')
ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k2.clus.col)+
  geom_bar()
# table(k2$cluster, status[idx, c('Sex')])

dx.def <- status[idx,]$most_general == 'bacterial' | status[idx,]$most_general == 'viral'
dx.grey <- status[idx,]$most_general == 'greyb' | status[idx,]$most_general == 'greyv'
dx.def

dim(pair1[dx.def,])
length(status[idx,]$most_general[dx.def])

# DEF DX PC1 PC2
ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$most_general[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=k2.clus.col)+
  ggtitle("Cluster Assignment of Bacterial and Viral Against first two PC")















############ LIMMA ############
e.set.f <- as.data.frame(e.set.f)
e.set.f$label <- as.character(status[idx,c('most_general')])
e.set.f$sex <- status[idx,c('Sex')]
e.set.f$age <- status[idx,c('Age..months.')]

dim(e.set.f)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

# rm(e.set, e.set.i, e.set.t)

X.sub <- X[,]
dim(X.sub)
X.sub.m <- apply(X.sub, 1, mean)
X.sub.sd <- apply(X.sub, 1, sd)

plot(X.sub.m, X.sub.sd)

### DESIGN MATRIX
design <- model.matrix(~label + sex + age + 0, data = e.set.f)
colnames(design)<- c("bct","greyb","greyv", 'vrl', 'sexM', 'age')

design[1:10,]
dim(design)
colSums(design)

contrast.matrix<- makeContrasts("bct-vrl", 'bct-greyv', levels=design)
contrast.matrix
# colnames(fit$coefficients)

fit <- lmFit(X, design)

hist(fit$Amean)
plotSA(fit)
abline(v=5)

keep <- fit$Amean > 5
sum(keep)
fit2<- contrasts.fit(fit, contrast.matrix)
dim(fit2)
fit2 <- eBayes(fit2[keep,], trend = TRUE)
dim(fit2)
plotSA(fit2)

lfc <- 1.5
pval <- 0.05

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'bct-vrl')
dim(results)
head(results)
summary(results)
vennDiagram(results, include = 'both')
vennCounts(results, include = 'both')


results.hits <- union(rownames(X[keep,])[results[,1] == 1]
                      ,rownames(X[keep,])[results[,1] == -1])
results.hits
length(results.hits)

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'bct-vrl', number = 77)
all.hits <- topTable(fit2, number=nrow(fit2), coef = 'bct-vrl')
top.hits
dim(top.hits)
dim(all.hits)

intersect(results.hits, rownames(top.hits))


ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  ggtitle("Volcano Plot of Log Fold Change Against -log10 P Value
    Cutoff - Fold Change:1, P Val:0.05")





# dim(results)
# intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))
# rownames(results)[intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))]
# 
# results.hits
# 
# intersect(results.hits, rownames(results)[intersect(which(results[,1] == 1 | results[,1] == -1), which(results[,2] == 1 | results[,2] == -1))])

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
### PCA
X.t <- t(X)
dim(X.t[,results.hits])
X.pca <- prcomp(X.t[,results.hits], scale = TRUE)

# summary(e.set.pca)

pair1 <- as.data.frame(X.pca$x[,1:2])
pair2 <- as.data.frame(X.pca$x[,3:4])

fviz_eig(X.pca)

ve <- X.pca$sdev^2
pve <- ve/sum(ve)*100
pve

fviz_nbclust(X.t[,results.hits], kmeans, method = "wss")
fviz_nbclust(X.t[,results.hits], kmeans, method = "silhouette")

# gap_stat <- clusGap(X.t[,results.hits], FUN = kmeans, nstart = 25,
                    # K.max = 20, B = 25)
# fviz_gap_stat(gap_stat)

### K2
k2 <- kmeans(X.t[,results.hits], centers = 2, nstart = 500)
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
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k2.clus.col)+
  geom_bar()
# table(k2$cluster, status[idx, c('Sex')])

dx.def <- status[idx,]$most_general == 'bacterial' | status[idx,]$most_general == 'viral'
dx.grey <- status[idx,]$most_general == 'greyb' | status[idx,]$most_general == 'greyv'
dx.def

dim(pair1[dx.def,])
length(status[idx,]$most_general[dx.def])

# DEF DX PC1 PC2
ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$most_general[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=k2.clus.col)+
  ggtitle("Cluster Assignment of Bacterial and Viral Against first two PC")

# DEF DX PC3 PC4
ggplot(pair2[dx.def,], aes(PC3, PC4)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$most_general[dx.def]), size=2) +
  xlab(paste0("PC3: (", round(pve[3],2), '%)') ) +
  ylab(paste0("PC4: (", round(pve[4],2), '%)') ) +
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=k2.clus.col)+
  ggtitle("Cluster Assignment of Bacterial and Viral Against first two PC")

# M F PC1 PC2
ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$Sex[dx.def]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=k2.clus.col)+
  ggtitle("Cluster Assignment of Male, Female Against first two PC")

# Grey Dx PC1 PC2
ggplot(pair1[dx.grey,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.grey], shape=status[idx,]$most_general[dx.grey]), size=2) +
  xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
  ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
  geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
  scale_color_manual(values=k2.clus.col)+
  ggtitle("Cluster Assignment of Probable Bacterial and Probable Viral Against first two PC")




ggplot(status[idx,], aes(most_general, Age..months., fill=k2$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k2.clus.col)+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")

p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k2$cluster[clean.idx])) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  xlab('') +
  scale_fill_manual(values=k2.clus.col)+
  ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k2$cluster[clean.idx]))+
  geom_boxplot()+
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_fill_manual(values=k2.clus.col)
gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)


filter.bct <- status[idx,]$most_general == 'bacterial'
filter.clus <- k2$cluster == 1
filter.comb <- filter.bct & filter.clus
status[idx,]$Diagnosis[filter.comb]

filter.clus <- k2$cluster == 2
filter.comb <- filter.bct & filter.clus
status[idx,]$Diagnosis[filter.comb]
# View(status[idx,][filter.comb, c('my_category_2', 'most_general',
#                                  'more_general', 'Age..months.', 'Sex', 'WBC',
#                                  'array.contemporary.CRP', 'Diagnosis')])

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
k3 <- kmeans(X.t[,results.hits], centers = 3, nstart = 25)
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
k5 <- kmeans(X.t[,results.hits], centers = 5, nstart = 25)
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
# k7 <- kmeans(X.t[,results.hits], centers = 7, nstart = 27)
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