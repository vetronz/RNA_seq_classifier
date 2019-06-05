

############# limma differential expression gram +ve gram -ve #############
dim(X.t)

X.g <- as.data.frame(X.t[status[idx,]$most_general == 'bacterial',])
dim(X.g)
class(X.g)
a <- status[idx,]$category == 'E' | status[idx,]$category == 'F'

sum(status[idx,]$category == 'F')

b <- status[idx,]$most_general == 'bacterial'
sum(a == b)

X.g$gram <- droplevels(status[idx,][a,]$category)


### DESIGN MATRIX
# site
design <- model.matrix(~gram + 0, data = X.g)

colnames(design)<- c('gram.pos', 'gram.neg')

design[1:10,]
dim(design)
colSums(design)

contrast.matrix<- makeContrasts("gram.pos-gram.neg", levels=design)
contrast.matrix
# colnames(fit$coefficients)


dim(t(X.g[-ncol(X.g)]))
dim(design)

fit <- lmFit(t(X.g[-ncol(X.g)]), design)

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

lfc <- 0.5
pval <- 0.1

results <- decideTests(fit2, method='global', p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'gram.pos-gram.neg')
dim(results)
head(results)
summary(results)
vennDiagram(results, include = 'both')
# vennCounts(results, include = 'both')

dim(X.g[-nrow(X.g)])
colnames(X.g[-ncol(X.g)])

length(colnames(X.g[-ncol(X.g)])[keep])
dim(results)
colnames(X.g[-ncol(X.g)])[keep][results == 1]
colnames(X.g[-ncol(X.g)])[keep][results == -1]
results.tot <- union(colnames(X.g[-ncol(X.g)])[keep][results == 1], colnames(X.g[-ncol(X.g)])[keep][results == -1])

top.hits <- topTable(fit2, p.value = pval, adjust.method = 'BH', lfc=lfc, coef = 'gram.pos-gram.neg', number = 19)
head(top.hits)
all.hits <- topTable(fit2, number=nrow(fit2), coef = 'gram.pos-gram.neg')
dim(top.hits)
dim(all.hits)

intersect(results.tot, rownames(top.hits))

ggplot(all.hits, aes(y=-log10(adj.P.Val), x=logFC)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(pval), linetype="longdash", colour="grey", size=1) +
  geom_vline(xintercept = lfc, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -(lfc), linetype="longdash", colour="#2C467A", size=1)+
  ggtitle("Volcano Plot of Log Fold Change Against -log10 P Value
          Cutoff - Fold Change:1, P Val:0.05")

dim(X.t)
dim(X.t[,results.tot])
colnames(X.t[,results.tot])
gram.hits<-colnames(X.t[,results.tot])
gram.hits
# View(X.t[,results.tot])
X.t[,results.tot][a,]



####### hierarchecal clustering #######
# Dissimilarity matrix
d <- dist(t(X.t[,results.tot][a,]), method = "euclidean")
d <- dist(X.t[,results.tot][a,], method = "euclidean")

dim(t(X.t[,results.tot][a,]))

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(d, method = "average" )
# hc1 <- hclust(d, method = "ward.D" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.7, hang = -1)
rect.hclust(hc1, k = 2, border = c(2,4))
rect.hclust(hc1, k = 3, border = c(2,4))

# outlier
dim(X.t[,results.tot][a,][rownames(X.t[,results.tot][a,]) != 'bacterialgpos_19_SMH',])
e.set.g <- X.t[,results.tot][a,][rownames(X.t[,results.tot][a,]) != 'bacterialgpos_19_SMH',]
status.g <- status[idx,][a,][status[idx,][a,]$my_category_2 != 'bacterialgpos_19_SMH',]
dim(e.set.g)
dim(status.g)


fviz_nbclust(e.set.g, FUN = hcut, method = "wss")
fviz_nbclust(e.set.g, FUN = hcut, method = "silhouette")
gap_stat <- clusGap(e.set.g, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

# d <- dist(t(e.set.g), method = "euclidean")
d <- dist(e.set.g, method = "euclidean")
hc1 <- hclust(d, method = "average" )
# hc1 <- hclust(d, method = "ward.D" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.7, hang = -1)
rect.hclust(hc1, k = 2, border = c(2,4))
rect.hclust(hc1, k = 3, border = c(2,4))

# Cut tree into 2 groups
sub_grp <- cutree(hc1, k = 3)
sub_grp
table(droplevels(status.g$category), sub_grp)
chisq.test(table(droplevels(status.g$category), sub_grp))

fviz_cluster(list(data = e.set.g, cluster = sub_grp))

# Hierarchical clustering using Complete Linkage
d1 <- as.dist(1-cor(t(e.set.g), method="pearson"))
d2 <- as.dist(1-cor(e.set.g, method="spearman"))


hc2 <- hclust(d1, method = "average" )
hc3 <- hclust(d2, method = "average" )

# Plot the obtained dendrogram
plot(hc2, cex = 0.7, hang = -1)
rect.hclust(hc2, k = 2, border = c(2,4))
rect.hclust(hc2, k = 3, border = c(2,4))

# Plot the obtained dendrogram
plot(hc3, cex = 0.7, hang = -1)
rect.hclust(hc3, k = 2, border = c(2,4))
rect.hclust(hc3, k = 3, border = c(2,4))

clustRows <- hclust(d1, method="average")
clustColumns <- hclust(d2, method="average")

# module.assign <- cutree(clustColumns, k=3)
# module.assign == 2
# module.assign[module.assign == 2]

module.assign <- cutree(clustRows, k=2)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- rainbow(length(unique(module.assign)), start=0.1, end=0.9) 
module.color <- module.color[as.vector(module.assign)] 
myheatcolors2 <- colorRampPalette(colors=c("blue", "yellow","red"))(100)
# produce a static heatmap of DEGs ----
#plot the hclust results as a heatmap
heatmap.2(e.set.g,
          Rowv=as.dendrogram(clustRows), 
          Colv=NA,
          RowSideColors=module.color,
          col=myheatcolors2, scale='row', labRow=NA,
          density.info="none", trace="none",  
          cexRow=1, cexCol=1, margins=c(8,20))


heatmap.2(e.set.g,
          Rowv=as.dendrogram(clustRows), 
          Colv=as.dendrogram(clustColumns),
          col=redgreen(100),
          scale="row",
          margins = c(7, 7),
          cexCol = 0.7,
          labRow = F,
          main = "Heatmap.2",
          trace = "none")


library(heatmaply) #for making interactive heatmaps using plotly
p<-heatmaply(e.set.g,
             colors = myheatcolors2,
             Rowv=as.dendrogram(clustRows),
             RowSideColors=module.color,
             # showticklabels=c(FALSE,FALSE),
             scale='row')
p
api_create(p, filename = "heatmap_micro")

p<-heatmaply(mtcars, k_col = 2, k_row = 3) %>% layout(margin = list(l = 130, b = 40))
api_create(p, filename = "heatmap_test")

library("illuminaHumanv4.db")

getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)

head(illumina)
nrow(illumina)

module.assign[module.assign == 2]
# 5720482 2570300 2100196  990768 3360343
trans <- c(3180392, 2570300, 2100196, 5720482, 7650358, 4780075, 7650433)
# trans <- c(7650358, 5720482, 2570300, 3180392, 5090754)
# trans <- c(3180392)
# trans <- c(2100196)

which(illumina$Array_Address_Id %in% trans)

probeID <- illumina$Probe_Id[which(illumina$Array_Address_Id %in% trans)]
probeID

x <- illuminaHumanv4CHR
x <- illuminaHumanv4NUID
x <- illuminaHumanv4ALIAS2PROBE
x <- illuminaHumanv4ENSEMBL
x <- illuminaHumanv4GENENAME
x <- illuminaHumanv4GO
x <- illuminaHumanv4MAP
x <- illuminaHumanv4REFSEQ
data.frame(trans, Gene=unlist(mget(x = probeID, envir = x)))




# fuzzy clustering allocation
k4.df %>%
  select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
  group_by(cluster, most_general) %>%
  summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))


which(k4.df$most_general == 'bacterial' & k4.df$cluster == 2)
which(k4.df$most_general == 'bacterial' & k4.df$cluster == 1)

View(k4.df[which(k4.df$most_general == 'bacterial' & k4.df$cluster == 2),])
View(k4.df[which(k4.df$most_general == 'bacterial' & k4.df$cluster == 1),])

k4$membership[which(status[idx,]$most_general == 'greyb'),]

which.max(k4$membership[which(status[idx,]$most_general == 'greyb'),][,2])

which(status[idx,]$most_general == 'greyb')[24]

status[idx,]$my_category_2[131]
pair3D[131,]
k4$cluster

fcm_centroids <- k4$centers
cor(t(fcm_centroids))


### K7
k7.pal <- c('#09f70d', '#2909f7',"#f70d09", '#5dddbb', '#fc37f5', '#ccf246', "#f76409")
set.seed(18)
k7 <- kmeans(X.t, centers = 7, nstart = 100)
k7$cluster <- as.factor(k7$cluster)

table(k7$cluster, droplevels(status[idx,]$most_general))
table(k7$cluster, droplevels(status[idx,]$more_general))

p<-ggplot(status[idx,], aes(most_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
p<-ggplot(status[idx,], aes(more_general, fill=k7$cluster)) +
  labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions.more)+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  geom_bar()
ggplotly(p)

# 2D PCA Dx Group
plot_ly(pair1, type="scatter", x = ~PC1, y = ~PC2, mode = "markers", color = ~k7$cluster, size = status[idx,]$Age..months.,
        colors=k7.pal, text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  layout(title = 'Cluster Assignment on PC 1-2', xaxis=x, yaxis=y)

# 3D PCA Dx Group
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster, size = status[idx,]$Age..months.,
        colors = c(k7.pal), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Cluster Assignment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

plot_ly(pair3D[clean.idx,], x = ~PC1, y = ~PC2, z = ~PC3, color = ~k7$cluster[clean.idx],
        colors = c(k7.pal), size = ~crp, text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2[clean.idx], '<br>Diagnosis: ',status[idx,]$Diagnosis[clean.idx])) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, CRP Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))


p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k7$cluster[clean.idx])) + geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('WBC Count') +
  xlab('') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k7.pal, 'Cluster')+
  ggtitle("Box Plot of Cluster WBC and CRP Count by Diagnosis")
p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k7$cluster[clean.idx]))+
  geom_boxplot(position=position_dodge(width=0.8)) +
  ylab('CRP Count') +
  xlab('Diagnosis') +
  scale_x_discrete(limits = positions[-5])+
  scale_fill_manual(values=k7.pal, 'Cluster')
p <- gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)

ggplot(status[idx,], aes(most_general, Age..months., fill=k7$cluster)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_fill_manual(values=k7.pal, 'Cluster')+
  ggtitle("Age (months) by Diagnostic Group, Split by Cluster")


# K7 analysis
k7.df <- status[idx,][, c('category', 'my_category_2', 'most_general', 'more_general', 'site',
                          'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
k7.df$cluster <- k4$cluster
dim(k7.df)
View(k7.df)

dim(k7.df[k7.df$most_general == 'bacterial',])
k7.df[k7.df$most_general == 'bacterial',]

ggplot(k7.df[k7.df$most_general == 'bacterial',], aes(cluster, fill=category)) +
  labs(title = "Barplot of Microbiology Results by Cluster", x = "Diagnosis", y = "Counts")+
  scale_fill_manual(values=c('#6c5ddd', '#dd5dbb'), 'Micro')+
  geom_bar()







###########################################################################################
### prediction
dim(X.t)

scale01 <- function(x){
  (x - min(x)) / (max(x) - min(x))
}

X.s<-data.frame(apply(X.t, 2, scale01))
dim(X.s)

X.s$bacterial <- status[idx,]$most_general == 'bacterial'
dim(X.s)

# train test split
train <- sample(seq(1, dim(X.s)[1]), dim(X.s)[1]*0.7, replace = FALSE)
test <- seq(1, dim(X.s)[1])[-train]

train.x <- X.s[train,]
test.x <- X.s[test,][-ncol(X.s)]
y <- ifelse(status[idx,]$most_general[test] == 'bacterial', TRUE, FALSE)
dim(train.x)
dim(test.x)
y

## fit model using `f`
# colnames(X.s)
# model <- neuralnet(bacterial ~ X1740360 + X450189, data = train.x, hidden=c(4,3), threshold=0.01)

# hi hi hi git
roc.a <- NULL
roc.t <- NULL
h.n <- 9
rep <- 10

for(j in 1:rep){
  for(i in 1 : h.n) {
    model <- neuralnet(bacterial ~ ., data = train.x, hidden=i, threshold=0.01)
    pred <- predict(model, test.x, type="class")
    roc.a[i] <- prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
      performance(measure = "auc") %>%
      .@y.values
  }
  roc.t <- append(roc.t, roc.a)
}


roc.t

df <- data.frame(matrix(unlist(roc.t), nrow=length(roc.t), byrow=T))
colnames(df) <- 'a'
df
as.data.frame(split(df, 1:h.n))
df <- as.data.frame(split(df, 1:h.n))
df.t <- as.data.frame(t(df))
df.t$h.n <- seq(1:nrow(df.t))

boxplot(df, use.cols=TRUE)

# # boxplot reshaping

# colnames(df.t) <- c('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6', 'rep7', 'rep8', 'rep9', 'rep10', 'rep11', 'rep12', 'h.n')
# df.t
# 
# dat.m <- melt(df.t,id.vars='h.n', measure.vars=c('rep1', 'rep2', 'rep3', 'rep4', 'rep5', 'rep6'))
# dat.m
# 
# ggplot(dat.m) + geom_boxplot(aes(x=h.n, y=value, color=variable))
# 
# # melt(df, )
# dat.m <- melt(dat,id.vars='ID', measure.vars=c('Freq','Freq.1','Freq.2'))
# 
# roc.m <- as.data.frame(apply(df, 2, mean))
# roc.m$sd <- apply(df, 2, sd)
# roc.m$hidden.nodes <- seq.int(nrow(roc.m))
# roc.m
# colnames(roc.m)[1] <- 'mean.ROC'
# 
# roc.m
# p<-ggplot(roc.m, aes(x=hidden.nodes, y=mean.ROC))+geom_jitter()
# p<-ggplot(roc.m, aes(x=hidden.nodes, y=mean.ROC))+geom_point()



model <- neuralnet(bacterial ~ ., data = train.x, hidden=3, threshold=0.01)
pred <- predict(model, test.x, type="class")

prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "tpr", x.measure = "fpr") %>%
  plot()

prediction(pred, status[idx,]$most_general[test] == 'bacterial') %>%
  performance(measure = "auc") %>%
  .@y.values

