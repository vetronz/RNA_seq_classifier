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
library(plotly)

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

Sys.setenv("plotly_username"="vetronz1992")
Sys.setenv("plotly_api_key"="Wtx9CzYqbl9iC8EzXp2B")

dx.cols.2 <- c("#ed0404", "#165bfc")
dx.cols <- c("#ed0404", "#fc5716", '#16fc31', '#384ebc', "#bc38ab")
dx.cols.f <- c("#ed0404", "#fc5716", '#16fc31','#38bc9d', '#55f1fc', '#0934f4', '#384ebc', "#bc38ab")
sex.cols <- c('#fc1676', '#16acfc')


positions <- c('bacterial', 'greyb', 'greyv', 'viral', 'HC')
positions.f <- c('bacterial', 'greyb', 'greyv', 'flu', 'RSV', 'adeno', 'viralother', 'HC')


site.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc")
cat.pal <- c("#ed0404", "#fc5716", '#d7fc35', '#35c7fc', '#16fc31', '#464647', "#165bfc", '#16fc31', '#464647', "#165bfc")

dim(e.set)

idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'HC'
sum(idx)


which(as.character(status[idx,]$array.contemporary.CRP) == 'na')
length(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'))

which(status[idx,]$WBC == 0)
clean<-union(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'), which(status[idx,]$WBC == 0))

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]
wbc <- status[idx,]$WBC[clean.idx]
crp <- as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))
cor(wbc, crp)
ggplot(status[idx,][clean.idx,], aes(x=WBC, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))) +
  geom_point(shape=1)+
  geom_smooth(method=lm)




# Diagnosis Breakdown
table(droplevels(status[idx,]$most_general))
table(droplevels(status[idx,]$more_general))

dx <- c('bacterial', 'greyb', 'greyv', 'viral', 'HC')
counts <- c(52, 42, 5, 92, 62)
data <- data.frame(dx, counts)
p <- plot_ly(data, x = ~dx, y = ~counts, type = 'bar', marker = list(color = dx.cols)) %>%
  layout(title = 'Barplot of Diagnostic Group Breakdown',
         yaxis = list(title = 'Count'),
         xaxis = list(title = 'Figure 1.2'),
         barmode = 'group')
p
# api_create(p, filename = "barplot_dx_breakdown")

dx <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother', 'HC')
counts <- c(52, 42, 5, 23, 23, 27, 19, 62)
data <- data.frame(dx, counts)
p <- plot_ly(data, x = ~dx, y = ~counts, type = 'bar', marker = list(color = dx.cols.f)) %>%
  layout(title = 'Barplot of Diagnostic Group Breakdown',
         yaxis = list(title = 'Count'),
         xaxis = list(title = 'Diagnosis'),
         barmode = 'group')
p
# api_create(p, filename = "barplot_dx_breakdown_full")

# sex
# sex <- c('M', 'F')
# bacterial <- c(22, 30)
# greyb <- c(24,18)
# greyv <- c(4,1)
# viral <- c(65,27)
# HC <- c(33,29)
# df <- data.frame(bacterial, greyb, greyv, viral, HC)
# df.2 <- mutate(df, sex = factor(c('M','F')))
# df.3 <- gather(df.2, dx, count, -sex)
# df.3
# p<-ggplot(df.3, aes(x = dx, y = count, fill = sex)) + 
#   geom_bar(position = "fill",stat = "identity")+
#   scale_fill_manual(values=sex.cols)+
#   labs(title = "Barplot of Gender Proportions Within Diagnostic Groups", x = "Diagnosis", y = "Proportion")
# ggplotly(p)
# api_create(p, filename = "prop_plot_sex")

dx.m <- c(22, 24, 4, 65, 33)
dx.f <- c(30, 18, 1, 27, 29)
df.3 <- data.frame(dx, dx.m, dx.f)

p <- plot_ly(df.3, x = ~dx, y = ~dx.m, type = 'bar', name = 'male', marker = list(color = '#16acfc')) %>%
  add_trace(y = ~dx.f, name = 'female', marker = list(color = '#fc1676')) %>%
  layout(title = 'Barplot of Diagnostic Group by Sex',
         yaxis = list(title = 'Count'),
         xaxis = list(title = 'Figure 1.2'),
         barmode = 'group')
p
api_create(p, filename = "prop_plot_sex")

table(status[idx,]$Sex)
chisq.test(table(status[idx,]$Sex))

addmargins(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))
chisq.test(table(droplevels(status[idx,]$most_general), status[idx,]$Sex))




### site
# p <- ggplot(status[idx,], aes(most_general, fill=site)) +
#   labs(title = "Barplot of Diagnostic Group by Recruitment Site", x = "Diagnosis", y = "Counts")+
#   scale_fill_manual(values=site.pal)+
#   geom_bar()
# ggplotly(p)

dx <- c('bacterial', 'greyb', 'greyv', 'viral', 'HC')
dx.chw <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,1])
dx.euc101 <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,2])
dx.fed <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,3])
dx.oxf <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,4])
dx.smh <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,5])
dx.sot <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,6])
dx.ucsd <- c(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))[,7])
site <- c(dx.chw, dx.euc101, dx.fed, dx.oxf, dx.smh, dx.sot, dx.ucsd)
site.df <- data.frame(dx.chw, dx.euc101, dx.fed, dx.oxf, dx.smh, dx.sot, dx.ucsd)
site.df <- site.df[c(rownames(site.df)[1:3], c("viral","HC")),] # swap order or viral and HC
site.df

text <- c(paste0('Total: ', 52), paste0('Total: ', 42), paste0('Total: ', 5), paste0('Total: ', 92), paste0('Total: ', 62))

p <- plot_ly(site.df, x = ~dx, y = ~dx.chw, type = 'bar', name = 'CHW', marker = list(color = '#ed0404'), text=text) %>%
  add_trace(y = ~dx.euc101, name = 'EUC101', marker = list(color = '##fc5716')) %>%
  add_trace(y = ~dx.fed, name = 'FED', marker = list(color = '##d7fc35')) %>%
  add_trace(y = ~dx.oxf, name = 'OXF', marker = list(color = '#35c7fc')) %>%
  add_trace(y = ~dx.smh, name = 'SMH', marker = list(color = '#16fc31')) %>%
  add_trace(y = ~dx.sot, name = 'SOT', marker = list(color = '#464647')) %>%
  add_trace(y = ~dx.ucsd, name = 'UCSD', marker = list(color = '#165bfc')) %>%
  layout(title = 'Barplot of Diagnostic Group by Recruitment Site',
         yaxis = list(title = 'Count'),
         xaxis = list(title = 'Figure 1.1'),
         barmode = 'stack')
p
api_create(p, filename = "barplot_recruitment")

table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site))
chisq.test(table(droplevels(status[idx,]$most_general), droplevels(status[idx,]$site)))



plot.code <- ifelse(status[idx,]$most_general == 'bacterial', 1,
       ifelse(status[idx,]$most_general == 'greyb', 2,
              ifelse(status[idx,]$most_general == 'greyv', 3,
                     ifelse(status[idx,]$most_general == 'viral', 4,
                            ifelse(status[idx,]$most_general == 'HC', 5, 0)))))


### Age Breakdown
p<-plot_ly(status[idx,][order(plot.code),], x = ~droplevels(most_general), y = ~Age..months.,
           color = status[idx,][order(plot.code),]$Sex, colors=c(sex.cols), type = "box",
           text= ~paste0('<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis, '<br>Age: ', Age..months.)) %>%
  layout(boxmode = "group",
         title = 'Box and Whisker Plot of Age (momths) by Diagnostic Group, Split by Gender', 
         xaxis = list(title = 'Diagnosis'),
         yaxis = list(title = 'Age'))
p
# api_create(p, filename = "box_whisker_age")

# check sums
status[idx,][order(plot.code),] %>%
  group_by(most_general, Sex) %>%
  summarise(age.m = median(Age..months.))

# stats
ddply(status[idx,], ~most_general, summarise, mean=mean(Age..months.))
a <- as.data.frame(droplevels(status[idx,]$most_general))
colnames(a) <- 'dx'
a$age <-status[idx,]$Age..months.
head(a)
summary(a)
b <- ddply(a,~dx)
anova.res <- aov(age ~ dx, data = b)
summary(anova.res)



### Inflamatory Marker Breakdown
v <- status[idx,][clean.idx,]$most_general == 'bacterial' | status[idx,][clean.idx,]$most_general == 'greyb' |
  status[idx,][clean.idx,]$most_general == 'greyv' | status[idx,][clean.idx,]$most_general == 'viral'
status[idx,][clean.idx,][v,]

p<-plot_ly(status[idx,][clean.idx,][v,], x = ~droplevels(status[idx,][clean.idx,][v,]$most_general), y = ~status[idx,][clean.idx,][v,]$WBC,
        color = status[idx,][clean.idx,][v,]$Sex, colors=c(sex.cols), type = "box",
        text= ~paste0('<br>label:', my_category_2, '<br>Age: ', Age..months., '<br>WBC: ', WBC, '<br>CRP: ', array.contemporary.CRP, '<br>Diagnosis: ',Diagnosis)) %>%
  layout(boxmode = "group",
         title = 'Box and Whisker Plot of WBC Count by Diagnostic Group, Split by Gender', 
         xaxis = list(title = 'Diagnosis'),
         yaxis = list(title = 'WBC Count'))
p
api_create(p, filename = "box_whisker_wbc")

# check sums
status[idx,][clean.idx,] %>%
  group_by(most_general) %>%
  summarise(crp.m = median(WBC))

p<-plot_ly(status[idx,][clean.idx,][v,], x = ~droplevels(status[idx,][clean.idx,][v,]$most_general), y = ~as.numeric((as.character(status[idx,][clean.idx,][v,]$array.contemporary.CRP))),
           color = status[idx,][clean.idx,][v,]$Sex, colors=c(sex.cols), type = "box",
           text= ~paste0('<br>label:', my_category_2, '<br>Age: ', Age..months., '<br>WBC: ', WBC, '<br>CRP: ', array.contemporary.CRP, '<br>Diagnosis: ',Diagnosis)) %>%
  layout(boxmode = "group",
         title = 'Box and Whisker Plot of CRP Count by Diagnostic Group, Split by Gender', 
         xaxis = list(title = 'Diagnosis'),
         yaxis = list(title = 'CRP Count'))
p
api_create(p, filename = "box_whisker_crp")


crp <- as.numeric(as.character(status[idx,][clean.idx,]$array.contemporary.CRP))
crp[status[idx,][clean.idx,]$Sex=='M' & status[idx,][clean.idx,]$most_general == 'bacterial']
crp[status[idx,][clean.idx,]$Sex=='F' & status[idx,][clean.idx,]$most_general == 'bacterial']

t.test(crp[status[idx,][clean.idx,]$Sex=='M' & status[idx,][clean.idx,]$most_general == 'bacterial'],
       crp[status[idx,][clean.idx,]$Sex=='F' & status[idx,][clean.idx,]$most_general == 'bacterial'])
# p-value = 0.6789
# no stat sig difference in crp between males and females in bct group

dim(status)
dim(e.set[,idx])
dim(status[idx,])









####### PCA #######
# full pca
full.pca <- prcomp(t(e.set[,idx]), scale=TRUE)

# filter pca
e.set.f <- as.data.frame(e.set[,idx])
e.set.f <- as.data.frame(t(e.set.f))
dim(e.set.f)

# mean/variance calculations
x_var <- apply(e.set[,idx], 1, var)
x_mean <- apply(e.set[,idx], 1, mean)
# hist(x_mean)
# hist(log2(x_mean))
plot(log2(x_mean), log2(x_var), pch='.')
abline(v=log2(5), col='red')

plot_ly(x = log2(x_mean), y=x_var, color = x_mean, type='scatter',
        mode = 'markers', marker = list(pacity= 0.5)) %>%
  layout(shapes=list(type='line', x0= log2(5), x1= log2(5), y0=0, y1=10, line=list(dash='dot', width=1)),
         title = 'Figure XYZ',
         xaxis = list(title = "X-Axis", showgrid = TRUE),
         yaxis = list(title = "Y-Axis", showgrid = TRUE))
  

X <- e.set[,idx][which(x_mean > 5),]
dim(X)

X.t <- t(X)
X.pca <- prcomp(X.t, scale = TRUE)


# set pca
# pca <- full.pca
pca <- X.pca
pair1 <- as.data.frame(pca$x[,1:2])
pair2 <- as.data.frame(pca$x[,3:4])
pair3D <- as.data.frame(pca$x[,1:3])

fviz_eig(pca)

ve <- pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]


# most_gen 2D
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, color = ~droplevels(status[idx,]$most_general), size = status[idx,]$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
         xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
         yaxis = list(title = paste0("PC1: (", round(pve[2],2), '%)')))
p
# api_create(p, filename = "2d_pca")

# most_gen 3D
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~droplevels(status[idx,]$most_general), size = status[idx,]$Age..months.,
             colors=c(dx.cols), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Diagnostic Groups by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))
p
# api_create(p, filename = "3d_pca")

# # category
# plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$category, size = status[idx,]$Age..months.,
#         colors=c(cat.pal), text= ~paste0('category: ', status[idx,]$category, '<br>age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
#   add_markers() %>%
#   layout(title = 'Category Groups by PCA 1-2-3, Age Size Mapping',
#          scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
#                       yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
#                       zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

# sex
p <- plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$Sex, size = ~status[idx,]$Age..months.,
             colors = c(sex.cols), text= ~paste0('age: ', status[idx,]$Age..months., '<br>WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Gender by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))

p
# api_create(p, filename = "sex_pca")


# site
plot_ly(pair3D, x = ~PC1, y = ~PC2, z = ~PC3, color = ~status[idx,]$site, size = ~status[idx,]$Age..months.,
        colors = c(site.pal), text= ~paste0('WBC: ', wbc, '<br>CRP: ',crp, '<br>label:',status[idx,]$my_category_2, '<br>Diagnosis: ',status[idx,]$Diagnosis)) %>%
  add_markers() %>%
  layout(title = 'Site Recruitment by PCA 1-2-3, Age Size Mapping',
         scene = list(xaxis = list(title = paste0("PC1: (", round(pve[1],2), '%)')),
                      yaxis = list(title = paste0("PC2: (", round(pve[2],2), '%)')),
                      zaxis = list(title = paste0("PC3: (", round(pve[3],2), '%)'))))




############################## CLUSTERING ##############################
### PCA

fviz_eig(pca)

ve <- pca$sdev^2
pve <- ve/sum(ve)*100
pve[1:5]

fviz_nbclust(X.t, kmeans, method = "wss")
fviz_nbclust(X.t, kmeans, method = "silhouette")

gap_stat <- clusGap(X.t, FUN = kmeans, nstart = 10,
                    K.max = 20, B = 25)
fviz_gap_stat(gap_stat)
#
# ### K2
# k2 <- kmeans(X.t[,results.hits], centers = 2, nstart = 500)
# str(k2)
# k2$cluster <- as.factor(k2$cluster)
# 
# 
# k2.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
#                                    'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
# k2.df$cluster <- k2$cluster[clean.idx]
# dim(k2.df)
# k2.df[1:5,(ncol(k2.df)-3): ncol(k2.df)]
# k2.df$array.contemporary.CRP <- as.numeric(as.character(k2.df$array.contemporary.CRP))
# 
# addmargins(table(k2$cluster, droplevels(status[idx,]$most_general)))
# addmargins(table(k2$cluster, droplevels(status[idx,]$more_general)))
# 
# k2.clus.col <- c("#ed0404", "#165bfc")
# positions.more <- c('bacterial', 'greyb', 'greyv', 'adeno', 'flu', 'RSV', 'viralother')
# ggplot(status[idx,], aes(more_general, fill=k2$cluster)) +
#   labs(title = "Barplot of Diagnostic Groups by Cluster", x = "Diagnosis", y = "Counts")+
#   scale_x_discrete(limits = positions.more)+
#   scale_fill_manual(values=k2.clus.col)+
#   geom_bar()
# # table(k2$cluster, status[idx, c('Sex')])
# 
# dx.def <- status[idx,]$most_general == 'bacterial' | status[idx,]$most_general == 'viral'
# dx.grey <- status[idx,]$most_general == 'greyb' | status[idx,]$most_general == 'greyv'
# dx.def
# 
# dim(pair1[dx.def,])
# length(status[idx,]$most_general[dx.def])
# 
# # DEF DX PC1 PC2
# ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$most_general[dx.def]), size=2) +
#   xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
#   ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
#   geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
#   scale_color_manual(values=k2.clus.col)+
#   ggtitle("Cluster Assignment of Bacterial and Viral Against first two PC")
# 
# # DEF DX PC3 PC4
# ggplot(pair2[dx.def,], aes(PC3, PC4)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$most_general[dx.def]), size=2) +
#   xlab(paste0("PC3: (", round(pve[3],2), '%)') ) +
#   ylab(paste0("PC4: (", round(pve[4],2), '%)') ) +
#   geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
#   scale_color_manual(values=k2.clus.col)+
#   ggtitle("Cluster Assignment of Bacterial and Viral Against first two PC")
# 
# # M F PC1 PC2
# ggplot(pair1[dx.def,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.def], shape=status[idx,]$Sex[dx.def]), size=2) +
#   xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
#   ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
#   geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
#   scale_color_manual(values=k2.clus.col)+
#   ggtitle("Cluster Assignment of Male, Female Against first two PC")
# 
# # Grey Dx PC1 PC2
# ggplot(pair1[dx.grey,], aes(PC1, PC2)) + geom_point(aes(color=k2$cluster[dx.grey], shape=status[idx,]$most_general[dx.grey]), size=2) +
#   xlab(paste0("PC1: (", round(pve[1],2), '%)') ) +
#   ylab(paste0("PC2: (", round(pve[2],2), '%)') ) +
#   geom_hline(yintercept = 0, linetype="longdash", colour="grey", size=1) +
#   geom_vline(xintercept = 0, linetype="longdash", colour="grey", size=1) +
#   scale_color_manual(values=k2.clus.col)+
#   ggtitle("Cluster Assignment of Probable Bacterial and Probable Viral Against first two PC")
# 
# 
# 
# 
# ggplot(status[idx,], aes(most_general, Age..months., fill=k2$cluster)) + geom_boxplot()+
#   scale_x_discrete(limits = positions)+
#   xlab('Diagnostic Group') +
#   ylab('Age') +
#   scale_fill_manual(values=k2.clus.col)+
#   ggtitle("Age (months) by Diagnostic Group, Split by Cluster")
# 
# p1.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k2$cluster[clean.idx])) + geom_boxplot() +
#   # scale_y_continuous(limits = c(0, 50))+
#   ylab('WBC Count') +
#   xlab('') +
#   scale_fill_manual(values=k2.clus.col)+
#   ggtitle("Box Plot of WBC and CRP Count by Diagnosis")
# p2.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
#   as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k2$cluster[clean.idx]))+
#   geom_boxplot()+
#   ylab('CRP Count') +
#   xlab('Diagnosis') +
#   scale_fill_manual(values=k2.clus.col)
# gridExtra::grid.arrange(p1.wbc, p2.crp, nrow = 2)
# 
# 
# filter.bct <- status[idx,]$most_general == 'bacterial'
# filter.clus <- k2$cluster == 1
# filter.comb <- filter.bct & filter.clus
# status[idx,]$Diagnosis[filter.comb]
# 
# filter.clus <- k2$cluster == 2
# filter.comb <- filter.bct & filter.clus
# status[idx,]$Diagnosis[filter.comb]
# # View(status[idx,][filter.comb, c('my_category_2', 'most_general',
# #                                  'more_general', 'Age..months.', 'Sex', 'WBC',
# #                                  'array.contemporary.CRP', 'Diagnosis')])
# 
# detach('package:plyr', unload = TRUE, character.only = TRUE)
# k2.df %>%
#   select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
#   group_by(cluster, most_general) %>%
#   summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))
# 
# k2.df$WBC[k2.df$cluster == 1 & k2.df$most_general == 'bacterial']
# k2.df$WBC[k2.df$cluster == 2 & k2.df$most_general == 'bacterial']
# t.test(k2.df$WBC[k2.df$cluster == 1 & k2.df$most_general == 'bacterial'], k2.df$WBC[k2.df$cluster == 2 & k2.df$most_general == 'bacterial'])
# # not enough power to detect sig diff within bct group given few bct assigned to class 2
# # no difference with crp
# 
# k2.df$WBC[k2.df$cluster == 1]
# k2.df$WBC[k2.df$cluster == 2]
# t.test(k2.df$WBC[k2.df$cluster == 1], k2.df$WBC[k2.df$cluster == 2])
# # t = 4.628, df = 138.86, p-value = 8.385e-06
# 
# k2.df$array.contemporary.CRP[k2.df$cluster == 1]
# k2.df$array.contemporary.CRP[k2.df$cluster == 2]
# t.test(k2.df$array.contemporary.CRP[k2.df$cluster == 1]
#        ,k2.df$array.contemporary.CRP[k2.df$cluster == 2])
# # t = 8.1468, df = 127.89, p-value = 2.92e-13
# 
# 
# 
# 
# ### K3
# k3 <- kmeans(X.t[,results.hits], centers = 3, nstart = 25)
# # str(k3)
# k3$cluster <- as.factor(k3$cluster)
# 
# addmargins(table(k3$cluster, droplevels(status[idx,]$most_general)))
# addmargins(table(k3$cluster, droplevels(status[idx,]$more_general)))
# 
# k3.df <- status[idx,][clean.idx, c('my_category_2', 'most_general', 'more_general',
#                                    'Age..months.', 'Sex', 'WBC', 'array.contemporary.CRP', 'Diagnosis')]
# k3.df$cluster <- k3$cluster[clean.idx]
# dim(k3.df)
# k3.df[1:8,(ncol(k2.df)-3): ncol(k2.df)]
# k3.df$array.contemporary.CRP <- as.numeric(as.character(k3.df$array.contemporary.CRP))
# 
# # had to un attach plyr to get working
# k3.df %>%
#   select(WBC, cluster, most_general, array.contemporary.CRP, Age..months., Sex) %>%
#   group_by(cluster, Sex) %>%
#   summarise(wbc.m = mean(WBC), crp.m = mean(array.contemporary.CRP), age.m = mean(Age..months.))
# 
# k3.clus.col <- c("#fcac16","#ed0404", "#165bfc")
# ggplot(status[idx,], aes(more_general, fill=k3$cluster)) +
#   labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
#   scale_x_discrete(limits = positions.more)+
#   scale_fill_manual(values=k3.clus.col)+
#   geom_bar()
# 
# ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k3$cluster, shape=status[idx,]$most_general), size=2) +
#   xlab("First Principal Component") +
#   ylab("Second Principal Component") +
#   scale_color_manual(values=k3.clus.col)+
#   ggtitle("Cluster Assignment of First Two Principal Components")
# 
# ggplot(status[idx,], aes(k3$cluster, Age..months., fill=Sex)) + geom_boxplot()+
#   xlab('Diagnostic Group') +
#   ylab('Age') +
#   ggtitle("Age (months) by Diagnostic Group, Split by Gender")
# 
# k3.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k3$cluster[clean.idx])) + geom_boxplot() +
#   # scale_y_continuous(limits = c(0, 50))+
#   ylab('WBC Count') +
#   scale_x_discrete(limits = positions)+
#   scale_fill_manual(values=k3.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# 
# k3.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
#                                                as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k3$cluster[clean.idx]))+
#   geom_boxplot()+
#   # scale_y_continuous(limits = c(0, 230))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_fill_manual(values=k3.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k3.wbc, k3.crp, nrow = 2)
# 
# k3.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, color=k3$cluster[clean.idx])) + geom_jitter(alpha = .9, width = 0.2) +
#   # scale_y_continuous(limits = c(0, 50))+
#   ylab('WBC Count') +
#   scale_x_discrete(limits = positions)+
#   scale_color_manual(values=k3.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# 
# k3.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
#                                                as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k3$cluster[clean.idx]))+
#   geom_jitter(alpha = .9, width = 0.2)+
#   # scale_y_continuous(limits = c(0, 230))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_color_manual(values=k3.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k3.wbc, k3.crp, nrow = 2)
# 
# 
# filter.bct <- status[idx,]$most_general == 'bacterial'
# # filter.bct <- status[idx,]$most_general == 'greyb'
# filter.clus <- k3$cluster == 2
# filter.clus <- k3$cluster == 1
# filter.clus <- k3$cluster == 3
# filter.clus <- k3$cluster == 1 | k3$cluster == 3
# 
# filter.comb <- filter.bct & filter.clus
# 
# filter.comb
# # View(status[idx,])
# status[idx,]$Diagnosis[filter.comb]
# View(status[idx,][filter.comb, c('my_category_2', 'most_general',
#                                  'more_general', 'Age..months.', 'Sex', 'WBC',
#                                  'array.contemporary.CRP', 'Diagnosis')])
# 
# 
# 
# 
# ### K5
# k5 <- kmeans(X.t[,results.hits], centers = 5, nstart = 25)
# k5$cluster <- as.factor(k5$cluster)
# 
# addmargins(table(k5$cluster, droplevels(status[idx,]$most_general)))
# addmargins(table(k5$cluster, droplevels(status[idx,]$more_general)))
# 
# k5.clus.col <- c("#fcac16","#165bfc", '#16fcd2', '#fc16af', "#ed0404")
# ggplot(status[idx,], aes(more_general, fill=k5$cluster)) +
#   labs(title = "Barplot of Diagnostic Group Breakdown by Gender", x = "Diagnosis", y = "Counts")+
#   scale_x_discrete(limits = positions.more)+
#   scale_fill_manual(values=k5.clus.col)+
#   geom_bar()
# 
# ggplot(pair1, aes(PC1, PC2)) + geom_point(aes(color=k5$cluster, shape=status[idx,]$most_general), size=2) +
#   xlab("First Principal Component") +
#   ylab("Second Principal Component") +
#   scale_color_manual(values=k5.clus.col)+
#   ggtitle("Cluster Assignment of First Two Principal Components")
# 
# k5.wbc <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, fill=k5$cluster[clean.idx])) + geom_boxplot() +
#   # scale_y_continuous(limits = c(0, 50))+
#   ylab('WBC Count') +
#   scale_x_discrete(limits = positions)+
#   scale_fill_manual(values=k5.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# k5.crp <- ggplot(status[idx,][clean.idx,], aes(most_general,
#   as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), fill=k5$cluster[clean.idx]))+
#   geom_boxplot()+
#   # scale_y_continuous(limits = c(0, 250))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_fill_manual(values=k5.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k5.wbc, k5.crp, nrow = 2)
# 
# k5.wbc <- ggplot(status[idx,][clean.idx,], aes(more_general, WBC, color=k5$cluster[clean.idx])) +
#   geom_jitter(alpha = .9, width = 0.2) +
#   # scale_y_continuous(limits = c(0, 50))+
#   ylab('WBC Count') +
#   # scale_x_discrete(limits = positions)+
#   scale_color_manual(values=k5.clus.col)+
#   ggtitle("Jitter Plot: WBC Count by Diagnosis")
# k5.crp <- ggplot(status[idx,][clean.idx,], aes(more_general,
#   as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])), color=k5$cluster[clean.idx]))+
#   geom_jitter(alpha = .9, width = 0.2)+
#   # scale_y_continuous(limits = c(0, 250))+
#   xlab('Cluster Assignment') +
#   ylab('CRP Count') +
#   scale_color_manual(values=k5.clus.col)+
#   ggtitle("Jitter Plot: CRP Count by Diagnosis")
# gridExtra::grid.arrange(k5.wbc, k5.crp, nrow = 2)


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