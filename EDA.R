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
  status['most_general'] == 'greyv' |
  status['most_general'] == 'OD'
sum(idx)

# e.set.f <- e.set[,idx]
# status.f <- status[idx,]
# dim(e.set.f)
# dim(status.f)

attributes(status)$names

levels(status[idx,]$array.contemporary.CRP)


which(as.character(status[idx,]$array.contemporary.CRP) == 'na')
length(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'))

which(status[idx,]$WBC == 0)
clean<-union(which(as.character(status[idx,]$array.contemporary.CRP) == 'na'), which(status[idx,]$WBC == 0))
length(a)

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]
clean.idx

status[idx,]$array.contemporary.CRP[clean.idx]
as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))

cor(status[idx,]$WBC[clean.idx], as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))


ggplot(status[idx,][clean.idx,], aes(x=WBC, as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx])))) +
  geom_point(shape=1)+
  geom_smooth(method=lm)




# Gener Breakdown
positions <- c('bacterial', 'greyb', 'OD', 'greyv', 'viral')
ggplot(status[idx,], aes(most_general, fill=most_general, alpha=Sex)) +
  scale_alpha_manual(values=c(0.6, 1)) +
  labs(title = "Diagnostic Groups", x = "Diagnosis", y = "Counts")+
  scale_x_discrete(limits = positions)+
  geom_bar()

ggplot(status[idx,], aes(most_general, group=Sex)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="day") +
  facet_grid(~Sex) +
  scale_x_discrete(limits = positions)+
  scale_y_continuous(labels = scales::percent)

table(status[idx,]$Sex)
chisq.test(table(status[idx,]$Sex))

tab.g <- table(droplevels(status[idx,]$most_general), status[idx,]$Sex)
tab.g
chisq.test(tab.g)


# Age Breakdown
ggplot(status[idx,], aes(status[idx,]$most_general, Age..months., fill=Sex)) + geom_boxplot()+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age') +
  scale_x_discrete(limits = positions)+
  ggtitle("Age (months) by Diagnostic Group, Split by Gender")

ggplot(status[idx,], aes(status[idx,]$most_general, Age..months., color=Sex)) + geom_jitter(alpha = .9, width = 0.2)+
  scale_x_discrete(limits = positions)+
  xlab('Diagnostic Group') +
  ylab('Age (months)') +
  ggtitle("Age (months) by Diagnostic Group, Split by Gender")

ddply(status[idx,], ~most_general, summarise, mean=mean(Age..months.))

a <- as.data.frame(droplevels(status[idx,]$most_general))
colnames(a) <- 'dx'
a$age <-status[idx,]$Age..months.
head(a)
summary(a)

b <- ddply(a,~dx)

anova.res <- aov(age ~ dx, data = b)
summary(anova.res)

dim(status[clean.idx,])




# Inflamatory Marker Breakdown
positions <- c('bacterial', 'greyb', 'greyv', 'viral')
p1 <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, color=Sex)) + geom_jitter(alpha = .8, width = .2) +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")

p2 <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,][clean.idx,]$array.contemporary.CRP)), color=Sex))+
  geom_jitter(alpha = .8, width = .2)+
  # scale_y_continuous(limits = c(0, 230))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_x_discrete(limits = positions)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")

gridExtra::grid.arrange(p1, p2, nrow = 2)


p3 <- ggplot(status[idx,][clean.idx,], aes(most_general, WBC, color=Sex)) + geom_boxplot() +
  # scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  scale_x_discrete(limits = positions)+
  ggtitle("Jitter Plot: WBC Count by Diagnosis")

p4 <- ggplot(status[idx,][clean.idx,], aes(most_general,
  as.numeric(as.character(status[idx,][clean.idx,]$array.contemporary.CRP)), color=Sex))+
  geom_boxplot()+
  # scale_y_continuous(limits = c(0, 230))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  scale_x_discrete(limits = positions)+
  ggtitle("Jitter Plot: CRP Count by Diagnosis")

gridExtra::grid.arrange(p3, p4, nrow = 2)



# WBC
ggplot(status[idx,], aes(Sex, WBC)) + geom_boxplot() +
  xlab('Cluster Assignment') +
  ylab('WBC Count') +
  ggtitle("Box Plot of WBC Count by Cluster Assignment")

ggplot(status[idx,], aes(WBC, fill=Sex)) + geom_histogram(bins = 20)
  xlab('WBC') +
  ylab('Density') +
  ggtitle("Histogram of WBC Count by Cluster Assignment")

ggplot(status[idx,], aes(WBC, fill=Sex)) + geom_density(alpha=.5)+
  xlab('WBC') +
  ylab('Density') +
  ggtitle("Density Plot of WBC Count by Cluster Assignment")

ps <- (seq(0,499) + 0.5)/500
qs <- quantile(status[idx,]$WBC, ps)
normalqs <- qnorm(ps, mean(status[idx,]$WBC), sd(status[idx,]$WBC))
plot(normalqs,qs,xlab="Normal Percentiles",ylab="WBC Percentiles")
abline(0,1)

# bin
q <- quantile(status[idx,]$WBC, seq(0,1,0.2))
wbc_bin <- cut(status[idx,]$WBC, breaks=q, include.lowest = T)

tab <- with(status[idx,], table(wbc_bin, Sex))
print(prop.table(tab, 2))
chisq.test(tab)

df1 <- melt(ddply(status[idx,],.(Sex),function(x){prop.table(table(wbc_bin))}),id.vars = 1)

ggplot(df1, aes(x = variable,y = value)) +
  facet_wrap(~Sex, nrow=2, ncol=1) +
  scale_y_continuous(labels=percent) +
  geom_bar(stat = "identity") +
  xlab('Binned WBC Levels') +
  ylab('Percentage Bin Assignment') +
  ggtitle("% Assignment to Binned WBC Levels by Cluster")


# CRP
ggplot(status[idx,], aes(Sex, as.numeric(array.contemporary.CRP))) + geom_boxplot() +
  xlab('Gender') +
  ylab('CRP') +
  ggtitle("Box Plot of CRP Count by Gender")

ggplot(status[idx,], aes(as.numeric(array.contemporary.CRP), fill=Sex)) + geom_histogram(bins = 10)
xlab('CRP') +
  ylab('Density') +
  ggtitle("Histogram of CRP by Gender")

ggplot(status[idx,], aes(as.numeric(array.contemporary.CRP), fill=Sex)) + geom_density(alpha=.5)+
  xlab('CRP') +
  ylab('Density') +
  ggtitle("Density Plot of CRP by Gender")


ps <- (seq(0,499) + 0.5)/500
qs <- quantile(as.numeric(status[idx,]$array.contemporary.CRP), ps)
normalqs <- qnorm(ps, mean(as.numeric(status[idx,]$array.contemporary.CRP)), sd(as.numeric(status[idx,]$array.contemporary.CRP)))
plot(normalqs,qs,xlab="Normal Percentiles",ylab="CRP Percentiles")
abline(0,1)

# bin
q <- quantile(as.numeric(status[idx,]$array.contemporary.CRP), seq(0,1,0.25))
crp_bin <- cut(as.numeric(status[idx,]$array.contemporary.CRP), breaks=q, include.lowest = T)

tab <- with(status[idx,], table(crp_bin, Sex))
print(prop.table(tab, 2))
chisq.test(tab)

df1 <- melt(ddply(status[idx,],.(Sex),function(x){prop.table(table(crp_bin))}),id.vars = 1)
df2 <- melt(ddply(status[idx,],.(WBC),function(x){prop.table(table(wbc_bin))}),id.vars = 1)
df1
df2

ggplot(df1, aes(x = variable,y = value, fill=Sex)) +
  facet_wrap(~Sex, nrow=2, ncol=1) +
  scale_y_continuous(labels=percent) +
  geom_bar(stat = "identity") +
  xlab('Binned CRP Levels') +
  ylab('Percentage Bin Assignment') +
  ggtitle("% Assignment to Binned CRP Levels by Cluster")

















# end