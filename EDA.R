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

which(status[idx,]$array.contemporary.CRP == 'na')
length(which(status[idx,]$array.contemporary.CRP == 'na'))

crp.idx <- seq(1:nrow(status[idx,]))[-(which(status[idx,]$array.contemporary.CRP == 'na'))]
crp.idx

as.numeric(status[idx,]$array.contemporary.CRP[crp.idx])
cor(status[idx,]$WBC[crp.idx], as.numeric(status[idx,]$array.contemporary.CRP[crp.idx]))


# Gener Breakdown
ggplot(status[idx,], aes(most_general, group=Sex)) + 
  geom_bar(aes(y = ..prop.., fill = factor(..x..)), stat="count") +
  geom_text(aes( label = scales::percent(..prop..),
                 y= ..prop.. ), stat= "count", vjust = -.5) +
  labs(y = "Percent", fill="day") +
  facet_grid(~Sex) +
  scale_y_continuous(labels = scales::percent)
# ggplot(status[idx,], aes(Sex, fill=Sex)) + geom_bar(position="dodge")+
#   facet_wrap(~most_general, nrow=3, ncol=2) +
#   xlab('Cluster') +
#   ylab('Counts') +
#   ggtitle("Bar Plot of Gender Split Between Clusters")

# ggplot(status[idx,], aes(most_general, fill=most_general)) + 
#   geom_bar()+
#   facet_grid(~Sex)

table(status[idx,]$Sex)
tab.g <- table(droplevels(status[idx,]$most_general), status[idx,]$Sex)
tab.g
chisq.test(tab.g)


# Diagnosis Breakdowns
p1 <- ggplot(status[idx,], aes(most_general, WBC, color=Sex)) + geom_jitter(alpha = .5, width = .15) +
  scale_y_continuous(limits = c(0, 50))+
  ylab('WBC Count') +
  ggtitle("Jitter Plot: WBC Count by Diagnosis")

p2 <- ggplot(status[idx,][crp.idx,], aes(most_general,
  as.numeric(status[idx,]$array.contemporary.CRP[crp.idx]), color=Sex))+
  geom_jitter(alpha = .5, width = .15)+
  scale_y_continuous(limits = c(0, 210))+
  xlab('Cluster Assignment') +
  ylab('CRP Count') +
  ggtitle("Jitter Plot: CRP Count by Diagnosis")
gridExtra::grid.arrange(p1, p2, nrow = 2)



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