
getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

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
clean

clean.idx <- seq(1:nrow(status[idx,]))[-(clean)]

status[idx,]$array.contemporary.CRP[clean.idx]
wbc <- status[idx,]$WBC[clean.idx]
crp <- as.numeric(as.character(status[idx,]$array.contemporary.CRP[clean.idx]))

