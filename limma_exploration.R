library(cluster)
library(factoextra)
library(gridExtra)
library(tidyverse)
library(plyr)
require(reshape) # for melt()
require(scales) # for percent
library(limma)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('e.set')
# saveRDS(e.set, file = 'e.set')
dim(e.set)


e.set[1:4,1:4]

?lmFit
fit <- lmFit(e.set)
# dim(fit)

# attributes(fit[1,1])

fit <- eBayes(fit)

topTable(fit)
attributes(fit)