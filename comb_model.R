library(Biobase)
library(GEOquery)
library(dplyr)
#getwd()
#setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()

gset.m <- getGEO('GSE72809')
gset.i <- getGEO('GSE42026')

# gset.m <-readRDS(file = "gset_GSE72809")
# gset.i <-readRDS(file = "gset_GSE42026")

gset.m.df <- exprs(gset.m[[1]])
gset.i.df <- exprs(gset.i[[1]])
n1 <- dim(gset.m.df)[2]
n2 <- dim(gset.i.df)[2]

gset.m.df[1:5,1]
gset.i.df[1:5,1]

attributes(gset.m[[1]])
# gset.m[[1]]$`category:ch1`
phenoData(gset.m[[1]])$'category:ch1'
class(gset.m[[1]])
gset.m[[1]]

class(gset.m[[1]])
#remove(gset.m, gset.i)

length(rownames(gset.m.df))
length(rownames(gset.i.df))

#length(colnames(gset.m.t.df))
#length(colnames(gset.i.t.df))

# transpose
gset.m.t.df <- as.data.frame(t(gset.m.df))
gset.i.t.df <- as.data.frame(t(gset.i.df))

#remove(gset.m.df, gset.i.df)
print('got to the transpose')

common <-intersect(colnames(gset.m.t.df), colnames(gset.i.t.df))
length(common) # 39426
common[1]

dim(gset.m.t.df[,common])
dim(gset.i.t.df[,common])

gset.m.c <- gset.m.t.df[,common]
gset.i.c <- gset.i.t.df[,common]

# standardize
gset.m.s <- as.data.frame(scale(gset.m.c))
gset.i.s <- as.data.frame(scale(gset.i.c))
#remove(gset.m.c, gset.i.c, gset.m.t.df, gset.i.t.df)

dim(gset.m.s)
dim(gset.i.s)

# labels
y1 <- phenoData(gset.m[[1]])$'category:ch1'
y2 <- gset.i[[1]]$'infecting pathogen:ch1'

y1.t <- ifelse(y1 == 'Definite Bacterial', 1, 0)
y2.t <- ifelse(y2 == 'gram positive bacterial infection', 1, 0)

gset.m.s$label <- y1.t
gset.i.s$label <- y2.t

gset.m.s[,(ncol(gset.m.s)-5):ncol(gset.m.s)]
gset.i.s[,(ncol(gset.i.s)-5):ncol(gset.i.s)]

# merge dataframes
gset.c <- rbind(gset.m.s, gset.i.s)
dim(gset.c)
#remove(gset.i.df, gset.m.df, gset.i.s, gset.m.s)

gset.c[1,1:5] # first row of the combined # GSM1872417
gset.m.s[1,1:5] # GSM1872417

gset.c[(n1+1),1:5] # GSM1030788
gset.i.s[1,1:5] # GSM1030788

# verify correct label assignment
gset.c <-readRDS(file = 'gset_GSE72809_GSE42026')
gset_f_GSE72809 <- readRDS(file = 'gset_f_GSE72809')
# GSMList(gse)[[1]]

dim(gset.c)
sample <- rownames(gset.c)[c(1,50,100)]
sample
g
GSMList(gset_f_GSE72809)[['GSM1872417']]

#remove(gset_f_GSE72809)

# save / load
getwd()
saveRDS(gset.c, file = 'gset_GSE72809_GSE42026')
getwd()
setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
gset.c <-readRDS(file = 'gset_GSE72809_GSE42026')


# verify labels
gset_f_GSE72809 <-readRDS(file = 'gset_f_GSE72809')
gset_f_GSE42026 <-readRDS(file = 'gset_f_GSE42026')

GSMList(gset_f_GSE72809)[c(1,53,145)] # def baf, def vrl, prob vrl
GSMList(gset_f_GSE42026)[c(1,19,52,71)] # def bac, control, influenzae, rsv

gset.m[[1]]$'category:ch1'[c(1,53,145)] # "Definite Bacterial" "Definite Viral"     "Probable Bacterial"
gset.i[[1]]$'infecting pathogen:ch1'[c(1,19,52,71)] # "gram positive bacterial infection" "none" "Influenza A H1N1/09" "RSV" 

gset.c$label[c(1,52,53,292,293)] # check transition points between bct and other

#remove(gset_f_GSE42026, gset_f_GSE72809)
#remove(gset.i, gset.m)


# PCA
# check var 1 and mean 0
apply(gset.c[,1:10], 2, var)
apply(gset.c[,1:10], 2, mean)

print('computing the cov')

gset.cov <- cov(gset.c)
saveRDS(gset.cov, file = 'gset.cov')

# to do
# pca plots to check no systemic skewing
# partition gset.c into bac and other set for differential gene expression analysis
# DGE
# construct classifier from SDG










# 2
