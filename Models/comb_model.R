library(Biobase)
library(GEOquery)
library(dplyr)
library(ggfortify)
getwd()
setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()

# gset.m <- getGEO('GSE72809')
# gset.i <- getGEO('GSE42026')

gset.m <-readRDS(file = "gset_GSE72809")
gset.i <-readRDS(file = "gset_GSE42026")

gset.m.df <- exprs(gset.m[[1]])
gset.i.df <- exprs(gset.i[[1]])

n1 <- dim(gset.m.df)[2]
n2 <- dim(gset.i.df)[2]

gset.m.df[1:5,1]
gset.i.df[1:5,1]

attributes(gset.m[[1]])
phenoData(gset.m[[1]])$'category:ch1'
class(gset.m[[1]])

#remove(gset.m, gset.i)
cat(paste('features gset.m:', length(rownames(gset.m.df))))
cat(paste('features gset.i:', length(rownames(gset.i.df))))

# transpose
gset.m.t.df <- as.data.frame(t(gset.m.df))
gset.i.t.df <- as.data.frame(t(gset.i.df))

#remove(gset.m.df, gset.i.df)

common <-intersect(colnames(gset.m.t.df), colnames(gset.i.t.df))
length(common) # 39426

dim(gset.m.t.df[,common])
dim(gset.i.t.df[,common])

gset.m.c <- gset.m.t.df[,common]
gset.i.c <- gset.i.t.df[,common]

gset.m.c[1:5,1:5]
gset.i.c[1:5,1:5]

# standardize
gset.m.s <- as.data.frame(scale(gset.m.c))
gset.i.s <- as.data.frame(scale(gset.i.c))
#remove(gset.m.c, gset.i.c, gset.m.t.df, gset.i.t.df)

dim(gset.m.s)
dim(gset.i.s)

# labels
y1 <- phenoData(gset.m[[1]])$'category:ch1'
y2 <- gset.i[[1]]$'infecting pathogen:ch1'

y1 <- ifelse(y1 == 'Definite Bacterial', 1, 0)
y2 <- ifelse(y2 == 'gram positive bacterial infection', 1, 0)

gset.m.s$label <- y1
gset.i.s$label <- y2

# gset.m.s$exp <- 1
# gset.i.s$exp <- 2

gset.m.s[1:5,(ncol(gset.m.s)-2):ncol(gset.m.s)]
gset.i.s[1:5,(ncol(gset.i.s)-2):ncol(gset.i.s)]

# merge dataframes
gset.c <- rbind(gset.m.s, gset.i.s)
dim(gset.c)
gset.c[1:5,(ncol(gset.c)-2):ncol(gset.c)]
#remove(gset.i.df, gset.m.df, gset.i.s, gset.m.s)

gset.c[1,1:4] # first row of the combined # GSM1872417
gset.m.s[1,1:4] # GSM1872417

gset.c[(n1+1),1:4] # GSM1030788
gset.i.s[1,1:4] # GSM1030788


# save / load
getwd()
# saveRDS(gset.c, file = 'gset_GSE72809_GSE42026')
# getwd()
# setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
# gset.c <-readRDS(file = 'gset_GSE72809_GSE42026')


## verify labels
# gset_f_GSE72809 <-readRDS(file = 'gset_f_GSE72809')
# gset_f_GSE42026 <-readRDS(file = 'gset_f_GSE42026')
# 
# GSMList(gset_f_GSE72809)[c(1,53,145)] # def baf, def vrl, prob vrl
# GSMList(gset_f_GSE42026)[c(1,19,52,71)] # def bac, control, influenzae, rsv
# 
# gset.m[[1]]$'category:ch1'[c(1,53,145)]
  # "Definite Bacterial" "Definite Viral" "Probable Bacterial"
# gset.i[[1]]$'infecting pathogen:ch1'[c(1,19,52,71)]
  # "gram positive bacterial infection" "none" "Influenza A H1N1/09" "RSV" 
# gset.c$label[c(1,52,53,292,293)] # check transition points between bct and other

#remove(gset_f_GSE42026, gset_f_GSE72809)
#remove(gset.i, gset.m)


# PCA
# check var 1 and mean 0
apply(gset.c[,1:10], 2, var)
apply(gset.c[,1:10], 2, mean)

gset.pca <- prcomp(gset.c, scale. = FALSE)

dim(gset.pca$rotation)
gset.pca$rotation[1:5,1:4]
gset.pca$x[1:5,1:4]

summary(gset.pca)

plot(gset.pca, type = 'l')
# biplot(gset.pca)
# autoplot(gset.pca)
autoplot(gset.pca, data = gset.c, colour = 'exp')
# autoplot(gset.pca, data = gset.c, label = TRUE, shape = FALSE, loading = TRUE, colour = 'label')




biplot(pca_result, scale = 0)
# pca_result <- prcomp(USArrests, scale = TRUE)

# gset.cov <- cov(gset.c)
# saveRDS(gset.cov, file = 'gset.cov')


# to do
# pca plots to check no systemic skewing
# partition gset.c into bac and other set for differential gene expression analysis
# DGE
# construct classifier from SDG






# 2
