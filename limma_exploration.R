library(tidyverse)
library(plyr)
library(limma)
library(DESeq2)
library(Biobase)

setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

rm(list=setdiff(ls(), 'all'))
load('esets.RData')
ls()

dim(e.set)
e.set.t <- t(e.set)
dim(e.set.t)

# unique(status$most_general)

# e.set.bv <- e.set.t[status[,'most_general'] == 'bacterial' | status[,'most_general'] == 'viral',]
e.set.f <- e.set.t[status[,'most_general'] == 'bacterial'
                   | status[,'most_general'] == 'viral'
                   | status[,'most_general'] == 'greyb'
                   | status[,'most_general'] == 'greyv'
                   ,]

label.f <- status[status['most_general'] == 'bacterial' |
                    status[,'most_general'] == 'viral' |
                    status[,'most_general'] == 'greyv' |
                    status[,'most_general'] == 'greyb',
                  c('most_general')]


label <- as.character(label.f)
X <- as.matrix(t(e.set.f))
dim(X)

e.set.f <- as.data.frame(e.set.f)
e.set.f$label <- label

dim(e.set.f)
length(label)
e.set.f[1:5, (ncol(e.set.f)-3):ncol(e.set.f)]

### DESIGN MATRIX
design <- model.matrix(~label, data = e.set.f)
bct.int <- ifelse(design[,2] == 1 | design[,3] == 1 | design[,4] == 1, 0, 1)
bct.int
design[,1] = bct.int
design[1:5,]
dim(design)

# design <- model.matrix(~ e.set.bv$label) # equaivalent

colSums(design)
colnames(design)<- c("bct","greyb","greyv", 'vrl')

contrast.matrix<- makeContrasts("bct-greyb", levels=design)
contrast.matrix
# colnames(fit$coefficients)

fit <- lmFit(X, design)
fit2<- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
topTable(fit2, adjust="BH", number=10)


p.threshold <- 0.05
result$pval.threshold <- as.logical(result$adj.P.Val < p.threshold)
head(result)

results <- decideTests(fit2, method="global", adjust.method="BH", p.value=0.05,lfc=0)
summary(results)
vennDiagram(results)


# expression vs variability
hist(fit2$Amean)
plotSA(fit2)
# rule of thumb is to keep genes with exp > 5-10 in > k arays
exprs.thresh <- 4
keep <- fit2$Amean > exprs.thresh

fit2 <- eBayes(fit2[keep,], trend=TRUE)
plotSA(fit2)

results <- decideTests(fit2, method="global", adjust.method="BH", p.value=0.05)
dim(results)



colnames(fit$coefficients)
tab.results <- topTable(fit, coef = colnames(fit$coefficients)[4], adjust.method = 'BH', number = ncol(e.set))
# tab.results <- topTable(fit2, adjust.method = 'BH', number = 5)
dim(tab.results)
head(tab.results)
class(tab.results)

p.threshold <- 0.05
tab.results$pval.threshold <- as.logical(tab.results$adj.P.Val < p.threshold)
head(tab.results)

# sig.trans.p <- sig.trans
sig.trans <- rownames(tab.results)[which(tab.results$pval.threshold)]
head(sig.trans)
length(sig.trans)

results <- decideTests(fit[,'labelviral'], method="global", adjust.method="BH", p.value=0.05,lfc=0)
summary(results)
vennDiagram(results)

length(intersect(sig.trans, sig.trans.p))





## multiple groups
set.seed(42)
sd<- 0.3*sqrt(4/rchisq(100,df=4))
sd
y2<- matrix(rnorm(100*9,sd=sd),100,9)
dim(y2)

rownames(y2)<- paste("Gene",1:100)

y2[1:2,4:6]<- y2[1:2,4:6] + 2
y2[1:2,7:9]<- y2[1:2,7:9] + 2

head(y2)
dim(y2)
class(y2)

y2.t <- as.data.frame(t(y2))
dim(y2.t)
y2.t[,98:101]

y2.t$label <- c('A','A','A','B','B','B','C','C','C')

f<- factor(c(rep("A",3),rep("B",3),rep("C",3)),levels=c("A","B","C"))
f


design <- model.matrix(~label, data = y2.t)
design[,1][4:9]=0
# design<- model.matrix(~0+f)

colnames(design)<- c("A","B","C")

# contrast.matrix<- makeContrasts("B-A", "C-A", "C-B", levels=design)
contrast.matrix<- makeContrasts("C-A", levels=design)
contrast.matrix

fit1<- lmFit(y2,design)
# colnames(fit1$coefficients)
fit2<- contrasts.fit(fit1, contrast.matrix)
fit2<- eBayes(fit2)
topTable(fit2, adjust="BH")








?makeContrasts




# end