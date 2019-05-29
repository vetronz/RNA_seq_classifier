source("http://www.bioconductor.org/biocLite.R")
biocLite("GEOquery")
biocLite("ggplot2")
biocLite("colorspace")
library("colorspace")
library("GEOquery")
library("Biobase")
library("GEOquery")
library("ggplot2")
library("limma")
install.packages("ggplot2", dependencies = TRUE)


#Download GSE file, put it in the current directory, and load it:
my.data <- getGEO('GSE68004', destdir=".")

str(my.data) #it's a list


my.dataset<-my.data[[1]] # and needs unlisting

#phenodata look at the phentypic data
pData(my.dataset)

#and look at the structure
str(pData(my.dataset))

#assign the phenotypic data into a new vector called phenos
phenos<-pData(my.dataset)
#change the class of the age column into numeric
class(phenos$"age (mos.):ch1")<-"numeric"

#and create a ggplot of the age ~ group
p<-ggplot(phenos, aes(phenos$"final condition:ch1", phenos$"age (mos.):ch1"))
p+ geom_boxplot()


#feature data - look at the probe illumina data
featureData(my.dataset)@varMetadata
#and their structure
str(featureData(my.dataset)@data)
# the probe_ids
featureData(my.dataset)@data$Probe_Id
# and the gene symbols
featureData(my.dataset)@data$Symbol



#boxplot of the expression values
boxplot(exprs(my.dataset[1:2400,1:45]))

#normalisation - bg correction - log2
lumi.N.Q <- lumiExpresso(my.dataset, normalize.param = list(method='quantile'))

#retrieve the expression values
eset<-exprs(lumi.N.Q)

#plot the normalised and logged boxplot
boxplot(eset[1:2400,1:45])

#look at the expression set structure
str(eset)

#Expression Set (exprSet) with 
#47323 genes
#162 samples

#its rownames
rownames(eset)[1:10]
#[1] "ILMN_1343291" "ILMN_1343295" "ILMN_1651199" "ILMN_1651209" "ILMN_1651210"
#[6] "ILMN_1651221" "ILMN_1651228" "ILMN_1651229" "ILMN_1651230" "ILMN_1651232"

#and column names
colnames(eset)
#[1] "GSM1660729" "GSM1660730" "GSM1660731" "GSM1660732" "GSM1660733" "GSM1660734"

#create a vector with the diagnostic groups
group<-pData(my.dataset)$"final condition:ch1"
#[1] "HAdV"    "HAdV"    "HAdV"    "HAdV"    "HAdV"    "HAdV"    "HAdV"    "HAdV"   



#look up a gene of interest - is it there??
"IFI44L" %in% featureData(my.dataset)@data$Symbol
match("IFI44L", featureData(my.dataset)@data$Symbol)

#and where is it?
eset[match("IFI44L", featureData(my.dataset)@data$Symbol),]

#create a dataframe with the expression value of the genes and the diagnostic groups
df<-data.frame(x=eset[match("SNCA", featureData(my.dataset)@data$Symbol),],grp=phenos$"final condition:ch1")
p<-ggplot(df, aes(df$grp, df$x))
p+ geom_boxplot()


#### An appropriate design matrix can be created 
design<-model.matrix(~0 + factor(group))
colnames(design) <-c("cKD","GAS","GAS_SF","Adeno","healthy","inKD")
###

#and a linear model fitted using
fit <- lmFit(eset, design)

#To make all pair-wise comparisons between the three groups the appropriate contrast matrix can be
contrast.matrix <- makeContrasts(Adeno-healthy, GAS-healthy, levels=design)

#bayesian fit
fit.cb <- eBayes(contrasts.fit(fit, contrast.matrix))

#look at the structure
str(fit.cb)
#and the top significantly differentially probes
head(topTable(fit.cb, coef=1, adjust="fdr",p.value = 0.01))

#save into txt files
write.table(topTable(fit.cb, coef=1, adjust="fdr",number=49000), file="adeno-HC.txt", sep="\t")
write.table(topTable(fit.cb, coef=2, adjust="fdr",number=49000), file="gas-HC.txt", sep="\t")

#load packages if needed
BiocManager::instal("HumanIDMapping")

#rename the probe ids and use gene symbols for downstream pathway analysis
nuIDs <- rownames(fit.cb$p.value)
mappedIDs <- IlluminaID2nuID(nuIDs, lib.mapping = "lumiHumanIDMapping")[,2] #geneid
#mappedIDs <- nuID2probeID(nuIDs, lib.mapping = "lumiHumanIDMapping") #probeid
genenames <- mappedIDs

rownames(fit.cb$p.value) <- genenames

#pathway analysis using tmod
biocLite("tmod")
library("tmod")

#look at the adeno vs hc genes
res.1 <- tmodLimmaTest(fit.cb, coef=1, genenames, tmodFunc = tmodCERNOtest)
head(tmodSummary(res.1))
str(res.1)

#plot the results
pdf(file="adeno.pdf", height=30, width=8)
tmodPanelPlot(res.1, text.cex = 0.8, plot.cex=1, pval.thr = 0.00001, filter.empty.rows = T, filter.unknown = F, col.labels = colnames(fit.cb$p.value))
dev.off()

#look at the GAS-hc genes and plot the results
res.2 <- tmodLimmaTest(fit.cb, coef=2, genenames, tmodFunc = tmodCERNOtest)

pdf(file="GAS1.pdf", height=30, width=8)
tmodPanelPlot(res.2, text.cex = 0.8, plot.cex=1, pval.thr = 0.00001, filter.empty.rows = T, filter.unknown = F, col.labels = colnames(fit.cb$lods))
dev.off()

#create object for comparative analysis between the two comparisons (adeno vs hc and GAS vs hc)
res.3 <- tmodLimmaTest(fit.cb, genenames, tmodFunc = tmodCERNOtest)

pie <- tmodLimmaDecideTests(fit.cb, genenames, pval.thr = 0.01,lfc.thr=0.5)

#and create plot
pdf(file="tmodpanelpielfc1.pdf", height=20, width=9)
tmodPanelPlot(res.3, pie=pie, pie.style = "pie", text.cex=0.7, filter.empty.rows = T, col.labels.style = "top", col.labels = colnames(fit.cb$p.value))
dev.off()