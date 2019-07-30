library(Biobase)
library(GEOquery)
# library(ggplot2)
# library(tidyr)
# # library(rsample)
# library(caret)
# # library(h2o)
# library(dplyr)
# library(tidyverse)  # data manipulation and visualization
# library(modelr)     # provides easy pipeline modeling functions
# library(broom)

## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("inSilicoMerging")
# errors, will need to downgrade R version if i want to use db merge
# possible as per:
# https://www.r-bloggers.com/installation-of-r-3-5-on-ubuntu-18-04-lts-and-tips-for-spatial-packages/
# but it will be a pain


# GET GEO RETURNS A LIST
# discovery: GEO accession: GSE72809
# gset.d <- getGEO('GSE72809')

# IRIS VALIDATION: GEO accession: GSE72810
# gset.v <- getGEO('GSE72810')


### SAVING
# saveRDS(gset.d, file = 'gset_GSE72809')
# saveRDS(gset.d, file = 'gset_GSE72810')

## LOADING
setwd('/home/patrick/Code/Gitworkspace/RNA_seq_classifier/Datasets')
gset <-readRDS(file = "gset_XYZ")


# extract the expression set using exprs method. note double [[]] to index into the list
gset.d.df <- exprs(gset.d[[1]])
gset.v.df <- exprs(gset.v[[1]])

dim(gset.d.df)
dim(gset.v.df)

gset.d.df[1:5,1:5]
gset.v.df[1:5,1:5]


a <- phenoData(gset.d[[1]])
attributes(a)
df <- attr(a, 'data')
dim(df)

# pData method is equivalent to passing data attr to phenoData
dim(pData(gset.d[[1]]))


class(gset.d[1])
class(gset.d[[1]])
gset.d[[1]]
phenoData(gset.d[1])

colnames(phenoData(gset[[1]]))
phenoData(gset[[1]])$'category:ch1'
sum(phenoData(gset[[1]])$'category:ch1' != 'Control') # correct number after removing outlier



