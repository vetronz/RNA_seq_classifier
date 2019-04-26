setwd('/Users/patrickhedley-miller/code/R/infxRNAseq')
#setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')
#getwd()


load('esets.RData')
rm(list=setdiff(ls(), "x"))

class(e.set)
e.set[1:5,1:5]
dim(e.set)
dim(status)

dim(e.set.i)
dim(status.iris)

class(status)
dim(status)

dim(status)
class(status)
status[1:4,1:4]
attributes(status)

status$Diagnosis
status$category

e.set.t <- t(e.set)
e.set.i.t <- t(e.set.i)

View(status.iris)

labels.i <- rownames(e.set.i.t)
labels.s <- status.iris$My_code
common <- match(label.i, labels.s)
common

as.character(status.iris$My_code[common])
as.character(status.iris$most_general[common])

e.set.i.t[1:10,1:3]

# my_code provides the index in the e.set dataframes
# we can use this to check you hvae pulled the correct labels from the status df
# thinks looking good preliminary
# next step is to essentially redo yesterdays work cheking pca for batch effect with the new data.