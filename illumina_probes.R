library('biomaRt')

getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
nrow(illumina)

getwd()
setwd('/home/patrick/Code/R')
# setwd('/Users/patrickhedley-miller/code/gitWorkspace/infxRNAseq')

# rm(list=setdiff(ls(), 'all'))
load('esets.RData')


# trans <- c(7650358, 5720482, 2570300, 3180392, 5090754)
# sum(illumina$Array_Address_Id %in% trans)

probe.id <- illumina$Probe_Id


probeID <- illumina$Probe_Id[which(illumina$Array_Address_Id %in% trans)]

df.1$Gene <- as.character(df.1$Gene)

conversion <- illuminaHumanv4ENSEMBL
# conversion <- illuminaHumanv4GENENAME
# x <- illuminaHumanv4CHR
# x <- illuminaHumanv4NUID
# x <- illuminaHumanv4ALIAS2PROBE
# x <- illuminaHumanv4GO
# x <- illuminaHumanv4MAP
# x <- illuminaHumanv4REFSEQ

con <- mget(x = probe.id, envir = conversion)

lapply(con[670:680], function(x){
  paste0('the probe is: ', x)
})

unlist(lapply(con, function(x){x[1]}))
unlist(lapply(con, function(x){x[2]}))
con[1:100]
illumina$ensemb <- as.character(con)


lapply(con, function(x){
  x[1]
}) %>% unlist
lapply(con, function(x){
  if(length(x)>=2){
    x[2]
  } else {
    NA
  }}) %>% unlist



get_element <- function(l, n){
  lapply(l, function(elem, n){
    if(length(elem)>=n){
      elem[n]
    } else {
      NA
    }}, n) %>% unlist  
}


View(illumina)

ENSG00000137959

df.1 <- data.frame(Gene=unlist(mget(x = probeID, envir = x)))
View(df.1)


# test 
# ILMN_1770772
# ENSG00000092009
df.1[which(rownames(df.1) == "ILMN_1770772"),]
rownames(df.1)[which(df.1$Gene == 'ENSG00000092009')]

rownames(df.1)[which(df.1$Gene == 'ENSG00000137959')]
# "ILMN_1835092" "ILMN_1723912"

df.1[which(rownames(df.1) == "ILMN_1835092"),]
df.1[which(rownames(df.1) == "ILMN_1723912"),]

# ILMN_1729749 3180392 ENSG00000138646
# Q9UII4 (HERC5_HUMAN)
# induced by bacterial lipopolysaccharides (LPS) found on outer membrane of gram -ve



getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
library(readxl)
# read_excel reads both xls and xlsx files
# read_excel("my-old-spreadsheet.xls")
a <- read_excel("Mega_iris_probe_comp.xlsx")

e.set[1:10,1:5]

a[1:10,]
View(a[1:10,])
dim(a)
class(a[[1]])

'ILMN_2650615' %in% a[[1]]
'ILMN_1767362' %in% a[[1]]

head(b)
b<-substr(b, 6, nchar(b))
b
gram.hits
gram.hits %in% b
# conda install -c bioconda bioconductor-biomart

colnames(e.set.f[1]) %in% b


















#
