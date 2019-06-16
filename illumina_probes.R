library("illuminaHumanv4.db")
getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
nrow(illumina)

module.assign[module.assign == 2]
# 5720482 2570300 2100196  990768 3360343
trans <- c(5720482, 2570300, 2100196,  990768, 3360343)
trans <- c(7650358, 5720482, 2570300, 3180392, 5090754)
trans <- c(830440)
which(illumina$Array_Address_Id %in% trans)

probeID <- illumina$Probe_Id[which(illumina$Array_Address_Id %in% trans)]

# 830440                     STAM binding protein like 1 ENSG00000138134

x <- illuminaHumanv4CHR
x <- illuminaHumanv4NUID
x <- illuminaHumanv4ALIAS2PROBE
x <- illuminaHumanv4ENSEMBL
x <- illuminaHumanv4GENENAME
x <- illuminaHumanv4GO
x <- illuminaHumanv4MAP
x <- illuminaHumanv4REFSEQ
data.frame(Gene=unlist(mget(x = probeID, envir = x)))



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
