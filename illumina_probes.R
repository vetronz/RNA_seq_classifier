getwd()
setwd('~/Downloads')
illumina <- read.table('HumanHT-12_V4_0_R2_15002873_B.txt', sep = '\t', skip = 8, stringsAsFactors = FALSE, fill = TRUE, header = TRUE)

head(illumina)
nrow(illumina)
length(unique(illumina$Array_Address_Id))

sort(table(illumina$Array_Address_Id),decreasing=TRUE, includeNA="always")[1:10]

which(table(illumina$Array_Address_Id, useNA="always") > 1)
which(is.na(illumina$Array_Address_Id))
illumina[113:116,]

head(illumina)

table(as.numeric(colnames(e.set.f)) %in% illumina$Array_Address_Id)

dim(illumina)
a
a %in% illumina$Array_Address_Id





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

library("illuminaHumanv4.db")
x <- illuminaHumanv4CHR


x <- illuminaHumanv4NUID
x <- illuminaHumanv4ALIAS2PROBE
x <- illuminaHumanv4ENSEMBL
x <- illuminaHumanv4GENENAME
x <- illuminaHumanv4GO
x <- illuminaHumanv4MAP
x <- illuminaHumanv4REFSEQ
data.frame(Gene=unlist(mget(x = probeID,envir = x)))



probeID=c("ILMN_1690170", "ILMN_2410826", "ILMN_1675640", "ILMN_1801246",
          "ILMN_1658247", "ILMN_1740938", "ILMN_1657871", "ILMN_1769520",
          "ILMN_1778401")

















#
