library('biomaRt')
library('xml2')

listMarts()

ensembl=useMart("ensembl")

datasets <- listDatasets(ensembl)
head(datasets)

ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# filters
filters = listFilters(ensembl)
filters

# attributes
attributes = listAttributes(ensembl)
attributes[1:5,]


###### illumina run ######
searchAttributes(mart = ensembl, pattern = "illum")
searchAttributes(mart = ensembl, pattern = "ensembl.*id")
# ensembl_gene_id
# ensembl_transcript_id

searchFilters(mart = ensembl, pattern = "illum")

getwd()
setwd('~/Documents/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
nrow(illumina)

ill.probs <- illumina$Probe_Id
# ill.probs <- illumina$Probe_Id[1:200]
# ill.probs <- illumina$Probe_Id[10000:19999]
# ill.probs <- illumina$Probe_Id[40000:47323]
ill.probs[1:5]
length(ill.probs)

df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id'), 
      filters = 'illumina_humanht_12_v4', 
      values = ill.probs, 
      mart = ensembl)

df.1
View(df.1)
dim(df.1)
length(unique(df.1$illumina_humanht_12_v4))
length(unique(df.1$hgnc_symbol))

#  NOTE. NOT SURE WHY THIS IS HAPPENING BUT ILL.PROBES WHICH YOU ARE PASSING TO GETBM FUNCTION
# DO NOT APPEAR TO BE RETURNED WHEN YOU LOOK FOR THE ILL PROBES IN THE data.frame
# NEED TO INVESTIGATE THIS UPON RETURN

# write.csv(df.1, file = "ill_hugo_47323.csv", row.names=TRUE)
# ill_hugo_10000
ill_hugo_10000 <- read.table('ill_hugo_10000.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)

head(ill_hugo_10000)

# getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id'), 
#               filters = 'ensembl_gene_id', 
#               values = a, 
#               mart = ensembl)
# 
# a <- c('ENSG00000236172')













