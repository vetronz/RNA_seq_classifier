library('biomaRt')
library('xml2')
library("illuminaHumanv4.db")
library('lumi')
library(lumiHumanIDMapping)

listMarts()
ensembl=useMart("ensembl")

datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

# filters
filters = listFilters(ensembl)
filters[1:5,]

# attributes
attributes = listAttributes(ensembl)
attributes[1:5,]


###### illumina run ######
# searchAttributes(mart = ensembl, pattern = "illum")
# searchAttributes(mart = ensembl, pattern = "ensembl.*id")

# searchFilters(mart = ensembl, pattern = "illum")
# searchFilters(mart = ensembl, pattern = "array")

getwd()
setwd('~/Documents/Masters/RNA_seq_classifier/Data')
illumina <- read.table('ill_probe.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
head(illumina)
length(unique(illumina$Probe_Id))
length(unique(illumina$Array_Address_Id))


setwd('/home/patrick/Code/R')
load('esets.RData')
setwd('~/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')

# idx <- status['most_general'] == 'bacterial' |
#   status['most_general'] == 'viral' |
#   status['most_general'] == 'greyb' |
#   status['most_general'] == 'greyv'|
#   status['most_general'] == 'greyu'
# sum(idx)



e.set.f <- as.data.frame(e.set)
# e.set.f <- as.data.frame(e.set[,idx])
dim(e.set.f)

e.set.f$Probe_Id <- illumina$Probe_Id[match(rownames(e.set.f),illumina$Array_Address_Id)]
e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)]


con <- IlluminaID2nuID(e.set.f$Probe_Id, lib.mapping='lumiHumanIDMapping', chipVersion = NULL)
con.df <- as.data.frame(con)
dim(con.df)
head(con.df)

# check the e.set.f probe id match that in the conversion df
sum(con.df$Probe_Id == e.set.f$Probe_Id)

# check the conversion df if gene and symbol are equivalent = FALSE
sum(as.character(con.df$ILMN_Gene) == as.character(con.df$Symbol))

# reads in the mixture2 file used by Pre to check against
mixture2 <- read.table('mixture2.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(mixture2)

sum(mixture2$Gene == as.character(con.df$ILMN_Gene))
sum(mixture2$Gene == as.character(con.df$Symbol))

length(intersect(mixture2$Gene, as.character(con.df$ILMN_Gene)))
length(intersect(mixture2$Gene, as.character(con.df$Symbol)))
# looks like gene names is closer than gene symbol to the Gene names used by Pre
# note that the sum check does not equal to the set check becuase intersect only counts unique overlap


# check the discrepency between mixture2 and the lumi gene names
which(mixture2$Gene != as.character(con.df$ILMN_Gene))
mixture2$Gene[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]
con.df$ILMN_Gene[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]
# looks like just a capitalization difference

ill.probes <- e.set.f$Probe_Id[which(mixture2$Gene != as.character(con.df$ILMN_Gene))]
ill.probes
length(ill.probes)
# ill.probes <- e.set.f$illumina_probe[1:15000]

# biomaRt check of the illumina probes over which there is discrepency
df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'illumina_humanht_12_v4', 
              values = ill.probes, 
              mart = ensembl)
View(df.1)

# looks like when we do the biomaRt conversion that the hgnc gene names are capitalized similarly to illumina gene names
# however we have some emptyy rows in our biomaRt conversion
# double checking these in our con.df conversion we can see they are not empty

df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == '']
match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id)
con.df[match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id),]

dim(e.set.f)
dim(con.df)
e.set.f$Gene <- con.df$ILMN_Gene
View(cbind(Gene=con.df$ILMN_Gene, e.set.f))
View(cbind(Gene=mixture2$Gene, e.set.f))

# rip the gene column from the con.df / mixture2 and stick it on my version as col1
e.set.f$Probe_Id <- NULL
e.set.r <- round(e.set.f,5)

e.set.out <- cbind(Gene=con.df$ILMN_Gene, e.set.r)
rownames(e.set.out) <- NULL
dim(e.set.out)
View(head(e.set.out))

# saving and reading 
# setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')
# write.table(mixture2, file = "mixture2_rep2.txt", sep='\t', row.names=FALSE, quote = FALSE) # worked!
write.table(e.set.out, file = "eset1.txt", sep='\t', row.names=FALSE, quote = FALSE) # worked!!!!!!!!!!




# sql
# mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19
# select distinct G.gene,N.value from ensGtp as G, ensemblToGeneName as N where G.transcript=N.name and G.gene in ("ENSG00000183742", "ENSG00000121410") ;











