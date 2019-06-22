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

idx <- status['most_general'] == 'bacterial' |
  status['most_general'] == 'viral' |
  status['most_general'] == 'greyb' |
  status['most_general'] == 'greyv'|
  status['most_general'] == 'greyu'
sum(idx)

e.set.f <- as.data.frame(e.set[,idx])
dim(e.set.f)

# rownames(e.set.f)[1:10]

e.set.f$illumina_probe <- illumina$Probe_Id[match(rownames(e.set.f),illumina$Array_Address_Id)]
e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)]


# testing
IlluminaID2nuID(e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)][5], lib.mapping='lumiHumanIDMapping', chipVersion = NULL)

# Search_Key    ILMN_Gene Accession     Symbol Probe_Id       Array_Address_Id nuID                
# 6450255 "NM_182762.2" "7A5"     "NM_182762.2" "7A5"  "ILMN_1762337" "0006450255"     "Ku8QhfS0n_hIOABXuE"
# 2570615 "NM_130786.2" "A1BG"    "NM_130786.2" "A1BG" "ILMN_2055271" "0002570615"     "fqPEquJRRlSVSfL.8A"
# 6370619 "NM_130786.2" "A1BG"    "NM_130786.2" "A1BG" "ILMN_1736007" "0006370619"     "ckiehnugOno9d7vf1Q"
# 2600039 "NM_138932.1" "A1CF"    "NM_138932.1" "A1CF" "ILMN_2383229" "0002600039"     "x57Vw5B5Fbt5JUnQkI"
# 2650615 "NM_138933.1" "A1CF"    "NM_014576.2" "A1CF" "ILMN_1806310" "0002650615"     "ritxUH.kuHlYqjozpE"

nuID2IlluminaID('Ku8QhfS0n_hIOABXuE')
nuID2probeID('Ku8QhfS0n_hIOABXuE')

nuID2EntrezID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')
nuID2RefSeqID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')
nuID2targetID('Ku8QhfS0n_hIOABXuE', lib.mapping='lumiHumanIDMapping')



# real
con <- IlluminaID2nuID(e.set.f$illumina_probe, lib.mapping='lumiHumanIDMapping', chipVersion = NULL)
con.df <- as.data.frame(con)
dim(con.df)
head(con.df)

# check the e.set.f probe id match that in the conversion df
sum(con.df$Probe_Id == e.set.f$illumina_probe)

# check the conversion df if gene and symbol are equivalent = FALSE
sum(as.character(con.df$ILMN_Gene) == as.character(con.df$Symbol))

sum(mix2$Gene == as.character(con.df$ILMN_Gene))
sum(mix2$Gene == as.character(con.df$Symbol))
# looks like the gene name is closer to that used by pre in the mixture2

which(mix2$Gene != as.character(con.df$ILMN_Gene))
mix2$Gene[which(mix2$Gene != as.character(con.df$ILMN_Gene))]
con.df$ILMN_Gene[which(mix2$Gene != as.character(con.df$ILMN_Gene))]

ill.probes <- e.set.f$illumina_probe[which(mix2$Gene != as.character(con.df$ILMN_Gene))]

ill.probes
length(ill.probes)
# ill.probes <- e.set.f$illumina_probe[1:15000]

df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol', 'ensembl_gene_id'), 
              filters = 'illumina_humanht_12_v4', 
              values = ill.probes, 
              mart = ensembl)
View(df.1)

# looks like when we do the biomaRt conversion that the hgnc gene names are capitalized just like our con.df ilmn_gene
# however we have some emptyy rows in our biomaRt conversion
# double checking these in our con.df conversion we can see they are not empty

df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == '']
match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id)
con.df[match(df.1$illumina_humanht_12_v4[df.1$hgnc_symbol == ''], con.df$Probe_Id),]

dim(e.set.f)
dim(con.df)
e.set.f$Gene <- con.df$ILMN_Gene

# visual check that e.set.f array, probe and genes match con.df values
con.df[1:5,]
e.set.f[1:5,(ncol(e.set.f)-4):ncol(e.set.f)]

e.set.f$illumina_probe <- NULL
e.set.f$Gene
colnames(e.set.f)
View(e.set.f)
dim(e.set.f)

cols.names <- colnames(e.set.f)
cols.names <-c('Gene', cols.names)
cols.reorder <- cols.names[-length(cols.names)]

setwd('/home/patrick/Documents/Masters/RNA_seq_classifier/Data/Ciber_sort')
mix_one <- e.set.f[, cols.reorder]
dim(mix_one)
# write.table(mix_one, file = "mix_one.txt", ,sep='\t', row.names=TRUE)
mix_test <- read.table('mix_one.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)

dim(mix_test)
View(mix_test)


# pre mixture

mix2 <- read.table('mixture2.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(mix2)
View(mix2)

LM22.ref <- read.table('LM22-ref-sample.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
dim(LM22.ref)
View(LM22.ref)


###############################################################################################
# ill.probs <- illumina$Probe_Id
# ill.probs <- illumina$Probe_Id[1:20000]
# ill.probs <- illumina$Probe_Id[10000:19999]
# ill.probs <- illumina$Probe_Id[40000:47323]
# ill.probs
# length(ill.probs)
# 
# df.1 <- getBM(attributes=c('illumina_humanht_12_v4', 'hgnc_symbol'), 
#       filters = 'illumina_humanht_12_v4', 
#       values = ill.probs, 
#       mart = ensembl)
# 
# df.1
# View(df.1)
# dim(df.1)
# length(unique(df.1$illumina_humanht_12_v4))
# length(unique(df.1$hgnc_symbol))



# ill_hugo_10000 <- read.table('ill_hugo_10000.csv', sep = ',', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
# dim(ill_hugo_10000)



mysql  --user=genome --host=genome-mysql.cse.ucsc.edu -A -D hg19
select distinct G.gene,N.value from ensGtp as G, ensemblToGeneName as N where G.transcript=N.name and G.gene in ("ENSG00000183742", "ENSG00000121410") ;












# 
# LM22_ref <- read.table('LM22-ref-sample.txt', sep = '\t', stringsAsFactors = FALSE, fill = FALSE, header = TRUE)
# View(LM22_ref)
# lm<-LM22_ref$Relabel
# length(lm)
# 
# df.3 <- getBM(attributes=c('hgnc_symbol', 'illumina_humanht_12_v4'), 
#               filters = 'hgnc_symbol', 
#               values = lm, 
#               mart = ensembl)
# dim(df.3)
# 
# View(df.3)
# 


